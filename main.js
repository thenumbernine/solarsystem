import {vec3, mat4, quat, glMatrix} from '/js/gl-matrix-3.4.1/index.js';
glMatrix.setMatrixArrayType(Array);	//use double rather than float precision with gl-matrix
import {quatZAxis} from '/js/gl-util.js';
import {Br, Hr, Div} from '/js/dom.js';
import {show, hide, assert, animate, mathClamp, mathDeg} from '/js/util.js';
import {ids, cfg, urlparams, floatToGLSL} from './globals.js';
import {SmallBody, SmallBodies} from './smallbodies.js';
import {Galaxy, Galaxies} from './galaxies.js';
import {PickObject} from './pickobject.js';
import {Mouse3D} from '/js/mouse3d.js';
import {StarSystem} from './starsystem.js';
import {Planet} from './planet.js';

import {overlayTexts} from './overlaytexts.js';
import {metersPerUnits} from './units.js';
import {ui} from './ui.js';
import {starSystemsExtra} from './starsystems.js';
import {calcTides} from './calctides.js';
import {milkyWay} from './milkyway.js';
import {starfield} from './starfield.js';
import {skyCube} from './skycube.js';
import {
	solarSystemBarycentricToPlanetGeodetic,
	planetGeodeticToSolarSystemBarycentric,
} from './vec.js';

const orbitPathIndexForType = {elliptic:0, hyperbolic:1, parabolic:2};	//used in the shader

// track ball motion variables
let orbitGeodeticLocation;
let orbitDistance;
const orbitOffset = [0,0,0];
const orbitZoomFactor = .0003;	// upon mousewheel

let mouseDir;
let smallBodies;
let galaxies;

let resetDistanceOnSelectingNewSystem = false;

let colorShader;

let pointObj;
let lineObj;
let quadObj;

let planetAmbientLight = .1;

function clearGeodeticLocation() {
	const glutil = ui.glutil;
	if (orbitGeodeticLocation === undefined) return;

	orbitDistance = cfg.orbitTargetDistance = cfg.orbitTarget.radius * .2;

	glutil.view.fovY = 90;
	glutil.updateProjection();

	//...and in one fell swoop, restore the original camera angle
	//TODO spread this out over a few frames
	quat.copy(glutil.view.angle, orbitGeodeticLocation.lastViewAngle);

	orbitGeodeticLocation = undefined;
}

function resetOrbitViewPos() {
	const glutil = ui.glutil;
	let viewAngleZAxis = vec3.create();
	quatZAxis(viewAngleZAxis, glutil.view.angle);
	vec3.scale(glutil.view.pos, viewAngleZAxis, orbitDistance + (cfg.orbitTarget.equatorialRadius || cfg.orbitTarget.radius || 0));
	vec3.add(glutil.view.pos, glutil.view.pos, orbitOffset);
}

cfg.refreshMeasureText = function() {
	ids.measureMin.innerText = cfg.orbitTarget.measureMin === undefined ? '' : (cfg.orbitTarget.measureMin.toExponential() + ' m/s^2');
	ids.measureMax.innerText = cfg.orbitTarget.measureMax === undefined ? '' : (cfg.orbitTarget.measureMax.toExponential() + ' m/s^2');
}

let gravityWellZScale = 1;
let gravityWellTargetZScale = 1;
let planetPointVisRatio = .001;
let showFPS = false;
{
	let delta = vec3.create();//[];
	let viewAngleInv = quat.create();
	let invRotMat = mat4.create();
	let viewPosInv = vec3.create();
	let viewfwd = vec3.create();

	//fps counter
	let frames = 0;
	let lastFPSTime = Date.now();

	//used for new gpu update of tide tex
	let updatePlanetStateBuffer = new Float32Array(1);

	cfg.drawScene = function(picking) {
		const glutil = ui.glutil;
		const gl = ui.gl;

		//should picking count towards fps? nah
		if (!picking) {
			frames++;
			let thisTime = Date.now();
			if (thisTime - lastFPSTime > 1000) {
				let fps = frames * 1000 / (thisTime - lastFPSTime);
				if (showFPS) {
					console.log('fps '+fps);
				}
				frames = 0;
				lastFPSTime = thisTime;
			}

			cfg.flatEarthCoeff += cfg.flatEarthConvCoeff * (cfg.targetFlatEarthCoeff - cfg.flatEarthCoeff);
		}

		if (cfg.overlayShowCurrentPosition &&
			cfg.orbitTarget instanceof Planet
		) {
			//update the min and max to reflect the current position
			let x = glutil.view.pos[0] + cfg.orbitTarget.pos[0];
			let y = glutil.view.pos[1] + cfg.orbitTarget.pos[1];
			let z = glutil.view.pos[2] + cfg.orbitTarget.pos[2];
			let t = cfg.calcMetricForce([x,y,z], cfg.orbitTarget);
			ids.measureMin.innerText = t === undefined ? '' : t.toExponential() + ' m/s^2';
			ids.measureMax.innerText = t === undefined ? '' : t.toExponential() + ' m/s^2';
		}

		//update the planet state texture
		// right now it's only used for tide calcs so
		//1) only update it if tide calcs are on
		//2) only include planets that are enabled for tide calcs
		// more discernment later when i make it general purpose
		let useOverlay = cfg.overlayShowOrbitTarget && cfg.displayMethod != 'None';
		if (useOverlay && !picking) {

			let targetSize = cfg.orbitStarSystem.planetStateTex.width * cfg.orbitStarSystem.planetStateTex.height * 4;

			if (updatePlanetStateBuffer.length != targetSize) {
				updatePlanetStateBuffer = new Float32Array(targetSize);
			}

			//gather pos and mass
			for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
				let planet = cfg.orbitStarSystem.planets[planetIndex];
				updatePlanetStateBuffer[0 + 4 * planetIndex] = planet.pos[0];
				updatePlanetStateBuffer[1 + 4 * planetIndex] = planet.pos[1];
				updatePlanetStateBuffer[2 + 4 * planetIndex] = planet.pos[2];

				if (!cfg.planetInfluences[planetIndex]) {
					//if we're not using this planet then set the mass to zero
					// works the same as not using it
					updatePlanetStateBuffer[3 + 4 * planetIndex] = 0;
				} else {
					updatePlanetStateBuffer[3 + 4 * planetIndex] =
						(planet.mass === undefined ? 0 : planet.mass)
						* 1e-18;	//shader mass precision bias
				}
			}

			//update orbit star system's planet state tex
			cfg.orbitStarSystem.planetStateTex.bind();

			gl.texSubImage2D(
				gl.TEXTURE_2D,
				0,
				0,
				0,
				cfg.orbitStarSystem.planetStateTex.width,
				cfg.orbitStarSystem.planetStateTex.height,
				gl.RGBA,
				gl.FLOAT,
				updatePlanetStateBuffer);

			cfg.orbitStarSystem.planetStateTex.unbind();
		}

		//get our calculations for orientation

		quat.conjugate(viewAngleInv, glutil.view.angle);
		mat4.fromQuat(invRotMat, viewAngleInv);

		mat4.copy(glutil.scene.mvMat, invRotMat);

		//TODO pull from matrix
		quatZAxis(viewfwd, glutil.view.angle);
		vec3.scale(viewfwd, viewfwd, -1);

		let tanFovY = Math.tan(glutil.view.fovY * Math.PI / 360);

		let distFromSolarSystemInM = vec3.length(glutil.view.pos);
		let distFromSolarSystemInLyr = distFromSolarSystemInM / metersPerUnits.lyr;
		let distFromSolarSystemInPc = distFromSolarSystemInM / metersPerUnits.pc;
		let distFromSolarSystemInMpc = distFromSolarSystemInM / metersPerUnits.Mpc;


		//draw from furthest to nearest, so the varying-scaled objects don't get conflicting depth information
		// be sure to clear the depth buffer between each scale

		//first comes the sky cube (with no depth information)

		skyCube.draw(
			picking,
			distFromSolarSystemInLyr
		);

		//next comes galaxies

		//next comes the milky way
		//render the milkyWay sceneObj in Mpc units rather than meters (because it's hitting the limit of fp accuracy)
		// set up Mpc scaled view
		vec3.scale(viewPosInv, glutil.view.pos, -1/cfg.interGalacticRenderScale);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);

		milkyWay.draw(
			tanFovY,
			picking,
			distFromSolarSystemInLyr,
			distFromSolarSystemInMpc
		);

		galaxies.draw(
			tanFovY,
			picking,
			distFromSolarSystemInMpc
		);

		//next comes local stars of the milky way

		starfield.draw(
			tanFovY,
			picking,
			viewPosInv,
			invRotMat,
			distFromSolarSystemInLyr,
			distFromSolarSystemInMpc
		);

		//last is planets

		gl.clear(gl.DEPTH_BUFFER_BIT)
		vec3.scale(viewPosInv, glutil.view.pos, -1);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);

		//flat earth
		const solarSystem = starSystemsExtra.solarSystem;
		const earth = solarSystem.planets[solarSystem.indexes.Earth];
assert(earth);
		cfg.flatEarthRelativeEarthPos = [];
		vec3.sub(cfg.flatEarthRelativeEarthPos, earth.pos, cfg.orbitTarget.pos);
		cfg.flatEarthRelativeEarthNorthDir = [0,0,1];
		quatZAxis(cfg.flatEarthRelativeEarthNorthDir, earth.angle);

		//update all our sceneObjs' flatEarth uniforms here
		lineObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
		lineObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
		lineObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;
		pointObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
		pointObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
		pointObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;

		//draw debug lines
		if (!picking) {
			for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
				let planet = cfg.orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;
				if (planet.pos === undefined) continue;	//comets don't have pos yet, but I'm working on that
				if (cfg.orbitTarget.pos === undefined) continue;

				/*
				if (planet.isComet) continue;
				if (planet.isAsteroid) continue;
				if (planet.parent !== solarSystem.planets[solarSystem.indexes.Sun] &&
					planet !== solarSystem.planets[solarSystem.indexes.Sun] &&
					planet !== solarSystem.planets[solarSystem.indexes.Moon]) continue;
				*/


				//draw lines to other planets
				if (cfg.showLinesToOtherPlanets &&
					cfg.orbitStarSystem == solarSystem &&
					cfg.orbitTarget !== planet &&
					(!planet.parent || planet.parent.index == 0))
				{

					//while here, update lines

					vec3.sub(delta, planet.pos, cfg.orbitTarget.pos);

					let dist = vec3.length(delta);
					vec3.scale(delta, delta, 1/dist);

					lineObj.attrs.vertex.buffer.data[0] = delta[0] * cfg.orbitTarget.radius;
					lineObj.attrs.vertex.buffer.data[1] = delta[1] * cfg.orbitTarget.radius;
					lineObj.attrs.vertex.buffer.data[2] = delta[2] * cfg.orbitTarget.radius;
					lineObj.attrs.vertex.buffer.data[3] = delta[0] * orbitDistance;
					lineObj.attrs.vertex.buffer.data[4] = delta[1] * orbitDistance;
					lineObj.attrs.vertex.buffer.data[5] = delta[2] * orbitDistance;

					lineObj.attrs.vertex.buffer.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show velocity vectors
				if (cfg.showVelocityVectors) {
					vec3.sub(delta, planet.pos, cfg.orbitTarget.pos);
					lineObj.attrs.vertex.buffer.data[0] = delta[0];
					lineObj.attrs.vertex.buffer.data[1] = delta[1];
					lineObj.attrs.vertex.buffer.data[2] = delta[2];
					lineObj.attrs.vertex.buffer.data[3] = delta[0] + planet.vel[0] * cfg.velocityVectorScale;
					lineObj.attrs.vertex.buffer.data[4] = delta[1] + planet.vel[1] * cfg.velocityVectorScale;
					lineObj.attrs.vertex.buffer.data[5] = delta[2] + planet.vel[2] * cfg.velocityVectorScale;
					lineObj.attrs.vertex.buffer.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show rotation axis
				if (cfg.showRotationAxis) {
					vec3.sub(delta, planet.pos, cfg.orbitTarget.pos);
					let axis = [0,0,1];
					quatZAxis(axis, planet.angle);
					lineObj.attrs.vertex.buffer.data[0] = delta[0] + axis[0] * 2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[1] = delta[1] + axis[1] * 2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[2] = delta[2] + axis[2] * 2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[3] = delta[0] + axis[0] * -2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[4] = delta[1] + axis[1] * -2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[5] = delta[2] + axis[2] * -2 * planet.radius;
					lineObj.attrs.vertex.buffer.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show orbit axis
				if (cfg.showOrbitAxis && planet.orbitAxis) {
					vec3.sub(delta, planet.pos, cfg.orbitTarget.pos);
					lineObj.attrs.vertex.buffer.data[0] = delta[0] + planet.orbitAxis[0] * 2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[1] = delta[1] + planet.orbitAxis[1] * 2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[2] = delta[2] + planet.orbitAxis[2] * 2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[3] = delta[0] + planet.orbitAxis[0] * -2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[4] = delta[1] + planet.orbitAxis[1] * -2 * planet.radius;
					lineObj.attrs.vertex.buffer.data[5] = delta[2] + planet.orbitAxis[2] * -2 * planet.radius;
					lineObj.attrs.vertex.buffer.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show ellipses major and minor axis
				if (cfg.showEllipseAxis &&
					planet.keplerianOrbitalElements &&
					planet.keplerianOrbitalElements.A &&
					planet.keplerianOrbitalElements.B)
				{
					if (planet.parent) {
						vec3.sub(delta, planet.parent.pos, cfg.orbitTarget.pos);
					} else {
						vec3.copy(delta, cfg.orbitTarget.pos);
					}
					lineObj.attrs.vertex.buffer.data[0] = delta[0] - planet.keplerianOrbitalElements.A[0];
					lineObj.attrs.vertex.buffer.data[1] = delta[1] - planet.keplerianOrbitalElements.A[1];
					lineObj.attrs.vertex.buffer.data[2] = delta[2] - planet.keplerianOrbitalElements.A[2];
					lineObj.attrs.vertex.buffer.data[3] = delta[0] + planet.keplerianOrbitalElements.A[0];
					lineObj.attrs.vertex.buffer.data[4] = delta[1] + planet.keplerianOrbitalElements.A[1];
					lineObj.attrs.vertex.buffer.data[5] = delta[2] + planet.keplerianOrbitalElements.A[2];
					lineObj.attrs.vertex.buffer.updateData();
					lineObj.draw({uniforms : { color : planet.color }});

					if (planet.parent) {
						vec3.sub(delta, planet.parent.pos, cfg.orbitTarget.pos);
					} else {
						vec3.copy(delta, cfg.orbitTarget.pos);
					}
					lineObj.attrs.vertex.buffer.data[0] = delta[0] - planet.keplerianOrbitalElements.B[0];
					lineObj.attrs.vertex.buffer.data[1] = delta[1] - planet.keplerianOrbitalElements.B[1];
					lineObj.attrs.vertex.buffer.data[2] = delta[2] - planet.keplerianOrbitalElements.B[2];
					lineObj.attrs.vertex.buffer.data[3] = delta[0] + planet.keplerianOrbitalElements.B[0];
					lineObj.attrs.vertex.buffer.data[4] = delta[1] + planet.keplerianOrbitalElements.B[1];
					lineObj.attrs.vertex.buffer.data[5] = delta[2] + planet.keplerianOrbitalElements.B[2];
					lineObj.attrs.vertex.buffer.updateData();
					lineObj.draw({uniforms : {
						color : [
							planet.color[0] * .3,
							planet.color[1] * .3,
							planet.color[2] * .3,
							1
						]
					}});
				}
			}
		}

		//update planet vis ratio
		//for picking do we need to? TODO make sure the view doesn't change between the last render and a pick
		for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
			let planet = cfg.orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			//update vis ratio
			let dx = planet.pos[0] - glutil.view.pos[0] - cfg.orbitTarget.pos[0];
			let dy = planet.pos[1] - glutil.view.pos[1] - cfg.orbitTarget.pos[1];
			let dz = planet.pos[2] - glutil.view.pos[2] - cfg.orbitTarget.pos[2];
			//approximated pixel width with a fov of 90 degrees
			planet.visRatio = cfg.planetScaleExaggeration * planet.radius / (Math.sqrt(dx * dx + dy * dy + dz * dz) * tanFovY);
		}

		//draw sphere planets
		for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
			let planet = cfg.orbitStarSystem.planets[planetIndex];

			if (planet.sceneObj) {

				let canSee = !planet.hide && planet.visRatio >= planetPointVisRatio;

				//if the planet is visible then
				// if there is no texture then
				//  load the texture and bind it
				if (canSee) {
					//if we have no image ... then load it
					if (!planet.img) {
						if (planet.imgURL) {
							planet.imgIsLoading = true;
							planet.img = new Image();
							planet.img.addEventListener('load', e => {
								//console.log('finished loading tex for planet '+planet.name);
								planet.imgIsLoading = undefined;
								delete planet.imgIsLoading;
							});
							planet.img.addEventListener('error', e => {
								console.log('failed to find texture for planet '+planet.name);
								planet.imgIsLoading = undefined;
							});
							//console.log('loading planet '+planet.name+' tex '+planet.imgURL);
							planet.img.src = planet.imgURL;
						}
					//if we do have an image ... see if it's done yet
					} else if (planet.imgIsLoading) {
					//if it's done ... use it
					} else if (!planet.tex) {
						//console.log('done loading planet '+planet.name+', now uploading to GPU');
						planet.tex = new glutil.Texture2D({
							flipY : true,
							data : planet.img,
							minFilter : gl.LINEAR_MIPMAP_LINEAR,
							magFilter : gl.LINEAR,
							generateMipmap : true,
						});
					}
				} else {
				//if the planet isn't visible then
				// if there is a texture then
				//  unload it
					if (planet.tex) {
						//console.log('unloading planet '+planet.name+' from the GPU');
						gl.deleteTexture(planet.tex.obj);
						planet.tex = undefined;
						delete planet.tex;
					}
				}

				if (canSee) {
					planet.updateSceneObj();

					//update scene object
					//don't forget one is shared among all planets
					vec3.sub(planet.sceneObj.uniforms.pos, planet.pos, cfg.orbitTarget.pos);
					quat.copy(planet.sceneObj.uniforms.angle, planet.angle);

					//only allow ring obj if there's a radii and a sphere
					//... this way i can just copy over the pos and angle
					if (planet.ringObj !== undefined) {
						vec3.sub(planet.ringObj.pos, planet.pos, cfg.orbitTarget.pos);
						quat.copy(planet.ringObj.angle, planet.sceneObj.uniforms.angle);
					}

					//calculate sun position for lighting
					for (let starIndex = 0; starIndex < cfg.orbitStarSystem.stars.length; ++starIndex) {
						let star = cfg.orbitStarSystem.stars[starIndex];
						let tmp = [];
						vec3.sub(tmp, star.pos, planet.pos);
						planet.sceneObj.uniforms.sunDir[0+3*starIndex] = tmp[0];
						planet.sceneObj.uniforms.sunDir[1+3*starIndex] = tmp[1];
						planet.sceneObj.uniforms.sunDir[2+3*starIndex] = tmp[2];
					}

					//webkit bug
					if (!picking) {
						planet.sceneObj.shader.use();
						gl.uniform3fv(
							gl.getUniformLocation(
								planet.sceneObj.shader.obj, 'sunDir[0]'
							), planet.sceneObj.uniforms.sunDir);
					}

					//update ellipsoid parameters
					planet.sceneObj.uniforms.equatorialRadius = planet.equatorialRadius !== undefined ? planet.equatorialRadius : planet.radius;
					planet.sceneObj.uniforms.inverseFlattening = planet.inverseFlattening !== undefined ? planet.inverseFlattening : .5;
					if (planet.ringObj !== undefined) {
						planet.sceneObj.uniforms.ringMinRadius = planet.ringRadiusRange[0];
						planet.sceneObj.uniforms.ringMaxRadius = planet.ringRadiusRange[1];
					}
					planet.sceneObj.uniforms.color = planet.color;

					calcTides.drawPlanet(planet);

					//TODO FIXME something is overriding this between the update calc and here
					// making me need to re-assign it ...
					planet.sceneObj.uniforms.forceMin = planet.forceMin;
					planet.sceneObj.uniforms.forceMax = planet.forceMax;

					planet.sceneObj.uniforms.ambient = planet.type == 'star' ? 1 : planetAmbientLight;
					planet.sceneObj.uniforms.scaleExaggeration = cfg.planetScaleExaggeration;

					//flat earth
					planet.sceneObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
					planet.sceneObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
					planet.sceneObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;


					if (picking) {
						cfg.pickObject.drawPlanet(planet.sceneObj, planet);
					} else {
						if (planet.type != 'barycenter') {
							planet.sceneObj.draw();
							overlayTexts.add(planet);
						}
					}

					//show latitude and longitude lines
					if (cfg.showLatAndLonLines && !picking) {
						vec3.copy(planetLatLonObj.pos, planet.sceneObj.uniforms.pos);
						quat.copy(planetLatLonObj.angle, planet.sceneObj.uniforms.angle);
						planetLatLonObj.uniforms.equatorialRadius = planet.equatorialRadius !== undefined ? planet.equatorialRadius : planet.radius;
						planetLatLonObj.uniforms.inverseFlattening = planet.inverseFlattening !== undefined ? planet.inverseFlattening : .5;
						planetLatLonObj.uniforms.scaleExaggeration = cfg.planetScaleExaggeration;
						planetLatLonObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
						planetLatLonObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
						planetLatLonObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;
						planetLatLonObj.draw();
					}
				}
			}
		}

		//draw rings and transparent objects last
		//disable depth writing so orbits are drawn in front and behind them
		//do this last so (without depth writing) other planets in the background don't show up in front of the rings
		if (!picking) {
			for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
				let planet = cfg.orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;

				if (planet.sceneObj &&
					planet.visRatio >= planetPointVisRatio &&
					planet.ringObj !== undefined)
				{
					gl.disable(gl.CULL_FACE);
					gl.depthMask(false);
					//to provide to the ring shader:
					//1) vector from planet to the sun (for fwd-vs-back scattering)
					//2) radius and eccentricity of planet (for self-shadows)
					//3) axis of planet rotation / normal of ring plane (for unlit side) ... can be derived from the angle, which is copied into the ring object

					//TODO cache rotation axis -- for here and for showing axis (like we already do for orbit axis)
					let axis = [];
					let sunDir = [];
					vec3.sub(sunDir, planet.pos, cfg.orbitStarSystem.planets[cfg.orbitStarSystem.indexes.Sun].pos);
					vec3.normalize(sunDir, sunDir);
					quatZAxis(axis, planet.ringObj.angle);
					let axisDotSun = vec3.dot(axis, sunDir);
					let lookingAtLitSide = vec3.dot(viewfwd, axis) * axisDotSun;
					lookingAtLitSide = (lookingAtLitSide < 0 ? -1 : 1) * Math.pow(Math.abs(lookingAtLitSide), 1/4);	//preserve sign, so raise to odd power
					planet.ringObj.uniforms.lookingAtLitSide = mathClamp(.5 + .5 * lookingAtLitSide, 0, 1);
					let viewDotSun = vec3.dot(viewfwd, sunDir);
					planet.ringObj.uniforms.backToFrontLitBlend = .5 - .5 * viewDotSun;	//clamp(sqrt(sqrt(dot( normalize(sun.pos - planet.pos), axis ) * dot( viewfwd, axis ))), 0., 1.) * -.5 + .5

					//have to recalculate these because the uniforms are shared, so they've been overwritten
					//...and the ring objs have to be drawn last for transparency reasons...

					//calculate sun position for lighting
					for (let starIndex = 0; starIndex < cfg.orbitStarSystem.stars.length; ++starIndex) {
						let star = cfg.orbitStarSystem.stars[starIndex];
						let tmp = [];
						vec3.sub(tmp, star.pos, planet.pos);
						planet.ringObj.uniforms.sunDir[0+3*starIndex] = tmp[0];
						planet.ringObj.uniforms.sunDir[1+3*starIndex] = tmp[1];
						planet.ringObj.uniforms.sunDir[2+3*starIndex] = tmp[2];
					}
					vec3.sub(planet.ringObj.uniforms.pos, planet.pos, cfg.orbitTarget.pos);
					quat.copy(planet.ringObj.uniforms.angle, planet.angle);
					planet.ringObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
					planet.ringObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
					planet.ringObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;

					//webkit bug
					planet.ringObj.shader.use();
					gl.uniform3fv(
						gl.getUniformLocation(
							planet.ringObj.shader.obj, 'sunDir[0]'
						), planet.ringObj.uniforms.sunDir);

					planet.ringObj.draw();

					gl.depthMask(true);
					gl.enable(gl.CULL_FACE);
				}
			}
		}

		//draw point planets
		// goes slow when comets are included
		//TODO make a buffer with a vertex per-planet rather than changing a single vertex
		for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
			let planet = cfg.orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			if (!planet.sceneObj || planet.visRatio < planetPointVisRatio) {
				if (picking) {
					//condition to skip moons ...
					if (planet.orbitVisRatio === undefined || planet.orbitVisRatio >= .03) {
						vec3.sub(pointObj.attrs.vertex.buffer.data, planet.pos, cfg.orbitTarget.pos);
						pointObj.attrs.vertex.buffer.updateData();
						cfg.pickObject.drawPoints({
							sceneObj : pointObj,
							targetCallback : planet,
							pointSize : 4,
							pointSizeScaleWithDist : false
						});
					}
				} else if (cfg.showPlanetsAsDistantPoints) {
					vec3.sub(pointObj.attrs.vertex.buffer.data, planet.pos, cfg.orbitTarget.pos);
					pointObj.attrs.vertex.buffer.updateData();
					pointObj.draw({
						uniforms : {
							color : planet.color,
							pointSize : 4
						}
					});
				} else {	//if (cfg.showStars) {
					//draw as a pixel with apparent magnitude
					//use the starfield shader
					vec3.sub(pointObj.attrs.vertex.buffer.data, planet.pos, cfg.orbitTarget.pos);
					pointObj.attrs.vertex.buffer.updateData();
					pointObj.attrs.velocity.buffer.data[0] = planet.vel[0];
					pointObj.attrs.velocity.buffer.data[1] = planet.vel[1];
					pointObj.attrs.velocity.buffer.data[2] = planet.vel[2];
					pointObj.attrs.velocity.buffer.updateData();
					pointObj.attrs.luminosity.buffer.data[0] = planet.luminosity || 10;
					pointObj.attrs.luminosity.buffer.updateData();
					pointObj.attrs.temperature.buffer.data[0] = planet.temperature || 20000;
					pointObj.attrs.temperature.buffer.updateData();
					pointObj.draw({
						shader : starfield.colorIndexShader,
						//disable depth test too?
						blend : [gl.SRC_ALPHA, gl.ONE],
						uniforms : {
							tempTex : 0,
							starTex : 1,
							starPointSizeBias : cfg.starPointSizeBias,
							starPointSizeScale : cfg.starPointSizeScale,
							starsPointAlpha : starsPointAlpha
							//old vars
							//pointSize : 1,
							//pointSizeMax : 5,
						},
						texs : [starfield.tempTex, starfield.starTex]
					});
				}
			}
		}

		//show small bodies in the solar system

		smallBodies.draw(
			tanFovY,
			picking,
			viewPosInv,
			invRotMat,
			distFromSolarSystemInM
		);

		//draw mouse-over highlight
		if (cfg.mouseOverTarget !== undefined &&
			!picking)
		{
			gl.disable(gl.DEPTH_TEST);
			pointObj.attrs.vertex.buffer.data[0] = cfg.mouseOverTarget.pos[0] - cfg.orbitTarget.pos[0];
			pointObj.attrs.vertex.buffer.data[1] = cfg.mouseOverTarget.pos[1] - cfg.orbitTarget.pos[1];
			pointObj.attrs.vertex.buffer.data[2] = cfg.mouseOverTarget.pos[2] - cfg.orbitTarget.pos[2];
			pointObj.attrs.vertex.buffer.updateData();
			pointObj.draw({
				uniforms : {
					color : [0,1,0,1],
					pointSize : 8
				},
				blend : [gl.SRC_ALPHA, gl.ONE]
			});
			gl.enable(gl.DEPTH_TEST);
		}

		//draw orbits
		if (cfg.showOrbits &&
			!picking)
		{
			window.orbitPathsDrawn = 0;
			for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
				let planet = cfg.orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;

				if (planet.renderOrbit) {

					let semiMajorAxis = planet.keplerianOrbitalElements.semiMajorAxis;
					let eccentricity = planet.keplerianOrbitalElements.eccentricity;
					let distPeriapsis = semiMajorAxis * (1 + eccentricity);	//largest distance from the parent planet

					//vector from view to parent planet
					let parentPlanet = planet.parent;
					vec3.sub(delta, parentPlanet.pos, cfg.orbitTarget.pos);
					vec3.sub(delta, delta, glutil.view.pos);
					let deltaLength = vec3.length(delta);

					planet.orbitVisRatio = distPeriapsis / (deltaLength * tanFovY);

					if (planet.orbitVisRatio > planetPointVisRatio) {
						const orbitPathSceneObj = starSystemsExtra.orbitPathSceneObj;
						//recenter around orbiting planet
						vec3.sub(orbitPathSceneObj.pos, planet.parent.pos, cfg.orbitTarget.pos);

						orbitPathSceneObj.uniforms.color = planet.color;
						orbitPathSceneObj.uniforms.A = planet.keplerianOrbitalElements.A;
						orbitPathSceneObj.uniforms.B = planet.keplerianOrbitalElements.B;
						orbitPathSceneObj.uniforms.eccentricity = planet.keplerianOrbitalElements.eccentricity;
						orbitPathSceneObj.uniforms.eccentricAnomaly = planet.keplerianOrbitalElements.eccentricAnomaly;
						orbitPathSceneObj.uniforms.orbitType = orbitPathIndexForType[planet.keplerianOrbitalElements.orbitType];
						orbitPathSceneObj.uniforms.fractionOffset = planet.keplerianOrbitalElements.fractionOffset || 0;
						orbitPathSceneObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
						orbitPathSceneObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
						orbitPathSceneObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;

						orbitPathSceneObj.draw();

						++orbitPathsDrawn;
					}

					if (planet.orbitVisRatio > .1) {
						if (planet.type != 'barycenter') {
							overlayTexts.add(planet);
						}
					}
				}
			}
		}

		if (cfg.showGravityWell &&
			!picking)
		{
			//do transformation math in double
			//TODO just give gl-matrix a type param in its init
			let glMvMat = [];
			let viewAngleInvd = [];
			quat.conjugate(viewAngleInvd, glutil.view.angle);
			quat.normalize(viewAngleInvd, viewAngleInvd);	//normalize in double precision
			mat4.fromQuat(glMvMat, viewAngleInvd);
			let viewPosInvd = [];
			vec3.negate(viewPosInvd, glutil.view.pos);
			mat4.translate(glMvMat, glMvMat, viewPosInvd);

			let orbitBasis = [];
			let orbitBasisInv = [];
			let mvMat = [];
			for (let planetIndex = 0; planetIndex < cfg.orbitStarSystem.planets.length; ++planetIndex) {
				let planet = cfg.orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;
				if (planet.radius === undefined) continue;
				//max radial dist is R * Math.pow(100, cfg.gravityWellRadialMaxLog100)
				if (planet.visRatio * Math.pow(100, cfg.gravityWellRadialMaxLog100) < .001) continue;

				mat4.identity(orbitBasis);
				mat4.identity(orbitBasisInv);
				for (let i = 0; i < 16; ++i) {
					mvMat[i] = glMvMat[i];
				}
				for (let i = 0; i < 3; ++i) {
					orbitBasis[i] = planet.orbitBasis[0][i];
					orbitBasis[4+i] = planet.orbitBasis[1][i];
					orbitBasis[8+i] = planet.orbitBasis[2][i];
				}
				//gravWellObjBasis = orbitBasis * gravityWellZScale * zOffsetByRadius * orbitBasis^-1 * planetPosBasis

				//calc this for non-sceneObj planets.  maybe I should store it as a member variable?
				let relPos = [
					planet.pos[0] - cfg.orbitTarget.pos[0],
					planet.pos[1] - cfg.orbitTarget.pos[1],
					planet.pos[2] - cfg.orbitTarget.pos[2]
				];

				//translate to planet center
				mat4.translate(mvMat, mvMat, relPos);

				//align gravity well with orbit plane
				mat4.multiply(mvMat, mvMat, orbitBasis);

				//align base gravity well with base of planet
				mat4.translate(mvMat, mvMat, [0,0,-planet.radius]);

				//scale for visual effect
				//too flat for earth, too sharp for sun
				if (cfg.gravityWellScaleFixed) {
					gravityWellTargetZScale = cfg.gravityWellScaleFixedValue;
				} else if (cfg.gravityWellScaleNormalized) {
					//causes larger wells to be much smaller and sharper ...
					//let gravityWellTargetZScale = 2000000 * Math.log(1 + z);
					//normalized visually per-planet.  scale is not 1-1 across planets
					gravityWellTargetZScale = 1 / cfg.orbitTarget.gravityWellScalar;
				}
				gravityWellZScale += .01 * (gravityWellTargetZScale - gravityWellZScale);
				if (gravityWellZScale !== gravityWellZScale) gravityWellZScale = 1;
				mat4.scale(mvMat, mvMat, [1,1,gravityWellZScale * planet.gravityWellScalar]);

				//apply orbit basis inverse
				mat4.transpose(orbitBasisInv, orbitBasis);
				mat4.multiply(mvMat, mvMat, orbitBasisInv);

				for (let i = 0; i < 16; ++i) {
					planet.gravWellObj.uniforms.mvMat[i] = mvMat[i];
				}

				planet.gravWellObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
				planet.gravWellObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
				planet.gravWellObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;

				planet.gravWellObj.draw();
			}
		}
	};
}

//unitsToTest should be largest to smallest
function unitsToStr(measureInBaseUnits, otherUnitsInBaseUnits, unitsToTest) {
	let bestUnit = 'm';
	let bestDist = measureInBaseUnits;
	for (let i = 0; i < unitsToTest.length; ++i) {
		let unit = unitsToTest[i];
		let ratio = measureInBaseUnits / otherUnitsInBaseUnits[unit];
		if (ratio > .1) {
			bestUnit = unit;
			bestDist = ratio;
			break;
		}
	}
	return bestDist.toFixed(4)+' '+bestUnit;
}

let distanceToStrUnits = ['Mpc', 'lyr', 'AU', 'km', 'm'];
function distanceToStr(distInMeters) {
	return unitsToStr(distInMeters, metersPerUnits, distanceToStrUnits);
}

let timeToStrUnits = ['years', 'days', 'm', 'h', 's'];
let daysPerUnits = {
	years : 365.4996816895717,
	days : 1,
	h : 1/24,
	m : 1/(60*24),
	s : 1/(60*60*24),
}
function timeToStr(timeInDays) {
	return unitsToStr(timeInDays, daysPerUnits, timeToStrUnits);
}

function refreshOrbitTargetDistanceText() {
	ids.orbitTargetDistance.innerText = distanceToStr(cfg.orbitTargetDistance);
}

cfg.refreshCurrentTimeText = function() {
	const astro = undefined;
	if (astro) {
		ids.currentTimeText.innerText = astro.julianToCalendar(cfg.julianDate);
	}
}

let primaryPlanetHorizonIDs = [10, 199, 299, 301, 399, 499, 599, 699, 799, 899, 999];

function calendarToJulian(d) {
	let year = d.getUTCFullYear();
	let month = d.getUTCMonth();
	let day = d.getUTCDate();
	let hour = d.getUTCHours();
	let min = d.getUTCMinutes();
	let sec = d.getUTCSeconds();
//http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
//testing against results in Wikipedia page: http://en.wikipedia.org/wiki/Julian_day
// (notice I couldn't recreate these results with the algorithm on the same page)
//for UTC date 2000 jan 1 12:00:00 this gives 2451545 - right on
//for UTC date 2013 jan 1 00:30:00 this gives 2456294.0208333335 - right on
	let oneBasedMonth = 1 + month;
	let astronomicalYear = year < 0 ? year + 1 : year;
	if (oneBasedMonth <= 2) {	//jan or feb
		--astronomicalYear;
		oneBasedMonth += 12;
	}
	let a = Math.floor(astronomicalYear / 100);
	let b = Math.floor(a / 4);
	let c = 2 - a + b;
	let e = Math.floor(365.25 * (astronomicalYear + 4716));
	let f = Math.floor(30.6001 * (oneBasedMonth + 1));
	let jdn = c + day + e + f - 1524.5;
	jdn += (hour + (min + sec / 60) / 60) / 24;
	return jdn;
}

function init() {
	ui.init();
	const gl = ui.gl;
	const glutil = ui.glutil;

	//post-WebGL init:

	const depthConstant = 1e-6;//2 / Math.log(1e+7 + 1);
	ui.ModifiedDepthShaderProgram = class extends glutil.Program {
		constructor(args) {
			// maximizing depth range: http://outerra.blogspot.com/2012/11/maximizing-depth-buffer-range-and.html
			//  but, true to the original design, I actually want more detail up front.  Maybe I'll map it to get more detail up front than it already has?
			args.vertexCode = `
uniform float zNear, zFar, depthConstant;
//float depthfunction(vec4 v) {
	//return (log(v.w + 1.) * depthConstant - 1.) * v.w;
	//return (2.0 * log(v.w / zNear) / log(zFar / zNear) - 1.) * v.w;
//}
#define depthfunction(x) ((x).z)

#define earthRadius	6.37101e+6	//TODO get this from solarSystem.planets[solarSystem.indexes.Earth].radius

uniform vec3 earthPos;
uniform vec3 earthNorthDir;
uniform float flatEarthCoeff;	//0. is original, 1. is fully flat-earth
vec4 flatEarthXForm(vec4 pos) {
	if (flatEarthCoeff == 0.) return pos;

	pos.xyz -= earthPos;

	float z = dot(pos.xyz, earthNorthDir);
	vec3 xy = pos.xyz - earthNorthDir * z;	//no 2D basis chosen mind you
	float r2 = length(xy);
	float r3 = length(pos);
	vec3 unitxy = xy / r2;
	float theta = atan(r2, z);

	// so we are mapping the spherical theta to the new cylindrical distance
	// but what about objects beyond the earth?  they need to be mapped above earth
	// and cannot exceed a certain distance... how to do this ... hyperbolic mapping?
	vec3 newPos = unitxy * theta;	// times planet radius ... len xyz always messes us up.

	//works fine for earth, but stars are too far out
	//newPos *= r3;

	//better for distant objects like orbits and stars
	newPos *= earthRadius;

	//what about the new z?
	// for vertices within the earth radius, they map to z=0
	// for vertices outside that radius, map to z=radius and map to xy =
	// well, for most objects away from earth, they will be at theta=pi/2 ... so how do i map them in more diverse locations?
	if (r3 > 1.1 * earthRadius) {
		z = earthRadius;
		newPos.xyz += z * earthNorthDir;
	}

	// past some threshold, wrap antarctica around the bottom.
	// the threshold angle is probably whatever the angle of the triangles touching the south pole on the sphere model.
	if (theta > 3.1) {
		newPos = pos.xyz;
	}

	pos.xyz = mix(pos.xyz, newPos.xyz, flatEarthCoeff);
	pos.xyz += earthPos;
	return pos;
}
	` + (args.vertexCode || '');
			if (args.uniforms === undefined) args.uniforms = {};
			args.uniforms.zNear = glutil.view.zNear;
			args.uniforms.zFar = glutil.view.zFar;
			args.uniforms.depthConstant = depthConstant;

			super(args);
		}
	}

	/*glutil.view.angle[0] = -0.4693271591372717;
	glutil.view.angle[1] = 0.7157221264895661;
	glutil.view.angle[2] = -0.4298784661116332;
	glutil.view.angle[3] = 0.28753844912098436;*/
	glutil.view.zNear = 1e+4;
	glutil.view.zFar = 1e+25;


	//now that GL is initialized, and before we are referencing planets,
	// build the solar sytsem

	starSystemsExtra.initSolarSystem();

	//once the planets are initialized, build the overlay side panel
	ui.initSidePanel();

	//cfg.julianDate = horizonsDynamicData.julianDate;	//get most recent recorded date from planets
	//TODO get current julian date from client
	// integrate forward by the small timestep between the two
	cfg.julianDate = calendarToJulian(new Date());

	cfg.initJulianDate = cfg.julianDate;
	cfg.refreshCurrentTimeText();

	const solarSystem = starSystemsExtra.solarSystem;
	solarSystem.copyPlanets(solarSystem.initPlanets, solarSystem.planets);

/*
	let texURLs = [];
	primaryPlanetHorizonIDs.forEach(horizonID => {
		let planet = solarSystem.planetForHorizonID[horizonID];
		texURLs.push('textures/'+planet.name.toLowerCase()+'.png');
	});
	texURLs.push('textures/jupiter-rings-color.png');
	texURLs.push('textures/saturn-rings-color.png');
	texURLs.push('textures/saturn-rings-transparency.png');
	texURLs.push('textures/saturn-rings-back-scattered.png');
	texURLs.push('textures/saturn-rings-forward-scattered.png');
	texURLs.push('textures/saturn-rings-unlit-side.png');
	//pick the cubemap textures once we know our max cubemap texture size -- 512 on firefox =( but 1024 on chrome =)
	texURLs.splice(texURLs.length, 0, skyCube.texURLs);
	console.log('loading '+texURLs.join(', '));
	preload(texURLs, () => {
		hide(ids.loadingDiv);
	}, percent =>{
		ids.loading.setAttribute('value', parseInt(100*percent));
	});
*/
	ids.infoPanel.addEventListener('touchmove', e => {
		if (!ui.showBodyInfo) {
			if (!e.target.classList.contains('scrollable')) {
				//only touch devices will scroll the div even when browsers know not to
				//so for touch scroll, explicitly tell it not to scroll
				e.preventDefault();
			}
		}
	});

	ids.toggleBodyInfo.style.transform = 'rotate(-90deg)';
	//ids.infoPanel.style.height = (ids.infoDiv.offset().top - ids.infoPanel.offset().top))+'px';
	//ids.infoPanel.style.bottom = '-25px';
	let currentAnimation;
	ids.toggleBodyInfo.addEventListener('click', e => {
		if (!ui.showBodyInfo) {
			ui.showBodyInfo = true;
			const infoDivTop = ids.infoPanel.offsetTop;
			const infoDivDestTop = ids.timeControlDiv.offsetTop + ids.timeControlDiv.offsetHeight;
			ids.infoPanel.style.height = (window.innerHeight - infoDivDestTop)+'px';

			//go from bottom-aligned to top-algined
			ids.infoPanel.style.bottom = 'auto';
			ids.infoPanel.style.top = infoDivTop+'px';
			ids.infoPanel.style.overflow = 'scroll';
			//interpolate up to the time controls ... or some fixed distance to it
			if (currentAnimation) currentAnimation.stop();
			currentAnimation = animate({
				duration : cfg.slideDuration,
				callback : frac => {
					const degrees = 180 * frac - 90;
					ids.toggleBodyInfo.style.transform = 'rotate('+degrees+'deg)';
					ids.infoPanel.style.top = (infoDivTop*(1-frac) + infoDivDestTop*frac)+'px';
				},
				done : () => {
					currentAnimation = undefined;
				},
			});
		} else {
			ui.showBodyInfo = false;
			ids.infoPanel.style.height = '104px';
			let infoDivBottom = window.innerHeight - (ids.infoPanel.offsetTop + ids.infoPanel.offsetHeight);
			let infoDivDestBottom = -25;

			//go from bottom-aligned to top-algined
			ids.infoPanel.style.top = 'auto';
			ids.infoPanel.style.bottom = infoDivBottom;
			ids.infoPanel.style.overflow = 'visible';
			//interpolate up to the time controls ... or some fixed distance to it
			if (currentAnimation) currentAnimation.stop();
			currentAnimation = animate({
				duration : cfg.slideDuration,
				callback : frac => {
					const degrees = 180 * (1 - frac) - 90;
					ids.toggleBodyInfo.style.transform = 'rotate('+degrees+'deg)';
					ids.infoPanel.style.bottom = (infoDivBottom*(1-frac) + infoDivDestBottom*frac)+'px';
				},
				done : () => {
					currentAnimation = undefined;
				},
			});
		}
	});

	show(ids.menu);
	show(ids.timeDiv);
	//show(ids.infoPanel);	//don't show here -- instead show upon first setOrbitTarget, when the css positioning gets fixed

	colorShader = new ui.ModifiedDepthShaderProgram({
		vertexCode : `
in vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;
void main() {
	vec4 vtx4 = vec4(vertex, 1.);
	vtx4 = flatEarthXForm(vtx4);
	vtx4 = mvMat * vtx4;
	gl_Position = projMat * vtx4;
	gl_PointSize = pointSize;
	gl_Position.z = depthfunction(gl_Position);
}
`,
		fragmentCode : `
uniform vec4 color;
out vec4 fragColor;
void main() {
	fragColor = color;
}
`,
		uniforms : {
			color : [1,1,1,1],
			pointSize : 4
		}
	});

	cfg.latLonShader = new ui.ModifiedDepthShaderProgram({
		vertexCode : `
in vec2 vertex;	//lat/lon pairs
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;
uniform float inverseFlattening;
uniform float scaleExaggeration;

` + cfg.geodeticPositionCode
+ cfg.quatRotateCode + `

void main() {
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	vec4 vtx4 = vec4(vtx3, 1.);
	vtx4 = flatEarthXForm(vtx4);
	vtx4 = mvMat * vtx4;
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
`,
		fragmentCode : `
uniform vec4 color;
out vec4 fragColor;
void main() {
	fragColor = color;
}
`,
		uniforms : {
			color : [1,1,1,1]
		}
	});

	calcTides.initGL();

	console.log('init overlay slider...');
	ids.overlaySlider.value = 100 * cfg.heatAlpha;
	ids.overlaySlider.addEventListener('input', e => {
		cfg.heatAlpha = ids.overlaySlider.value / 100;
		calcTides.updateHeatMapAlpha(cfg.heatAlpha);
	});
	calcTides.updateHeatMapAlpha(cfg.heatAlpha);

	pointObj = new glutil.SceneObject({
		mode : gl.POINTS,
		shader : colorShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 3,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			}),
			velocity : new glutil.ArrayBuffer({
				dim : 3,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			}),
			luminosity : new glutil.ArrayBuffer({
				dim : 1,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			}),
			temperature : new glutil.ArrayBuffer({
				dim : 1,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			})
		},
		uniforms : {
			pointSize : 4
		},
		parent : null
	});

	lineObj = new glutil.SceneObject({
		mode : gl.LINES,
		shader : colorShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				count : 2,
				usage : gl.DYNAMIC_DRAW
			})
		},
		uniforms : {
			color : [1,1,1,1]
		},
		parent : null,
		static : true
	});

	quadObj = new glutil.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 2,
				data : [-1,-1, 1,-1, -1,1, 1,1]
			}),
			texCoord : new glutil.ArrayBuffer({
				dim : 2,
				data : [0,0, 1,0, 0,1, 1,1]
			})
		},
		parent : null,
		static : true
	});


	//init our planet shaders

	console.log('solarSystemExtra.initPlanetSceneObjs...');
	starSystemsExtra.initPlanetSceneObjs();

	console.log('milkyWay.init...');
	milkyWay.init();

	//init stars now that shaders are made
	console.log('starfield.init...');
	starfield.init();

	//and small bodies in the solar system
	console.log('smallBodies.init...');
	smallBodies = new SmallBodies();

	//and galaxies
	console.log('galaxies.init...');
	galaxies = new Galaxies();

	//initialize here after webgl canvas is up
	cfg.pickObject = new PickObject();

	//only do this after the orbit shader is made
	starSystemsExtra.initExoplanets();

	skyCube.init();

	initScene();
}

cfg.setOrbitTarget = function(newTarget) {
	let selectingNewSystem = false;
	if (newTarget === undefined) {
		newTarget = cfg.orbitStarSystem.stars[0];
	}
	if (newTarget instanceof StarSystem) {
		let targetSystem = newTarget;
		let i = 0;
		for (; i < targetSystem.planets.length; ++i) {
			if (targetSystem.planets[i].type !== 'barycenter') {
				newTarget = targetSystem.planets[i];
				selectingNewSystem = true;
				break;
			}
		}
		if (i == targetSystem.planets.length) {
			//assume we get something?
			//will this get reached for the named planets? if so, give them a star!
			console.log("failed to find a planet in the system to select!");
		}
	}
	if (newTarget instanceof Planet) {
		//if we're changing star systems...
		if (newTarget.starSystem !== cfg.orbitStarSystem) {
			selectingNewSystem = true;
			cfg.orbitStarSystem = newTarget.starSystem;
		}
	}

	//query the server for this small body's orbit data
	if (newTarget instanceof SmallBody) {
		newTarget.onSelect();
	}

	if (selectingNewSystem) {
		if (resetDistanceOnSelectingNewSystem) {
			cfg.orbitTargetDistance = Math.max(100000, newTarget.radius || newTarget.equatorialRadius || 0);
			for (let i = 0; i < cfg.orbitStarSystem.planets.length; ++i) {
				let planet = cfg.orbitStarSystem.planets[i];

				//zoom out past orbit for planets
				//TODO don't do this for stars orbiting the milky way
				if (planet.sourceData !== undefined && planet.sourceData.semiMajorAxis !== undefined) {
					cfg.orbitTargetDistance = Math.max(cfg.orbitTargetDistance, planet.sourceData.semiMajorAxis);
				}
				if (planet.keplerianOrbitalElements !== undefined && planet.keplerianOrbitalElements.semiMajorAxis !== undefined) {
					cfg.orbitTargetDistance = Math.max(cfg.orbitTargetDistance, planet.keplerianOrbitalElements.semiMajorAxis);
				}

				refreshOrbitTargetDistanceText();
			}
		}
		cfg.orbitStarSystem.updatePlanetsPos();
	}

	if (cfg.orbitTarget !== undefined) vec3.sub(orbitOffset, cfg.orbitTarget.pos, newTarget.pos);
	cfg.orbitTarget = newTarget

	//reset orbit distance so it doesn't lock to surface immediately
	if (orbitGeodeticLocation !== undefined) {
		//TODO this is what kicks the next frame out of fixed-view mode,
		// and the code that does the kicking out resets the distance as if the user had zoomed out (i.e. relatively close)
		// so somehow impose this distance over the other distance
		orbitDistance = cfg.orbitTargetDistance = cfg.orbitTarget.radius * 10;

		//and re-orient the camera (because locking on surface had it flipped around / will flip it back around)
		quat.copy(orbitGeodeticLocation.lastViewAngle, glutil.view.angle);
	}

	//refresh info div

	ids.orbitTargetText.innerText = cfg.orbitTarget.name;
	ids.infoDiv.innerHTML = '';
	ids.infoDiv.appendChild(Hr());
	ids.infoDiv.appendChild(Div({innerText:'Type: '+cfg.orbitTarget.type}));
	if (cfg.orbitTarget.mass !== undefined) {
		ids.infoDiv.appendChild(Div({innerText:'Mass: '+cfg.orbitTarget.mass+' kg'}));
	}
	if (cfg.orbitTarget.equatorialRadius !== undefined) {
		ids.infoDiv.appendChild(Div({innerText:'Equatorial Radius: '+distanceToStr(cfg.orbitTarget.equatorialRadius)}));
	} else if (cfg.orbitTarget.radius !== undefined) {
		ids.infoDiv.appendChild(Div({innerText:'Radius: '+distanceToStr(cfg.orbitTarget.radius)}));
	}
	if (cfg.orbitTarget.inverseFlattening !== undefined) {
		ids.infoDiv.appendChild(Div({innerText:'Inverse Flattening: '+cfg.orbitTarget.inverseFlattening}));
	}
	if (cfg.orbitTarget.rotationPeriod !== undefined) {
		ids.infoDiv.appendChild(Div({innerText:'Rotation Period: '+timeToStr(cfg.orbitTarget.rotationPeriod)}));
	}
	if (cfg.orbitTarget.ringRadiusRange !== undefined) {
		ids.infoDiv.appendChild(Div({innerText:'Ring Min Radius: '+distanceToStr(cfg.orbitTarget.ringRadiusRange[0])}));
		ids.infoDiv.appendChild(Div({innerText:'Ring Max Radius: '+distanceToStr(cfg.orbitTarget.ringRadiusRange[1])}));
	}
	if (cfg.orbitTarget.keplerianOrbitalElements) {
		if (cfg.orbitTarget.keplerianOrbitalElements.semiMajorAxis) {
			ids.infoDiv.appendChild(Div({innerText:'Semi-Major Axis: '+distanceToStr(cfg.orbitTarget.keplerianOrbitalElements.semiMajorAxis)}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.orbitType) {
			ids.infoDiv.appendChild(Div({innerText:'Orbit Type: '+cfg.orbitTarget.keplerianOrbitalElements.orbitType}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.eccentricity) {
			ids.infoDiv.appendChild(Div({innerText:'Eccentricity: '+cfg.orbitTarget.keplerianOrbitalElements.eccentricity}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.eccentricAnomaly) {
			ids.infoDiv.appendChild(Div({innerHTML:'Eccentric Anomaly: '+mathDeg(cfg.orbitTarget.keplerianOrbitalElements.eccentricAnomaly)+'&deg;'}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.longitudeOfAscendingNode) {
			ids.infoDiv.appendChild(Div({innerHTML:'Longitude of Ascending Node: '+mathDeg(cfg.orbitTarget.keplerianOrbitalElements.longitudeOfAscendingNode)+'&deg;'}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.argumentOfPeriapsis) {
			ids.infoDiv.appendChild(Div({innerHTML:'Argument of Pericenter: '+mathDeg(cfg.orbitTarget.keplerianOrbitalElements.argumentOfPeriapsis)+'&deg;'}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.inclination) {
			ids.infoDiv.appendChild(Div({innerHTML:'Inclination: '+mathDeg(cfg.orbitTarget.keplerianOrbitalElements.inclination)+'&deg;'}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.timeOfPeriapsisCrossing) {
			ids.infoDiv.appendChild(Div({innerHTML:'Time of Periapsis Crossing: '+timeToStr(cfg.orbitTarget.keplerianOrbitalElements.timeOfPeriapsisCrossing)}));
		}
		if (cfg.orbitTarget.keplerianOrbitalElements.orbitalPeriod) {
			ids.infoDiv.appendChild(Div({innerText:'Orbital Period: '+timeToStr(cfg.orbitTarget.keplerianOrbitalElements.orbitalPeriod)}));
		}
	}
	ids.infoDiv.appendChild(Br());
	ids.infoDiv.appendChild(Br());

	setTimeout(() => {
		//refresh offset if the infoPanel is hidden
		ids.infoPanel.style.bottom = '-25px';
		if (!ui.showBodyInfo) {
			ids.infoPanel.style.height = '104px';
		} else {
			const infoDivDestTop = ids.timeControlDiv.offsetTop + ids.timeControlDiv.offsetHeight;
			ids.infoPanel.style.height = (window.innerHeight - infoDivDestTop)+'px';
		}
		setTimeout(() => {
			show(ids.infoPanel);
		}, 1);
	}, 1);
}

//TODO use glutil.mouseDir?
/*
function mouseRay() {
	let viewX = glutil.view.pos[0];
	let viewY = glutil.view.pos[1];
	let viewZ = glutil.view.pos[2];
	let viewFwdX = -2 * (glutil.view.angle[0] * glutil.view.angle[2] + glutil.view.angle[3] * glutil.view.angle[1]);
	let viewFwdY = -2 * (glutil.view.angle[1] * glutil.view.angle[2] - glutil.view.angle[3] * glutil.view.angle[0]);
	let viewFwdZ = -(1 - 2 * (glutil.view.angle[0] * glutil.view.angle[0] + glutil.view.angle[1] * glutil.view.angle[1]));
	let viewRightX = 1 - 2 * (glutil.view.angle[1] * glutil.view.angle[1] + glutil.view.angle[2] * glutil.view.angle[2]);
	let viewRightY = 2 * (glutil.view.angle[0] * glutil.view.angle[1] + glutil.view.angle[2] * glutil.view.angle[3]);
	let viewRightZ = 2 * (glutil.view.angle[0] * glutil.view.angle[2] - glutil.view.angle[3] * glutil.view.angle[1]);
	let viewUpX = 2 * (glutil.view.angle[0] * glutil.view.angle[1] - glutil.view.angle[3] * glutil.view.angle[2]);
	let viewUpY = 1 - 2 * (glutil.view.angle[0] * glutil.view.angle[0] + glutil.view.angle[2] * glutil.view.angle[2]);
	let viewUpZ = 2 * (glutil.view.angle[1] * glutil.view.angle[2] + glutil.view.angle[3] * glutil.view.angle[0]);
	let aspectRatio = ui.canvas.width / ui.canvas.height;
	let mxf = cfg.mouse.xf * 2 - 1;
	let myf = 1 - cfg.mouse.yf * 2;
	let tanFovY = Math.tan(glutil.view.fovY * Math.PI / 360);
	let mouseDirX = viewFwdX + tanFovY * (viewRightX * mxf * aspectRatio + viewUpX * myf);
	let mouseDirY = viewFwdY + tanFovY * (viewRightY * mxf * aspectRatio + viewUpY * myf);
	let mouseDirZ = viewFwdZ + tanFovY * (viewRightZ * mxf * aspectRatio + viewUpZ * myf);
	let mouseDirLength = Math.sqrt(mouseDirX * mouseDirX + mouseDirY * mouseDirY + mouseDirZ * mouseDirZ);
	return [mouseDirX/mouseDirLength, mouseDirY/mouseDirLength, mouseDirZ/mouseDirLength];
}
*/

function initScene() {
	const gl = ui.gl;
	//gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL);
	gl.clearColor(0,0,0,0);

	//assign our initial orbiting solar system
	const solarSystem = starSystemsExtra.solarSystem;
	cfg.orbitStarSystem = solarSystem;

	let trackPlanetName = 'Earth';
	if (urlparams.get('target') !== null) {
		trackPlanetName = urlparams.get('target');
	}

	const trackPlanet = solarSystem.planets[solarSystem.indexes[trackPlanetName]];
	cfg.setOrbitTarget(trackPlanet);
	cfg.orbitTargetDistance = 2. * cfg.orbitTarget.radius;
	refreshOrbitTargetDistanceText();
	orbitDistance = cfg.orbitTargetDistance;

	let tmpQ = quat.create();
	cfg.mouse = new Mouse3D({
		pressObj : ui.canvas,
		move : function(dx,dy) {
			let rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			let tanHalfFov = Math.tan(glutil.view.fovY * Math.PI / 360);
			if (orbitGeodeticLocation !== undefined) {
				quat.setAxisAngle(tmpQ, [tanHalfFov * dy, tanHalfFov * dx, 0], rotAngle);
			} else {
				quat.setAxisAngle(tmpQ, [-dy, -dx, 0], rotAngle);
			}
			quat.multiply(glutil.view.angle, glutil.view.angle, tmpQ);
			quat.normalize(glutil.view.angle, glutil.view.angle);
		},
		zoom : function(zoomChange) {
			let scale = Math.exp(-orbitZoomFactor * zoomChange);

			if (orbitGeodeticLocation !== undefined) {
				glutil.view.fovY *= scale;
				glutil.updateProjection();
			} else {
				cfg.orbitTargetDistance = Math.min(cfg.orbitTargetDistance * scale, 4000 * metersPerUnits.Mpc);
			}

			refreshOrbitTargetDistanceText();
		},
		passiveMove : function() {
			//seems to not always update, but running this every frame is much slower
			cfg.pickObject.pick(false);
		},
		click : function() {
			if (cfg.mouse.isDragging) return;
			cfg.pickObject.pick(true);
		}
	});

	update();
}

function update() {
	const gl = ui.gl;

	/*
	some new thoughts:
	- only enable mouseover highlighting if orbit major axis (largest radius) exceeds point vis threshold
	*/
	overlayTexts.updateBegin();
	if (cfg.mouseOverTarget) overlayTexts.add(cfg.mouseOverTarget);

	//TODO instead of always showing the orbit target, show what it collapses into (moon -> earth -> solar system -> milky way)
	// but this would mean unifying all the overlay text show/hide conditions ...
	if (cfg.orbitTarget && cfg.orbitTarget.name && cfg.orbitTarget.pos) overlayTexts.add(cfg.orbitTarget);

	/* converage angle on target planet * /

	let delta = vec3.create();
	vec3.scale(delta, glutil.view.pos, -1);
	let fwd = vec3.create();
	quatZAxis(fwd, glutil.view.angle);
	vec3.scale(fwd, fwd, -1);
	let axis = vec3.create();
	vec3.cross(axis, fwd, delta);
	vec3.scale(axis, axis, 1/vec3.length(delta));	//divide out length of delta so axis length is sin(theta)
	let sinTheta = vec3.length(axis);
	let theta = 0;
	if (sinTheta > 1e-3) {
		vec3.scale(axis, axis, 1/sinTheta);	//normalize axis
		theta = Math.asin(mathClamp(sinTheta, -1,1));
		let q = quat.create();
		let convergeAngleCoeff = .5;
		quat.setAxisAngle(q, axis, theta * convergeAngleCoeff);
		quat.mul(glutil.view.angle, glutil.view.angle, q);
	}
	/**/

	// track ball orbit
	//assumes z points away from the planet

	//while fixed, zoom affects fov, so if our fov surpasses the max then zoom us back out
	if (glutil.view.fovY > 120) {
		orbitDistance = cfg.orbitTargetDistance = cfg.orbitTarget.radius * .2;
	}

	//if we are close enough to the planet then rotate with it
	let fixViewToSurface = cfg.orbitTargetDistance < cfg.orbitTarget.radius * .1 &&
		!(cfg.orbitTarget instanceof Galaxy);	//don't allow perspetive from galaxy "surfaces"
	if (fixViewToSurface && orbitGeodeticLocation === undefined) {
		let pos = [];
		vec3.add(pos, cfg.orbitTarget.pos, glutil.view.pos);

		orbitGeodeticLocation = solarSystemBarycentricToPlanetGeodetic(cfg.orbitTarget, pos);
		orbitGeodeticLocation.lastViewAngle = quat.create();
		quat.copy(orbitGeodeticLocation.lastViewAngle, glutil.view.angle);

		glutil.view.fovY = 120;

		//...and in one fell swoop, turn the camera around
		//TODO spread this out over a few frames
		let rot = [];
		quat.identity(rot);
		quat.rotateY(rot, rot, Math.PI);
		quat.multiply(glutil.view.angle, glutil.view.angle, rot);
	} else if (!fixViewToSurface && orbitGeodeticLocation !== undefined) {
		clearGeodeticLocation();
	}

	if (orbitGeodeticLocation !== undefined) {
		//fix the height at the surface
		//TODO spread this out over a few frames
		let height = (cfg.orbitTarget.equatorialRadius || cfg.orbitTarget.radius || 1000) * .01;//= orbitGeodeticLocation.height;
		orbitDistance = cfg.orbitTargetDistance = height;

		planetGeodeticToSolarSystemBarycentric(
			glutil.view.pos,
			cfg.orbitTarget,
			orbitGeodeticLocation.lat,
			orbitGeodeticLocation.lon,
			height);
		vec3.sub(glutil.view.pos, glutil.view.pos, cfg.orbitTarget.pos);
	} else {
		resetOrbitViewPos();
	}

	let orbitConvergeCoeff = .9;
	vec3.scale(orbitOffset, orbitOffset, orbitConvergeCoeff);
	{
		let logDist = Math.log(orbitDistance);
		let logTarget = Math.log(cfg.orbitTargetDistance);
		let coeff = .05;
		let newLogDist = (1 - coeff) * logDist + coeff * logTarget;
		orbitDistance = Math.exp(newLogDist);
	}

	if (!cfg.integrationPaused) {
		cfg.julianDate += cfg.integrateTimeStep;
		cfg.refreshCurrentTimeText();
	}

	if (cfg.julianDate !== cfg.lastJulianDate) {

		//recompute position by delta since init planets
		cfg.orbitStarSystem.updatePlanetsPos();

		//recompute angle based on sidereal rotation period
		for (let i = 0; i < cfg.orbitStarSystem.planets.length; ++i) {
			let planet = cfg.orbitStarSystem.planets[i];
			if (planet.rotationPeriod) {
				let angle = ((cfg.julianDate % planet.rotationPeriod) / planet.rotationPeriod) * 2 * Math.PI;

				//I really don't like this variable.  I don't think it should be used.
				if (planet.rotationOffset !== undefined) {
					angle += planet.rotationOffset;
				}

				let zrot = [0,0,0,1];
				quat.rotateZ(zrot, zrot, angle);
				quat.multiply(planet.angle, planet.tiltAngle, zrot);
			} else {
				quat.copy(planet.angle, planet.tiltAngle);
			}
		}

		//if we're not fixed on the surface but we are close enough to orbit then spin with the planet
		//TODO this is messing up on mercury, venus, pluto ... retrograde planets + mercury ... ?
		if (cfg.orbitTargetDistance < cfg.orbitTarget.radius * 10) {
			let orbitTargetAxis = [0,0,1];
			quatZAxis(orbitTargetAxis, cfg.orbitTarget.angle);

			let deltaJulianDate = cfg.julianDate - cfg.lastJulianDate;

			let deltaAngle = quat.create();
			quat.setAxisAngle(deltaAngle, orbitTargetAxis, deltaJulianDate * (cfg.orbitTarget.rotationPeriod || 0) * 2 * Math.PI);

			quat.multiply(glutil.view.angle, deltaAngle, glutil.view.angle);
			quat.normalize(glutil.view.angle, glutil.view.angle);

			//reset angle (to keep targets from offsetting when the ffwd/rewind buttons are pushed)
			resetOrbitViewPos();
		}

		//if we are zoomed in then keep the view following the mouse over target
		if (orbitGeodeticLocation !== undefined &&
			cfg.mouseOverTarget !== undefined)
		{
			let viewFwd = [];
			quatZAxis(viewFwd, glutil.view.angle);
			vec3.scale(viewFwd, viewFwd, -1);
			vec3.normalize(viewFwd, viewFwd);

			let destViewFwd = [];
			vec3.sub(destViewFwd, cfg.mouseOverTarget.pos, glutil.view.pos);
			vec3.normalize(destViewFwd, destViewFwd);

			let axis = [];
			vec3.cross(axis, viewFwd, destViewFwd);
			let sinTheta = vec3.length(axis);
			if (sinTheta > 1e-5) {
				let theta = Math.asin(mathClamp(sinTheta, -1, 1));
				let rot = [];
				quat.setAxisAngle(rot, axis, theta);
				quat.multiply(glutil.view.angle, rot, glutil.view.angle);
				quat.normalize(glutil.view.angle, glutil.view.angle);
			}
		}
	}

	cfg.lastJulianDate = cfg.julianDate;

	glutil.scene.setupMatrices();
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	if (window.showPickScene) {
		cfg.pickObject.pick(false, true);
	} else {
		cfg.drawScene(false);
	}
	glutil.clearAlpha();

	requestAnimationFrame(update);

	overlayTexts.updateEnd();
}

init();
