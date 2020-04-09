var orbitPathResolution = 500;
var ringResolution = 200;
var orbitPathIndexForType = {elliptic:0, hyperbolic:1, parabolic:2};	//used in the shader
var orbitPathShader;

var julianDate = 0;
var lastJulianDate = 0;
var initJulianDate = 0;

// track ball motion variables
var mouseOverTarget;
var orbitStarSystem;	//only do surface calculations for what star system we are in
var orbitTarget;
var orbitGeodeticLocation;
var orbitDistance;
var orbitOffset = [0,0,0];
var orbitTargetDistance;
var orbitZoomFactor = .0003;	// upon mousewheel

//if we're orbiting at 1AU then we can only click things at 1000 AU
var ratioOfOrbitDistanceToAllowSelection = 10000;

var mouse;
var mouseDir;

var resetDistanceOnSelectingNewSystem = false;

var gl;
var canvas;

var colorShader;
var latLonShader;

var pointObj;

var planetAmbientLight = .1;

var displayMethods = [
	'None',
	'Tangent Tidal',
	'Normal Tidal',
	'Total Tidal',
	'Tangent Gravitational',
	'Normal Gravitational',
	'Total Gravitational',
	'Tangent Total',
	'Normal Total',
	'Total'
];
var displayMethod = 'None';
var planetInfluences = [];

var showLinesToOtherPlanets = false;
var showVelocityVectors = false;
var velocityVectorScale = 30;
var showRotationAxis = false;
var showOrbitAxis = false;
var showLatAndLonLines = false;
var showGravityWell = false;
var showPlanetsAsDistantPoints = true;
var showOrbits = true;
var planetScaleExaggeration = 1;

var gravityWellScaleNormalized = true;
var gravityWellScaleFixed = false;
var gravityWellScaleFixedValue = 2000;
var gravityWellRadialMinLog100 = -1;
var gravityWellRadialMaxLog100 = 2;

var overlayShowOrbitTarget = true;
var overlayShowCurrentPosition = false;

var heatAlpha = .5;

var integrationPaused = true;
var defaultIntegrateTimeStep = 1/(24*60);
var integrateTimeStep = defaultIntegrateTimeStep;


var Galaxy = makeClass({
	init : function(args) {
		for (k in args) {
			this[k] = args[k];
		}
	}
});


function clearGeodeticLocation() {
	if (orbitGeodeticLocation === undefined) return;

	orbitDistance = orbitTargetDistance = orbitTarget.radius * .2;

	glutil.view.fovY = 90;
	glutil.updateProjection();

	//...and in one fell swoop, restore the original camera angle
	//TODO spread this out over a few frames
	quat.copy(glutil.view.angle, orbitGeodeticLocation.lastViewAngle);
	
	orbitGeodeticLocation = undefined;
}

function resetOrbitViewPos() {
	var viewAngleZAxis = vec3.create();
	vec3.quatZAxis(viewAngleZAxis, glutil.view.angle);
	vec3.scale(glutil.view.pos, viewAngleZAxis, orbitDistance + (orbitTarget.equatorialRadius || orbitTarget.radius || 0));
	vec3.add(glutil.view.pos, glutil.view.pos, orbitOffset);
}

function refreshMeasureText() {
	$('#measureMin').text(orbitTarget.measureMin === undefined ? '' : (orbitTarget.measureMin.toExponential() + ' m/s^2'));
	$('#measureMax').text(orbitTarget.measureMax === undefined ? '' : (orbitTarget.measureMax.toExponential() + ' m/s^2'));
}

var drawScene;
var gravityWellZScale = 1;
var gravityWellTargetZScale = 1;
var planetPointVisRatio = .001;
var showFPS = false;
(function(){
	var delta = vec3.create();//[];
	var viewAngleInv = quat.create();
	var invRotMat = mat4.create();
	var viewPosInv = vec3.create();
	var viewfwd = vec3.create();
	
	//fps counter
	var frames = 0;
	var lastFPSTime = Date.now();

	//used for new gpu update of tide tex
	var updatePlanetStateBuffer = new Float32Array(1);

	drawScene = function(picking) {

		//should picking count towards fps? nah
		if (!picking) {
			frames++;
			var thisTime = Date.now();
			if (thisTime - lastFPSTime > 1000) {
				var fps = frames * 1000 / (thisTime - lastFPSTime);
				if (showFPS) {
					console.log('fps '+fps);
				}
				frames = 0;
				lastFPSTime = thisTime;	
			}
		}

		if (overlayShowCurrentPosition && orbitTarget.isa(Planet)) {
			//update the min and max to reflect the current position
			var x = glutil.view.pos[0] + orbitTarget.pos[0];
			var y = glutil.view.pos[1] + orbitTarget.pos[1];
			var z = glutil.view.pos[2] + orbitTarget.pos[2];
			var t = calcMetricForce([x,y,z], orbitTarget);
			$('#measureMin').text(t === undefined ? '' : t.toExponential() + ' m/s^2');
			$('#measureMax').text(t === undefined ? '' : t.toExponential() + ' m/s^2');
		}
		
		//update the planet state texture
		// right now it's only used for tide calcs so
		//1) only update it if tide calcs are on
		//2) only include planets that are enabled for tide calcs
		// more discernment later when i make it general purpose
		var useOverlay = overlayShowOrbitTarget && displayMethod != 'None';
		if (useOverlay && !picking) {
		
			var targetSize = orbitStarSystem.planetStateTex.width * orbitStarSystem.planetStateTex.height * 4;
			
			if (updatePlanetStateBuffer.length != targetSize) {
				updatePlanetStateBuffer = new Float32Array(targetSize);
			}

			//gather pos and mass
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
				updatePlanetStateBuffer[0 + 4 * planetIndex] = planet.pos[0];
				updatePlanetStateBuffer[1 + 4 * planetIndex] = planet.pos[1];
				updatePlanetStateBuffer[2 + 4 * planetIndex] = planet.pos[2];
				
				if (!planetInfluences[planetIndex]) {
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
			orbitStarSystem.planetStateTex.bind();
			
			gl.texSubImage2D(
				gl.TEXTURE_2D,
				0,
				0,
				0,
				orbitStarSystem.planetStateTex.width,
				orbitStarSystem.planetStateTex.height,
				gl.RGBA,
				gl.FLOAT,
				updatePlanetStateBuffer);
			
			orbitStarSystem.planetStateTex.unbind();
		}

		//get our calculations for orientation
		
		quat.conjugate(viewAngleInv, glutil.view.angle);
		mat4.fromQuat(invRotMat, viewAngleInv);
		
		mat4.copy(glutil.scene.mvMat, invRotMat);

		//TODO pull from matrix
		vec3.quatZAxis(viewfwd, glutil.view.angle);
		vec3.scale(viewfwd, viewfwd, -1);
					
		var tanFovY = Math.tan(glutil.view.fovY * Math.PI / 360);
		
		var distFromSolarSystemInM = vec3.length(glutil.view.pos);
		var distFromSolarSystemInLyr = distFromSolarSystemInM / metersPerUnits.lyr;
		var distFromSolarSystemInPc = distFromSolarSystemInM / metersPerUnits.pc;
		var distFromSolarSystemInMpc = distFromSolarSystemInM / metersPerUnits.Mpc;

		
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
		vec3.scale(viewPosInv, glutil.view.pos, -1/interGalacticRenderScale);
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

		//draw debug lines 
		if (!picking) {
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;
				if (planet.pos === undefined) continue;	//comets don't have pos yet, but I'm working on that
				if (orbitTarget.pos === undefined) continue;
			
				/*
				if (planet.isComet) continue;
				if (planet.isAsteroid) continue;
				if (planet.parent !== solarSystem.planets[solarSystem.indexes.Sun] &&
					planet !== solarSystem.planets[solarSystem.indexes.Sun] &&
					planet !== solarSystem.planets[solarSystem.indexes.Moon]) continue;
				*/

				//draw lines to other planets
				if (showLinesToOtherPlanets &&
					orbitStarSystem == solarSystem &&
					orbitTarget !== planet &&
					(!planet.parent || planet.parent.index == 0))
				{
					//while here, update lines

					vec3.sub(delta, planet.pos, orbitTarget.pos);

					var dist = vec3.length(delta);
					vec3.scale(delta, delta, 1/dist);

					lineObj.attrs.vertex.data[0] = delta[0] * orbitTarget.radius;
					lineObj.attrs.vertex.data[1] = delta[1] * orbitTarget.radius;
					lineObj.attrs.vertex.data[2] = delta[2] * orbitTarget.radius;
					lineObj.attrs.vertex.data[3] = delta[0] * orbitDistance;
					lineObj.attrs.vertex.data[4] = delta[1] * orbitDistance;
					lineObj.attrs.vertex.data[5] = delta[2] * orbitDistance;

					lineObj.attrs.vertex.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show velocity vectors
				if (showVelocityVectors) {
					vec3.sub(delta, planet.pos, orbitTarget.pos);
					lineObj.attrs.vertex.data[0] = delta[0];
					lineObj.attrs.vertex.data[1] = delta[1];
					lineObj.attrs.vertex.data[2] = delta[2];
					lineObj.attrs.vertex.data[3] = delta[0] + planet.vel[0] * velocityVectorScale;
					lineObj.attrs.vertex.data[4] = delta[1] + planet.vel[1] * velocityVectorScale;
					lineObj.attrs.vertex.data[5] = delta[2] + planet.vel[2] * velocityVectorScale;
					lineObj.attrs.vertex.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show rotation axis
				if (showRotationAxis) {
					vec3.sub(delta, planet.pos, orbitTarget.pos);
					var axis = [0,0,1];
					vec3.quatZAxis(axis, planet.angle);
					lineObj.attrs.vertex.data[0] = delta[0] + axis[0] * 2 * planet.radius;
					lineObj.attrs.vertex.data[1] = delta[1] + axis[1] * 2 * planet.radius;
					lineObj.attrs.vertex.data[2] = delta[2] + axis[2] * 2 * planet.radius;
					lineObj.attrs.vertex.data[3] = delta[0] + axis[0] * -2 * planet.radius;
					lineObj.attrs.vertex.data[4] = delta[1] + axis[1] * -2 * planet.radius;
					lineObj.attrs.vertex.data[5] = delta[2] + axis[2] * -2 * planet.radius;
					lineObj.attrs.vertex.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}

				//show orbit axis
				if (showOrbitAxis && planet.orbitAxis) {
					vec3.sub(delta, planet.pos, orbitTarget.pos);
					lineObj.attrs.vertex.data[0] = delta[0] + planet.orbitAxis[0] * 2 * planet.radius;
					lineObj.attrs.vertex.data[1] = delta[1] + planet.orbitAxis[1] * 2 * planet.radius;
					lineObj.attrs.vertex.data[2] = delta[2] + planet.orbitAxis[2] * 2 * planet.radius;
					lineObj.attrs.vertex.data[3] = delta[0] + planet.orbitAxis[0] * -2 * planet.radius;
					lineObj.attrs.vertex.data[4] = delta[1] + planet.orbitAxis[1] * -2 * planet.radius;
					lineObj.attrs.vertex.data[5] = delta[2] + planet.orbitAxis[2] * -2 * planet.radius;
					lineObj.attrs.vertex.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}
			}
		}

		//update planet vis ratio
		//for picking do we need to? TODO make sure the view doesn't change between the last render and a pick
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			//update vis ratio
			var dx = planet.pos[0] - glutil.view.pos[0] - orbitTarget.pos[0];
			var dy = planet.pos[1] - glutil.view.pos[1] - orbitTarget.pos[1];
			var dz = planet.pos[2] - glutil.view.pos[2] - orbitTarget.pos[2];
			//approximated pixel width with a fov of 90 degrees
			planet.visRatio = planetScaleExaggeration * planet.radius / (Math.sqrt(dx * dx + dy * dy + dz * dz) * tanFovY);
		}

		//draw sphere planets
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
		(function() {
			var planet = orbitStarSystem.planets[planetIndex];

			if (planet.sceneObj) {
			
				var canSee = !planet.hide && planet.visRatio >= planetPointVisRatio;

				//if the planet is visible then
				// if there is no texture then
				//  load the texture and bind it
				if (canSee) {
					//if we have no image ... then load it 
					if (!planet.img) {
						if (planet.imgURL) {
							planet.imgIsLoading = true;
							planet.img = new Image();
							planet.img.onload = function() {
								//console.log('finished loading tex for planet '+planet.name);
								planet.imgIsLoading = undefined;
								delete planet.imgIsLoading;
							};
							planet.img.onerror = function() {
								console.log('failed to find texture for planet '+planet.name);
								planet.imgIsLoading = undefined;
							};
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
							generateMipmap : true
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
					vec3.sub(planet.sceneObj.uniforms.pos, planet.pos, orbitTarget.pos);
					quat.copy(planet.sceneObj.uniforms.angle, planet.angle);

					//only allow ring obj if there's a radii and a sphere
					//... this way i can just copy over the pos and angle
					if (planet.ringObj !== undefined) {
						vec3.sub(planet.ringObj.pos, planet.pos, orbitTarget.pos);
						quat.copy(planet.ringObj.angle, planet.sceneObj.uniforms.angle);
					}

					//calculate sun position for lighting
					for (var starIndex = 0; starIndex < orbitStarSystem.stars.length; ++starIndex) {
						var star = orbitStarSystem.stars[starIndex];
						var tmp = [];
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
					planet.sceneObj.uniforms.scaleExaggeration = planetScaleExaggeration;	
					
					if (picking) {
						pickObject.drawPlanet(planet.sceneObj, planet);
					} else {
						if (planet.type != 'barycenter') {
							planet.sceneObj.draw();
							overlayTexts.add(planet);
						}
					}

					//show latitude and longitude lines
					if (showLatAndLonLines && !picking) {
						vec3.copy(planetLatLonObj.pos, planet.sceneObj.uniforms.pos);
						quat.copy(planetLatLonObj.angle, planet.sceneObj.uniforms.angle);
						planetLatLonObj.uniforms.equatorialRadius = planet.equatorialRadius !== undefined ? planet.equatorialRadius : planet.radius;
						planetLatLonObj.uniforms.inverseFlattening = planet.inverseFlattening !== undefined ? planet.inverseFlattening : .5;
						planetLatLonObj.uniforms.scaleExaggeration = planetScaleExaggeration;
						planetLatLonObj.draw();
					}
				}
			}
		})();
		}

		//draw rings and transparent objects last
		//disable depth writing so orbits are drawn in front and behind them
		//do this last so (without depth writing) other planets in the background don't show up in front of the rings
		if (!picking) {
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
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
					var axis = [];
					var sunDir = [];
					vec3.sub(sunDir, planet.pos, orbitStarSystem.planets[orbitStarSystem.indexes.Sun].pos);
					vec3.normalize(sunDir, sunDir);
					vec3.quatZAxis(axis, planet.ringObj.angle);
					var axisDotSun = vec3.dot(axis, sunDir);
					var lookingAtLitSide = vec3.dot(viewfwd, axis) * axisDotSun;
					lookingAtLitSide = (lookingAtLitSide < 0 ? -1 : 1) * Math.pow(Math.abs(lookingAtLitSide), 1/4);	//preserve sign, so raise to odd power
					planet.ringObj.uniforms.lookingAtLitSide = Math.clamp(.5 + .5 * lookingAtLitSide, 0, 1);
					var viewDotSun = vec3.dot(viewfwd, sunDir);
					planet.ringObj.uniforms.backToFrontLitBlend = .5 - .5 * viewDotSun;	//clamp(sqrt(sqrt(dot( normalize(sun.pos - planet.pos), axis ) * dot( viewfwd, axis ))), 0., 1.) * -.5 + .5

					//have to recalculate these because the uniforms are shared, so they've been overwritten
					//...and the ring objs have to be drawn last for transparency reasons...

					//calculate sun position for lighting
					for (var starIndex = 0; starIndex < orbitStarSystem.stars.length; ++starIndex) {
						var star = orbitStarSystem.stars[starIndex];
						var tmp = [];
						vec3.sub(tmp, star.pos, planet.pos);
						planet.ringObj.uniforms.sunDir[0+3*starIndex] = tmp[0];
						planet.ringObj.uniforms.sunDir[1+3*starIndex] = tmp[1];
						planet.ringObj.uniforms.sunDir[2+3*starIndex] = tmp[2];
					}
					vec3.sub(planet.ringObj.uniforms.pos, planet.pos, orbitTarget.pos);
					quat.copy(planet.ringObj.uniforms.angle, planet.angle);

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
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			if (!planet.sceneObj || planet.visRatio < planetPointVisRatio) {
				if (picking) {
					//condition to skip moons ...
					if (planet.orbitVisRatio === undefined || planet.orbitVisRatio >= .03) { 
						vec3.sub(pointObj.attrs.vertex.data, planet.pos, orbitTarget.pos);
						pointObj.attrs.vertex.updateData();
						pickObject.drawPoints({
							sceneObj : pointObj,
							targetCallback : planet,
							pointSize : 4,
							pointSizeScaleWithDist : false
						});
					}
				} else if (showPlanetsAsDistantPoints) {
					vec3.sub(pointObj.attrs.vertex.data, planet.pos, orbitTarget.pos);
					pointObj.attrs.vertex.updateData();
					pointObj.draw({
						uniforms : {
							color : planet.color,
							pointSize : 4
						}
					});
				} else {	//if (showStars) {
					//draw as a pixel with apparent magnitude
					//use the starfield shader
					vec3.sub(pointObj.attrs.vertex.data, planet.pos, orbitTarget.pos);
					pointObj.attrs.vertex.updateData();
					pointObj.attrs.velocity.data[0] = planet.vel[0];
					pointObj.attrs.velocity.data[1] = planet.vel[1];
					pointObj.attrs.velocity.data[2] = planet.vel[2];
					pointObj.attrs.velocity.updateData();
					pointObj.attrs.absoluteMagnitude.data[0] = planet.absoluteMagnitude || 10;
					pointObj.attrs.absoluteMagnitude.updateData();
					pointObj.attrs.colorIndex.data[0] = planet.colorIndex || 0;
					pointObj.attrs.colorIndex.updateData();
					pointObj.draw({
						shader : starfield.colorIndexShader,
						//disable depth test too?
						blend : [gl.SRC_ALPHA, gl.ONE],
						uniforms : {
							colorIndexTex : 0,
							starTex : 1,
							pointSize : 1,
							pointSizeMax : 5,
							visibleMagnitudeBias : starsVisibleMagnitudeBias
						},
						texs : [starfield.colorIndexTex, starfield.starTex]
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
		if (mouseOverTarget !== undefined &&
			!picking)
		{
			gl.disable(gl.DEPTH_TEST);
			pointObj.attrs.vertex.data[0] = mouseOverTarget.pos[0] - orbitTarget.pos[0];
			pointObj.attrs.vertex.data[1] = mouseOverTarget.pos[1] - orbitTarget.pos[1];
			pointObj.attrs.vertex.data[2] = mouseOverTarget.pos[2] - orbitTarget.pos[2];
			pointObj.attrs.vertex.updateData();
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
		if (showOrbits &&
			!picking)
		{
			window.orbitPathsDrawn = 0;
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;

				if (planet.renderOrbit) {

					var semiMajorAxis = planet.keplerianOrbitalElements.semiMajorAxis;
					var eccentricity = planet.keplerianOrbitalElements.eccentricity;
					var distPeriapsis = semiMajorAxis * (1 + eccentricity);	//largest distance from the parent planet

					//vector from view to parent planet
					var parentPlanet = planet.parent;
					vec3.sub(delta, parentPlanet.pos, orbitTarget.pos);
					vec3.sub(delta, delta, glutil.view.pos);
					var deltaLength = vec3.length(delta);

					planet.orbitVisRatio = distPeriapsis / (deltaLength * tanFovY);

					if (planet.orbitVisRatio > planetPointVisRatio) {
						//recenter around orbiting planet
						vec3.sub(orbitPathSceneObj.pos, planet.parent.pos, orbitTarget.pos);

						orbitPathSceneObj.uniforms.color = planet.color;
						orbitPathSceneObj.uniforms.A = planet.keplerianOrbitalElements.A;
						orbitPathSceneObj.uniforms.B = planet.keplerianOrbitalElements.B;
						orbitPathSceneObj.uniforms.eccentricity = planet.keplerianOrbitalElements.eccentricity;
						orbitPathSceneObj.uniforms.eccentricAnomaly = planet.keplerianOrbitalElements.eccentricAnomaly;
						orbitPathSceneObj.uniforms.orbitType = orbitPathIndexForType[planet.keplerianOrbitalElements.orbitType];
						orbitPathSceneObj.uniforms.fractionOffset = planet.keplerianOrbitalElements.fractionOffset || 0;

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

		if (showGravityWell &&
			!picking)
		{
			//do transformation math in double
			//TODO just give gl-matrix a type param in its init
			var glMvMat = [];
			var viewAngleInvd = [];
			quat.conjugate(viewAngleInvd, glutil.view.angle);
			quat.normalize(viewAngleInvd, viewAngleInvd);	//normalize in double precision
			mat4.fromQuat(glMvMat, viewAngleInvd);
			var viewPosInvd = [];
			vec3.negate(viewPosInvd, glutil.view.pos);
			mat4.translate(glMvMat, glMvMat, viewPosInvd);

			var orbitBasis = [];
			var orbitBasisInv = [];
			var mvMat = [];
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;
				if (planet.radius === undefined) continue;
				//max radial dist is R * Math.pow(100, gravityWellRadialMaxLog100)
				if (planet.visRatio * Math.pow(100, gravityWellRadialMaxLog100) < .001) continue;

				mat4.identity(orbitBasis);
				mat4.identity(orbitBasisInv);
				for (var i = 0; i < 16; ++i) {
					mvMat[i] = glMvMat[i];
				}
				for (var i = 0; i < 3; ++i) {
					orbitBasis[i] = planet.orbitBasis[0][i];
					orbitBasis[4+i] = planet.orbitBasis[1][i];
					orbitBasis[8+i] = planet.orbitBasis[2][i];
				}
				//gravWellObjBasis = orbitBasis * gravityWellZScale * zOffsetByRadius * orbitBasis^-1 * planetPosBasis

				//calc this for non-sceneObj planets.  maybe I should store it as a member variable?
				var relPos = [
					planet.pos[0] - orbitTarget.pos[0],
					planet.pos[1] - orbitTarget.pos[1],
					planet.pos[2] - orbitTarget.pos[2]
				];

				//translate to planet center
				mat4.translate(mvMat, mvMat, relPos);

				//align gravity well with orbit plane
				mat4.multiply(mvMat, mvMat, orbitBasis);

				//align base gravity well with base of planet
				mat4.translate(mvMat, mvMat, [0,0,-planet.radius]);

				//scale for visual effect
				//too flat for earth, too sharp for sun
				if (gravityWellScaleFixed) {
					gravityWellTargetZScale = gravityWellScaleFixedValue;
				} else if (gravityWellScaleNormalized) {
					//causes larger wells to be much smaller and sharper ...
					//var gravityWellTargetZScale = 2000000 * Math.log(1 + z);
					//normalized visually per-planet.  scale is not 1-1 across planets
					gravityWellTargetZScale = 1 / orbitTarget.gravityWellScalar;
				}
				gravityWellZScale += .01 * (gravityWellTargetZScale - gravityWellZScale);
				if (gravityWellZScale !== gravityWellZScale) gravityWellZScale = 1;
				mat4.scale(mvMat, mvMat, [1,1,gravityWellZScale * planet.gravityWellScalar]);

				//apply orbit basis inverse
				mat4.transpose(orbitBasisInv, orbitBasis);
				mat4.multiply(mvMat, mvMat, orbitBasisInv);

				for (var i = 0; i < 16; ++i) {
					planet.gravWellObj.uniforms.mvMat[i] = mvMat[i];
				}

				planet.gravWellObj.draw();
			}
		}
	};
})();

//unitsToTest should be largest to smallest
function unitsToStr(measureInBaseUnits, otherUnitsInBaseUnits, unitsToTest) {
	var bestUnit = 'm';
	var bestDist = measureInBaseUnits;
	for (var i = 0; i < unitsToTest.length; ++i) {
		var unit = unitsToTest[i];
		var ratio = measureInBaseUnits / otherUnitsInBaseUnits[unit];
		if (ratio > .1) {
			bestUnit = unit;
			bestDist = ratio;
			break;
		}
	}
	return bestDist.toFixed(4)+' '+bestUnit;
}

var distanceToStrUnits = ['Mpc', 'lyr', 'AU', 'km', 'm'];
function distanceToStr(distInMeters) {
	return unitsToStr(distInMeters, metersPerUnits, distanceToStrUnits);
}

var timeToStrUnits = ['years', 'days', 'm', 'h', 's'];
var daysPerUnits = {
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
	$('#orbitTargetDistance').text(distanceToStr(orbitTargetDistance));
}

function refreshCurrentTimeText() {
if (astro) {
	$('#currentTimeText').text(astro.julianToCalendar(julianDate));
}
}

var primaryPlanetHorizonIDs = [10, 199, 299, 301, 399, 499, 599, 699, 799, 899, 999];

var ModifiedDepthShaderProgram;


function floatToGLSL(x) {
	x = ''+x;
	if (x.indexOf('.') == -1 && x.indexOf('e') == -1) x += '.';
	return x;
}

function mat3ToGLSL(m) {
	var s = 'mat3(';
	for (var i = 0; i < 9; ++i) {
		if (i > 0) s += ',';
		s += floatToGLSL(m[i]);
	}
	s += ')';
	return s;
}

//GLSL code generator
function unravelForLoop(varname, start, end, code) {
	var lines = [];
	for (var i = start; i <= end; ++i) {
		lines.push('#define '+varname+' '+i);
		lines.push(code);
		lines.push('#undef '+varname);
	}
	return lines.join('\n')+'\n';
};

var geodeticPositionCode = mlstr(function(){/*

//takes in a vec2 of lat/lon
//returns the ellipsoid coordinates in the planet's frame of reference
//I could move the uniforms here, but then the function would be imposing on the shader
#define M_PI 3.141592653589793115997963468544185161590576171875
vec3 geodeticPosition(vec2 latLon) {
	float phi = latLon.x * M_PI / 180.;
	float lambda = latLon.y * M_PI / 180.;
	float cosPhi = cos(phi);
	float sinPhi = sin(phi);
	float eccentricitySquared = (2. * inverseFlattening - 1.) / (inverseFlattening * inverseFlattening);
	float sinPhiSquared = sinPhi * sinPhi;
	float N = equatorialRadius / sqrt(1. - eccentricitySquared * sinPhiSquared);
	const float height = 0.;
	float NPlusH = N + height;	//plus height, if you want?  no one is using height at the moment.  heightmaps someday...
	return vec3(
		NPlusH * cosPhi * cos(lambda),
		NPlusH * cosPhi * sin(lambda),
		(N * (1. - eccentricitySquared) + height) * sinPhi);
}
*/});

var quatRotateCode = mlstr(function(){/*
vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
*/});

//used by milkyway and skycube
//which really both do the same thing: display the milky way in the foreground or background
var coordinateSystemCode = 
'const mat3 eclipticalToGalactic = '+mat3ToGLSL(eclipticalToGalacticTransform)+';\n'
+'const mat3 equatorialToEcliptical = '+mat3ToGLSL(equatorialToEclipticalTransform)+';\n';



$(document).ready(init);

function calendarToJulian(d) {
	var year = d.getUTCFullYear();
	var month = d.getUTCMonth();
	var day = d.getUTCDate();
	var hour = d.getUTCHours();
	var min = d.getUTCMinutes();
	var sec = d.getUTCSeconds();
//http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
//testing against results in Wikipedia page: http://en.wikipedia.org/wiki/Julian_day
// (notice I couldn't recreate these results with the algorithm on the same page)
//for UTC date 2000 jan 1 12:00:00 this gives 2451545 - right on
//for UTC date 2013 jan 1 00:30:00 this gives 2456294.0208333335 - right on
	var oneBasedMonth = 1 + month;
	var astronomicalYear = year < 0 ? year + 1 : year;
	if (oneBasedMonth <= 2) {	//jan or feb
		--astronomicalYear;
		oneBasedMonth += 12;
	}
	var a = Math.floor(astronomicalYear / 100);
	var b = Math.floor(a / 4);
	var c = 2 - a + b;
	var e = Math.floor(365.25 * (astronomicalYear + 4716));
	var f = Math.floor(30.6001 * (oneBasedMonth + 1));
	var jdn = c + day + e + f - 1524.5;
	jdn += (hour + (min + sec / 60) / 60) / 24;
	return jdn;
}

function init() {
	ui.init();

	//post-WebGL init:
	
	var depthConstant = 1e-6;//2 / Math.log(1e+7 + 1);
	ModifiedDepthShaderProgram = makeClass({
		super : glutil.ShaderProgram,
		init : function(args) {

			// maximizing depth range: http://outerra.blogspot.com/2012/11/maximizing-depth-buffer-range-and.html
			//  but, true to the original design, I actually want more detail up front.  Maybe I'll map it to get more detail up front than it already has?
			args.vertexCode = mlstr(function(){/*
uniform float zNear, zFar, depthConstant;
//float depthfunction(vec4 v) {
	//return (log(v.w + 1.) * depthConstant - 1.) * v.w;
	//return (2.0 * log(v.w / zNear) / log(zFar / zNear) - 1.) * v.w;
//}
#define depthfunction(x) ((x).z)
	*/}) + (args.vertexCode || '');
			if (args.uniforms === undefined) args.uniforms = {};
			args.uniforms.zNear = glutil.view.zNear;
			args.uniforms.zFar = glutil.view.zFar;
			args.uniforms.depthConstant = depthConstant;
			args.vertexPrecision = 'best';
			args.fragmentPrecision = 'best';
			ModifiedDepthShaderProgram.super.call(this, args);
		}
	});

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

	//julianDate = horizonsDynamicData.julianDate;	//get most recent recorded date from planets
	//TODO get current julian date from client
	// integrate forward by the small timestep between the two
	julianDate = calendarToJulian(new Date());

	initJulianDate = julianDate;
	refreshCurrentTimeText();

	solarSystem.copyPlanets(solarSystem.initPlanets, solarSystem.planets);

/*
	var texURLs = [];
	$.each(primaryPlanetHorizonIDs, function(_,horizonID) {
		var planet = solarSystem.planetForHorizonID[horizonID];
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
	$(texURLs).preload(function(){
		$('#loadingDiv').hide();
	}, function(percent){
		$('#loading').attr('value', parseInt(100*percent));
	});
*/
	$('#infoPanel').bind('touchmove', function(e) {
		if (!showBodyInfo) {
			if (!e.target.classList.contains('scrollable')) {
				//only touch devices will scroll the div even when browsers know not to
				//so for touch scroll, explicitly tell it not to scroll
				e.preventDefault();
			}
		}
	});

	setCSSRotation($('#toggleBodyInfo'), -90);
	//$('#infoPanel').css('height', $('#infoDiv').offset().top - $('#infoPanel').offset().top);
	//$('#infoPanel').css('bottom', -25);
	$('#toggleBodyInfo').click(function() {
		if (!showBodyInfo) {
			showBodyInfo = true;
			var infoDivTop = $('#infoPanel').offset().top;
			var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
			$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);

			$('#infoPanel')
				//go from bottom-aligned to top-algined
				.css('bottom', 'auto')
				.css('top', infoDivTop)
				.css('overflow', 'scroll')
				//interpolate up to the time controls ... or some fixed distance to it
				.stop(true)
				.animate(
					{
						top : infoDivDestTop+'px'
					},
					{
						duration : slideDuration,
						step : function(now, fx) {
							var frac = 1 - (now - infoDivDestTop) / (infoDivTop - infoDivDestTop);
							var degrees = 180 * frac - 90;
							setCSSRotation($('#toggleBodyInfo'), degrees);
						}
					}
				);
		} else {
			showBodyInfo = false;
			$('#infoPanel').css('height', '104px');
			var infoDivBottom = window.innerHeight - ($('#infoPanel').offset().top + $('#infoPanel').height());
			var infoDivDestBottom = -25;

			$('#infoPanel')
				//go from bottom-aligned to top-algined
				.css('top', 'auto')
				.css('bottom', infoDivBottom)
				.css('overflow', 'visible')
				//interpolate up to the time controls ... or some fixed distance to it
				.stop(true)
				.animate(
					{
						bottom : infoDivDestBottom+'px'
					},
					{
						duration : slideDuration,
						step : function(now, fx) {
							var frac = (now - infoDivDestBottom) / (infoDivBottom - infoDivDestBottom);
							var degrees = 180 * frac - 90;
							setCSSRotation($('#toggleBodyInfo'), degrees);
						}
					}
				);

		}
	});

	$('#menu').show();
	$('#timeDiv').show();
	//$('#infoPanel').show();	//don't show here -- instead show upon first setOrbitTarget, when the css positioning gets fixed

	colorShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex, 1.);
	gl_Position = projMat * vtx4;
	gl_PointSize = pointSize;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
void main() {
	gl_FragColor = color;
}
*/}),
		uniforms : {
			color : [1,1,1,1],
			pointSize : 4
		}
	});

	latLonShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;	//lat/lon pairs
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;
uniform float inverseFlattening;
uniform float scaleExaggeration;

*/}) + geodeticPositionCode + quatRotateCode + mlstr(function(){/*

void main() {
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
		*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
void main() {
	gl_FragColor = color;
}
		*/}),
		uniforms : {
			color : [1,1,1,1]
		}
	});

	calcTides.initGL();

	console.log('init overlay slider...');
	$('#overlaySlider').slider({
		range : 'max',
		width : '200px',
		min : 0,
		max : 100,
		value : 100 * heatAlpha,
		slide : function(event, ui) {
			heatAlpha = ui.value / 100;
			calcTides.updateHeatMapAlpha(heatAlpha);
		}
	});
	calcTides.updateHeatMapAlpha(heatAlpha);

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
			absoluteMagnitude : new glutil.ArrayBuffer({
				dim : 1,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			}),
			colorIndex : new glutil.ArrayBuffer({
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
	galaxies.init();

	//initialize here after webgl canvas is up	
	pickObject = new PickObject();
	
	//only do this after the orbit shader is made
	starSystemsExtra.initExoplanets();

	skyCube.init();

	initScene();
}

function setOrbitTarget(newTarget) {
	var selectingNewSystem = false;
	if (newTarget === undefined) {
		newTarget = orbitStarSystem.stars[0];
	}
	if (newTarget.isa && newTarget.isa(StarSystem)) {
		var targetSystem = newTarget;
		var i = 0;
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
	if (newTarget.isa && newTarget.isa(Planet)) {
		//if we're changing star systems...
		if (newTarget.starSystem !== orbitStarSystem) {
			selectingNewSystem = true;
			orbitStarSystem = newTarget.starSystem;
		}
	}

	//query the server for this small body's orbit data
	if (newTarget.isa && newTarget.isa(SmallBody)) {
		newTarget.onSelect();
	}

	if (selectingNewSystem) {
		if (resetDistanceOnSelectingNewSystem) {
			orbitTargetDistance = Math.max(100000, newTarget.radius || newTarget.equatorialRadius || 0);
			for (var i = 0; i < orbitStarSystem.planets.length; ++i) {
				var planet = orbitStarSystem.planets[i];
		
				//zoom out past orbit for planets
				//TODO don't do this for stars orbiting the milky way
				if (planet.sourceData !== undefined && planet.sourceData.semiMajorAxis !== undefined) {
					orbitTargetDistance = Math.max(orbitTargetDistance, planet.sourceData.semiMajorAxis);
				}
				if (planet.keplerianOrbitalElements !== undefined && planet.keplerianOrbitalElements.semiMajorAxis !== undefined) {
					orbitTargetDistance = Math.max(orbitTargetDistance, planet.keplerianOrbitalElements.semiMajorAxis);
				}
				
				refreshOrbitTargetDistanceText();
			}
		}
		orbitStarSystem.updatePlanetsPos();
	}

	if (orbitTarget !== undefined) vec3.sub(orbitOffset, orbitTarget.pos, newTarget.pos);
	orbitTarget = newTarget

	//reset orbit distance so it doesn't lock to surface immediately
	if (orbitGeodeticLocation !== undefined) {
		//TODO this is what kicks the next frame out of fixed-view mode,
		// and the code that does the kicking out resets the distance as if the user had zoomed out (i.e. relatively close)
		// so somehow impose this distance over the other distance
		orbitDistance = orbitTargetDistance = orbitTarget.radius * 10;

		//and re-orient the camera (because locking on surface had it flipped around / will flip it back around)
		quat.copy(orbitGeodeticLocation.lastViewAngle, glutil.view.angle); 
	}

	//refresh info div

	$('#orbitTargetText').text(orbitTarget.name);
	$('#infoDiv').empty();
	$('#infoDiv').append($('<hr>'));
	$('#infoDiv').append($('<div>', {text:'Type: '+orbitTarget.type}));
	if (orbitTarget.mass !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Mass: '+orbitTarget.mass+' kg'}));
	}
	if (orbitTarget.equatorialRadius !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Equatorial Radius: '+distanceToStr(orbitTarget.equatorialRadius)}));
	} else if (orbitTarget.radius !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Radius: '+distanceToStr(orbitTarget.radius)}));
	}
	if (orbitTarget.inverseFlattening !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Inverse Flattening: '+orbitTarget.inverseFlattening}));
	}
	if (orbitTarget.rotationPeriod !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Rotation Period: '+timeToStr(orbitTarget.rotationPeriod)}));
	}
	if (orbitTarget.ringRadiusRange !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Ring Min Radius: '+distanceToStr(orbitTarget.ringRadiusRange[0])}));
		$('#infoDiv').append($('<div>', {text:'Ring Max Radius: '+distanceToStr(orbitTarget.ringRadiusRange[1])}));
	}
	if (orbitTarget.keplerianOrbitalElements) {
		if (orbitTarget.keplerianOrbitalElements.semiMajorAxis) {
			$('#infoDiv').append($('<div>', {text:'Semi-Major Axis: '+distanceToStr(orbitTarget.keplerianOrbitalElements.semiMajorAxis)}));
		}
		if (orbitTarget.keplerianOrbitalElements.orbitType) {
			$('#infoDiv').append($('<div>', {text:'Orbit Type: '+orbitTarget.keplerianOrbitalElements.orbitType}));
		}
		if (orbitTarget.keplerianOrbitalElements.eccentricity) {
			$('#infoDiv').append($('<div>', {text:'Eccentricity: '+orbitTarget.keplerianOrbitalElements.eccentricity}));
		}
		if (orbitTarget.keplerianOrbitalElements.eccentricAnomaly) {
			$('#infoDiv').append($('<div>', {html:'Eccentric Anomaly: '+Math.deg(orbitTarget.keplerianOrbitalElements.eccentricAnomaly)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.longitudeOfAscendingNode) {
			$('#infoDiv').append($('<div>', {html:'Longitude of Ascending Node: '+Math.deg(orbitTarget.keplerianOrbitalElements.longitudeOfAscendingNode)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.argumentOfPeriapsis) {
			$('#infoDiv').append($('<div>', {html:'Argument of Pericenter: '+Math.deg(orbitTarget.keplerianOrbitalElements.argumentOfPeriapsis)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.inclination) {
			$('#infoDiv').append($('<div>', {html:'Inclination: '+Math.deg(orbitTarget.keplerianOrbitalElements.inclination)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.timeOfPeriapsisCrossing) {
			$('#infoDiv').append($('<div>', {html:'Time of Periapsis Crossing: '+timeToStr(orbitTarget.keplerianOrbitalElements.timeOfPeriapsisCrossing)}));
		}
		if (orbitTarget.keplerianOrbitalElements.orbitalPeriod) {
			$('#infoDiv').append($('<div>', {text:'Orbital Period: '+timeToStr(orbitTarget.keplerianOrbitalElements.orbitalPeriod)}));
		}
	}
	$('#infoDiv').append($('<br>'));
	$('#infoDiv').append($('<br>'));

	setTimeout(function() {
		//refresh offset if the infoPanel is hidden
		$('#infoPanel').css('bottom', -25);
		if (!showBodyInfo) {
			$('#infoPanel').css('height', '104px');
		} else {
			var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
			$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);
		}
		setTimeout(function() {
			$('#infoPanel').show();
		}, 1);
	}, 1);
}

//TODO use glutil.mouseDir?
/*
function mouseRay() {
	var viewX = glutil.view.pos[0];
	var viewY = glutil.view.pos[1];
	var viewZ = glutil.view.pos[2];
	var viewFwdX = -2 * (glutil.view.angle[0] * glutil.view.angle[2] + glutil.view.angle[3] * glutil.view.angle[1]);
	var viewFwdY = -2 * (glutil.view.angle[1] * glutil.view.angle[2] - glutil.view.angle[3] * glutil.view.angle[0]);
	var viewFwdZ = -(1 - 2 * (glutil.view.angle[0] * glutil.view.angle[0] + glutil.view.angle[1] * glutil.view.angle[1]));
	var viewRightX = 1 - 2 * (glutil.view.angle[1] * glutil.view.angle[1] + glutil.view.angle[2] * glutil.view.angle[2]);
	var viewRightY = 2 * (glutil.view.angle[0] * glutil.view.angle[1] + glutil.view.angle[2] * glutil.view.angle[3]);
	var viewRightZ = 2 * (glutil.view.angle[0] * glutil.view.angle[2] - glutil.view.angle[3] * glutil.view.angle[1]);
	var viewUpX = 2 * (glutil.view.angle[0] * glutil.view.angle[1] - glutil.view.angle[3] * glutil.view.angle[2]);
	var viewUpY = 1 - 2 * (glutil.view.angle[0] * glutil.view.angle[0] + glutil.view.angle[2] * glutil.view.angle[2]);
	var viewUpZ = 2 * (glutil.view.angle[1] * glutil.view.angle[2] + glutil.view.angle[3] * glutil.view.angle[0]);
	var aspectRatio = canvas.width / canvas.height;
	var mxf = mouse.xf * 2 - 1;
	var myf = 1 - mouse.yf * 2;
	var tanFovY = Math.tan(glutil.view.fovY * Math.PI / 360);
	var mouseDirX = viewFwdX + tanFovY * (viewRightX * mxf * aspectRatio + viewUpX * myf);
	var mouseDirY = viewFwdY + tanFovY * (viewRightY * mxf * aspectRatio + viewUpY * myf);
	var mouseDirZ = viewFwdZ + tanFovY * (viewRightZ * mxf * aspectRatio + viewUpZ * myf);
	var mouseDirLength = Math.sqrt(mouseDirX * mouseDirX + mouseDirY * mouseDirY + mouseDirZ * mouseDirZ);
	return [mouseDirX/mouseDirLength, mouseDirY/mouseDirLength, mouseDirZ/mouseDirLength];
}
*/

function initScene() {
	//gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL);
	gl.clearColor(0,0,0,0);

	//assign our initial orbiting solar system
	orbitStarSystem = solarSystem;

	var trackPlanetName = 'Earth';
	if ($.url().param('target') !== undefined) {
		trackPlanetName = $.url().param('target');
	}

	var trackPlanet = solarSystem.planets[solarSystem.indexes[trackPlanetName]];
	setOrbitTarget(trackPlanet);
	orbitTargetDistance = 2. * orbitTarget.radius;
	refreshOrbitTargetDistanceText();
	orbitDistance = orbitTargetDistance;

	var tmpQ = quat.create();
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			var tanHalfFov = Math.tan(glutil.view.fovY * Math.PI / 360);
			if (orbitGeodeticLocation !== undefined) {
				quat.setAxisAngle(tmpQ, [tanHalfFov * dy, tanHalfFov * dx, 0], rotAngle);
			} else {
				quat.setAxisAngle(tmpQ, [-dy, -dx, 0], rotAngle);
			}
			quat.multiply(glutil.view.angle, glutil.view.angle, tmpQ);
			quat.normalize(glutil.view.angle, glutil.view.angle);
		},
		zoom : function(zoomChange) {
			var scale = Math.exp(-orbitZoomFactor * zoomChange);
		
			if (orbitGeodeticLocation !== undefined) {
				glutil.view.fovY *= scale;
				glutil.updateProjection();
			} else {
				orbitTargetDistance = Math.min(orbitTargetDistance * scale, 4000 * metersPerUnits.Mpc);
			}
			
			refreshOrbitTargetDistanceText();
		},
		passiveMove : function() {
			//seems to not always update, but running this every frame is much slower
			pickObject.pick(false);
		},
		click : function() {
			if (mouse.isDragging) return;
			pickObject.pick(true);
		}
	});

	update();
}

function update() {

	/*
	some new thoughts:
	- only enable mouseover highlighting if orbit major axis (largest radius) exceeds point vis threshold
	*/
	overlayTexts.updateBegin();
	if (mouseOverTarget) overlayTexts.add(mouseOverTarget);
	
	//TODO instead of always showing the orbit target, show what it collapses into (moon -> earth -> solar system -> milky way)
	// but this would mean unifying all the overlay text show/hide conditions ...
	if (orbitTarget && orbitTarget.name && orbitTarget.pos) overlayTexts.add(orbitTarget);
	
	/* converage angle on target planet * /
	
	var delta = vec3.create();
	vec3.scale(delta, glutil.view.pos, -1);
	var fwd = vec3.create();
	vec3.quatZAxis(fwd, glutil.view.angle);
	vec3.scale(fwd, fwd, -1);
	var axis = vec3.create();
	vec3.cross(axis, fwd, delta);
	vec3.scale(axis, axis, 1/vec3.length(delta));	//divide out length of delta so axis length is sin(theta)
	var sinTheta = vec3.length(axis);
	var theta = 0;
	if (sinTheta > 1e-3) {
		vec3.scale(axis, axis, 1/sinTheta);	//normalize axis
		theta = Math.asin(Math.clamp(sinTheta, -1,1));
		var q = quat.create();
		var convergeAngleCoeff = .5;
		quat.setAxisAngle(q, axis, theta * convergeAngleCoeff);
		quat.mul(glutil.view.angle, glutil.view.angle, q);
	}
	/**/
	
	// track ball orbit
	//assumes z points away from the planet

	//while fixed, zoom affects fov, so if our fov surpasses the max then zoom us back out
	if (glutil.view.fovY > 120) {
		orbitDistance = orbitTargetDistance = orbitTarget.radius * .2;
	}

	//if we are close enough to the planet then rotate with it
	var fixViewToSurface = orbitTargetDistance < orbitTarget.radius * .1 &&
		(!orbitTarget.isa || !orbitTarget.isa(Galaxy));	//don't allow perspetive from galaxy "surfaces"
	if (fixViewToSurface && orbitGeodeticLocation === undefined) {
		var pos = [];
		vec3.add(pos, orbitTarget.pos, glutil.view.pos);
		
		orbitGeodeticLocation = solarSystemBarycentricToPlanetGeodetic(orbitTarget, pos);
		orbitGeodeticLocation.lastViewAngle = quat.create();
		quat.copy(orbitGeodeticLocation.lastViewAngle, glutil.view.angle);

		glutil.view.fovY = 120;
		
		//...and in one fell swoop, turn the camera around
		//TODO spread this out over a few frames
		var rot = [];
		quat.identity(rot);
		quat.rotateY(rot, rot, Math.PI);
		quat.multiply(glutil.view.angle, glutil.view.angle, rot);
	} else if (!fixViewToSurface && orbitGeodeticLocation !== undefined) {
		clearGeodeticLocation();
	}

	if (orbitGeodeticLocation !== undefined) {
		//fix the height at the surface
		//TODO spread this out over a few frames
		var height = (orbitTarget.equatorialRadius || orbitTarget.radius || 1000) * .01;//= orbitGeodeticLocation.height;
		orbitDistance = orbitTargetDistance = height;
		
		planetGeodeticToSolarSystemBarycentric(
			glutil.view.pos,
			orbitTarget,
			orbitGeodeticLocation.lat,
			orbitGeodeticLocation.lon,
			height);
		vec3.sub(glutil.view.pos, glutil.view.pos, orbitTarget.pos);
	} else {
		resetOrbitViewPos();
	}

	var orbitConvergeCoeff = .9;
	vec3.scale(orbitOffset, orbitOffset, orbitConvergeCoeff);
	{
		var logDist = Math.log(orbitDistance);
		var logTarget = Math.log(orbitTargetDistance);
		var coeff = .05;
		var newLogDist = (1 - coeff) * logDist + coeff * logTarget;
		orbitDistance = Math.exp(newLogDist);
	}

	if (!integrationPaused) {
		julianDate += integrateTimeStep;
		refreshCurrentTimeText();
	}
		
	if (julianDate !== lastJulianDate) {
		
		//recompute position by delta since init planets
		orbitStarSystem.updatePlanetsPos();

		//recompute angle based on sidereal rotation period
		for (var i = 0; i < orbitStarSystem.planets.length; ++i) {
			var planet = orbitStarSystem.planets[i];
			if (planet.rotationPeriod) {
				var angle = ((julianDate % planet.rotationPeriod) / planet.rotationPeriod) * 2 * Math.PI;
			
				//I really don't like this variable.  I don't think it should be used.
				if (planet.rotationOffset !== undefined) {
					angle += planet.rotationOffset;
				}
				
				var zrot = [0,0,0,1];
				quat.rotateZ(zrot, zrot, angle);
				quat.multiply(planet.angle, planet.tiltAngle, zrot);
			} else {
				quat.copy(planet.angle, planet.tiltAngle);
			}
		}

		//if we're not fixed on the surface but we are close enough to orbit then spin with the planet
		//TODO this is messing up on mercury, venus, pluto ... retrograde planets + mercury ... ?
		if (orbitTargetDistance < orbitTarget.radius * 10) {
			var orbitTargetAxis = [0,0,1];
			vec3.quatZAxis(orbitTargetAxis, orbitTarget.angle);
	
			var deltaJulianDate = julianDate - lastJulianDate;
			
			var deltaAngle = quat.create();
			quat.setAxisAngle(deltaAngle, orbitTargetAxis, deltaJulianDate * (orbitTarget.rotationPeriod || 0) * 2 * Math.PI);
		
			quat.multiply(glutil.view.angle, deltaAngle, glutil.view.angle); 
			quat.normalize(glutil.view.angle, glutil.view.angle);
			
			//reset angle (to keep targets from offsetting when the ffwd/rewind buttons are pushed)
			resetOrbitViewPos();
		}

		//if we are zoomed in then keep the view following the mouse over target
		if (orbitGeodeticLocation !== undefined &&
			mouseOverTarget !== undefined)
		{
			var viewFwd = [];
			vec3.quatZAxis(viewFwd, glutil.view.angle);
			vec3.scale(viewFwd, viewFwd, -1);
			vec3.normalize(viewFwd, viewFwd);

			var destViewFwd = [];
			vec3.sub(destViewFwd, mouseOverTarget.pos, glutil.view.pos); 
			vec3.normalize(destViewFwd, destViewFwd);
			
			var axis = [];
			vec3.cross(axis, viewFwd, destViewFwd);
			var sinTheta = vec3.length(axis);
			if (sinTheta > 1e-5) {
				var theta = Math.asin(Math.clamp(sinTheta, -1, 1));
				var rot = [];
				quat.setAxisAngle(rot, axis, theta);
				quat.multiply(glutil.view.angle, rot, glutil.view.angle);
				quat.normalize(glutil.view.angle, glutil.view.angle);
			}
		}
	}

	lastJulianDate = julianDate;

	glutil.scene.setupMatrices();
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	if (window.showPickScene) {
		pickObject.pick(false, true);
	} else {
		drawScene(false);
	}
	glutil.clearAlpha();

	requestAnimFrame(update);
	
	overlayTexts.updateEnd();
}
