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
						planet.sceneObj.draw();
						overlayTexts.add(planet);
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
						overlayTexts.add(planet);
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

var latitudeMin = -90;
var latitudeMax = 90;
var latitudeStep = 5;
var longitudeMin = -180;
var longitudeMax = 180;
var longitudeStep = 5;
var latitudeDivisions = Math.floor((latitudeMax-latitudeMin)/latitudeStep);
var longitudeDivisions = Math.floor((longitudeMax-longitudeMin)/longitudeStep);
var tideTexWidth = 128;
var tideTexHeight = 128;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;

	//fix side panel heights
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('height', window.innerHeight);
	});

	//fix info panel height
	if (showBodyInfo) {
		var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
		$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);
	}

	glutil.resize();

	//TODO fix mouseline to work with this
	//glutil.view.fovY = Math.min(1, canvas.height / canvas.width) * 90;

	/*
	var aspectRatio = canvas.width / canvas.height;
	var nearHeight = Math.tan(glutil.view.fovY * Math.PI / 360);
	var nearWidth = aspectRatio * nearHeight;
	*/

	/** /
	//setup infinite projection matrix
	//http://www.terathon.com/gdc07_lengyel.pdf
	//http://www.gamasutra.com/view/feature/131351/the_mechanics_of_robust_stencil_.php?page=2
	var epsilon = 0;
	mat4.identity(glutil.scene.projMat);
	glutil.scene.projMat[0] = glutil.view.zNear / nearWidth;
	glutil.scene.projMat[5] = glutil.view.zNear / nearHeight;
	glutil.scene.projMat[10] = epsilon-1;
	glutil.scene.projMat[11] = (epsilon - 2) * glutil.view.zNear;
	glutil.scene.projMat[14] = -1;
	glutil.scene.projMat[15] = 0;
	/**/

	/* or we could just linearly ramp the depth over the range we want, since the rest is going to full anyways * /
	//(-z - zNear) / (zFar - zNear) * 2 - 1
	//z * (-2 / (zFar - zNear)) + ( - 2 * zNear / (zFar - zNear) - 1)
	var resZNear = 1e+4;
	var resZFar = 1e+8;
	mat4.identity(glutil.scene.projMat);
	glutil.scene.projMat[0] = glutil.view.zNear / nearWidth;
	glutil.scene.projMat[5] = glutil.view.zNear / nearHeight;
	glutil.scene.projMat[10] = -2 * (resZFar - resZNear);
	glutil.scene.projMat[11] = -(1 + 2 * resZNear / (resZFar - resZNear));
	glutil.scene.projMat[14] = -1;
	glutil.scene.projMat[15] = 0;
	/**/
}

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

$(document).ready(init1);

var slideDuration = 500;
var slideWidth = 300;
var currentOpenSidePanelID = undefined;
var showBodyInfo = false;
var allSidePanelIDs = [];

function showSidePanel(sidePanelID) {
	if (currentOpenSidePanelID !== undefined) {
		hideSidePanel(currentOpenSidePanelID, true);
	}
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('z-index', 1);
	});
	var sidePanel = $('#'+sidePanelID);
	sidePanel.css('z-index', 2);
	sidePanel.show();
	$('#menu').animate(
		{
			left : slideWidth
		}, {
			duration : slideDuration,
			step : function(now, fx) {
				var degrees = now / slideWidth * 180;
				setCSSRotation($(this), degrees);
			}
		}
	);
	sidePanel.animate({left:0}, {duration:slideDuration});
	currentOpenSidePanelID = sidePanelID;
}

function hideSidePanel(sidePanelID, dontMoveOpenButton) {
	if (sidePanelID === currentOpenSidePanelID) currentOpenSidePanelID = undefined;
	var sidePanel = $('#'+sidePanelID);
	if (!dontMoveOpenButton) {
		$('#menu').animate(
			{
				left : 0
			},
			{
				duration : slideDuration,
				step : function(now, fx) {
					var degrees = (now / slideWidth) * 180;
					setCSSRotation($(this), degrees);
				}
			}
		);
	}
	sidePanel.animate({left:-slideWidth}, {
		duration:slideDuration,
		complete:function() {
			sidePanel.hide();
		}
	});
}

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

function getSmallBodyNameFromRow(row) {
	var name = row.name;
	if (row.idNumber) {
		name = row.idNumber+'/'+name;
	}
	return name;
}

function removeSmallBody(row) {
	var name = getSmallBodyNameFromRow(row);

	//only add if it's already there
	if (solarSystem.indexes[name] === undefined) return;
		
	var index = solarSystem.indexes[name];
	var planet = solarSystem.planets[index];
	solarSystem.initPlanets.splice(index, 1);
	solarSystem.planets.splice(index, 1);
	//now remap indexes
	for (var i = index; i < solarSystem.planets.length; ++i) {
		solarSystem.planets[index].index = i;
	}
	if (orbitTarget === planet) {
		setOrbitTarget(solarSystem.planets[solarSystem.indexes.Sun]);
	}
	//TODO destruct WebGL geometry?  or is it gc'd automatically?
	//now rebuild indexes
	solarSystem.indexes = {};
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		solarSystem.indexes[solarSystem.planets[i].name] = i;
	}
}

function addSmallBody(row) {
	var name = getSmallBodyNameFromRow(row);

	//only add if it's not there
	if (solarSystem.indexes[name] !== undefined) return solarSystem.planets[solarSystem.indexes[name]];

	//add the row to the bodies

	var index = solarSystem.planets.length;
	var planet = mergeInto(new Planet(), {
		name : name,
		isComet : row.bodyType == 'comet',
		isAsteroid : row.bodyType == 'numbered asteroid' || row.bodyType == 'unnumbered asteroid',
		sourceData : row,
		parent : solarSystem.planets[solarSystem.indexes.Sun],
		starSystem : solarSystem,
		index : index
	});

	solarSystem.planets.push(planet);

	//add copy to initPlanets for when we get reset() working
	solarSystem.initPlanets[index] = planet.clone();

	solarSystem.indexes[planet.name] = index;

	planet.initColorSchRadiusAngle();
	planet.initSceneLatLonLineObjs();
	calcKeplerianOrbitalElements(planet, false);
	recomputePlanetAlongOrbit(planet);

	return planet;
}

function init1() {
	allSidePanelIDs = [
		'mainSidePanel',
		'displayOptionsSidePanel',
		'overlaySidePanel',
		'solarSystemSidePanel',
		'smallBodiesSidePanel',
		'starSystemsSidePanel',
		'controlsSidePanel'
	];
	mainSidePanel = $('#mainSidePanel');

	//keep track of what menu is open
	//if none are open, open the main menu
	//if any are ... if it's not main
	$('#menu').click(function() {
		if (currentOpenSidePanelID === undefined) {
			showSidePanel('mainSidePanel');
		} else if (currentOpenSidePanelID == 'mainSidePanel') {
			hideSidePanel(currentOpenSidePanelID);
		} else {
			showSidePanel('mainSidePanel');
		}
	});

	$.each([
		{buttonID:'mainButtonDisplayOptions', divID:'displayOptionsSidePanel'},
		{buttonID:'mainButtonOverlay', divID:'overlaySidePanel'},
		{buttonID:'mainButtonSolarSystem', divID:'solarSystemSidePanel'},
		{buttonID:'mainButtonSmallBodies', divID:'smallBodiesSidePanel'},
		{buttonID:'mainButtonStarSystems', divID:'starSystemsSidePanel'},
		{buttonID:'mainButtonControls', divID:'controlsSidePanel'},
	], function(i, info) {
		$('#'+info.buttonID).click(function() {
			showSidePanel(info.divID);
		});
	});

	$('#reset').click(function() {
		integrationPaused = true;
		integrateTimeStep = defaultIntegrateTimeStep;
		julianDate = initJulianDate;
		refreshCurrentTimeText();
	});
	$('#play').click(function() {
		integrationPaused = false;
		integrateTimeStep = Math.abs(integrateTimeStep);
	});
	$('#reverse').click(function() {
		integrationPaused = false;
		integrateTimeStep = -Math.abs(integrateTimeStep);
	});
	$('#pause').click(function() {
		integrationPaused = true;
	});
	//fast/slow vs ffwd/rewind?
	$('#ffwd').click(function() {
		integrateTimeStep *= 2;
	});
	$('#rewind').click(function() {
		integrateTimeStep /= 2;
	});

	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute',
			background : 'red'
		}
	}).prependTo(document.body).get(0);

	$(canvas).disableSelection();

	try {
		glutil = new GLUtil({canvas:canvas});
		gl = glutil.context;
	} catch (e) {
		$('#menu').remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}

	$(document).keydown(function(e) {
		switch (e.keyCode) {
		case 32:	//space
			integrationPaused = !integrationPaused;
			break;
		case 38:	//up
			integrateTimeStep *= 2;
			break;
		case 40:	//down
			integrateTimeStep /= 2;
			break;
		case 39:	//right
			integrationPaused = false;
			integrateTimeStep = Math.abs(integrateTimeStep);
			break;
		case 37:	//left
			integrationPaused = false;
			integrateTimeStep = -Math.abs(integrateTimeStep);
			break;
		}
	});

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

	// overlay side panel


	var overlaySidePanelMetric = $('#overlaySidePanelMetric');
	$.each(metricInfos, function(metricIndex, metricInfo) {
		var radio = $('<input>', {
			type : 'radio',
			name : 'calculationMetric',
			value : metricIndex,
			click : function() {
				metric = new metricInfo.classObj();
				calcTides.invalidateForces();
			}
		})
			.attr('name', 'metric')
			.appendTo(overlaySidePanelMetric);
		if (metric.isa(metricInfo.classObj)) radio.attr('checked', 'checked');
		$('<span>', {text:metricInfo.name}).appendTo(overlaySidePanelMetric);
		$('<br>').appendTo(overlaySidePanelMetric);
	});

	var overlaySidePanelContents = $('#overlaySidePanelContents');
	$('<span>', {text:'Overlay:'}).appendTo(overlaySidePanelContents);
	$('<br>').appendTo(overlaySidePanelContents);
	$.each(displayMethods, function(displayMethodIndex,thisDisplayMethod) {
		var radio = $('<input>', {
			type : 'radio',
			name : 'displayMethods',
			value : displayMethodIndex,
			click : function() {
				displayMethod = thisDisplayMethod;
				calcTides.invalidateForces();
			}
		})
			.attr('name', 'display')
			.appendTo(overlaySidePanelContents);
		if (thisDisplayMethod == displayMethod) radio.attr('checked', 'checked');
		$('<span>', {text:thisDisplayMethod}).appendTo(overlaySidePanelContents);
		$('<br>').appendTo(overlaySidePanelContents);
	});
	$('<br>').appendTo(overlaySidePanelContents);
	$('<span>', {text:'Influencing Planets:'}).appendTo(overlaySidePanelContents);
	$('<br>').appendTo(overlaySidePanelContents);

	//add radio buttons hierarchically ...
	var overlayControlsForPlanets = {};

	var HierarchicalCheckboxControl = makeClass({
		/*
		args:
			title		<- used to identify this checkbox
			change		<- callback upon changing the checkbox value
			clickTitle	<- callback upon clicking the title
			isChecked
			... and anything else that can be referenced through this.args
		*/
		init : function(args) {
			this.args = args;
			this.childControls = [];
			this.div = $('<div>', {
				css : {paddingLeft:'5px'}
			});

			var thiz = this;
			this.checkbox = $('<input>', {
				type : 'checkbox',
				change : function() {
					//refresh all parent controls' moon checkboxes -- to whiteout or grey them
					for (var c = thiz.parentControls; c; c = c.parentControls) {
						c.recomputeMoonCheckbox();
					}

					args.change.call(thiz);
				}
			})
				.prop('checked', args.isChecked)
				.appendTo(this.div);

			$('<span>', {
				text : args.title,
				css : {
					textDecoration : 'underline',
					cursor : 'pointer',
				},
				click : function(e) {
					args.clickTitle.call(thiz);
				}
			}).appendTo(this.div);

			this.toggleChildDiv = $('<span>', {
				css : {
					paddingLeft : '10px',
					cursor : 'pointer'
				}
			}).appendTo(this.div);

			this.moonCheckbox = $('<input>', {
				type : 'checkbox',
				change : function() {
					//select / deselect all children
					thiz.setAllChildren($(this).prop('checked'));
				}
			})
				.prop('checked', 1)
				.appendTo(this.toggleChildDiv);

			$('<span>', {
				css : {
					cursor : 'pointer'
				},
				text : '...',
				click : function() {
					if (thiz.childDiv.css('display') == 'none') {
						thiz.childDiv.show();
					} else {
						thiz.childDiv.hide();
					}
				}
			}).appendTo(this.toggleChildDiv);

			$('<br>').appendTo(this.div);

			this.childDiv = $('<div>').appendTo(this.div);

		},
		addChild : function(childControl) {
			childControl.div.appendTo(this.childDiv);
			childControl.parentControls = this;
			this.childControls.push(childControl);
		},
		setAllChildren : function(checked) {
			for (var i = 0; i < this.childControls.length; ++i) {
				var ch = this.childControls[i].checkbox;
				if (checked) {
					if (!ch.prop('checked')) { ch.prop('checked', 1); ch.trigger('change'); }
				} else {
					if (ch.prop('checked')) { ch.prop('checked', 0); ch.trigger('change'); }
				}
				this.childControls[i].setAllChildren(checked);
			}
		},
		recomputeMoonCheckbox : function() {
			var numChecked = 0;
			var total = 0;
			for (var i = 0; i < this.childControls.length; ++i) {
				++total;
				//check the child
				if (this.childControls[i].checkbox.prop('checked')) {
					++numChecked;
				}
				//check the child's children if they exist
				if (this.childControls[i].childControls.length > 0) {
					if (this.childControls[i].moonCheckbox.prop('checked')) {
						if (this.childControls[i].moonCheckbox.prop('indeterminate')) {
							numChecked += .5;
						} else {
							++numChecked;
						}
					}
					++total;
				}
			}
			if (numChecked == 0) {
				this.moonCheckbox.prop('checked', 0)
					.prop('indeterminate', 0);
			} else if (numChecked == total) {
				this.moonCheckbox.prop('checked', 1)
					.prop('indeterminate', 0);
			} else {
				this.moonCheckbox.prop('checked', 1)
					.prop('indeterminate', 1);
			}
		}
	});

	//TODO add star systems
	$.each(solarSystem.planets, function(planetIndex, planet) {
		//if any other planet doesn't have recorded mass then skip it
		if (planet.mass === undefined) return;

		var parentPlanet = planet.parent;
		if (parentPlanet !== undefined) {
			if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
		}

		var controls = new HierarchicalCheckboxControl({
			title : planet.name,
			isChecked : true,
			change : function() {
				planetInfluences[this.args.planetIndex] = this.checkbox.is(':checked');
				calcTides.invalidateForces();
			},
			clickTitle : function() {
				setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
			},
			planetIndex : planetIndex
		});

		//add to parent or to the side panel
		if (parentPlanet === undefined) {
			controls.div.appendTo(overlaySidePanelContents);
		} else {
			overlayControlsForPlanets[parentPlanet.index].addChild(controls);
		}

		planetInfluences[planetIndex] = true;

		overlayControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
	});

	for (var planetIndex in overlayControlsForPlanets) {
		var planetIndex = +planetIndex;
		var controls = overlayControlsForPlanets[planetIndex];

		controls.recomputeMoonCheckbox();

		//if a planet didn't get any children, hide its 'toggle child div'
		if (controls.childDiv.children().length == 0) {
			controls.toggleChildDiv.hide();
			continue;
		}

		//if a planet got children and it's not the sun or earth then hide by default
		if (planetIndex !== solarSystem.indexes.Sun &&
			planetIndex !== solarSystem.indexes.Earth)
		{
			controls.childDiv.hide();
		}
	}

	$('<br>').appendTo(overlaySidePanelContents);


	// display options side panel

	var radioGroups = [
		['gravityWellScaleNormalized', 'gravityWellScaleFixed'],
		['overlayShowOrbitTarget', 'overlayShowCurrentPosition']
	];

	$.each([
		'showLinesToOtherPlanets',
		'showVelocityVectors',
		'showRotationAxis',
		'showOrbitAxis',
		'showLatAndLonLines',
		'showGravityWell',
		'showNames',				//overlays.js
		'showOrbits',
		'showSmallBodies',			//smallbodies.js
		'allowSelectSmallBodies',	//smallbodies.js
		'showStars',				//in starfield.js
		'allowSelectStars',			//in starfield.js
		'showGalaxies',				//in galaxies.js
		'allowSelectGalaxies',
		'showPlanetsAsDistantPoints',
		//radio
		'gravityWellScaleNormalized',
		'gravityWellScaleFixed',
		//radio
		'overlayShowOrbitTarget',
		'overlayShowCurrentPosition'
	], function(_, toggle) {
		var checkbox = $('#'+toggle);
		if (window[toggle]) checkbox.attr('checked', 'checked');
		checkbox.change(function() {
			window[toggle] = checkbox.is(':checked');
			
			var found = false;
			for (var i = 0; i < radioGroups.length; ++i) {
				var group = radioGroups[i];
				for (var j = 0; j < group.length; ++j) {
					if (group[j] == toggle) {
						for (var k = 0; k < group.length; ++k) {
							if (k == j) continue;
							window[group[k]] = false;
						}
						found = true;
						break;
					}
				}
				if (found) break;
			}
		});
	});

	$.each([
		'gravityWellScaleFixedValue',
		'starsVisibleMagnitudeBias',	//in starfield.js
		'planetScaleExaggeration'
	], function(_, toggle) {
		(function(){
			var textfield = $('#'+toggle);
			textfield.val(window[toggle]);
			textfield.change(function() {
				window[toggle] = textfield.val();
			});
		})();
	});


	// celestial bodies side panel


	var solarSystemSidePanel = $('#solarSystemSidePanel');

	//add radio buttons hierarchically ...
	var celestialBodiesControlsForPlanets = {};

	var cometParent;

	$.each(solarSystem.planets, function(planetIndex,planet) {

		var parentPlanet = planet.parent;
		if (parentPlanet !== undefined) {
			if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
		}

		var controls = new HierarchicalCheckboxControl({
			title : planet.name,
			isChecked : !planet.hide,
			change : function() {
				solarSystem.planets[this.args.planetIndex].hide = !this.checkbox.is(':checked');
			},
			clickTitle : function() {
				setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
			},
			planetIndex : planetIndex
		});

		if (parentPlanet === undefined) {
			controls.div.appendTo($('#celestialBodiesVisiblePlanets'));
		} else {
			celestialBodiesControlsForPlanets[parentPlanet.index].addChild(controls);
		}

		celestialBodiesControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
	});

	if (cometParent) cometParent.recomputeMoonCheckbox();
	for (var planetIndex in celestialBodiesControlsForPlanets) {
		var planetIndex = +planetIndex;
		var controls = celestialBodiesControlsForPlanets[planetIndex];

		controls.recomputeMoonCheckbox();

		//if a planet didn't get any children, hide its 'toggle child div'
		if (controls.childDiv.children().length == 0) {
			controls.toggleChildDiv.hide();
			continue;
		}

		//if a planet got children and it's not the sun or earth then hide by default
		if (planetIndex !== solarSystem.indexes.Sun &&
			planetIndex !== solarSystem.indexes.Earth)
		{
			controls.childDiv.hide();
		}
	}

	$('<br>').appendTo($('#celestialBodiesVisiblePlanets'));

	//these are added to the end of the result
	//they should get greyed upon new query (search, prev, next click)
	//they should be regenerated upon new content
	//they should be re-enabled upon error
	var nextButton = undefined;
	var prevButton = undefined;

	var searchResults = [];

	var searchLastID = 0;	//uid to inc so if we search twice, the first can be invalidated

	//dataSource = 'remote' for remote queries, 'local' for the currently-selected planets
	var processSearch = function(pageIndex, dataSource) {
		var button = $('#celestialBodiesSearch');
		var searchText = $('#celestialBodiesSearchText');
		var searchStr = searchText.val();
		button.prop('disabled', 1);
		searchText.val('searching...');
		searchText.prop('disabled', 1);

		if (prevButton) prevButton.prop('disabled', 1);
		if (nextButton) nextButton.prop('disabled', 1);

		var searchID = ++searchLastID;

		var processResults = function(results) {
			if (searchID < searchLastID-1) return;	//someone else has searched

			searchText.val(searchStr);
			searchText.prop('disabled', 0);
			button.prop('disabled', 0);

			$('#celestialBodiesSearchToggleAll').prop('checked', 0);	//disable too?

			var pageSize = 20;	//fixed atm
			var pageMax = Math.floor((results.count-1) / pageSize);

			searchResults = [];

			var resultsDiv = $('#celestialBodiesSearchResults');
			resultsDiv.empty();
			$.each(results.rows, function(i,row) {
				var rowDiv = $('<div>');
				rowDiv.appendTo(resultsDiv);

				var name = row.name;
				if (row.idNumber) {
					name = row.idNumber+'/'+name;
				}

				var titleSpan;
				var checkbox = $('<input>', {
					type : 'checkbox',
					change : function() {
						if (!$(this).is(':checked')) {	//uncheck checkbox => remove planet
							removeSmallBody(row);
							titleSpan.css({textDecoration:'', cursor:''});
						} else {	//check checkbox => add planet
							addSmallBody(row);
							titleSpan.css({textDecoration:'underline', cursor:'pointer'});
						}
					}
				})
					.prop('checked', solarSystem.indexes[name] !== undefined)
					.appendTo(rowDiv);

				titleSpan = $('<span>', {
					text : name,
					click : function() {
						var targetPlanet = solarSystem.planets[solarSystem.indexes[name]];
						if (targetPlanet !== undefined) setOrbitTarget(targetPlanet);
					}
				}).appendTo(rowDiv);
				//TODO put an 'add' button next to each
				//on clicking it, add the body to the planet list
				//and repopulate the 'extra' div

				searchResults.push({
					checkbox : checkbox,
					div : rowDiv,
					data : row,
					name : name
				});

			});

			//and add prev/next/pages if there is
			if (pageMax > 0) {
				var changePage = function(dir) {
					//TODO remove or grey out results as well?
					if (nextButton) nextButton.prop('disabled', 1);
					if (prevButton) prevButton.prop('disabled', 1);
					processSearch(pageIndex+dir, dataSource);
				};
				if (pageIndex > 0) {
					prevButton = $('<button>', {
						text : 'Prev',
						click : function() {
							changePage(-1);
						}
					}).appendTo($('#celestialBodiesSearchResults'));
				}
				if (pageIndex < pageMax) {
					nextButton = $('<button>', {
						text : 'Next',
						click : function() {
							changePage(1);
						}
					}).appendTo($('#celestialBodiesSearchResults'));
				}
				$('<span>', {text:(pageIndex+1)+' of '+pageMax}).appendTo($('#celestialBodiesSearchResults'));
			}
		};

		if (dataSource == 'remote') {
			$.ajax({
				url : '/solarsystem/jpl-ssd-smallbody/search.lua',
				dataType : 'json',
				data : {
					comet : $('#celestialBodiesSearchComets').prop('checked')?1:0,
					numbered : $('#celestialBodiesSearchNumbered').prop('checked')?1:0,
					unnumbered : $('#celestialBodiesSearchUnnumbered').prop('checked')?1:0,
					text : searchStr,
					page : pageIndex+1	//1-based
				},
				cache : false
			}).error(function() {
				searchText.val(searchStr);
				searchText.prop('disabled', 0);
				button.prop('disabled', 0);
				//TODO animate background color of search text
				if (prevButton) prevButton.prop('disabled', 0);
				if (nextButton) nextButton.prop('disabled', 0);

				var warning = $('<div>', {text:'Connection Failed!', css:{color:'red'}});
				$('#celestialBodiesSearchWarning').after(warning);
				setTimeout(function() {
					//after five seconds, fade away
					warning.animate({
						height : 0,
						opacity : 0
					}, {
						duration : 500,
						complete : function() {
							warning.remove();
						}
					});
				}, 3000);
			}).done(processResults);
		} else if (dataSource == 'local') {
			var rows = [];
			var searchingComets = $('#celestialBodiesSearchComets').prop('checked');
			var searchingNumbered = $('#celestialBodiesSearchNumbered').prop('checked');
			var searchingUnnumbered = $('#celestialBodiesSearchUnnumbered').prop('checked');
			for (var i = 0; i < solarSystem.planets.length; ++i) {
				var planet = solarSystem.planets[i];
				var row = planet.sourceData;
				if (row) {
					if ((row.bodyType == 'comet' && searchingComets) ||
						(row.bodyType == 'numbered asteroid' && searchingNumbered) ||
						(row.bodyType == 'unnumbered asteroid' && searchingUnnumbered))
					{
						rows.push(row);
					}
				}
			}
			
/*hack for debugging* /
rows = [];
if (false) {
   //test case for hyperbolic
	rows.push({
		"perihelionDistance" : 209232954020.56,
		"inclination" : 2.2519519080902,
		"timeOfPerihelionPassage" : 2456955.82823,
		"argumentOfPeriapsis" : 0.042602963442406,
		"bodyType" : "comet",
		"orbitSolutionReference" : "JPL 101",
		"epoch" : 56691,
		"eccentricity" : 1.00074241,
		"idNumber" : "C",
		"name" : "2013 A1 (Siding Spring)",
		"longitudeOfAscendingNode" : 5.2530270564044
	});
//target C1/2013 Siding Spring data: on 2014-11-15 10:30:00
//position m: 	64031628815.774, -196629902235.24, 57392197050.865
//velocity m/day: -1916173862.297, -455287182.34414, 2315832701.4279
	
	//test case for elliptical
	rows.push({
		"perihelionDistance" : 87661077532.81,
		"inclination" : 2.8320181936429,
		"timeOfPerihelionPassage" : 2446467.39532,
		"argumentOfPeriapsis" : 1.9431185149437,
		"bodyType" : "comet",
		"orbitSolutionReference" : "JPL J863/77",
		"epoch" : 49400,
		"eccentricity" : 0.96714291,
		"idNumber" : "1P",
		"name" : "Halley",
		"longitudeOfAscendingNode" : 1.0196227452785
	});
}
if (true) {	//recent asteroid passing by
	rows.push({
		"meanAnomaly":6.1790425090414,
		"longitudeOfAscendingNode":2.2116873367796,
		"inclination":0.41440347267775,
		"name":"2004 BL86",
		"eccentricity":0.40307315,
		"idNumber":"357439",
		"bodyType":"numbered asteroid",
		"orbitSolutionReference":"JPL 34",
		"absoluteMagnitude":19.1,
		"argumentOfPeriapsis":5.4324242142291,
		"magnitudeSlopeParameter":0.15,
		"semiMajorAxis":224726235521.07,
		"epoch":57000
	});
}
/**/
			
			processResults({rows:rows, count:rows.length});
		} else {
			throw "got an unknown data source request " + dataSource;
		}
	};

	$('#celestialBodiesSearchText').keydown(function(e){
		if (e.keyCode == 13) {
			$('#celestialBodiesSearch').trigger('click');
		}
	});

	//change a check box, immediately update search results
	$.each([
		$('#celestialBodiesSearchComets'),
		$('#celestialBodiesSearchNumbered'),
		$('#celestialBodiesSearchUnnumbered'),
		$('#celestialBodiesSearchVisible')
	], function(_,checkbox) {
		checkbox.change(function() {
			$('#celestialBodiesSearch').trigger('click');
		});
	});

	$('#celestialBodiesSearchToggleAll').click(function() {
		var checked = $(this).is(':checked');
		$.each(searchResults, function(i,result) {
			if (result.checkbox.prop('checked') != checked) {
				result.checkbox.trigger('click');
			}
		});
	});

	$('#celestialBodiesSearch').click(function() {
		processSearch(0,
			$('#celestialBodiesSearchVisible').prop('checked') ? 'local' : 'remote'
		);
	});
	$('#celestialBodiesSearch').trigger('click');	//fire one off


	// rest of the init


	$(window).resize(resize);
	resize();

	//julianDate = horizonsDynamicData.julianDate;	//get most recent recorded date from planets
	//TODO get current julian date from client
	// integrate forward by the small timestep between the two
	julianDate = calendarToJulian(new Date());

	initJulianDate = julianDate;
	refreshCurrentTimeText();

	solarSystem.copyPlanets(solarSystem.initPlanets, solarSystem.planets);

	init2();
}

//call this once we get back our star data
function initStarsControls() {
	//TODO maybe some pages or something for this, like the asteroid/smallbody search?
	for (var i = 0; i < starSystems.length; ++i) {
		(function(){
			var starSystem = starSystems[i];
			$('<div>', {
				css : {
					textDecoration : 'underline',
					cursor : 'pointer',
					paddingLeft : '10px'
				},
				click : function() {
					setOrbitTarget(starSystem);
				},
				text : starSystem.name
			}).appendTo($('#starSystemContents'));
		})();
	}
}

function init2() {
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

	init3();
}

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

function initExoplanets() {
	var processResults = function(results) {
		//process results
		$.each(results.systems, function(i,systemInfo) {
			var systemName = assertExists(systemInfo, 'name');

			var rightAscension = systemInfo.rightAscension;
			if (rightAscension === undefined) {
				console.log('failed to find right ascension for system '+systemName);
				return;
			}
			var declination = systemInfo.declination;
			if (declination === undefined) {
				console.log('failed to find declination for system '+systemName);
				return;
			}
			var cosRA = Math.cos(rightAscension);
			var sinRA = Math.sin(rightAscension);
			var cosDec = Math.cos(declination);
			var sinDec = Math.sin(declination);
			//convert to coordinates
			var pos = [];
			pos[0] = cosRA * cosDec;
			pos[1] = sinRA * cosDec;
			pos[2] = sinDec;
			//rotate for earth's tilt
			var epsilon = Math.rad(23 + 1/60*(26 + 1/60*(21.4119)));
			var cosEps = Math.cos(epsilon);
			var sinEps = Math.sin(epsilon);
			var yn = cosEps * pos[1] + sinEps * pos[2];
			pos[2] = -sinEps * pos[1] + cosEps * pos[2];
			pos[1] = yn;
			//distance
			var distance = systemInfo.distance;
			if (distance === undefined) {
//				console.log('failed to find distance for system '+systemName);
				return;
			}
			pos[0] *= distance;
			pos[1] *= distance;
			pos[2] *= distance;


			var starSystem = new StarSystem();
			starSystem.name = systemName;
			starSystem.sourceData = systemInfo;
			vec3.copy(starSystem.pos, pos);

			//TODO absoluate magnitude of the collective system (sum of all parts?)

			var minAbsMag = undefined;

			$.each(assertExists(systemInfo, 'bodies'), function(j, bodyInfo) {

				var name = assertExists(bodyInfo, 'name');
				var radius = bodyInfo.radius;
				if (radius === undefined) {
					if (bodyInfo.type !== 'barycenter') {
						//console.log('no radius for body '+name);
						//if planets don't have radii they can be estimated by the planet's mass or distance from sun or both?
						radius = 7e+7;	// use a jupiter radius
					}
				}

				var mass = bodyInfo.mass;
				if (mass === undefined) {
					if (bodyInfo.density && bodyInfo.radius) {
						mass = 4/3 * Math.PI * bodyInfo.density * radius * radius * radius;
					}
				}
				if (mass === undefined) {
					//this prevents us from an orbit ..
					//console.log('no mass for body '+name);
					mass = 2e+27;	//use a jupiter mass
				}

				var body = mergeInto(new Planet(), {
					type : assertExists(bodyInfo, 'type'),
					name : name,
					mass : mass,
					radius : radius,	//for planets/suns this is the radius of the planet.  does it exist for barycenters?
					//visualMagnitude : bodyInfo.visualMagnitude,	// TODO convert to absolute magnitude for star shader?
					//spectralType : bodyInfo.spectralType,			// or this one?
					//temperature : bodyInfo.temperature,			// or this one?
					sourceData : bodyInfo,
					parent : bodyInfo.parent,
					starSystem : starSystem,
					hide : bodyInfo.type === 'barycenter',
					isExoplanet : true
				});

				//hacks for orbit
				bodyInfo.longitudeOfAscendingNode = bodyInfo.longitudeOfAscendingNode || 0;
				bodyInfo.inclination = bodyInfo.inclination || 0;
				bodyInfo.eccentricity = bodyInfo.eccentricity || 0;
				bodyInfo.meanAnomaly = Math.PI * 2 * (j + Math.random()) / systemInfo.bodies.length;

				starSystem.planets.push(body);
				if (body.type == 'star') starSystem.stars.push(body);
				
				if (bodyInfo.visualMagnitude !== undefined) {
					//absolute magnitude based on visual magnitude
					body.magnitude = bodyInfo.visualMagnitude - 5 * (Math.log10(distance / metersPerUnits.pc) - 1);
					//...???
					minAbsMag = minAbsMag === undefined ? body.magnitude : Math.min(minAbsMag, body.magnitude);
				}
			});
			if (minAbsMag !== undefined) {
				starSystem.magnitude = minAbsMag;
			}

			starSystem.doneBuildingPlanets();

			for (var i = 0; i < starSystem.planets.length; ++i) {
				var body = starSystem.planets[i];

				//further hacks for orbit, now that the parent pointer has been established
				if (body.sourceData.semiMajorAxis === undefined) {
					if (body.parent &&
						body.parent.type === 'barycenter' &&
						//this is what I don't like about the open exoplanet database:
						//sometimes 'separation' or 'semimajoraxis' is the distance to the parent
						// sometimes it's the distance to the children
						(body.parent.sourceData.separation !== undefined ||
						body.parent.sourceData.semiMajorAxis !== undefined))
					{
						body.sourceData.semiMajorAxis = body.parent.sourceData.semiMajorAxis ||
							body.parent.sourceData.separation || 0;
						//TODO also remove semiMajorAxis from barycenter so it doesn't get offset from the planet's center?
						// also -- do this all in exoplanet file preprocessing?
					}
				}
				//longitude of periapsis = longitude of ascending node + argument of periapsis
				//argument of periapsis = longitude of periapsis - longitude of ascending node
				if (body.sourceData.longitudeOfPeriapsis === undefined) {
					body.sourceData.longitudeOfPeriapsis = 0;
				}
				body.sourceData.argumentOfPeriapsis = body.sourceData.longitudeOfPeriapsis - body.sourceData.longitudeOfAscendingNode;
				if (body.sourceData.argumentOfPeriapsis == 0) {
					body.sourceData.argumentOfPeriapsis = Math.PI * 2 * i / starSystem.planets.length;
				}

				vec3.add(body.pos, body.pos, starSystem.pos);

				body.initColorSchRadiusAngle();
				body.initSceneLatLonLineObjs();

				if (body.pos[0] !== body.pos[0]) {
					console.log('system '+starSystem.name+' planet '+body.name+' has bad pos');
				}

				calcKeplerianOrbitalElements(body, false);
			}

			starSystem.initPlanets = starSystem.clonePlanets();
			starSystemForNames[starSystem.name] = starSystem;
			starSystem.index = starSystems.length;
			starSystems.push(starSystem);
		});

		//now that we've built all our star system data ... add it to the star field
		starfield.addStarSystems();
	};

	//just over a meg, so I might as well ajax it
	var exoplanetURL = 'exoplanet/openExoplanetCatalog.json';
	$.ajax({
		url : exoplanetURL,
		dataType : 'json'
	}).error(function() {
		console.log('failed to get exoplanets from '+exoplanetURL+' , trying again...');
		setTimeout(function() {
			initExoplanets();
		}, 5000);
	}).done(function(data) {
		processResults(data);
	});
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

function init3() {

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

	calcTides.init3();

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
	//initialize here after webgl canvas is up	
	pickObject = new PickObject();

	starSystemsExtra.initSolarSystemImgUrls();
				
	initOrbitPaths();

	//ring texture for Jupiter
	//http://www.celestiamotherlode.net/catalog/jupiter.php
	(function(){
		var jupiterRingShader = new ModifiedDepthShaderProgram({
			//vertex code matches Saturn
			vertexCode :
quatRotateCode
+ mlstr(function(){/*
#define M_PI 3.1415926535897931
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float ringMinRadius;
uniform float ringMaxRadius;
uniform vec3 pos;
uniform vec4 angle;
varying vec2 texCoordv;
varying vec3 worldPosv;	//varying in world coordinates

void main() {
	texCoordv = vertex;

	//vertex is the uv of the texcoord wrapped around the annulus
	//u coord is radial, v coord is angular
	float rad = mix(ringMinRadius, ringMaxRadius, vertex.x);
	float theta = 2. * M_PI * vertex.y;
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);

	//rotate the planet-local-frame vector into modelview space
	vec4 modelPos = vec4(cosTheta * rad, sinTheta * rad, 0., 1.);

	//rotate the planet-local-frame vector into modelview space
	vec4 eyePos = mvMat * modelPos;
	worldPosv = quatRotate(angle, modelPos.xyz) + pos;

	gl_Position = projMat * eyePos;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
#define NUM_STARS 1
varying vec2 texCoordv;
varying vec3 worldPosv;

uniform sampler2D colorTex;

//vector from planet to sun
uniform vec3 sunDir[NUM_STARS];
uniform vec3 pos;

uniform float planetRadius;

float sphereIntersect(vec3 startPos, vec3 dir, vec3 spherePos, float sphereRadius) {
	vec3 sphereToStart = startPos - spherePos;
	float a = dot(dir, dir);
	float b = 2. * dot(sphereToStart, dir);
	float c = dot(sphereToStart, sphereToStart) - sphereRadius * sphereRadius;
	float discr = b*b - 4.*a*c;
	if (discr < 0.) return 1.;	//full trace, no intersection
	float sqrtDiscr = sqrt(discr);
	float invDenom = .5 / a;
	float t1 = (-b - sqrtDiscr) * invDenom;
	float t2 = (-b + sqrtDiscr) * invDenom;
	if (t1 >= 0. && t1 <= 1.) return t1;
	if (t2 >= 0. && t2 <= 1.) return t2;
	return 1.;
}

void main() {
	float luminance = 1.;
	//notice that I have to scale down the meters here for shader accuracy to work
	luminance = min(luminance, step(1., sphereIntersect(worldPosv * 1e-8, sunDir[0] * 1e-8, pos * 1e-8, planetRadius * 1e-8)));

	gl_FragColor = texture2D(colorTex, texCoordv);
	gl_FragColor.rgb *= sqrt(luminance);
}
*/})
		});

		var planet = solarSystem.planets[solarSystem.indexes.Jupiter];

		var texSrcInfo = [
			{field:'ringColorTex', url:'textures/jupiter-rings-color.png', format:gl.RGBA},
		];
		var numLoaded = 0;
		var onLoadRingTex = function(url) {
			++numLoaded;
			if (numLoaded < texSrcInfo.length) return;
			if (numLoaded > texSrcInfo.length) throw "already created the rings!";

			//done! create the ring object
			var vertexes = [];
			for (var i = 0; i < ringResolution; ++i) {
				var f = i / (ringResolution - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}

			var texs = [];
			for (var i = 0; i < texSrcInfo.length; ++i) {
				texs[i] = planet[texSrcInfo[i].field];
			}

			//and a ring object
			planet.ringObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				shader : jupiterRingShader,
				attrs : {
					vertex : new glutil.ArrayBuffer({
						dim : 2,
						data : vertexes
					})
				},
				uniforms : {
					colorTex : 0,
					ringMinRadius : planet.ringRadiusRange[0],
					ringMaxRadius : planet.ringRadiusRange[1],
					planetRadius : planet.radius,
					sunDir : [0,0,0, 0,0,0, 0,0,0, 0,0,0],
					pos : [0,0,0],
					angle : [0,0,0,1]
				},
				blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
				texs : texs,
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null
			});
		};
		$.each(texSrcInfo, function(i,info) {
			planet[info.field] = new glutil.Texture2D({
				url : info.url,
				format : info.format,
				internalFormat : info.format,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				wrap : {
					s :  gl.CLAMP_TO_EDGE,
					t : gl.REPEAT
				},
				generateMipmap : true,
				onload : onLoadRingTex
			});
		});
	})();


	//Saturn's rings
	(function(){
		var saturnRingShader = new ModifiedDepthShaderProgram({
			vertexCode : 
quatRotateCode
+ mlstr(function(){/*
#define M_PI 3.1415926535897931
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float ringMinRadius;
uniform float ringMaxRadius;

//these are actually baked into the mvMat, but I need them separate for converting coordinates to world coordinates (rather than eye coordinates)
//pos is in the scene, which is offset by the oribitting planet
//techinically the offset doesn't matter since i'm using it to offset vertices from the planet, then sphere-intersect-testing rays from those vertices against the planet
uniform vec3 pos;
uniform vec4 angle;

varying vec2 texCoordv;
varying vec3 worldPosv;	//varying in world coordinates

void main() {
	texCoordv = vertex;

	//vertex is the uv of the texcoord wrapped around the annulus
	//u coord is radial, v coord is angular
	float rad = mix(ringMinRadius, ringMaxRadius, vertex.x);
	float theta = 2. * M_PI * vertex.y;
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);

	vec4 modelPos = vec4(cosTheta * rad, sinTheta * rad, 0., 1.);

	//rotate the planet-local-frame vector into modelview space
	vec4 eyePos = mvMat * modelPos;
	worldPosv = quatRotate(angle, modelPos.xyz) + pos;

	gl_Position = projMat * eyePos;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
#define NUM_STARS 1
//described here:
//http://www.mmedia.is/~bjj/data/s_rings/index.html

varying vec2 texCoordv;
varying vec3 worldPosv;

uniform sampler2D colorTex;
uniform sampler2D backScatteredTex;	//from the sun looking at the rings
uniform sampler2D forwardScatteredTex;	//looking at the rings towards the sun
uniform sampler2D transparencyTex;
uniform sampler2D unlitSideTex;

uniform float lookingAtLitSide;
uniform float backToFrontLitBlend;

//vector from planet to sun
uniform vec3 sunDir[NUM_STARS];
uniform vec3 pos;

uniform float planetRadius;

float sphereIntersect(vec3 startPos, vec3 dir, vec3 spherePos, float sphereRadius) {
	vec3 sphereToStart = startPos - spherePos;
	float a = dot(dir, dir);
	float b = 2. * dot(sphereToStart, dir);
	float c = dot(sphereToStart, sphereToStart) - sphereRadius * sphereRadius;
	float discr = b*b - 4.*a*c;
	if (discr < 0.) return 1.;	//full trace, no intersection
	float sqrtDiscr = sqrt(discr);
	float invDenom = .5 / a;
	float t1 = (-b - sqrtDiscr) * invDenom;
	float t2 = (-b + sqrtDiscr) * invDenom;
	if (t1 >= 0. && t1 <= 1.) return t1;
	if (t2 >= 0. && t2 <= 1.) return t2;
	return 1.;
}

void main() {
	vec3 color = texture2D(colorTex, texCoordv).rgb;
	gl_FragColor.rgb = color;

	float backScattered = texture2D(backScatteredTex, texCoordv).r;
	float forwardScattered = texture2D(forwardScatteredTex, texCoordv).r;
	float litSide = mix(backScattered, forwardScattered, backToFrontLitBlend);

	float unlitSide = texture2D(unlitSideTex, texCoordv).r;

	//shadow from the sun
	//raytrace from worldPosv along sunDir (non-normalized)
	//if it intersects with a sphere at pos with radius planetRadius then we're in shadow
	//if you want to do a soft shadow based on the sun's radius and distance (length of sunDir) then you can

	float luminance = 1.;
	//notice that I have to scale down the meters here for shader accuracy to work
	luminance = min(luminance, step(1., sphereIntersect(worldPosv * 1e-8, sunDir[0] * 1e-8, pos * 1e-8, planetRadius * 1e-8)));
	luminance *= mix(unlitSide, litSide, lookingAtLitSide);

	gl_FragColor.rgb *= sqrt(luminance);

	float transparency = texture2D(transparencyTex, texCoordv).r;
	gl_FragColor.a = transparency;
}
*/})
		});

		var planet = solarSystem.planets[solarSystem.indexes.Saturn];

		var texSrcInfo = [
			{field:'ringColorTex', url:'textures/saturn-rings-color.png', format:gl.RGBA},
			{field:'ringBackScatteredTex', url:'textures/saturn-rings-back-scattered.png', format:gl.LUMINANCE},
			{field:'ringForwardScatteredTex', url:'textures/saturn-rings-forward-scattered.png', format:gl.LUMINANCE},
			{field:'ringTransparencyTex', url:'textures/saturn-rings-transparency.png', format:gl.LUMINANCE},
			{field:'ringUnlitSideTex', url:'textures/saturn-rings-unlit-side.png', format:gl.LUMINANCE}
		];
		var numLoaded = 0;
		var onLoadRingTex = function(url) {
			++numLoaded;
			if (numLoaded < texSrcInfo.length) return;
			if (numLoaded > texSrcInfo.length) throw "already created the rings!";

			//done! create the ring object
			var vertexes = [];
			for (var i = 0; i < ringResolution; ++i) {
				var f = i / (ringResolution - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}

			var texs = [];
			for (var i = 0; i < texSrcInfo.length; ++i) {
				texs[i] = planet[texSrcInfo[i].field];
			}

			//and a ring object
			planet.ringObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				shader : saturnRingShader,
				attrs : {
					vertex : new glutil.ArrayBuffer({
						dim : 2,
						data : vertexes
					})
				},
				uniforms : {
					colorTex : 0,
					backScatteredTex : 1,
					forwardScatteredTex : 2,
					transparencyTex : 3,
					unlitSideTex : 4,
					ringMinRadius : planet.ringRadiusRange[0],
					ringMaxRadius : planet.ringRadiusRange[1],
					planetRadius : planet.radius,
					lookingAtLitSide : 1,
					backToFrontLitBlend : 0,
					sunDir : [0,0,0, 0,0,0, 0,0,0, 0,0,0],
					pos : [0,0,0],
					angle : [0,0,0,1]
				},
				blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
				texs : texs,
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null
			});
		};
		$.each(texSrcInfo, function(i,info) {
			planet[info.field] = new glutil.Texture2D({
				url : info.url,
				format : info.format,
				internalFormat : info.format,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				wrap : {
					s :  gl.CLAMP_TO_EDGE,
					t : gl.REPEAT
				},
				generateMipmap : true,
				onload : onLoadRingTex
			});
		});
	})();
}

function calcKeplerianOrbitalElements(body, useVectorState) {
	var planetIndex = body.index;

	// based on body position and velocity, find plane of orbit
	var parentBody = body.parent;
	if (parentBody === undefined) {
		//console.log(body.name+' has no orbit parent');
		body.orbitAxis = [0,0,1];
		body.orbitBasis = [[1,0,0],[0,1,0],[0,0,1]];
		return;
	}

	//calculated variables:
	var semiMajorAxis = undefined;
	var eccentricity = undefined;
	var eccentricAnomaly = undefined;
	var orbitType = undefined;
	var meanAnomaly = undefined;
	var A, B;

	//used by planets to offset reconstructed orbit coordinates to exact position of body
	var checkPosToPosX = 0;
	var checkPosToPosY = 0;
	var checkPosToPosZ = 0;

	//if we're not using vector state then we're using keplerian orbital elements
	if (!useVectorState) {
		/*
		we have:
			eccentricity
			semiMajorAxis
			pericenterDistance
		we compute:
		*/

		eccentricity = assertExists(body.sourceData, 'eccentricity');

		var parabolicEccentricityEpsilon = 1e-7;
		if (Math.abs(eccentricity - 1) < parabolicEccentricityEpsilon) {
			orbitType = 'parabolic';
		} else if (eccentricity > 1) {
			orbitType = 'hyperbolic';
		} else {
			orbitType = 'elliptic';
		}

		var pericenterDistance;
		if (body.isComet) {
			pericenterDistance = assert(body.sourceData.perihelionDistance);
			if (orbitType !== 'parabolic') {
				semiMajorAxis = pericenterDistance / (1 - eccentricity);
			}	//otherwise, if it is parabolic, we don't get the semi-major axis ...
		} else if (body.isAsteroid) {
			semiMajorAxis = assert(body.sourceData.semiMajorAxis);
			pericenterDistance = semiMajorAxis * (1 - eccentricity);
		} else if (body.isExoplanet) {
			semiMajorAxis = body.sourceData.semiMajorAxis || 0;
		}

		var gravitationalParameter = gravitationalConstant * parentBody.mass;	//assuming the comet mass is negligible, since the comet mass is not provided
		var semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;
		
		//orbital period is only defined for circular and elliptical orbits (not parabolic or hyperbolic)
		var orbitalPeriod = undefined;
		if (orbitType === 'elliptic') {
			orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed / gravitationalParameter) / (60*60*24);	//julian day
		}

		var longitudeOfAscendingNode = assertExists(body.sourceData, 'longitudeOfAscendingNode');
		var cosAscending = Math.cos(longitudeOfAscendingNode);
		var sinAscending = Math.sin(longitudeOfAscendingNode);

		var argumentOfPeriapsis = assertExists(body.sourceData, 'argumentOfPeriapsis');
		var cosPericenter = Math.cos(argumentOfPeriapsis);
		var sinPericenter = Math.sin(argumentOfPeriapsis);

		var inclination = assertExists(body.sourceData, 'inclination');
		var cosInclination = Math.cos(inclination);
		var sinInclination = Math.sin(inclination);

		var oneMinusEccentricitySquared = 1 - eccentricity * eccentricity;
		//magnitude of A is a 
		A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
			 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
			 semiMajorAxis * sinPericenter * sinInclination];
		//magnitude of B is a * sqrt(|1 - e^2|)
		B = [semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * -(cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * cosPericenter * sinInclination];
		//inner product: A dot B = 0

		var timeOfPeriapsisCrossing;
		if (body.isComet) {
			timeOfPeriapsisCrossing = body.sourceData.timeOfPerihelionPassage;	//julian day
		}

		if (orbitType === 'parabolic') {
			eccentricAnomaly = Math.tan(argumentOfPeriapsis / 2);
			meanAnomaly = eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3; 
		} else if (orbitType === 'hyperbolic') {
			assert(timeOfPeriapsisCrossing !== undefined);	//only comets are hyperbolic, and all comets have timeOfPeriapsisCrossing defined
			meanAnomaly = Math.sqrt(-gravitationalParameter / semiMajorAxisCubed) * timeOfPeriapsisCrossing * 60*60*24;	//in seconds
		} else if (orbitType === 'elliptic') {
			//in theory I can say 
			//eccentricAnomaly = Math.acos((eccentricity + Math.cos(argumentOfPeriapsis)) / (1 + eccentricity * Math.cos(argumentOfPeriapsis)));
			// ... but Math.acos has a limited range ...

			if (body.isComet) {
				timeOfPeriapsisCrossing = body.sourceData.timeOfPerihelionPassage;	//julian day
				var timeSinceLastPeriapsisCrossing = initJulianDate - timeOfPeriapsisCrossing;
				meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * Math.PI / orbitalPeriod;
			} else if (body.isAsteroid) {
				meanAnomaly = body.sourceData.meanAnomalyAtEpoch + 2 * Math.PI / orbitalPeriod * (initJulianDate - body.sourceData.epoch);
			} else if (body.isExoplanet) {
				meanAnomaly = body.sourceData.meanAnomaly !== undefined ? body.sourceData.meanAnomaly : 0;
			} else {
			//	throw 'here';
			}
		} else {
			throw 'here';
		}
		
		//solve Newton Rhapson
		//for elliptical orbits:
		//	f(E) = M - E + e sin E = 0
		//	f'(E) = -1 + e cos E
		//for parabolic oribts:
		//	f(E) = M - E - E^3 / 3
		//	f'(E) = -1 - E^2
		//for hyperbolic orbits:
		//	f(E) = M - e sinh(E) - E
		//	f'(E) = -1 - e cosh(E)
		eccentricAnomaly = meanAnomaly;
		for (var i = 0; i < 10; ++i) {
			var func, deriv;
			if (orbitType === 'parabolic') {	//parabolic
				func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3;
				deriv = -1 - eccentricAnomaly * eccentricAnomaly;
			} else if (orbitType === 'elliptic') { 	//elliptical
				func = meanAnomaly - eccentricAnomaly + eccentricity * Math.sin(eccentricAnomaly);
				deriv = -1 + eccentricity * Math.cos(eccentricAnomaly);	//has zeroes ...
			} else if (orbitType === 'hyperbolic') {	//hyperbolic
				func = meanAnomaly + eccentricAnomaly - eccentricity  * Math.sinh(eccentricAnomaly);
				deriv = 1 - eccentricity * Math.cosh(eccentricAnomaly);
			} else {
				throw 'here';
			}

			var delta = func / deriv;
			if (Math.abs(delta) < 1e-15) break;
			eccentricAnomaly -= delta;
		}

		var sinEccentricAnomaly = Math.sin(eccentricAnomaly);
		var cosEccentricAnomaly = Math.cos(eccentricAnomaly);
		
		//parabolas and hyperbolas don't define orbitalPeriod
		// so no need to recalculate it
		if (orbitalPeriod !== undefined && meanAnomaly !== undefined) {
			timeOfPeriapsisCrossing = meanAnomaly * orbitalPeriod / (2 * Math.PI); //if it is a comet then we're just reversing the calculation above ...
		}

/*
r = radial distance
a = semi major axis
E = eccentric anomaly
nu = true anomaly

r^2 = a^2 (cos(E) - e)^2 + a^2 |1 - e^2| sin(E)^2
r^2 = a^2 (cos(E)^2 - 2 e cos(E) + e^2 + |1 - e^2| sin(E)^2)

=== for eccentric orbits ... e < 1 <=> e^2 < 1 <=> 1 - e^2 > 0
r^2 = a^2 (cos(E)^2 - 2 e cos(E) + e^2 + (1 - e^2) sin(E)^2)
r^2 = a^2 (cos(E)^2 - 2 e cos(E) + e^2 + sin(E)^2 - e^2 sin(E)^2)
r^2 = a^2 (1 - 2 e cos(E) + e^2 - e^2 sin(E)^2)
r^2 = a^2 (1 - 2 e cos(E) + e^2 - e^2 (1 - cos(E)^2))
r^2 = a^2 (1 - 2 e cos(E) + e^2 cos(E)^2)
r^2 = a^2 (1 - e cos(E))^2
r = a (1 - e cos(E))		true according to http://en.wikipedia.org/wiki/Eccentric_anomaly

r = a (1 - e^2) / (1 + e cos(nu)) is stated by http://www.bogan.ca/orbits/kepler/orbteqtn.html
therefore should be true:
(1 - e cos(E)) = (1 - e^2) / (1 + e cos(nu))
1 + e cos(nu) = (1 - e^2) / (1 - e cos(E))
e cos(nu) = (1 - e^2 - 1 + e cos(E)) / (1 - e cos(E))
cos(nu) = (cos(E) - e) / (1 - e cos(E))
...which checks out with Wikipedia: http://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly

r = a (1 - e^2) / (1 + e cos(nu))
r/a = (1 - e^2) / (1 + e cos(nu))
a/r = (1 + e cos(nu)) / (1 - e^2)
2 a/r - 1 = (1 + e cos(nu)) / (1 - e^2) - (1 - e^2)/(1 - e^2)
2 a/r - r/r = e (e + cos(nu)) / (1 - e^2)
(2 a - r)/(r a) = [e (e + cos(nu))] / [a (1 - e^2)]
v^2 = mu (2/r - 1/a) = mu (e^2 + e cos(nu)) / (a (1 - e^2))
the first of these is true: v^2 = mu (2/r - 1/a) according to both bogan.ca and wikipedia 

the second: v^2 = mu e (e + cos(nu)) / (a (1 - e^2)) should match v^2 = mu / (a (1 - e^2)) (1 - 2 cos(nu) + e^2)
is true only for e^2 + e cos(nu) = 1 - 2 cos(nu) + e^2
	... e cos(nu) = 1 - 2 cos(nu)
	... (2 + e) cos(nu) = 1
	... cos(nu) = 1/(e + 2)	...which doesn't look true ... so maybe that bogan.ca second velocity equation v^2 = (mu/p)(1 - 2 cos(nu) + e^2) is wrong?

=== for hyperbolic orbits ... e > 1 <=> e^2 > 1 <=> 1 - e^2 < 0
r^2 = a^2 (cos(E) 
	... we should get to 

*/
		var dt_dE;
		if (orbitType == 'parabolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 + eccentricAnomaly * eccentricAnomaly);
		} else if (orbitType == 'elliptic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 - eccentricity * Math.cos(eccentricAnomaly));
		} else if (orbitType == 'hyperbolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (eccentricity * Math.cosh(eccentricAnomaly) - 1);
		}
		var dE_dt = 1/dt_dE;
		//finally using http://en.wikipedia.org/wiki/Kepler_orbit like I should've in the first place ...
		var coeffA, coeffB;
		var coeffDerivA, coeffDerivB;
		if (orbitType == 'parabolic') {
			//...?
		} else if (orbitType == 'elliptic') { 
			coeffA = Math.cos(eccentricAnomaly) - eccentricity;
			coeffB = Math.sin(eccentricAnomaly);
			coeffDerivA = -Math.sin(eccentricAnomaly) * dE_dt;
			coeffDerivB = Math.cos(eccentricAnomaly) * dE_dt;
		} else if (orbitType == 'hyperbolic') {
			coeffA = eccentricity - Math.cosh(eccentricAnomaly);
			coeffB = Math.sinh(eccentricAnomaly);
			coeffDerivA = -Math.sinh(eccentricAnomaly) * dE_dt;
			coeffDerivB = Math.cosh(eccentricAnomaly) * dE_dt;
		}
		
		var posX = A[0] * coeffA + B[0] * coeffB;
		var posY = A[1] * coeffA + B[1] * coeffB;
		var posZ = A[2] * coeffA + B[2] * coeffB;
		//v^2 = (a^2 sin^2(E) + a^2 |1 - e^2| cos^2(E)) * mu/a^3 / (1 - e cos(E))
		//v^2 = (sin^2(E) + |1 - e^2| cos^2(E)) * mu/(a (1 - e cos(E)))
		
		//v^2 should be = mu/(a(1-e^2)) * (1 + e^2 - 2 cos(nu))	<- for nu = true anomaly
		//v^2 should also be = mu (2/r - 1/a) = mu (2a - r) / (r a)
		//... then (2 a - r) / (r a) should = (1 - 2 cos(nu) + e^2) / (a (1 - e^2))
		//...	2 a/r - 1 should = (1 - 2 cos(nu) + e^2) / ((1 + e) (1 - e))
		//...	2 a/r should = (1-e^2)/(1-e^2) + (1 - 2 cos(nu) + e^2) / (1-e^2)
		//...	2 a/r should = (2 - 2 cos(nu)) / (1-e^2)
		//...	a/r should = (1 - cos(nu)) / (1 - e^2)
		//...	r/a should = (1 - e^2) / (1 - cos(nu))
		//...	r should = a (1 - e^2) / (1 - cos(nu))
		var velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
		var velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
		var velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;
		body.pos[0] = posX + parentBody.pos[0];
		body.pos[1] = posY + parentBody.pos[1];
		body.pos[2] = posZ + parentBody.pos[2];
		body.vel[0] = velX + parentBody.vel[0];
		body.vel[1] = velY + parentBody.vel[1];
		body.vel[2] = velZ + parentBody.vel[2];
		vec3.copy(solarSystem.initPlanets[planetIndex].pos, body.pos);
		vec3.copy(solarSystem.initPlanets[planetIndex].vel, body.vel);

		body.keplerianOrbitalElements = {
			semiMajorAxis : semiMajorAxis,
			eccentricity : eccentricity,
			eccentricAnomaly : eccentricAnomaly,
			longitudeOfAscendingNode : longitudeOfAscendingNode,
			argumentOfPeriapsis : argumentOfPeriapsis,
			inclination : inclination,
			timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
			meanAnomaly : meanAnomaly,
			orbitType : orbitType,
			orbitalPeriod : orbitalPeriod,	//only exists for elliptical orbits
			A : A,
			B : B
		};

	//using vector state (pos & vel) as provided by Horizons
	} else {

		orbitType = 'elliptic';	//better be ...

		//consider position relative to orbiting parent
		// should I be doing the same thing with the velocity?  probably...
		var posX = body.pos[0] - parentBody.pos[0];
		var posY = body.pos[1] - parentBody.pos[1];
		var posZ = body.pos[2] - parentBody.pos[2];

		//convert from m/day to m/s to coincide with the units of our gravitational constant
		var velX = (body.vel[0] - parentBody.vel[0]) / (60 * 60 * 24);
		var velY = (body.vel[1] - parentBody.vel[1]) / (60 * 60 * 24);
		var velZ = (body.vel[2] - parentBody.vel[2]) / (60 * 60 * 24);

		var posDotVel = posX * velX + posY * velY + posZ * velZ;	//m^2/s

		var angularMomentumX = posY * velZ - posZ * velY; //m^2/s
		var angularMomentumY = posZ * velX - posX * velZ;
		var angularMomentumZ = posX * velY - posY * velX;
		var angularMomentumMagSq = angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY + angularMomentumZ * angularMomentumZ;		//m^4/s^2
		var angularMomentumMag = Math.sqrt(angularMomentumMagSq);
		if (angularMomentumMag < 1e-9) {
			body.orbitAxis = [0,0,1];
		} else {
			var axisX = angularMomentumX / angularMomentumMag;
			var axisY = angularMomentumY / angularMomentumMag;
			var axisZ = angularMomentumZ / angularMomentumMag;
			body.orbitAxis = [axisX, axisY, axisZ];
		}

		var basisX = [0,0,0];
		var basisY = [0,0,0];
		var basisZ = body.orbitAxis;
		//TODO use the inclination and longitudeOfAscendingNode
		calcBasis(basisX, basisY, basisZ);
		//a[j][i] = a_ij, so our indexing is backwards, but our storage is column-major
		body.orbitBasis = [basisX, basisY, basisZ];

		//now decompose the relative position in the coordinates of the orbit basis
		//i've eliminated all but one of the rotation degrees of freedom ...

		//http://www.mathworks.com/matlabcentral/fileexchange/31333-orbital-elements-from-positionvelocity-vectors/content/vec2orbElem.m
		//http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

		var velSq = velX * velX + velY * velY + velZ * velZ;		//(m/s)^2
		var distanceToParent = Math.sqrt(posX * posX + posY * posY + posZ * posZ);		//m
		var gravitationalParameter = gravitationalConstant * ((body.mass || 0) + parentBody.mass);	//m^3 / (kg s^2) * kg = m^3 / s^2
		var specificOrbitalEnergy  = .5 * velSq - gravitationalParameter / distanceToParent;		//m^2 / s^2 - m^3 / s^2 / m = m^2/s^2, supposed to be negative for elliptical orbits
		semiMajorAxis = -.5 * gravitationalParameter / specificOrbitalEnergy;		//m^3/s^2 / (m^2/s^2) = m
		var semiLatusRectum = angularMomentumMagSq / gravitationalParameter;			//m^4/s^2 / (m^3/s^2) = m
		eccentricity = Math.sqrt(1 - semiLatusRectum / semiMajorAxis);				//unitless (assuming elliptical orbit)

		var cosEccentricAnomaly = (1 - distanceToParent / semiMajorAxis) / eccentricity;						//unitless
		var sinEccentricAnomaly = posDotVel / (eccentricity * Math.sqrt(gravitationalParameter * semiMajorAxis));	//m^2/s / sqrt(m^3/s^2 * m) = m^2/s / sqrt(m^4/s^2) = m^2/s / (m^2/s) = unitless
		eccentricAnomaly = Math.atan2(sinEccentricAnomaly, cosEccentricAnomaly);	//radians (unitless)

		var sinInclination = Math.sqrt(angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY) / angularMomentumMag;	//unitless
		var cosInclination = angularMomentumZ / angularMomentumMag;	//unitless
		var inclination = Math.atan2(sinInclination, cosInclination);

		var sinPericenter = ((velX * angularMomentumY - velY * angularMomentumX) / gravitationalParameter - posZ / distanceToParent) / (eccentricity * sinInclination);
		var cosPericenter = (angularMomentumMag * velZ / gravitationalParameter - (angularMomentumX * posY - angularMomentumY * posX) / (angularMomentumMag * distanceToParent)) / (eccentricity * sinInclination);
		var argumentOfPeriapsis = Math.atan(sinPericenter, cosPericenter);

		var cosAscending = -angularMomentumY / (angularMomentumMag * sinInclination);
		var sinAscending = angularMomentumX / (angularMomentumMag * sinInclination);
		var longitudeOfAscendingNode = Math.atan2(sinAscending, cosAscending);

		var semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;	//m^3
		var orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed  / gravitationalParameter) / (60*60*24);	//julian day
//override for Earth
if (body.orbitalPeriod !== undefined) orbitalPeriod = body.orbitalPeriod;
		
		var timeOfPeriapsisCrossing = -(eccentricAnomaly - eccentricity * sinEccentricAnomaly) / Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (60*60*24);	//julian day

		A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
				 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
				 semiMajorAxis * sinPericenter * sinInclination];
		B = [-semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
				 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
				 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * cosPericenter * sinInclination];

		/*
		to convert back:
		pos[i] = A[i] * (cosEccentricAnomaly - eccentricity) + B[i] * sinEccentricAnomaly
		rDot[i] = (-A[i] * sinEccentricAnomaly + B[i] * cosEccentricAnomaly) * Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (1 - eccentricity * cosEccentricAnomaly)
		*/
		var checkPosX = A[0] * (cosEccentricAnomaly - eccentricity) + B[0] * sinEccentricAnomaly;
		var checkPosY = A[1] * (cosEccentricAnomaly - eccentricity) + B[1] * sinEccentricAnomaly;
		var checkPosZ = A[2] * (cosEccentricAnomaly - eccentricity) + B[2] * sinEccentricAnomaly;

		checkPosToPosX = checkPosX - posX;
		checkPosToPosY = checkPosY - posY;
		checkPosToPosZ = checkPosZ - posZ;
		var checkPosToPosDist = Math.sqrt(checkPosToPosX * checkPosToPosX + checkPosToPosY * checkPosToPosY + checkPosToPosZ * checkPosToPosZ);
		var checkPosError = checkPosToPosDist / distanceToParent;
		if (checkPosError === checkPosError) {
			if (checkPosError > 1e-5) {	//only report significant error
				console.log(body.name+' error of reconstructed position '+ checkPosError);
			}
		} else {	//NaN? debug!

		/*
			for (k in body) {
				var v = body[k];
				if (k != 'name' && typeof(v) != 'function') {
					console.log(body.name, k, v);
				}
			}
		*/
			console.log(body.name+' has no orbit info.  mass: '+body.mass+' radius: '+body.radius);
		}
				
		meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * Math.PI / orbitalPeriod;

		body.keplerianOrbitalElements = {
			relVelSq : velSq,
			gravitationalParameter : gravitationalParameter,
			specificOrbitalEnergy : specificOrbitalEnergy,
			distanceToParent : distanceToParent,
			semiMajorAxis : semiMajorAxis,
			semiLatusRectum : semiLatusRectum,
			inclination : inclination,
			argumentOfPeriapsis : argumentOfPeriapsis,
			longitudeOfAscendingNode : longitudeOfAscendingNode,
			timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
			meanAnomaly : meanAnomaly,
			orbitalPeriod : orbitalPeriod,
			//the following are used for the orbit path shader:
			orbitType : orbitType,
			eccentricity : eccentricity,
			eccentricAnomaly : eccentricAnomaly,
			A : A,
			B : B,
			//this is used for drawing but not in the shader
			fractionOffset : 0
		};

		//not NaN, we successfully reconstructed the position
		if (checkPosError !== checkPosError) {
			return;
		}
	}

	//for elliptic orbits,
	// we can accumulate & store the largest semi major axis of all children
	//but for parabolic/hyperbolic ...
	// there is no bound ... so maybe those objects should be stored/rendered separately?
	if (parentBody !== undefined && orbitType == 'elliptic') {
		parentBody.maxSemiMajorAxisOfEllipticalSatellites = Math.max(semiMajorAxis, parentBody.maxSemiMajorAxisOfEllipticalSatellites);
	}

	body.renderOrbit = true;

	/*
	integration can be simulated along Keplerian orbits using the A and B vectors ...
		time advanced = simulated date - data snapshot date
		orbit fraction advanced = time advanced / orbital period
		body position and velocity can then be computed using A and B's evaluation and derivatives
		then the alpha of the orbit needs to be updated ...
	*/
}

function recomputePlanetAlongOrbit(planet) {
	var timeAdvanced = julianDate - initJulianDate;
	if (!planet.parent) return;
	var ke = planet.keplerianOrbitalElements;
	var orbitType = ke.orbitType;

	var meanMotion = undefined;
	if (ke.orbitalPeriod !== undefined) {	//elliptical
		meanMotion = 2 * Math.PI / ke.orbitalPeriod;
	} else if (ke.timeOfPeriapsisCrossing !== undefined) {	//hyperbolic ... shouldn't be using mean motion?
		meanMotion = ke.meanAnomaly / (julianDate - ke.timeOfPeriapsisCrossing);
	} else {
		throw 'here';
	}

	//TODO don't use meanMotion for hyperbolic orbits
	var fractionOffset = timeAdvanced * meanMotion / (2 * Math.PI); 
	var theta = timeAdvanced * meanMotion;
	var pathEccentricAnomaly = ke.eccentricAnomaly + theta;
	var A = ke.A;
	var B = ke.B;

	//matches above
	var dt_dE;
	var semiMajorAxisCubed = ke.semiMajorAxis * ke.semiMajorAxis * ke.semiMajorAxis;
	if (orbitType == 'parabolic') {
		dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (1 + pathEccentricAnomaly * pathEccentricAnomaly);
	} else if (orbitType == 'elliptic') {
		dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (1 - ke.eccentricity * Math.cos(pathEccentricAnomaly));
	} else if (orbitType == 'hyperbolic') {
		dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (ke.eccentricity * Math.cosh(pathEccentricAnomaly) - 1);
	}
	var dE_dt = 1/dt_dE;
	var coeffA, coeffB;
	var coeffDerivA, coeffDerivB;
	if (orbitType == 'parabolic') {
		//...?
	} else if (orbitType == 'elliptic') { 
		coeffA = Math.cos(pathEccentricAnomaly) - ke.eccentricity;
		coeffB = Math.sin(pathEccentricAnomaly);
		coeffDerivA = -Math.sin(pathEccentricAnomaly) * dE_dt;
		coeffDerivB = Math.cos(pathEccentricAnomaly) * dE_dt;
	} else if (orbitType == 'hyperbolic') {
		coeffA = ke.eccentricity - Math.cosh(pathEccentricAnomaly);
		coeffB = Math.sinh(pathEccentricAnomaly);
		coeffDerivA = -Math.sinh(pathEccentricAnomaly) * dE_dt;
		coeffDerivB = Math.cosh(pathEccentricAnomaly) * dE_dt;
	}
	var posX = A[0] * coeffA + B[0] * coeffB;
	var posY = A[1] * coeffA + B[1] * coeffB;
	var posZ = A[2] * coeffA + B[2] * coeffB;
	var velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
	var velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
	var velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;
	
	planet.pos[0] = posX + planet.parent.pos[0];
	planet.pos[1] = posY + planet.parent.pos[1];
	planet.pos[2] = posZ + planet.parent.pos[2];
	planet.vel[0] = velX + planet.parent.vel[0];
	planet.vel[1] = velY + planet.parent.vel[1];
	planet.vel[2] = velZ + planet.parent.vel[2];
	
	planet.keplerianOrbitalElements.fractionOffset = fractionOffset;
}

function recomputePlanetsAlongOrbit() {
	var starSystem = orbitStarSystem;
	for (var i = 0; i < starSystem.planets.length; ++i) {	//or do it for all systems?
		var planet = starSystem.planets[i];
		recomputePlanetAlongOrbit(planet);
	}
}

function initOrbitPaths() {
	var gravWellShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec4 vertex;
varying float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex.xyz, 1.);
	alpha = vertex.w;
	gl_Position = projMat * vtx4;
	gl_PointSize = 4.;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
varying float alpha;
void main() {
	gl_FragColor = vec4(color.xyz, color.w * alpha);
}
*/})
	});

	//for rendering the orbital path
	var orbitPathShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec4 vertex;
varying float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 A, B;
uniform float eccentricity;
uniform float eccentricAnomaly;
uniform int orbitType;

float cosh(float x) { return .5 * (exp(x) + exp(-x)); }
float sinh(float x) { return .5 * (exp(x) - exp(-x)); }

#define TWO_PI 6.283185307179586231995926937088
void main() {
	float frac = vertex.x;
	float theta = frac * TWO_PI;
	float pathEccentricAnomaly = eccentricAnomaly + theta;

	float coeffA, coeffB;
	if (orbitType == 0) {	//elliptic
		coeffA = cos(pathEccentricAnomaly) - eccentricity;
		coeffB = sin(pathEccentricAnomaly);
	} else {	//if (orbitType == 1) {	//hyperbolic
		coeffA = eccentricity - cosh(pathEccentricAnomaly);
		coeffB = sinh(pathEccentricAnomaly);
	}
	vec3 pos = A * coeffA + B * coeffB;	
	
	alpha = frac; 

	gl_Position = projMat * (mvMat * vec4(pos, 1.));
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
uniform float fractionOffset;
varying float alpha;
void main() {
	float a = mod(alpha - fractionOffset, 1.);
	a = fract(a + 1.);
	a = a * .75 + .25;
	gl_FragColor = vec4(color.xyz, color.w * a);
}
*/})
	});

	//iterate around the eccentric anomaly to reconstruct the path
	var vertexes = [];
	for (var i = 0; i < orbitPathResolution; ++i) {
		vertexes.push(i / (orbitPathResolution - 1));
	}
	
	orbitPathSceneObj = new glutil.SceneObject({
		mode : gl.LINE_STRIP,
		shader : orbitPathShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({dim : 1, data : vertexes})
		},
		blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
		pos : [0,0,0],
		angle : [0,0,0,1],
		parent : null
	});


	//shader for recomputing planet positions
	//associated with texture per-planet that stores position (and maybe velocity)

	//TODO update these when we integrate!
	var calcOrbitPathStartTime = Date.now();
	//init solar system planets based on Horizons data
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		calcKeplerianOrbitalElements(solarSystem.planets[i], true);
	}
	var calcOrbitPathEndTime = Date.now();

	var setPlanetTiltAngleToMoonOrbitPlane = function(planetName, moonName) {
		var planet = solarSystem.planets[solarSystem.indexes[planetName]];
		var moon = solarSystem.planets[solarSystem.indexes[moonName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.longitudeOfAscendingNode);
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.inclination);
		quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
	};

	//accepts degrees
	//TODO start at orbit axis plane rather than earth's (ie J2000) orbit axis plane
	var setPlanetTiltAngleToFixedValue = function(planetName, inclination, tiltDirection) {
		if (tiltDirection === undefined) tiltDirection = 0;
		var planet = solarSystem.planets[solarSystem.indexes[planetName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, Math.rad(tiltDirection));
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, Math.rad(inclination));
		quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
	};

	setPlanetTiltAngleToFixedValue('Mercury', 2.11/60);		//TODO tilt from mercury orbit plane.  until then it's off
	setPlanetTiltAngleToFixedValue('Venus', 177.3);			//TODO tilt from venus orbit plane.  until then, this measures 175 degrees.
	setPlanetTiltAngleToFixedValue('Earth', 23 + 1/60*(26 + 1/60*(21.4119)), 180);
	setPlanetTiltAngleToMoonOrbitPlane('Mars', 'Phobos');		//ours: 25.79, exact: 25.19
	setPlanetTiltAngleToMoonOrbitPlane('Jupiter', 'Metis');		//ours: 3.12, exact: 3.13
	setPlanetTiltAngleToMoonOrbitPlane('Saturn', 'Atlas');		//ours: 26.75, exact: 26.73
	setPlanetTiltAngleToMoonOrbitPlane('Uranus', 'Cordelia');	//ours: 97.71, exact: 97.77
	setPlanetTiltAngleToMoonOrbitPlane('Neptune', 'Galatea');	//ours: 28.365, exact: 28.32
	setPlanetTiltAngleToMoonOrbitPlane('Pluto', 'Charon');		//ours: 119, exact: 119

	//now set all moon tilt angles to the orbit axis 
	var Sun = solarSystem.planets[0];
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		var planet = solarSystem.planets[i];
		if (planet.parent !== undefined &&
			planet.parent.parent !== undefined &&
			planet.parent.parent == Sun)
		{
			quat.rotateZ(planet.tiltAngle, planet.tiltAngle, planet.keplerianOrbitalElements.longitudeOfAscendingNode);
			quat.rotateX(planet.tiltAngle, planet.tiltAngle, planet.keplerianOrbitalElements.inclination);
			quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
		}
	}

	//looks like low grav wells run into fp accuracy issues
	//how about extracting the depth and storing normalized values?
	var calcGravWellStartTime = Date.now();
	for (var planetIndex = 0; planetIndex < solarSystem.planets.length; ++planetIndex) {
		var planet = solarSystem.planets[planetIndex];

		/*
		gravity wells

		calculate spacetime embedding radius
		  for radial distance r, radius R, Schwarzschild radius Rs
		 inside the planet  (r <= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R (r/R)^2 ))
		 outside the planet (r >= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R)) + sqrt(4Rs(r - Rs)) - sqrt(4Rs(R - Rs))
		*/
		var R = planet.radius;
		var Rs = planet.schwarzschildRadius;
		var R_Rs = R / Rs;
		var Rs_R = Rs / R;
		var R_sqrt_R_Rs = R * Math.sqrt(R_Rs);

		//extract the normalized scalar to post-multiply after transformation (so small fields will not lose accuracy)
		planet.gravityWellScalar = Math.sqrt(R / Rs) - Math.sqrt(R / Rs - 1);

		//populate first vertex
		var gravWellVtxs = [0,0,0,1];

		var gravWellIndexes = [];

		var rimax = 60;		//was 200 r and 60 th, but I added a lot of planets.  need to occlude these based on distance/angle ...
		var thimax = 60;

		for (var ri = 1; ri < rimax; ++ri) {
			var r = R * Math.pow(100, ri / rimax * (gravityWellRadialMaxLog100 - gravityWellRadialMinLog100) + gravityWellRadialMinLog100);
			//max radial dist is R * Math.pow(100, gravityWellRadialMaxLog100)

			var z;
			if (r <= R) {
				var r_R = r / R;
				z = R_sqrt_R_Rs * (1 - Math.sqrt(1 - Rs_R * r_R * r_R));
			} else {
				z = R_sqrt_R_Rs * (1 - Math.sqrt(1 - Rs_R))
					+ Math.sqrt(4 * Rs * (r - Rs))
					- Math.sqrt(4 * Rs * (R - Rs));
			}
			z /= planet.gravityWellScalar;

			for (var thi = 0; thi < thimax; ++thi) {
				var th = 2 * Math.PI * thi / thimax;

				var x = r * Math.cos(th);
				var y = r * Math.sin(th);

				//it would be great to include influences of other planets' gravity wells ...
				// but that would require us performing our calculation at some point in 3D space -- for the sake of the radial coordinate
				// I could approximate that if I knew the approximate plane of orbit of most all the planets ...
				// from there, take samples based on radial distance within that plane between planets ...
				// I could remap planes on a per-planet basis depending on which you are orbiting ...
				// or I could always try for something 3D ... exxhagerated pinch lattice vectors ...
				// to do this it might be best to get the chebyshev interval calculations working, or somehow calculate the ellipses of rotation of each planet
				//TODO also: recalculate the gravity well mesh when the planets change, to watch it in realtime

				gravWellVtxs.push(x * planet.orbitBasis[0][0] + y * planet.orbitBasis[1][0] + z * planet.orbitBasis[2][0]);
				gravWellVtxs.push(x * planet.orbitBasis[0][1] + y * planet.orbitBasis[1][1] + z * planet.orbitBasis[2][1]);
				gravWellVtxs.push(x * planet.orbitBasis[0][2] + y * planet.orbitBasis[1][2] + z * planet.orbitBasis[2][2]);
				var tau = ri/(rimax-1);
				gravWellVtxs.push(1 - tau * tau * tau * tau * tau * tau * tau * tau);

				if (ri == 1) {
					gravWellIndexes.push(0);
					gravWellIndexes.push(1 + thi + thimax * (ri-1));
				}
				if (ri < rimax-1) {
					gravWellIndexes.push(1 + thi + thimax * (ri-1));
					gravWellIndexes.push(1 + thi + thimax * ri);	//plus one radial
				}
				gravWellIndexes.push(1 + thi + thimax * (ri-1));
				gravWellIndexes.push(1 + ((thi+1)%thimax) + thimax * (ri-1));	//plus one tangential
			}
		}
		planet.gravWellObj = new glutil.SceneObject({
			mode : gl.LINES,
			indexes : new glutil.ElementArrayBuffer({
				data : gravWellIndexes
			}),
			shader : gravWellShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({
					dim : 4,
					data : gravWellVtxs
				})
			},
			uniforms : {
				color : [1,1,1,.2]
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			parent : null
		});
		//scenegraph is a mess
		//static objects have mvMat pointed to glutil.scene.mvMat
		//but dynamic objects only give control over position and angle
		//!static implies creating the object's matrices, but !static also implies overwriting them every draw call...
		planet.localMat = mat4.create();
		planet.gravWellObj.mvMat = mat4.create();
		planet.gravWellObj.uniforms.mvMat = planet.gravWellObj.mvMat;
	}
	var calcGravWellEndTime = Date.now();

	//on my machine this was 769ms
	console.log('calc orbit path time ',calcOrbitPathEndTime-calcOrbitPathStartTime,'ms');

	//on my machine this was 386ms
	console.log('calc grav well time ',calcGravWellEndTime-calcGravWellStartTime,'ms');

	//only do this after the orbit shader is made
	initExoplanets();

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
		recomputePlanetsAlongOrbit();
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
				orbitTargetDistance *= scale;
			}
			
			refreshOrbitTargetDistanceText();
		},
		passiveMove : function() {
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
		recomputePlanetsAlongOrbit();

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
