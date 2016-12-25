/*
all stars of the milky way
TODO OOP this, and make one per galaxy (which we have observed stars within) 
*/

var showStars = true;

//TODO merge with starSystems[] ... keep the StarField for point rendering of all StarSystems (or make it a Galaxy object, honestly, that's where thignsn are going)
// and remove StarInField ... make that just StarSystem (even for zero-planet systems)
//
//only instanciate these for the named stars.  87 in all.

var stars = new function() {
	this.maxDistInLyr = 5000;
	this.renderScale = 1e+10;
	
	this.init = function() {
		var thiz = this;
		var xhr = new XMLHttpRequest();
		xhr.open('GET', 'hyg/stardata.f32', true);
		xhr.responseType = 'arraybuffer';
		/* if we want a progress bar ...
		xhr.onprogress = function(e) {
			if (e.total) {
				progress.attr('value', parseInt(e.loaded / e.total * 100));
			}
		};
		*/
		var numElem = 8;
		xhr.onload = function(e) {
			var arrayBuffer = this.response;
			var data = new DataView(arrayBuffer);

			//units are in parsecs
			//don't forget velocity is not being rescaled (i'm not using it at the moment)
			var floatBuffer = new Float32Array(data.byteLength / Float32Array.BYTES_PER_ELEMENT);
			var len = floatBuffer.length;
			for (var j = 0; j < len; ++j) {
				var x = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);
				if (j % numElem < 3) {
					if (Math.abs(x) > 1e+20) {
						console.log('star '+Math.floor(j/numElem)+' has coordinate that position exceeds fp resolution');
					}
					//convert parsecs to meters
					x *= metersPerUnits.pc / thiz.renderScale;
				}
				floatBuffer[j] = x;
			}

			//now that we have the float buffer ...
			thiz.buffer = new glutil.ArrayBuffer({data : floatBuffer, dim : numElem});
			thiz.sceneObj = new glutil.SceneObject({
				mode : gl.POINTS,
				shader : colorIndexShader,
				texs : [colorIndexTex, starTex],
				attrs : {
					vertex : new glutil.Attribute({buffer : thiz.buffer, size : 3, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 0}),	//xyz abs-mag
					velocity : new glutil.Attribute({buffer : thiz.buffer, size : 3, stride : numElem * Float32Array.BYTES_PER_ELEM, offset : 3 * Float32Array.BYTES_PER_ELEMENT}),	//velocity
					absoluteMagnitude : new glutil.Attribute({buffer : thiz.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEM, offset : 6 * Float32Array.BYTES_PER_ELEMENT}),
					colorIndex : new glutil.Attribute({buffer : thiz.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 7 * Float32Array.BYTES_PER_ELEMENT})	//color-index
				},
				uniforms : {
					pointSize : 1,
					color : [1,1,1,1],
					visibleMagnitudeBias : starsVisibleMagnitudeBias,
					colorIndexTex : 0,
					starTex : 1
				},
				blend : [gl.SRC_ALPHA, gl.ONE],
				pos : [0,0,0],
				parent : null,
				static : false
			});

			//assign after all prototype buffer stuff is written, so StarField can call Star can use it during ctor

			//now that we've built all our star system data ... add it to the star field
			if (starSystems.length > 1) thiz.addStarSystems();

		};
		xhr.send();
	};

	//this is called after the exoplanets load
	this.addStarSystems = function() {
		if (stars.buffer === undefined) return;
console.log('adding star systems to star fields and vice versa');	
		assert(starSystems.length > 1);

		//add buffer points

		var array = stars.buffer.data;	//8 fields: x y z vx vy vz absmag colorIndex
		array = Array.prototype.slice.call(array);	//to js array

		for (var i = 0; i < starSystems.length; ++i) {
			var starSystem = starSystems[i];
			starSystem.starfieldIndex = array.length / 5;
			array.push(starSystem.pos[0] / this.renderScale);
			array.push(starSystem.pos[1] / this.renderScale);
			array.push(starSystem.pos[2] / this.renderScale);
			array.push(0);
			array.push(0);
			array.push(0);
			array.push(5);	//abs mag
			array.push(0);	//color index
		}

		stars.buffer.setData(array, gl.STATIC_DRAW, true);

		//then add named stars to the starSystem array
		// if they're not there ... and I'm pretty sure they're all not there ...

		for (var i = 0; i < namedStars.length; ++i) {
			var args = namedStars[i];

			if (({
				Sun:1,
				"Kapteyn's Star":1,
				"Fomalhaut":1,
			})[args.name]) continue;	//hmm... there are some double-ups ...
			
			var starSystem = new StarSystem();

			starSystem.name = args.name;

			//index in the stars.buffer, stride of stars.buffer.dim
			//preserved since the original 1000 or however many from the HYG database are not moved
			starSystem.starfieldIndex = args.index;

			//note stars.buffer holds x y z abs-mag color-index
			starSystem.pos = [
				stars.buffer.data[0 + starSystem.starfieldIndex * stars.buffer.dim] * this.renderScale,
				stars.buffer.data[1 + starSystem.starfieldIndex * stars.buffer.dim] * this.renderScale,
				stars.buffer.data[2 + starSystem.starfieldIndex * stars.buffer.dim] * this.renderScale
			];

			//starSystem.sourceData = ...

			//now add a single star to the new star system ...

			var body = mergeInto(new Planet(), {
				type : 'star',
				name : args.name,
				//mass : undefined,
				//radius : undefined,
				//sourceData : {},
				//parent : undefined,
				starSytsem : starSystem,
				hide : false,
				isExoplanet : true
			});
			vec3.add(body.pos, body.pos, starSystem.pos);
			initPlanetColorSchRadiusAngle(body);
			initPlanetSceneLatLonLineObjs(body);
			calcKeplerianOrbitalElements(body, false);
			
			starSystem.planets.push(body);
			starSystem.stars.push(body);
			
			starSystem.doneBuildingPlanets();
		
			starSystem.initPlanets = starSystem.clonePlanets();
			starSystemForNames[starSystem.name] = starSystem;
			starSystem.index = starSystems.length;
			starSystems.push(starSystem);

		}

		initStarsControls();
	};

	this.draw = function(
		tanFovY,
		picking,
		viewPosInv,
		invRotMat,
		distFromSolarSystemInLyr,
		distFromSolarSystemInMpc
	) {
		var thiz = this;

		if (!showStars) return;
			
		gl.clear(gl.DEPTH_BUFFER_BIT)
		vec3.scale(viewPosInv, glutil.view.pos, -1/this.renderScale);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);
		
		if (distFromSolarSystemInLyr < this.maxDistInLyr &&
			stars.sceneObj !== undefined && 
			orbitTarget !== undefined &&
			orbitTarget.pos !== undefined)
		{
			var pointSize = 1000 * canvas.width * Math.sqrt(distFromSolarSystemInMpc) * metersPerUnits.pc / this.renderScale / tanFovY;
			if (!picking) {
				gl.disable(gl.DEPTH_TEST);
				stars.sceneObj.uniforms.visibleMagnitudeBias = starsVisibleMagnitudeBias;
				stars.sceneObj.pos[0] = -orbitTarget.pos[0] / this.renderScale;
				stars.sceneObj.pos[1] = -orbitTarget.pos[1] / this.renderScale;
				stars.sceneObj.pos[2] = -orbitTarget.pos[2] / this.renderScale;
				stars.sceneObj.uniforms.pointSize = pointSize;
				stars.sceneObj.draw();
				gl.enable(gl.DEPTH_TEST);
			} else {
				//I've got to call 'draw' to have the SceneObject matrixes calculated correctly
				//that means I've got to push/pop glutil.scene.projMat and load it with the pick projMat
				//but I really can't use attrs or uniforms because GLUtil right now merges *only* and I need it to replace ...
				if (allowSelectStars) {
				
					/* TODO render via buffer?
					//also TODO:
					//if (list === starfield && target == orbitStarSystem) continue;
					//...and then there's the fact that the old ray code only checked the exoplanets
					// whereas this is rendering all hyg stars ...
					//so I'll turn it off for now
					pickObject.drawPoints({
						sceneObj : stars.sceneObj,
						targetCallback : function(index) {
							for (var i = 0; i < namedStars.length; ++i) {
								if (namedStars[i].index == index) {	//namedStars[i].index is the index in the hyg 120,000 star buffer
									//if (orbitStarSystem == starSystems[whatever starsystem index matches]) return
									return starSystems[i];
								}
							}
						},
						pointSize : pointSize,
						pointSizeScaleWithDist : true,
						//defined in shader
						pointSizeMin : .25,
						pointSizeMax : 5
					});
					*/
					/* until then ... manually? */
					$.each(starSystems, function(i,starSystem) {
						if (starSystem !== orbitStarSystem) {
							pointObj.attrs.vertex.data[0] = (starSystem.pos[0] - orbitTarget.pos[0]) / thiz.renderScale;
							pointObj.attrs.vertex.data[1] = (starSystem.pos[1] - orbitTarget.pos[1]) / thiz.renderScale;
							pointObj.attrs.vertex.data[2] = (starSystem.pos[2] - orbitTarget.pos[2]) / thiz.renderScale;
							pointObj.attrs.vertex.updateData();
							pickObject.drawPoints({
								sceneObj : pointObj,
								targetCallback : function() { return starSystem; },
								pointSize : pointSize,
								pointSizeScaleWithDist : true,
								pointSizeMin : .25,
								pointSizeMax : 5
							});
						}
					});
					/**/
				}
			}
		}

		//add overlay text
		if (!picking) {
			for (var starSystemIndex = 0; starSystemIndex < starSystems.length; ++starSystemIndex) {
				var starSystem = starSystems[starSystemIndex];
				
				var withinSolarSystem = false;
				for (var planetIndex = 0; planetIndex < starSystem.planets.length; ++planetIndex) {
					if (starSystem.planets[planetIndex].orbitVisRatio > 1) {
						withinSolarSystem = true;
						break;
					}
				}
				if (withinSolarSystem) continue;
				
				//if the star system is close enough
				//TODO use luminance
				var deltaX = starSystem.pos[0] - orbitTarget.pos[0];
				var deltaY = starSystem.pos[1] - orbitTarget.pos[1];
				var deltaZ = starSystem.pos[2] - orbitTarget.pos[2];
				var dist = Math.sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
				//TODO if the dist is lower than the max radius of the star system then don't draw it
				var ratio = dist / orbitTargetDistance
				if (0.1 < ratio && ratio < 10) {
					overlayTexts.add(starSystem);
				}
			}
		}
	};
};
