/*
all stars of the milky way
TODO OOP this, and make one per galaxy (which we have observed stars within) 
*/

var showStars = true;
var starPointSizeScale = 3;
var starPointSizeBias = -3;

var starPointAlpha = 1;
var allowSelectStars = true;

var bubbleStartFadeDistInLyr = .25;
var bubbleStopFadeDistInLyr = 1.25;

//TODO merge with starSystems[] ... keep the StarField for point rendering of all StarSystems (or make it a Galaxy object, honestly, that's where thignsn are going)
// and remove StarInField ... make that just StarSystem (even for zero-planet systems)
//
//only instanciate these for the named stars.  87 in all.

var drawConstellationColorScalar = 20;
var drawConstellationPointSizeMax = 10;

//hardwired from the color tex
var colorTempMin = 1000;
var colorTempMax = 40000;

var LSun = 3.828e+26; 	// Watts
var L0 = 3.0128e+28;	// Watts
var LSunOverL0 = LSun / L0;

//maps from constellation index to list of vertex indexes in the point cloud
//var indexesForConstellations = [];

//maps from point cloud index to constellation index 
var constellationForIndex = [];

var starfield = new function() {
	this.maxDistInLyr = 5000;
	//this.renderScale = 1e+5;	//znear=1e+4, zfar=1e+25, something that puts all stars between?
	//this.renderScale = 1e+10;	//worked wrt all the float precision of glsl to render points, but the shader Pc dist calcs are ~ machine precision
	//this.renderScale = 1e+11;	//draws
	//this.renderScale = 1e+12;	//draws
	//this.renderScale = 1e+13;	//draws
	//this.renderScale = 1e+14;	//draws
	//this.renderScale = 1e+15;	//draws some, but znear starts clipping
	//this.renderScale = 1e+16;	//doesn't draw
	//this.renderScale = 1e+17;	// gives us Pc/scale ~ 3.24 (on the order of 1) ... but it doesn't show. why not?  znear zfar?  but we have a custom depth function?
	this.renderScale = metersPerUnits.pc;	//1:1 Pc/scale
	
	this.init = function() {
		var thiz = this;
		
		var bubbleTexWidth = 64;
		var bubbleTexData = new Uint8Array(bubbleTexWidth * bubbleTexWidth * 3);
		for (var j = 0; j < bubbleTexWidth; ++j) {
			var y = 2 * ((j+.5) / bubbleTexWidth - .5);
			var ay = Math.abs(y);
			for (var i = 0; i < bubbleTexWidth; ++i) {
				var x = 2 * ((i+.5) / bubbleTexWidth - .5);
				var ax = Math.abs(x);
				var r = Math.sqrt(x*x + y*y);
				var lum = Math.pow(r,10) * (1 - r);
				bubbleTexData[0+3*(i+j*bubbleTexWidth)] = 255*Math.clamp(lum,0,1);
				bubbleTexData[1+3*(i+j*bubbleTexWidth)] = 255*Math.clamp(lum,0,1);
				bubbleTexData[2+3*(i+j*bubbleTexWidth)] = 255*Math.clamp(lum,0,1);
			}
		}
		this.bubbleTex = new glutil.Texture2D({
			width : bubbleTexWidth,
			height : bubbleTexWidth,
			internalFormat : gl.RGB,
			format : gl.RGB,
			type : gl.UNSIGNED_BYTE,
			magFilter : gl.LINEAR,
			minFilter : gl.LINEAR_MIPMAP_LINEAR,
			data : bubbleTexData,
			//alignment : 1,
			generateMipmap : true
		});
		
		console.log('init star tex...');
		var starTexSize = 64;
		var starTexData = new Uint8Array(starTexSize * starTexSize * 3);
		for (var j = 0; j < starTexSize; ++j) {
			for (var i = 0; i < starTexSize; ++i) {
				var u = ((i+.5) / starTexSize) * 2 - 1;
				var v = ((j+.5) / starTexSize) * 2 - 1;
				/* some kind of 4-point star tex * /				
				u *= .5;
				v *= .5;
				var av = Math.abs(v);
				var au = Math.abs(u);
				var rL2 = Math.sqrt(Math.pow(au,2) + Math.pow(av,2));
				var rL_2 = Math.pow(au,1/2) + Math.pow(av,1/2);
				var r = rL2 + rL_2;
				var sigma = 1;
				var rs = r / sigma;
				var l = Math.exp(-rs*rs);
				/**/
				
				/* ring tex */
				var r = Math.sqrt(u*u + v*v)
				var l = Math.exp(-50 * Math.pow(r - .75, 2))
				/**/
				
				l = Math.floor(255 * Math.clamp(l, 0, 1));
				starTexData[0+3*(i+j*starTexSize)] = l;
				starTexData[1+3*(i+j*starTexSize)] = l;
				starTexData[2+3*(i+j*starTexSize)] = l;
			}
		}
		this.starTex = new glutil.Texture2D({
			width : starTexSize,
			height : starTexSize,
			internalFormat : gl.RGB,
			format : gl.RGB,
			type : gl.UNSIGNED_BYTE,
			magFilter : gl.LINEAR,
			minFilter : gl.LINEAR_MIPMAP_LINEAR,
			data : starTexData,
			//alignment : 1,
			generateMipmap : true
		});
		/*{
			var level = 1;
			var lastStarTexData = starTexData;
			for (var w = starTexSize>>1; w; w>>=1, ++level) {
				var starTexData = new Uint8Array(w * w * 4);
				var lastW = w<<1;
				for (var j = 0; j < w; ++j) {
					for (var i = 0; i < w; ++i) {
						for (var k = 0; k < 4; ++k) {
							starTexData[k+4*(i+w*j)] = Math.max(
								lastStarTexData[k+4*((2*i+0)+2*w*(2*j+0))],
								lastStarTexData[k+4*((2*i+1)+2*w*(2*j+0))],
								lastStarTexData[k+4*((2*i+0)+2*w*(2*j+1))],
								lastStarTexData[k+4*((2*i+1)+2*w*(2*j+1))]);
						}
					}
				}
				this.starTex.bind();
				this.starTex.setImage({
					level : level,
					internalFormat : gl.RGBA,
					format : gl.RGBA,
					type : gl.UNSIGNED_BYTE,
					width : w,
					height : w,
					data : starTexData
				});
				lastStarTexData = starTexData;
			}
		}
		*/	
		
		var loadColorTempTex = function(done) {
			console.log('init color index tex...');
			//going by http://stackoverflow.com/questions/21977786/star-b-v-color-index-to-apparent-rgb-color
			//though this will be helpful too: http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
			var tempImg = new Image();
			tempImg.onload = function() {
				thiz.tempTex = new glutil.Texture2D({
					data : tempImg,
					magFilter : gl.LINEAR,
					minFilter : gl.NEAREST,
					wrap : {
						s : gl.CLAMP_TO_EDGE,
						t : gl.CLAMP_TO_EDGE
					}
				});
				if (done) done();
			};
			tempImg.onerror = function() {
				console.log('failed to find color index texture');
			};
			tempImg.src = 'colorForTemp.png';
		};	

		var initColorIndexShader = function(done) {
			console.log('init color index shader...');
			//currently only used by starfield shader
			//considering use with planet shader
			thiz.colorIndexShader = new ModifiedDepthShaderProgram({
				vertexCode :
'#define M_1_LOG_10 '+floatToGLSL(1/Math.log(10))+'\n'+
'#define colorTempMin '+floatToGLSL(colorTempMin)+'\n'+
'#define colorTempMax '+floatToGLSL(colorTempMax)+'\n'+
'#define LSunOverL0 '+floatToGLSL(LSunOverL0)+'\n'+
'#define log10UnitRatio '+floatToGLSL(Math.log10(thiz.renderScale / metersPerUnits.pc))+'\n'+
					mlstr(function(){/*
attribute vec3 vertex;
attribute vec3 velocity;
attribute float luminosity;		// in solar luminosity units
attribute float temperature;	// in K

varying float lumv;
varying vec3 tempcolor;
//varying float discardv;

uniform mat4 mvMat;
uniform mat4 projMat;
uniform float starPointSizeScale;
uniform float starPointSizeBias;
uniform sampler2D tempTex;

//uniform float sliceRMin, sliceRMax;
//uniform float sliceLumMin, sliceLumMax;

void main() {
//	discardv = 0.;

//	if (luminosity < sliceLumMin || luminosity > sliceLumMax) {
//		discardv = 1.;
//		return;
//	}

	vec4 vtx4 = vec4(vertex, 1.);
	vec4 vmv = mvMat * vtx4;
	gl_Position = projMat * vmv;
//	gl_Position.z = depthfunction(gl_Position);

	//how to calculate this in fragment space ...
	// coordinates are in Pc
	//distance based on the eye position
	//this is only in Pc if our scale is 1:1 with Pc
	//I'm just keeping the names the same with the offline/visualize-stars.lua in the hopes that things will break the least when porting over
	float distInPcSq = dot(vmv.xyz, vmv.xyz);
	
	//log(distInPc^2) = 2 log(distInPc)
	//so log10(distInPc) = .5 log10(distInPc^2)
	//so log10(distInPc) = .5 / log(10) * log(distInPc^2)
	float log10DistInPc = (.5 * M_1_LOG_10) * log(distInPcSq) 
//		+ M_1_LOG_10 * log10UnitRatio
	;

	float LStarOverLSun = luminosity;
	float LStarOverL0 = LSunOverL0 * LStarOverLSun;
	float absoluteMagnitude = (-2.5 * M_1_LOG_10) * log(LStarOverL0);	// abs magn
	
	//apparent magnitude:
	//M = absolute magnitude
	//m = apparent magnitude
	//d = distance in parsecs
	//m = M - 5 + 5 * log10(d)
	float apparentMagnitude = absoluteMagnitude - 5. + 5. * log10DistInPc;

	//ok now on to the point size ...
	//and magnitude ...
	//hmm ... 
	//HDR will have to factor into here somehow ...
	gl_PointSize = (6.5 - apparentMagnitude) * starPointSizeScale + starPointSizeBias;

	lumv = 1.;
	
	// if the point size is < .5 then just make the star dimmer instead
	const float pointSizeMin = .5;
	float dimmer = gl_PointSize - pointSizeMin;
	if (dimmer < 0.) {
		gl_PointSize = pointSizeMin;
		lumv *= pow(2., dimmer);
	}
	
	//pointSizeMax TODO param
	if (gl_PointSize > 50.) gl_PointSize = 50.;

	//calculate color
	float tempfrac = (temperature - colorTempMin) / (colorTempMax - colorTempMin);
	tempcolor = texture2D(tempTex, vec2(tempfrac, .5)).rgb;
}
*/}),
				fragmentCode : mlstr(function(){/*
varying float lumv;
varying vec3 tempcolor;
//varying float discardv;

uniform float starPointAlpha;
uniform sampler2D starTex;

void main() {
	//can you discard vertices?  or only in geometry/tesselation shaders?
	//if (discardv > 0.) discard;

	float lumf = lumv;

	gl_FragColor = vec4(tempcolor * lumf * starPointAlpha, 1.) 
		* texture2D(starTex, gl_PointCoord);
}
*/})
			});
		
			if (done) done();
		};

		var loadStarField = function(done) {
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
			
			//unnamed is index 0
//			for (var k = 0; k < constellations.length; ++k) {
//				indexesForConstellations[k] = []; 
//			}
			
			var numElem = 9;
			xhr.onload = function(e) {
				var arrayBuffer = this.response;
				var data = new DataView(arrayBuffer);

				//units are in parsecs
				//don't forget velocity is not being rescaled (i'm not using it at the moment)
				var floatBuffer = new Float32Array(data.byteLength / Float32Array.BYTES_PER_ELEMENT);
				var len = floatBuffer.length;
var numClose = 0;			

				//stats on the data to verify
				var xMin = Infinity;
				var xMax = -Infinity;
				var yMin = Infinity;
				var yMax = -Infinity;
				var zMin = Infinity;
				var zMax = -Infinity;
				var rMin = Infinity;
				var rMax = -Infinity;
				var lumMin = Infinity;
				var lumMax = -Infinity;
				var tempMin = Infinity;
				var tempMax = -Infinity;

console.log("rescaling data by "+(metersPerUnits.pc / thiz.renderScale));
				for (var j = 0; j < len; ++j) {
					var x = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);
					if (j % numElem < 3) {
						if (Math.abs(x) > 1e+20) {
							console.log('star '+Math.floor(j/numElem)+' has coordinate that position exceeds fp resolution');
						}
						//convert parsecs to meters
						//and downscale by the starfield render scale
						x *= metersPerUnits.pc / thiz.renderScale;
					}
					floatBuffer[j] = x;
				
					//process every row:
					if (j % numElem == numElem - 1) {
						var x = floatBuffer[j - numElem + 1];
						var y = floatBuffer[j - numElem + 2];
						var z = floatBuffer[j - numElem + 3];
						var lum = floatBuffer[j - numElem + 7];
						var temp = floatBuffer[j - numElem + 8];
						var r = Math.sqrt(x*x + y*y + z*z);
						if (r < metersPerUnits.pc / thiz.renderScale) {
							numClose++;
						}
						if (x < xMin) xMin = x;
						if (x > xMax) xMax = x;
						if (y < yMin) yMin = y;
						if (y > yMax) yMax = y;
						if (z < zMin) zMin = z;
						if (z > zMax) zMax = z;
						if (r < rMin) rMin = r;
						if (r > rMax) rMax = r;
						if (lum < lumMin) lumMin = lum;
						if (lum > lumMax) lumMax = lum;
						if (temp < tempMin) tempMin = temp;
						if (temp > tempMax) tempMax = temp;

//						//constellation
//						var constellationIndex = floatBuffer[j - numElem + 9];
//						var vertexIndex = (j-8)/numElem;
//						indexesForConstellations[constellationIndex].push(vertexIndex);
//						constellationForIndex[vertexIndex] = constellationIndex;
					}
				}

//				for (var j = 0; j < constellations.length; ++j) {
//					indexesForConstellations[j] = new glutil.ElementArrayBuffer({
//						data : indexesForConstellations[j],
//						keep : true
//					});
//				}

console.log('star x range', xMin, xMax);
console.log('star y range', yMin, yMax);
console.log('star z range', zMin, zMax);
console.log('star r range', rMin, rMax);
console.log('star lum range', lumMin, lumMax);
console.log('star temp range', tempMin, tempMax);
console.log('num stars total:', floatBuffer.length / numElem);
console.log('num stars within 1pc:', numClose);
				
				//now that we have the float buffer ...
				thiz.buffer = new glutil.ArrayBuffer({data : floatBuffer, dim : numElem});
				thiz.sceneObj = new glutil.SceneObject({
					mode : gl.POINTS,
					shader : thiz.colorIndexShader,
					attrs : {
						vertex : new glutil.Attribute({buffer : thiz.buffer, size : 3, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 0}),	//xyz in Parsecs
						velocity : new glutil.Attribute({buffer : thiz.buffer, size : 3, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 3 * Float32Array.BYTES_PER_ELEMENT}),	//velocity 
						luminosity : new glutil.Attribute({buffer : thiz.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 6 * Float32Array.BYTES_PER_ELEMENT}),	//luminosity in LSun units
						temperature : new glutil.Attribute({buffer : thiz.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 7 * Float32Array.BYTES_PER_ELEMENT})	//temperature in K
					},
					uniforms : {
						tempTex : 0,
						starTex : 1,
						starPointSizeBias : starPointSizeBias,
						starPointSizeScale : starPointSizeScale,
						starPointAlpha : starPointAlpha 
					},
					texs : [
						assert(thiz.tempTex),
						assert(thiz.starTex)
					],
					blend : [gl.SRC_ALPHA, gl.ONE],
					pos : [0,0,0],
					parent : null,
					static : false
				});

				//assign after all prototype buffer stuff is written, so StarField can call Star can use it during ctor

				//now that we've built all our star system data ... add it to the star field
				//TODO combine these datasets offline, since there is some overlap, and that's causing duplicate stars
				if (starSystems.length > 1) thiz.addStarSystems();

				if (done) done();
			

				//ugly ugly code
				// see if the starsystems are already oloaded
				thiz.addStarSystems();
			};
			xhr.send();
		};
	

		loadColorTempTex(function() {
			initColorIndexShader(function() {
				loadStarField();
			});
		});
	};

	//this is called after the exoplanets load
	this.addStarSystems = function() {
		if (this.alreadyAddedStarSystems) return;	// already processed
		if (!starSystemsHasGotResults) return;		// starsystems not ready yet?
		if (this.buffer === undefined) return;		// we're not ready yet?

/* add stars from open exoplanet data */
		{
			console.log('adding star systems to star fields and vice versa');	
			assert(starSystems.length > 1);

			//add buffer points

			var array = this.buffer.data;	//9 fields: x y z vx vy vz luminostiy temperature constellation
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
				array.push(5);	//luminosity
				array.push(0);	//temperature
				array.push(0);	//constellation
			}

			this.buffer.setData(array, gl.STATIC_DRAW, true);
		}
/**/

		//then add named stars to the starSystem array

		for (var i = 0; i < namedStars.length; ++i) {
			var args = namedStars[i];

			if (({
				Sun:1,
				"Kapteyn's Star":1,
				"Fomalhaut":1,
			})[args.name]) continue;	//hmm... there are some double-ups ...
			
			var starSystem = new StarSystem();

			starSystem.name = args.name;

			//index in the this.buffer, stride of this.buffer.dim
			//preserved since the original 1000 or however many from the HYG database are not moved
			starSystem.starfieldIndex = args.index;

			starSystem.constellationIndex = constellationForIndex[starSystem.starfieldIndex];

			//note this.buffer holds x y z abs-mag color-index
			starSystem.pos = [
				this.buffer.data[0 + starSystem.starfieldIndex * this.buffer.dim] * this.renderScale,
				this.buffer.data[1 + starSystem.starfieldIndex * this.buffer.dim] * this.renderScale,
				this.buffer.data[2 + starSystem.starfieldIndex * this.buffer.dim] * this.renderScale
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
			body.initColorSchRadiusAngle();
			body.initSceneLatLonLineObjs();
			body.getKOEFromSourceData();
			
			starSystem.planets.push(body);
			starSystem.stars.push(body);
			
			starSystem.doneBuildingPlanets();
		
			starSystem.initPlanets = starSystem.clonePlanets();
			starSystemForNames[starSystem.name] = starSystem;
			starSystem.index = starSystems.length;
			starSystems.push(starSystem);

		}

		starSystemsExtra.initStarsControls();

		this.alreadyAddedStarSystems = true;
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
			this.sceneObj !== undefined && 
			orbitTarget !== undefined &&
			orbitTarget.pos !== undefined)
		{
			var pointSize = 
				canvas.width 
				/ 100
				* metersPerUnits.pc 
				/ this.renderScale 
				/ tanFovY;

			this.sceneObj.uniforms.starPointSizeBias = starPointSizeBias;
			this.sceneObj.uniforms.starPointSizeScale = starPointSizeScale;
			this.sceneObj.uniforms.starPointAlpha = starPointAlpha;
			this.sceneObj.pos[0] = -orbitTarget.pos[0] / this.renderScale;
			this.sceneObj.pos[1] = -orbitTarget.pos[1] / this.renderScale;
			this.sceneObj.pos[2] = -orbitTarget.pos[2] / this.renderScale;

			if (!picking) {
				
				gl.disable(gl.DEPTH_TEST);

var pushZNear = glutil.view.zNear;
var pushZFar = glutil.view.zFar;
glutil.view.zNear = 1e-3;
glutil.view.zFar = 1e+6;
glutil.updateProjection();

				this.sceneObj.draw();

glutil.view.zNear = pushZNear;
glutil.view.zFar = pushZFar;
glutil.updateProjection();

				//only draw bubbles around stars once we're out of the star system
				//then fade them into display
				if (distFromSolarSystemInLyr > bubbleStartFadeDistInLyr) {
					//var colorScale = starfield.colorScale || .5;//(distFromSolarSystemInLyr - bubbleStartFadeDistInLyr) / (bubbleStopFadeDistInLyr - bubbleStartFadeDistInLyr);
					//colorScale = Math.clamp(colorScale, 0, 1);

					var pointSize = 
						canvas.width 
						/ 10
						* metersPerUnits.pc 
						/ this.renderScale 
						/ tanFovY;
					
					this.sceneObj.draw({
						uniforms : {
							starPointSizeBias : starPointSizeBias,
							starPointSizeScale : starPointSizeScale,
							starPointAlpha : starPointAlpha  
						},
						texs : [this.tempTex, this.bubbleTex]
					});
				}


				// draw constellations

//				for (var k = 0; k < constellations.length; ++k) {
//					if (displayConstellations[k]) {
//						this.sceneObj.geometry.indexes = indexesForConstellations[k];
//						this.sceneObj.draw({
//							uniforms : {
//								starPointSizeBias : starPointSizeBias,
//								starPointSizeScale : starPointSizeScale,
//								colorScale : drawConstellationColorScalar,
//							},
//							texs : [this.tempTex, this.bubbleTex]
//						});
//					
///**/
//						//and draw some bboxes around it
//						var minmax = ['min', 'max'];
//						var sunPos = solarSystem.planets[solarSystem.indexes.Sun].pos;	//I store the data wrt the sun's position
//						for (var v1 = 0; v1 < 8; ++v1) {
//							for (var edge = 0; edge < 3; ++edge) {
//								var v2 = v1 ^ (1 << edge);
//								if (v1 > v2) continue;
//
//								var con = constellations[k];
//
//								var ra1 = con.ra[ minmax[v1&1] ];
//								var ra2 = con.ra[ minmax[v2&1] ];
//								var dec1 = con.dec[ minmax[(v1>>1)&1] ];
//								var dec2 = con.dec[ minmax[(v2>>1)&1] ];
//								var dist1 = con.dist[ minmax[(v1>>2)&1] ];
//								var dist2 = con.dist[ minmax[(v2>>2)&1] ];
//
//								dist1 *= metersPerUnits.pc / this.renderScale;
//								dist2 *= metersPerUnits.pc / this.renderScale;
//								lineObj.attrs.vertex.data[0] = dist1 * Math.cos(ra1) * Math.cos(dec1) + (sunPos[0] - orbitTarget.pos[0]) / this.renderScale;
//								lineObj.attrs.vertex.data[1] = dist1 * Math.sin(ra1) * Math.cos(dec1) + (sunPos[1] - orbitTarget.pos[1]) / this.renderScale;
//								lineObj.attrs.vertex.data[2] = dist1 * Math.sin(dec1)                 + (sunPos[2] - orbitTarget.pos[2]) / this.renderScale;
//								lineObj.attrs.vertex.data[3] = dist2 * Math.cos(ra2) * Math.cos(dec2) + (sunPos[0] - orbitTarget.pos[0]) / this.renderScale;
//								lineObj.attrs.vertex.data[4] = dist2 * Math.sin(ra2) * Math.cos(dec2) + (sunPos[1] - orbitTarget.pos[1]) / this.renderScale;
//								lineObj.attrs.vertex.data[5] = dist2 * Math.sin(dec2)                 + (sunPos[2] - orbitTarget.pos[2]) / this.renderScale;
//					
//								lineObj.attrs.vertex.updateData();
//								lineObj.draw({uniforms : { color : [1,1,1,1] }});
//							}
//						}
///**/					
//					}
//				}
				this.sceneObj.geometry.indexes = undefined;



				gl.enable(gl.DEPTH_TEST);
			
			} else {
				//I've got to call 'draw' to have the SceneObject matrixes calculated correctly
				//that means I've got to push/pop glutil.scene.projMat and load it with the pick projMat
				//but I really can't use attrs or uniforms because GLUtil right now merges *only* and I need it to replace ...
				if (allowSelectStars) {
				
					/* TODO render via buffer? * /
					//also TODO:
					//if (list === this && target == orbitStarSystem) continue;
					//...and then there's the fact that the old ray code only checked the exoplanets
					// whereas this is rendering all hyg stars ...
					//so I'll turn it off for now
					pickObject.drawPoints({
						sceneObj : this.sceneObj,
						targetCallback : function(index) {
							for (var i = 0; i < namedStars.length; ++i) {
								if (namedStars[i].index == index) {	//namedStars[i].index is the index in the hyg 120,000 star buffer
									//if (orbitStarSystem == starSystems[whatever starsystem index matches]) return
									return starSystems[i];
								}
							}
						},
						pointSizeScaleWithDist : true,
						//defined in shader
						pointSizeMin : 0,
					});
					/**/
					/* until then ... manually? * /
					$.each(starSystems, function(i,starSystem) {
						if (starSystem !== orbitStarSystem) {
							if (starSystem.constellationIndex === undefined ||
								displayConstellations[starSystem.constellationIndex]
							) {
								pointObj.attrs.vertex.data[0] = (starSystem.pos[0] - orbitTarget.pos[0]) / thiz.renderScale;
								pointObj.attrs.vertex.data[1] = (starSystem.pos[1] - orbitTarget.pos[1]) / thiz.renderScale;
								pointObj.attrs.vertex.data[2] = (starSystem.pos[2] - orbitTarget.pos[2]) / thiz.renderScale;
								pointObj.attrs.vertex.updateData();
								pickObject.drawPoints({
									sceneObj : pointObj,
									targetCallback : function() { return starSystem; },
									pointSizeScaleWithDist : true,
									pointSizeMin : 0,
								});
							}
						}
					});
					/**/
					/* only visible constellations ... 
					but TODO we need a starSytsem object * /
					for (var k = 0; k < constellations.length; ++k) {
						if (displayConstellations[k]) {
							for (var j = 0; j < indexesForConstellations[k].count; ++j) {
								var i = indexesForConstellations[k].data[j];
								var x = thiz.buffer.data[0 + thiz.buffer.dim * i];
								var y = thiz.buffer.data[1 + thiz.buffer.dim * i];
								var z = thiz.buffer.data[2 + thiz.buffer.dim * i];
								pointObj.attrs.vertex.data[0] = x - orbitTarget.pos[0] / thiz.renderScale;
								pointObj.attrs.vertex.data[1] = y - orbitTarget.pos[1] / thiz.renderScale;
								pointObj.attrs.vertex.data[2] = z - orbitTarget.pos[2] / thiz.renderScale;
								pointObj.attrs.vertex.updateData();
								pickObject.drawPoints({
									sceneObj : pointObj,
									targetCallback : function() { return starSystem; },
									pointSizeScaleWithDist : true,
									pointSizeMin : 0,
									pointSizeMax : 5
								});
							}
						}
					}
					/**/
				}
			}
		}

		/*
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
		*/
	};
};
