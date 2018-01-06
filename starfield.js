/*
all stars of the milky way
TODO OOP this, and make one per galaxy (which we have observed stars within) 
*/

var showStars = true;
var starsVisibleMagnitudeBias = 0;
var allowSelectStars = true;

var bubbleStartFadeDistInLyr = .25;
var bubbleStopFadeDistInLyr = 1.25;

//TODO merge with starSystems[] ... keep the StarField for point rendering of all StarSystems (or make it a Galaxy object, honestly, that's where thignsn are going)
// and remove StarInField ... make that just StarSystem (even for zero-planet systems)
//
//only instanciate these for the named stars.  87 in all.

var starfield = new function() {
	this.maxDistInLyr = 5000;
	this.renderScale = 1e+10;
	
	this.init = function() {
		
		var colorIndexMin = 2.;
		var colorIndexMax = -.4;
		
		console.log('init color index tex...');
		//going by http://stackoverflow.com/questions/21977786/star-b-v-color-index-to-apparent-rgb-color
		//though this will be helpful too: http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
		var colorIndexTexWidth = 1024;
		this.colorIndexTex = new glutil.Texture2D({
			width : colorIndexTexWidth,
			height : 1,
			internalFormat : gl.RGBA,
			format : gl.RGBA,	//TODO function callback support for format == gl.RGB
			type : gl.UNSIGNED_BYTE,
			magFilter : gl.LINEAR,
			minFilter : gl.NEAREST,
			wrap : gl.CLAMP_TO_EDGE,
			data : function(texX, texY) {
//console.log('texX',texX);
				var frac = (texX + .5) / colorIndexTexWidth;
//console.log('frac',frac);
				var bv = frac * (colorIndexMax - colorIndexMin) + colorIndexMin;
//console.log('bv',bv);

//from here on out is the OP of the above article
if (false) {

				var t = 4600 * ((1 / ((.92 * bv) + 1.7)) + (1 / ((.92 * bv) + .62)));
//console.log('t',t);

				var x, y = 0;

				if (t>=1667 && t<=4000) {
					x = ((-0.2661239 * Math.pow(10,9)) / Math.pow(t,3)) + ((-0.2343580 * Math.pow(10,6)) / Math.pow(t,2)) + ((0.8776956 * Math.pow(10,3)) / t) + 0.179910;
				} else if (t > 4000 && t <= 25000) {
					x = ((-3.0258469 * Math.pow(10,9)) / Math.pow(t,3)) + ((2.1070379 * Math.pow(10,6)) / Math.pow(t,2)) + ((0.2226347 * Math.pow(10,3)) / t) + 0.240390;
				}
//console.log('x',x);

				if (t >= 1667 && t <= 2222) {
					y = -1.1063814 * Math.pow(x,3) - 1.34811020 * Math.pow(x,2) + 2.18555832 * x - 0.20219683;
				} else if (t > 2222 && t <= 4000) {
					y = -0.9549476 * Math.pow(x,3) - 1.37418593 * Math.pow(x,2) + 2.09137015 * x - 0.16748867;
				} else if (t > 4000 && t <= 25000) {
					y = 3.0817580 * Math.pow(x,3) - 5.87338670 * Math.pow(x,2) + 3.75112997 * x - 0.37001483;
				}
//console.log('y',y);

				//the rest is found at http://en.wikipedia.org/wiki/SRGB

				// xyY to XYZ, Y = 1
				var Y = (y == 0)? 0 : 1;
				var X = (y == 0)? 0 : (x * Y) / y;
				var Z = (y == 0)? 0 : ((1 - x - y) * Y) / y;
//console.log('X',X);
//console.log('Y',Y);
//console.log('Z',Z);

				var R_linear = 3.2406 * X - 1.5372 * Y - .4986 * Z;
				var G_linear = -.9689 * X + 1.8758 * Y + .0415 * Z;
				var B_linear = .0557 * X - .2040 * Y + 1.0570 * Z;
//console.log('R_linear',R_linear);
//console.log('G_linear',G_linear);
//console.log('B_linear',B_linear);

				var srgbGammaAdjust = function(C_linear) {
					var a = .055;
					if (C_linear <= .0031308) return 12.92 * C_linear;
					return (1 + a) * Math.pow(C_linear, 1/0.5) - a;
				};

				var R_srgb = srgbGammaAdjust(R_linear);
				var G_srgb = srgbGammaAdjust(G_linear);
				var B_srgb = srgbGammaAdjust(B_linear);
//console.log('R_srgb',R_srgb);
//console.log('G_srgb',G_srgb);
//console.log('B_srgb',B_srgb);

				var result = [
					Math.clamp(R_srgb, 0, 1),
					Math.clamp(G_srgb, 0, 1),
					Math.clamp(B_srgb, 0, 1),
					1];

				return result;

}	//...and this is the reply ...
if (true) {

				var t;  r=0.0; g=0.0; b=0.0; if (t<-0.4) t=-0.4; if (t> 2.0) t= 2.0;
				if ((bv>=-0.40)&&(bv<0.00)) { t=(bv+0.40)/(0.00+0.40); r=0.61+(0.11*t)+(0.1*t*t); }
				else if ((bv>= 0.00)&&(bv<0.40)) { t=(bv-0.00)/(0.40-0.00); r=0.83+(0.17*t)		  ; }
				else if ((bv>= 0.40)&&(bv<2.10)) { t=(bv-0.40)/(2.10-0.40); r=1.00				   ; }
				if ((bv>=-0.40)&&(bv<0.00)) { t=(bv+0.40)/(0.00+0.40); g=0.70+(0.07*t)+(0.1*t*t); }
				else if ((bv>= 0.00)&&(bv<0.40)) { t=(bv-0.00)/(0.40-0.00); g=0.87+(0.11*t)		  ; }
				else if ((bv>= 0.40)&&(bv<1.60)) { t=(bv-0.40)/(1.60-0.40); g=0.98-(0.16*t)		  ; }
				else if ((bv>= 1.60)&&(bv<2.00)) { t=(bv-1.60)/(2.00-1.60); g=0.82		 -(0.5*t*t); }
				if ((bv>=-0.40)&&(bv<0.40)) { t=(bv+0.40)/(0.40+0.40); b=1.00				   ; }
				else if ((bv>= 0.40)&&(bv<1.50)) { t=(bv-0.40)/(1.50-0.40); b=1.00-(0.47*t)+(0.1*t*t); }
				else if ((bv>= 1.50)&&(bv<1.94)) { t=(bv-1.50)/(1.94-1.50); b=0.63		 -(0.6*t*t); }

				return [r,g,b,1];
}
			}
		});

		console.log('init color index shader...');
		//currently only used by starfield shader
		//considering use with planet shader
		this.colorIndexShader = new ModifiedDepthShaderProgram({
			vertexCode :
'#define M_LOG_10 '+floatToGLSL(Math.log(10))+'\n'+
'#define COLOR_INDEX_MIN '+floatToGLSL(colorIndexMin)+'\n'+
'#define COLOR_INDEX_MAX '+floatToGLSL(colorIndexMax)+'\n'+
			mlstr(function(){/*
attribute vec3 vertex;
attribute vec3 velocity;
attribute float absoluteMagnitude;
attribute float colorIndex;

uniform mat4 mvMat;
uniform mat4 projMat;
uniform float visibleMagnitudeBias;
uniform float pointSize;
uniform float pointSizeMax;
uniform sampler2D colorIndexTex;

varying float alpha;
varying vec3 color;

void main() {
	vec4 vtx4 = mvMat * vec4(vertex, 1.);

	//calculate apparent magnitude, convert to alpha
	float distanceInParsecs = length(vtx4.xyz) / ( */}) 
				+ floatToGLSL(metersPerUnits.pc / starfield.renderScale) 
				+ mlstr(function(){/*);	//convert distance to parsecs
	
	//https://en.wikipedia.org/wiki/Apparent_magnitude
	float apparentMagnitude = absoluteMagnitude - 5. * (1. - log(distanceInParsecs) / M_LOG_10);
	
	//not sure where I got this one from ...
	float log100alpha = -.2*(apparentMagnitude - visibleMagnitudeBias);
	alpha = pow(100., log100alpha);
	
	alpha = clamp(alpha, 0., 1.);

	//calculate color
	color = texture2D(colorIndexTex, vec2((colorIndex - COLOR_INDEX_MIN) / (COLOR_INDEX_MAX - COLOR_INDEX_MIN), .5)).rgb;

	gl_Position = projMat * vtx4;
	
	//TODO point sprite / point spread function?
	gl_PointSize = pointSize / gl_Position.w;
	gl_PointSize = min(gl_PointSize, pointSizeMax);
	
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
uniform sampler2D starTex;
uniform float colorScale;
varying vec3 color;
varying float alpha;
void main() {
	gl_FragColor = colorScale * vec4(color, alpha) * texture2D(starTex, gl_PointCoord);
}
*/})
		});
		
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
		var starTexWidth = 64;
		var starTexData = new Uint8Array(starTexWidth * starTexWidth * 3);
		for (var j = 0; j < starTexWidth; ++j) {
			var y = (j+.5) / starTexWidth - .5;
			var ay = Math.abs(y);
			for (var i = 0; i < starTexWidth; ++i) {
				var x = (i+.5) / starTexWidth - .5;
				var ax = Math.abs(x);
				var rL2 = Math.sqrt(Math.pow(ax,2) + Math.pow(ay,2));
				var rL_2 = Math.pow(ax,1/2) + Math.pow(ay,1/2);
				var r = rL2 + rL_2;
				var sigma = 1;
				var rs = r / sigma;
				var lum = Math.exp(-rs*rs);
				starTexData[0+3*(i+j*starTexWidth)] = 255*Math.clamp(lum,0,1);
				starTexData[1+3*(i+j*starTexWidth)] = 255*Math.clamp(lum,0,1);
				starTexData[2+3*(i+j*starTexWidth)] = 255*Math.clamp(lum,0,1);
			}
		}
		this.starTex = new glutil.Texture2D({
			width : starTexWidth,
			height : starTexWidth,
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
			for (var w = starTexWidth>>1; w; w>>=1, ++level) {
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
var numClose = 0;			
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
				if (j%numElem==2) {
var x = floatBuffer[j-2];
var y = floatBuffer[j-1];
var z = floatBuffer[j];
var dist = Math.sqrt(x*x + y*y + z*z);
if (dist < metersPerUnits.pc / thiz.renderScale) {
	numClose++;
}
				}
			}
console.log('num stars within 1pc:', numClose);
			//now that we have the float buffer ...
			thiz.buffer = new glutil.ArrayBuffer({data : floatBuffer, dim : numElem});
			thiz.sceneObj = new glutil.SceneObject({
				mode : gl.POINTS,
				attrs : {
					vertex : new glutil.Attribute({buffer : thiz.buffer, size : 3, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 0}),	//xyz abs-mag
					velocity : new glutil.Attribute({buffer : thiz.buffer, size : 3, stride : numElem * Float32Array.BYTES_PER_ELEM, offset : 3 * Float32Array.BYTES_PER_ELEMENT}),	//velocity
					absoluteMagnitude : new glutil.Attribute({buffer : thiz.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEM, offset : 6 * Float32Array.BYTES_PER_ELEMENT}),
					colorIndex : new glutil.Attribute({buffer : thiz.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 7 * Float32Array.BYTES_PER_ELEMENT})	//color-index
				},
				uniforms : {
					colorIndexTex : 0,
					starTex : 1,
					pointSize : 1,
					pointSizeMax : 5,
					visibleMagnitudeBias : starsVisibleMagnitudeBias,
					colorScale : 1
				},
				shader : thiz.colorIndexShader,
				texs : [thiz.colorIndexTex, thiz.starTex],
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
		if (this.buffer === undefined) return;
console.log('adding star systems to star fields and vice versa');	
		assert(starSystems.length > 1);

		//add buffer points

		var array = this.buffer.data;	//8 fields: x y z vx vy vz absmag colorIndex
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

		this.buffer.setData(array, gl.STATIC_DRAW, true);

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

			//index in the this.buffer, stride of this.buffer.dim
			//preserved since the original 1000 or however many from the HYG database are not moved
			starSystem.starfieldIndex = args.index;

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
				
			this.sceneObj.uniforms.visibleMagnitudeBias = starsVisibleMagnitudeBias;
			this.sceneObj.pos[0] = -orbitTarget.pos[0] / this.renderScale;
			this.sceneObj.pos[1] = -orbitTarget.pos[1] / this.renderScale;
			this.sceneObj.pos[2] = -orbitTarget.pos[2] / this.renderScale;
			this.sceneObj.uniforms.pointSize = pointSize;
			
			if (!picking) {
				
				gl.disable(gl.DEPTH_TEST);
				
				this.sceneObj.draw();

				//only draw bubbles around stars once we're out of the star system
				//then fade them into display
				if (distFromSolarSystemInLyr > bubbleStartFadeDistInLyr) {
					var alpha = (distFromSolarSystemInLyr - bubbleStartFadeDistInLyr) / (bubbleStopFadeDistInLyr - bubbleStartFadeDistInLyr);
					alpha = Math.clamp(alpha, 0, 1);

					var pointSize = 
						canvas.width 
						/ 10
						* metersPerUnits.pc 
						/ this.renderScale 
						/ tanFovY;
					this.sceneObj.draw({
						uniforms : {
							pointSize : pointSize,
							pointSizeMax : 1000,
							visibleMagnitudeBias : 10,	//TODO just use a different shader
							colorScale : alpha * 2
						},
						texs : [this.colorIndexTex, this.bubbleTex]
					});
				}

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
						pointSize : pointSize,
						pointSizeScaleWithDist : true,
						//defined in shader
						pointSizeMin : 0,
						pointSizeMax : 5
					});
					/**/
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
								pointSizeMin : 0,
								pointSizeMax : 5
							});
						}
					});
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
