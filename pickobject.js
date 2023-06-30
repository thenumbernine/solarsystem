/*
in absence of coroutines i'm going to make a search object + callback
for planets, stars, etc this will resolve the callback in one frame
for galaxies it'll take ... 20 frames
this is really the job of a separate thread, but I haven't taken the time to learn how to navigate the excessive BS wrapped around multithreading javascript

here's my thought: if the galaxies can't be selected / seen then resolve the search in one frame
otherwise wait til the search is done
*/
import {mat4} from '/js/gl-matrix-3.4.1/index.js';
import {assert, assertExists} from '/js/util.js';
import {cfg} from './globals.js';
import {ui} from './ui.js';
import {Planet} from './planet.js';

class PickObject {
	constructor() {
		const gl = ui.gl;
		const ModifiedDepthShaderProgram = ui.ModifiedDepthShaderProgram;
		this.fboTexWidth = 32;
		this.fboTexHeight = 32;
		this.fboTex = new glutil.Texture2D({
			internalFormat : gl.RGBA32F,
			format : gl.RGBA,
			type : gl.FLOAT,
			width : this.fboTexWidth,
			height : this.fboTexHeight,
			magFilter : gl.NEAREST,
			minFilter : gl.NEAREST,
			wrap : {
				s : gl.CLAMP_TO_EDGE,
				t : gl.CLAMP_TO_EDGE,
			},
		});

		//do we need our own fbo?  this one needs depth info, so maybe
		this.fbo = new glutil.Framebuffer({
			width : this.fboTexWidth,
			height : this.fboTexHeight,
			useDepth : true,
		});

		this.fbo.bind();
		this.fbo.setColorAttachment(this.fboTex);
		this.fbo.check();
		this.fbo.unbind();

		//this is the shader to use with point clouds
		//point sets pass 'vertexID' as a sequential list of numbers
		this.pickPointShader = new ModifiedDepthShaderProgram({
			vertexCode : `
in vec3 vertex;

//there is no gl_VertexID in webgl ... so I'll just have to allocate 500k of sequential values
//ch0 holds the first 11 bytes, ch1 holds the next 11
in float vertexIDCh0;
in float vertexIDCh1;

uniform mat4 mvMat, projMat;
uniform float pointSize;
uniform float pointSizeMin;
uniform float pointSizeMax;
uniform bool pointSizeScaleWithDist;

out vec3 vertexIDv;

void main() {
	vertexIDv = vec3(
		mod(vertexIDCh0, 256.), //first 8 bits of ch0
		mod(floor(vertexIDCh0 / 256.), 8.) + 8. * mod(vertexIDCh1, 32.),	//next 3 of ch0 + first 5 of ch1
		floor(vertexIDCh1 / 32.));	//last 6 of ch1

	gl_Position = projMat * (mvMat * flatEarthXForm(vec4(vertex, 1.)));

	gl_PointSize = pointSize;
	if (pointSizeScaleWithDist) gl_PointSize /= gl_Position.w;

	//if a point is too small then discard it by throwing it offscreen
	//if (gl_PointSize < .5) gl_Position = vec4(100., 100., 100., -2.);

	gl_PointSize = clamp(gl_PointSize, pointSizeMin, pointSizeMax);

	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode : `
uniform vec3 startID;
in vec3 vertexIDv;
out vec4 fragColor;
void main() {
	vec4 v = vec4(startID + vertexIDv, 0.);
	vec3 carry = floor(v.xyz / 256.);
	v.xyz = v.xyz - 256. * carry;
	v.yzw = v.yzw + carry;
	fragColor = vec4(v.xyz / 256., 1.);
}
`,
			uniforms : {
				startID : [0,0,0],
				pointSize : 1,
				pointSizeMin : 0,
				pointSizeMax : 1,
				pointSizeScaleWithDist : true
			}
		});

		this.pickPlanetShader = new ModifiedDepthShaderProgram({
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
+ cfg.quatRotateCode
+ `
void main() {
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	gl_Position = projMat * (mvMat * flatEarthXForm(vec4(vtx3, 1.)));
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode : `
uniform vec3 id;
out vec4 fragColor;
void main() {
	fragColor = vec4(id / 256., 1.);
}
`,
			uniforms : {
				id : [0,0,0]
			}
		});

		this.pickProjMat = mat4.create();

		let maxPickID = 600000;	//max point size
		//http://stackoverflow.com/questions/27874983/webgl-how-to-use-integer-attributes-in-glsl/27884245#27884245
		let bitsPerChannel = 11;
		let vertexIDCh0 = new Float32Array(maxPickID);
		let vertexIDCh1 = new Float32Array(maxPickID);
		//as big as our largest point buffer
		for (let i = 0; i < maxPickID; ++i) {
			let mask = (1 << (bitsPerChannel - 1)) - 1;
			vertexIDCh0[i] = i & mask;
			vertexIDCh1[i] = (i >> bitsPerChannel) & mask;
		}
		this.vertexIDCh0Buffer = new glutil.ArrayBuffer({data : vertexIDCh0, dim : 1});
		this.vertexIDCh1Buffer = new glutil.ArrayBuffer({data : vertexIDCh1, dim : 1});

		//because i'm lazy ...
		//and it's more flexible to non-LInf metrics
		this.pixelOrder = [];
		for (let j = 0; j < this.fboTexHeight; ++j) {
			let y = j - this.fboTexHeight/2;
			for (let i = 0; i < this.fboTexWidth; ++i) {
				let x = i - this.fboTexWidth/2;
				this.pixelOrder.push([i,j,x*x + y*y]);
			}
		}
		this.pixelOrder.sort(function(a,b) {
			return a[2] - b[2];
		});

		this.callbacks = [];
	}

	pick(doChoose, skipProjection) {
		const gl = ui.gl;
		const canvas = ui.canvas;
		let viewport = gl.getParameter(gl.VIEWPORT);

		//pick window size
		let sizeX = this.fboTexWidth;
		let sizeY = this.fboTexHeight;
		let x = cfg.mouse.lastX;
		let y = canvas.height - cfg.mouse.lastY - 1;

		//mesa3d gluPickMatrix code: https://www.opengl.org/discussion_boards/showthread.php/184308-gluPickMatrix-Implementation
		//does glmatrix apply matrix operations lhs or rhs?  rhs I hope ..
		mat4.identity(this.pickProjMat);
if (!skipProjection) {
		mat4.translate(this.pickProjMat, this.pickProjMat, [(canvas.width - 2 * x) / sizeX, (canvas.height - 2 * y) / sizeY, 0]);
		mat4.scale(this.pickProjMat, this.pickProjMat, [canvas.width / sizeX, canvas.height / sizeY, 1]);
}
		mat4.multiply(this.pickProjMat, this.pickProjMat, glutil.scene.projMat);

		this.startPickID = 0xbf0000;
		this.pickID = this.startPickID;
		this.callbacks.length = 0;

		let foundIndex = 0;
if (!skipProjection) {
		gl.viewport(0, 0, sizeX, sizeY);
}
		let thiz = this;
		let fboCallback = function() {
			gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

			// run the render loop
			// set the scene to 'picking'
			// change the projection matrix to be a pick-matrix
			cfg.drawScene(true);

if (!skipProjection) {
			let pixels = new Float32Array(sizeX * sizeY * 4);
			gl.readPixels(0, 0, sizeX, sizeY, gl.RGBA, gl.FLOAT, pixels);
			for (let e = 0; e < thiz.pixelOrder.length; ++e) {
				let p = thiz.pixelOrder[e];
				let i = p[0];
				let j = p[1];
				let r = 256 * pixels[0 + 4 * (i + sizeX * j)];
				let g = 256 * pixels[1 + 4 * (i + sizeX * j)];
				let b = 256 * pixels[2 + 4 * (i + sizeX * j)];
				foundIndex = r | (g << 8) | (b << 16);
				if (foundIndex) break;
			}
}
		};
if (skipProjection) {
	fboCallback();
} else {
	this.fbo.draw({callback : fboCallback});
}

		let body = undefined;
		for (let i = 0; i < this.callbacks.length; ++i) {
			let cb = this.callbacks[i];
			if (foundIndex >= cb.start && foundIndex < cb.end) {
				if (cb.callbackObj instanceof Planet) {
					body = cb.callbackObj;
				} else {
					body = cb.callbackObj(foundIndex - cb.start);
				}
				break;
			}
		}

		cfg.mouseOverTarget = undefined;
		if (body !== undefined) {
			cfg.mouseOverTarget = body;
			if (body !== cfg.orbitTarget && doChoose) {
				setOrbitTarget(body);
				refreshMeasureText();
			}
		}

		gl.viewport.apply(gl, viewport);
	}

	registerCallback(callbackObj, count) {
		let i = this.pickID;
		this.callbacks.push({
			start : this.pickID,
			end : this.pickID + count,
			callbackObj : callbackObj
		});
		this.pickID += count;
		return [i & 255, (i >> 8) & 255, (i >> 16) & 255];
	}

	drawPoints(args) {
		const glutil = ui.glutil;
		let sceneObj = assertExists(args, 'sceneObj');
		let callbackObj = assertExists(args, 'targetCallback');
		let pointSize = assertExists(args, 'pointSize');
		let pointSizeMin = args.pointSizeMin !== undefined ? args.pointSizeMin : -Infinity;
		let pointSizeMax = args.pointSizeMax !== undefined ? args.pointSizeMax : Infinity;
		let pointSizeScaleWithDist = !!args.pointSizeScaleWithDist;
		let vertexAttr = sceneObj.attrs.vertex;
		let vertexBuffer = undefined;
		if (vertexAttr instanceof glutil.Attribute) {
			vertexBuffer = vertexAttr.buffer;
		} else if (vertexAttr instanceof glutil.ArrayBuffer) {
			vertexBuffer = vertexAttr;
		}

		let count = vertexBuffer.count || (vertexBuffer.data.length / vertexBuffer.dim);
		this.pickPointShader.use();
		this.pickPointShader.setAttrs({
			vertex : assert(sceneObj.attrs.vertex),
			vertexIDCh0 : this.vertexIDCh0Buffer,
			vertexIDCh1 : this.vertexIDCh1Buffer
		});
		sceneObj.setupMatrices();
		this.pickPointShader.setUniforms({
			projMat : this.pickProjMat,
			mvMat : sceneObj.uniforms.mvMat,
			startID : this.registerCallback(callbackObj, count),
			pointSize : pointSize,
			pointSizeMin : pointSizeMin,
			pointSizeMax : pointSizeMax,
			pointSizeScaleWithDist : pointSizeScaleWithDist
		});
		sceneObj.geometry.draw();
		this.pickPointShader.useNone();
	}

	drawPlanet(sceneObj, callbackObj) {
		this.pickPlanetShader.use();
		this.pickPlanetShader.setAttrs({
			vertex : sceneObj.attrs.vertex
		});
		sceneObj.setupMatrices();
		this.pickPlanetShader.setUniforms({
			projMat : this.pickProjMat,
			mvMat : sceneObj.uniforms.mvMat,
			pos : sceneObj.uniforms.pos,
			angle : sceneObj.uniforms.angle,
			equatorialRadius : sceneObj.uniforms.equatorialRadius,
			inverseFlattening : sceneObj.uniforms.inverseFlattening,
			scaleExaggeration : sceneObj.uniforms.scaleExaggeration,
			flatEarthCoeff : cfg.flatEarthCoeff,
			earthPos : cfg.flatEarthRelativeEarthPos,
			earthNorthDir : cfg.flatEarthRelativeEarthNorthDir,
			id : this.registerCallback(callbackObj, 1)
		});
		sceneObj.geometry.draw();
		this.pickPlanetShader.useNone();
	}
}

export {PickObject};
