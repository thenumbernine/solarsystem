import {mathClamp} from '/js/util.js';
import {cfg} from './globals.js';
import {ui} from './ui.js';

let skyCubeObj;
let skyCubeMaxBrightness = .3;
let skyCubeFadeOutStartDistInLyr = 100;
let skyCubeFadeOutEndDistInLyr = 200;

let skyCube = new function() {

	this.init = function() {
		const thiz = this;
		const glutil = ui.glutil;
		const gl = ui.gl;

		//looks like doing these in realtime will mean toning the detail down a bit ...
		//if the image is too big, how do we downsample the skymap without lagging the whole browser?  1) browsers start using LuaJIT (not going to happen, stupid JS) 2) provide pre-computed sampled down versions.

		let glMaxCubeMapTextureSize = gl.getParameter(gl.MAX_CUBE_MAP_TEXTURE_SIZE);

		this.texURLPrefixes = [
			'textures/sky-visible-cube-xp-',
			'textures/sky-visible-cube-xn-',
			'textures/sky-visible-cube-yp-',
			'textures/sky-visible-cube-yn-',
			'textures/sky-visible-cube-zp-',
			'textures/sky-visible-cube-zn-'
		];

		this.texURLs = this.texURLPrefixes.map(function(filename) {
			return filename + Math.min(1024, glMaxCubeMapTextureSize) +'.png';
		});

		new glutil.TextureCube({
			flipY : true,
			generateMipmap : true,
			magFilter : gl.LINEAR,
			minFilter : gl.LINEAR_MIPMAP_LINEAR,
			wrap : {
				s : gl.CLAMP_TO_EDGE,
				t : gl.CLAMP_TO_EDGE
			},
			urls : this.texURLs,
			done : function() {
				//init the "sky" cubemap (the galaxy background) once the texture for it loads
				thiz.texLoaded(this);
			}
		});
	};

	this.texLoaded = function(skyTex) {
		const glutil = ui.glutil;
		const gl = ui.gl;
		const ModifiedDepthShaderProgram = ui.ModifiedDepthShaderProgram;
		let cubeShader = new ModifiedDepthShaderProgram({
			vertexCode : `
in vec3 vertex;
out vec3 vertexv;
uniform mat4 projMat;

void main() {
	vertexv = vertex;
	gl_Position = projMat * vec4(vertex, 1.);
}
`,
			fragmentCode : `
in vec3 vertexv;
uniform samplerCube skyTex;
uniform float brightness;
//uniform vec4 angle;
uniform vec4 viewAngle;

` + cfg.coordinateSystemCode + cfg.quatRotateCode + `

out vec4 fragColor;
void main() {
	vec3 dir = vertexv;
	dir = quatRotate(viewAngle, dir);
	dir = eclipticalToGalactic * equatorialToEcliptical * dir;
	fragColor.rgb = brightness * texture(skyTex, dir).rgb;
	fragColor.a = 1.;
}
`,
			uniforms : {
				skyTex : 0
			}
		});

		let cubeVtxArray = new Float32Array(3*8);
		for (let i = 0; i < 8; i++) {
			cubeVtxArray[0+3*i] = glutil.view.zNear*10*(2*(i&1)-1);
			cubeVtxArray[1+3*i] = glutil.view.zNear*10*(2*((i>>1)&1)-1);
			cubeVtxArray[2+3*i] = glutil.view.zNear*10*(2*((i>>2)&1)-1);
		}


		let cubeIndexBuf = new glutil.ElementArrayBuffer({
			data : [
				5,7,3,3,1,5,		// <- each value has the x,y,z in the 0,1,2 bits (off = 0, on = 1)
				6,4,0,0,2,6,
				2,3,7,7,6,2,
				4,5,1,1,0,4,
				6,7,5,5,4,6,
				0,1,3,3,2,0
			]
		});

		//my galaxy texture is centered at x+ and lies in the xy plane
		skyCubeObj = new glutil.SceneObject({
			mode : gl.TRIANGLES,
			indexes : cubeIndexBuf,
			shader : cubeShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({
					data : cubeVtxArray
				})
			},
			uniforms : {
				viewAngle : glutil.view.angle
			},
			texs : [skyTex],
			parent : null
		});
	};

	this.draw = function(
		picking,
		distFromSolarSystemInLyr
	) {
		if (!skyCubeObj) return;
		if (picking) return;
		if (distFromSolarSystemInLyr >= skyCubeFadeOutEndDistInLyr) return;

		const gl = ui.gl;
		const brightness = skyCubeMaxBrightness * (1 - mathClamp((distFromSolarSystemInLyr - skyCubeFadeOutStartDistInLyr) / (skyCubeFadeOutEndDistInLyr - skyCubeFadeOutStartDistInLyr), 0, 1));

		gl.disable(gl.DEPTH_TEST);
		skyCubeObj.uniforms.brightness = brightness;
		skyCubeObj.draw();
		gl.enable(gl.DEPTH_TEST);
	};
};

export {skyCube};
