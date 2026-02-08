import {mathClamp} from '/js/util.js';
import {cfg, floatToGLSL} from './globals.js';
import {ui} from './ui.js';
import {metersPerUnits} from './units.js';
import {galaxyCenterInEquatorialCoordsInMpc} from './vec.js';

let milkyWay = new function(){
	this.fadeMinDistInLyr = 50;
	this.fadeMaxDistInLyr = 1000;

	this.init = function(
	) {
		const gl = ui.gl;
		const ModifiedDepthShaderProgram = ui.ModifiedDepthShaderProgram;
		let thiz = this;
		let img = new Image();
		img.addEventListener('load', e => {
			thiz.sceneObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				attrs : {
					vertex : new glutil.ArrayBuffer({dim : 2, data : [-.5, -.5, .5, -.5, -.5, .5, .5, .5]}),
					texCoord : new glutil.ArrayBuffer({dim : 2, data : [0, 0, 1, 0, 0, 1, 1, 1]})
				},
				shader : new ModifiedDepthShaderProgram({
					vertexCode :
`#define SCALE_OF_MILKY_WAY_SPRITE ` + floatToGLSL(160000 * metersPerUnits.lyr / cfg.interGalacticRenderScale) + `
in vec2 vertex;
in vec2 texCoord;
uniform mat4 mvMat;
uniform mat4 projMat;
out vec2 texCoordv;

` + cfg.coordinateSystemCode + `

void main() {
	//the image is rotated 90 degrees ...
	texCoordv = vec2(-texCoord.y, texCoord.x);
	vec3 vtx3 = vec3(vertex.x, vertex.y, 0.);
    vec3 modelPos = transpose(eclipticalToGalactic * equatorialToEcliptical) * (SCALE_OF_MILKY_WAY_SPRITE * vtx3);
    gl_Position = projMat * (mvMat * flatEarthXForm(vec4(modelPos, 1.)));
	gl_Position.z = depthfunction(gl_Position);
}
`,
					fragmentCode : `
uniform sampler2D tex;
uniform float fadeInAlpha;
in vec2 texCoordv;
out vec4 fragColor;
void main() {
	fragColor = texture(tex, texCoordv);
	fragColor.a *= fadeInAlpha;
}
`,
				}),
				uniforms : {
					tex : 0
				},
				texs : [
					new glutil.Texture2D({
						flipY : true,
						data : img,
						minFilter : gl.LINEAR_MIPMAP_LINEAR,
						magFilter : gl.LINEAR,
						generateMipmap : true
					})
				],
				blend : [gl.SRC_ALPHA, gl.ONE],
				parent : null,
				pos : [0,0,0],
				static : false
			});
			console.log("created milky way")
		});
		img.addEventListener('error', e => {
			console.log('failed to find texture for milky way');
		});
		img.src = 'textures/milkyway.png';
	};

	this.draw = function(
		tanFovY,
		picking,
		distFromSolarSystemInLyr,
		distFromSolarSystemInMpc
	) {
		if (!this.sceneObj) return;
		if (picking) return;
		//draw milky way if we're far enough out
		if (distFromSolarSystemInLyr <= this.fadeMinDistInLyr) return;

		if (cfg.orbitTarget !== undefined && cfg.orbitTarget.pos !== undefined) {
			this.sceneObj.pos[0] = galaxyCenterInEquatorialCoordsInMpc[0] * (metersPerUnits.Mpc / cfg.interGalacticRenderScale) - cfg.orbitTarget.pos[0] / cfg.interGalacticRenderScale;
			this.sceneObj.pos[1] = galaxyCenterInEquatorialCoordsInMpc[1] * (metersPerUnits.Mpc / cfg.interGalacticRenderScale) - cfg.orbitTarget.pos[1] / cfg.interGalacticRenderScale;
			this.sceneObj.pos[2] = galaxyCenterInEquatorialCoordsInMpc[2] * (metersPerUnits.Mpc / cfg.interGalacticRenderScale) - cfg.orbitTarget.pos[2] / cfg.interGalacticRenderScale;
		}

		//apply milky way local transforms to mpc mv mat
		const gl = ui.gl;
		gl.disable(gl.CULL_FACE);
		gl.depthMask(false);
		this.sceneObj.uniforms.fadeInAlpha = mathClamp((distFromSolarSystemInLyr - this.fadeMinDistInLyr) / (this.fadeMaxDistInLyr - this.fadeMinDistInLyr), 0, 1);
		this.sceneObj.draw();
		gl.enable(gl.CULL_FACE);
		gl.depthMask(true);
	};
};

export {milkyWay};
