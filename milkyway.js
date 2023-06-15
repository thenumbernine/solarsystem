let milkyWay = new function(){
	this.fadeMinDistInLyr = 50;
	this.fadeMaxDistInLyr = 1000;
	
	this.init = function(
	) {
		let thiz = this;
		let img = new Image();
		img.onload = function() {
			thiz.sceneObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				attrs : {
					vertex : new glutil.ArrayBuffer({dim : 2, data : [-.5, -.5, .5, -.5, -.5, .5, .5, .5]}),
					texCoord : new glutil.ArrayBuffer({dim : 2, data : [0, 0, 1, 0, 0, 1, 1, 1]})
				},
				shader : new ModifiedDepthShaderProgram({
					vertexCode :
`#define SCALE_OF_MILKY_WAY_SPRITE ` + floatToGLSL(160000 * metersPerUnits.lyr / interGalacticRenderScale) + `
in vec2 vertex;
in vec2 texCoord;
uniform mat4 mvMat;
uniform mat4 projMat;
out vec2 texCoordv;

` + coordinateSystemCode + `

mat3 transpose(mat3 m) {
	return mat3(m[0][0], m[1][0], m[2][0],
				m[0][1], m[1][1], m[2][1],
				m[0][2], m[1][2], m[2][2]);
}

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
		};
		img.onerror = function() {
			console.log('failed to find texture for milky way');
		};
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
			
		if (orbitTarget !== undefined && orbitTarget.pos !== undefined) {
			this.sceneObj.pos[0] = galaxyCenterInEquatorialCoordsInMpc[0] * (metersPerUnits.Mpc / interGalacticRenderScale) - orbitTarget.pos[0] / interGalacticRenderScale;
			this.sceneObj.pos[1] = galaxyCenterInEquatorialCoordsInMpc[1] * (metersPerUnits.Mpc / interGalacticRenderScale) - orbitTarget.pos[1] / interGalacticRenderScale;
			this.sceneObj.pos[2] = galaxyCenterInEquatorialCoordsInMpc[2] * (metersPerUnits.Mpc / interGalacticRenderScale) - orbitTarget.pos[2] / interGalacticRenderScale;
		}
	
		//apply milky way local transforms to mpc mv mat
		gl.disable(gl.CULL_FACE);
		gl.depthMask(false);
		this.sceneObj.uniforms.fadeInAlpha = Math.clamp((distFromSolarSystemInLyr - this.fadeMinDistInLyr) / (this.fadeMaxDistInLyr - this.fadeMinDistInLyr), 0, 1);
		this.sceneObj.draw();
		gl.enable(gl.CULL_FACE);
		gl.depthMask(true);
	};
};
