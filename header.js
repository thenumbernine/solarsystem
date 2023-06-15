//TODO put this in js/gl-util-kernel.js

let glutil;
GLUtil.prototype.oninit.push(function() {
	glutil = this;
	class KernelShader extends glutil.ShaderProgram {
		/*
		args:
			code : the fragment code
			varying : name of varying variable.  default 'pos'
			vertexCode : (optional) vertex code
			uniforms : { uniformName : uniformType }
					: { uniformName : [uniformType, initialValue] }
			texs : [texName]
				: [{texName : texType}]
			precision : (optional) mediump (default), highp, etc
		*/
		constructor(args) {
			let varyingVar = args.varying !== undefined ? args.varying : 'pos';
			
			let varyingCodePrefix = 'varying vec2 '+varyingVar+';\n';

			let fragmentCodePrefix = '';
			let uniforms = {};
			if (args.uniforms !== undefined) {
				$.each(args.uniforms, function(uniformName, uniformType) {
					if ($.isArray(uniformType)) {
						//save initial value
						uniforms[uniformName] = uniformType[1];
						uniformType = uniformType[0];
					}
					fragmentCodePrefix += 'uniform '+uniformType+' '+uniformName+';\n';
				});
			}
			if (args.texs !== undefined) {
				for (let i = 0; i < args.texs.length; ++i) {
					let v = args.texs[i];
					let name, vartype;
					if (typeof(v) == 'string') {
						name = v;
						vartype = 'sampler2D';
					} else {
						name = v[0];
						vartype = v[1];
					}
					fragmentCodePrefix += 'uniform '+vartype+' '+name+';\n';
					uniforms[name] = i;
				}
			}

			if (!glutil.KernelShader.prototype.kernelVertexShader) {
				glutil.KernelShader.prototype.kernelVertexShader = new glutil.VertexShader({
					code : `#version 300 es
precision `+glutil.vertexBestPrec+` float;
` + varyingCodePrefix.replace(/varying/g, 'out') + `
in vec2 vertex;
in vec2 texCoord;
void main() {
	` + varyingVar + ` = texCoord; 
	gl_Position = vec4(vertex, 0., 1.);
}
`
				});	
			}

			args.vertexShader = glutil.KernelShader.prototype.kernelVertexShader;
			args.fragmentCode = 
varyingCodePrefix.replace(/varying/g, 'in') 
+ fragmentCodePrefix 
+ args.code;
			delete args.code;
			args.uniforms = uniforms;	
			glutil.KernelShader.super.call(this, args);
		}
	}
	glutil.KernelShader = KernelShader;
});

function setCSSRotation(obj, degrees) {
	obj.css('-webkit-transform','rotate('+degrees+'deg)');
	obj.css('-moz-transform','rotate('+degrees+'deg)');
	obj.css('transform','rotate('+degrees+'deg)');
}
