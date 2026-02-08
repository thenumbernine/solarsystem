//TODO put this in js/gl-util-kernel.js

let glutil;
GLUtil.prototype.oninit.push(function() {
	glutil = this;
	class KernelShader extends glutil.Program {
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
			const varyingVar = args.varying !== undefined ? args.varying : 'pos';
			const varyingCodePrefix = 'varying vec2 '+varyingVar+';\n';
throw "TODO use unitQuad, so varyingVar = vertex.xy; gl_Position = vec4(vertex.xy * 2. - 1.), and pos varies [0,1]";
			const vertexCode =
varyingCodePrefix.replace(/varying/g, 'out')
+ `
in vec2 vertex;
in vec2 texCoord;
void main() {
	` + varyingVar + ` = texCoord;
	gl_Position = vec4(vertex, 0., 1.);
}
`;


			let fragmentCodePrefix = '';
			let uniforms = {};
			if (args.uniforms !== undefined) {
				Object.entries(args.uniforms).forEach(entry => {
					let [uniformName, uniformType] = entry;
					if (Array.isArray(uniformType)) {
						//save initial value
						uniforms[uniformName] = uniformType[1];
						uniformType = uniformType[0];
					}
					fragmentCodePrefix += 'uniform '+uniformType+' '+uniformName+';\n';
				});
			}
			if (args.texs !== undefined) {
				args.texs.forEach((v, i) => {
					let name, vartype;
					if (typeof(v) == 'string') {
						[name, vartype] = [v, 'sampler2D'];
					} else {
						[name, vartype] = v;
					}
					fragmentCodePrefix += 'uniform '+vartype+' '+name+';\n';
					uniforms[name] = i;
				});
			}

			super({
				vertexPrecision : args.precision,
				vertexCode : args.vertexCode !== undefined ? args.vertexCode : vertexCode,
				fragmentPrecision : args.precision,
				fragmentCode :
					varyingCodePrefix.replace(/varying/g, 'in')
					+ fragmentCodePrefix
					+ (args.code !== undefined ? args.code : ''),
				uniforms : uniforms,
			});
		}
	}
	glutil.KernelShader = KernelShader;
});
