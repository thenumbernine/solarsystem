//tried an experiment of doing surface calculations on the GPU
//it ran a lot faster than doing them in CPU for JS ... but the floating point accuracy was too low to get any good results back, even with double precision functions
//I might try worker threads later...
var CALCULATE_TIDES_WITH_GPU = false;

var SHOW_ALL_SMALL_BODIES_AT_ONCE = true;
var SHOW_ALL_SMALL_BODIES_WITH_DENSITY = false;

//TODO put this in js/gl-util-kernel.js

var glutil;
GLUtil.prototype.oninit.push(function() {
	glutil = this;
	glutil.KernelShader = makeClass({
		super : glutil.ShaderProgram,
		
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
		init : function(args) {
			var varyingVar = args.varying !== undefined ? args.varying : 'pos';
			
			var varyingCodePrefix = 'varying vec2 '+varyingVar+';\n';

			var fragmentCodePrefix = '';
			var uniforms = {};
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
				for (var i = 0; i < args.texs.length; ++i) {
					var v = args.texs[i];
					var name, vartype;
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
					code : 
						glutil.vertexPrecision + 
						varyingCodePrefix +
						mlstr(function(){/*
attribute vec2 vertex;
attribute vec2 texCoord;
void main() {
	*/}) + varyingVar + mlstr(function(){/* = texCoord; 
	gl_Position = vec4(vertex, 0., 1.);
}
*/})
				});	
			}

			args.vertexShader = glutil.KernelShader.prototype.kernelVertexShader;
			args.fragmentCode = glutil.fragmentPrecision + varyingCodePrefix + fragmentCodePrefix + args.code;
			delete args.code;
			args.uniforms = uniforms;	
			glutil.KernelShader.super.call(this, args);
		}
	});
});

function setCSSRotation(obj, degrees) {
	obj.css('-webkit-transform','rotate('+degrees+'deg)');
	obj.css('-moz-transform','rotate('+degrees+'deg)');
	obj.css('transform','rotate('+degrees+'deg)');
}
