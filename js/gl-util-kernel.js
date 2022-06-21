if (!GLUtil) throw "require gl-util.js before gl-util-kernel.js";

GLUtil.prototype.oninit.push(function(gl) {
	var glutil = this;

	if (!this.unitQuad) throw "require gl-util-unitquad.js before gl-util-kernel.js";

	var KernelShader = makeClass({
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
			
			
			var varyingCode = [
'varying vec2 '+varyingVar+';'].join('\n');
			var vertexCode = [
varyingCode,
'attribute vec2 vertex;',
'void main() {',
'	'+varyingVar+' = vertex.xy;',
'	gl_Position = vec4(vertex.xy * 2. - 1., 0., 1.);',
'}'].join('\n');
			var precision = 'mediump';
			if (args.precision !== undefined) precision = args.precision;
			var fragmentCodePrefix = 'precision '+precision+' float;\n' + varyingCode;
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
			var code;
			if (args.code !== undefined) code = args.code;
			if (args.codeID !== undefined) code = $('#'+args.codeID).text();
			KernelShader.super.call(this, {
				vertexCodeID : args.vertexCodeID,
				vertexCode : args.vertexCode !== undefined ? args.vertexCode : vertexCode,
				vertexPrecision : args.vertexPrecision,
				fragmentCode : fragmentCodePrefix + code,
				fragmentPrecision : args.fragmentPrecision,
				uniforms : uniforms
			});
		}
	});
	glutil.KernelShader = KernelShader;
});
