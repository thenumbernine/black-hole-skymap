GeodesicFBORenderer = makeClass({
	//1024x1024 is just too big for my current card to handle
	lightTexWidth : 512,
	lightTexHeight : 512,
	lightPosVelChannels : [
		{pos_or_vel:'pos'},
		{pos_or_vel:'vel'}
	],
	
	init : function(glutil) {
		this.glutil = glutil;
		if (!this.glutil.context.getExtension('OES_texture_float')) {
			throw 'This requires OES_texture_float';
		}
	},
	
	initScene : function(skyTex) {
		var thiz = this;
		$.each(thiz.lightPosVelChannels, function(_,channel) {
			//working on how to organize this
			channel.uniforms = [
				'blackHoleMass', 
				'warpBubbleThickness', 
				'warpBubbleRadius', 
				'warpBubbleVelocity', 
				'objectDist', 
				'deltaLambda'
			];
			
			channel.texs = [];
			channel.fbos = [];
			for (var history = 0; history < 2; ++history) {
				var texs = [];
				channel.texs[history] = texs;
				var fbos = [];
				channel.fbos[history] = fbos;
				for (var side = 0; side < 6; ++side) {
					var tex = new thiz.glutil.Texture2D({
						internalFormat : gl.RGBA,
						format : gl.RGBA,
						type : gl.FLOAT,
						width : thiz.lightTexWidth,
						height : thiz.lightTexHeight,
						magFilter : gl.NEAREST,
						minFilter : gl.NEAREST,
						wrap : {
							s : gl.CLAMP_TO_EDGE,
							t : gl.CLAMP_TO_EDGE
						}
					});
					texs[side] = tex;
					
					var fbo = new thiz.glutil.Framebuffer();
					gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.obj);
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex.obj, 0);
					gl.bindFramebuffer(gl.FRAMEBUFFER, null);
					fbos[side] = fbo;
				}
			}
			
			channel.shaders = {};
					
			$.each(objectTypes, function(_,objType) {
				var pos_or_vel = channel.pos_or_vel;
				
				channel.shaders[objType] = {};
				
				var shaderTypes = ['reset', 'iterate'];
				
				$.each(shaderTypes, function(_,shaderType) {

					var vsh = [];
					var fsh = [shaderCommonCode];
					if (shaderType == 'reset') {
						
						//new
						vsh.push(mlstr(function(){/*
attribute vec2 vertex;
varying vec3 vtxv;
uniform mat3 rotation;
void main() {
	vtxv = rotation * vec3(vertex * 2. - 1., 1.);
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/}));
							fsh.push(mlstr(function(){/*
uniform float objectDist;
*/}));
						if (pos_or_vel == 'pos') {
							fsh.push(mlstr(function(){/*
varying vec3 vtxv;
void main() {
	gl_FragColor = vec4(0.);
	gl_FragColor.x -= objectDist;
}
*/}));
						} else if (pos_or_vel == 'vel') {
							if (objType == 'Black Hole') {
								fsh.push(mlstr(function(){/*
uniform float blackHoleMass;

// g_tt = -(1 - m/(2r))^2 / (1 + m/(2r))^2
float g_tt(vec3 pos, vec3 vel) {
	float r = length(pos);
	float _m_2r = blackHoleMass / (2. * r);
	float _1_plus_m_2r = 1. + _m_2r;
	float _1_minus_m_2r = 1. - _m_2r;
	float ratio = _1_minus_m_2r / _1_plus_m_2r;
	return -ratio * ratio; 
}

// g_ti = 0
vec3 g_ti(vec3 pos, vec3 vel) { 
	return vec3(0., 0., 0.); 
}

// g_ij = delta_ij (1 + m/(2r))^4 
mat3 g_ij(vec3 pos, vec3 vel) {
	float r = length(pos);
	float _m_2r = blackHoleMass / (2. * r);
	float _1_plus_m_2r = 1. + _m_2r;
	float _1_plus_m_2r_sq = _1_plus_m_2r * _1_plus_m_2r;
	float _1_plus_m_2r_toTheFourth = _1_plus_m_2r_sq * _1_plus_m_2r_sq;
	return mat3(_1_plus_m_2r_toTheFourth);
}

*/}));						
							} else if (objType == 'Alcubierre Warp Drive Bubble') {
								fsh.push(mlstr(function(){/*
uniform float warpBubbleThickness;
uniform float warpBubbleVelocity;
uniform float warpBubbleRadius;

float Alcubierre_f(vec3 pos) {
	float rs = length(pos);
	float sigmaFront = warpBubbleThickness * (rs + warpBubbleRadius);
	float sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	float sigmaBack = warpBubbleThickness * (rs - warpBubbleRadius);
	float tanhSigmaCenter = tanh(sigmaCenter);
	float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanhSigmaCenter);
	return f;
}

float g_tt(vec3 pos, vec3 vel) {
	float f = Alcubierre_f(pos);
	float vf = f * warpBubbleVelocity;
	float vf2 = vf * vf;
	return -(1. - vf2); 
}

vec3 g_ti(vec3 pos, vec3 vel) {
	float f = Alcubierre_f(pos);
	float vf = f * warpBubbleVelocity;
	float g_tx = -vf;
	return vec3(g_tx, 0., 0.);
}

mat3 g_ij(vec3 pos, vec3 vel) {
	return mat3(1.);
}

*/}));
							}	

							fsh.push(mlstr(function(){/*

//when initializing our metric:
//g_ab v^a v^b = 0 for our metric g
//g_tt (v^t)^2 + 2 g_ti v^t v^i + g_ij v^i v^j = 0
//(v^t)^2 + (2 v^i g_ti / g_tt) v^t + g_ij v^i v^j / g_tt = 0
//v_t = v^i g_ti / (-g_tt) +- sqrt( (v^i g_ti / g_tt)^2 + g_ij v^i v^j / (-g_tt) )
// I'll go with the plus
float reset_w(vec3 pos, vec3 vel) {
	float _g_tt = g_tt(pos, vel);
	float negInv_g_tt = -1. / _g_tt;
	vec3 _g_ti = g_ti(pos, vel);
	float v_t = dot(vel, _g_ti);	//v_t = g_ti v^i
	float a = v_t * negInv_g_tt;
	float bsq = a * a + dot(vel, g_ij(pos, vel) * vel) * negInv_g_tt;
	if (bsq < 0.) bsq = 0.;
	return a + sqrt(bsq);
}

varying vec3 vtxv;
void main() {
	vec3 pos = vec3(-objectDist, 0., 0.);
	vec3 vel = normalize(vtxv);
	gl_FragColor.xyz = vel;
	gl_FragColor.w = reset_w(pos, vel);
}
*/}));						
						}
					} else if (shaderType == 'iterate') {
						vsh.push(mlstr(function(){/*
varying vec2 uv;
varying vec3 vtxv;
attribute vec2 vertex;
uniform mat3 rotation;
void main() {
	uv = vertex;
	vtxv = rotation * vec3(vertex * 2. - 1., 1.);
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/}));
						fsh.push(mlstr(function(){/*
varying vec2 uv;
varying vec3 vtxv;
uniform float deltaLambda;
uniform float objectDist;
uniform sampler2D lightPosTex;
uniform sampler2D lightVelTex;
*/}));
						if (pos_or_vel == 'vel') {
							if (objType == 'Black Hole') {
								fsh.push(mlstr(function(){/*
uniform float blackHoleMass;

vec4 accel(vec4 pos, vec4 vel) {
	float r = length(pos.xyz);
	float posDotVel = dot(pos.xyz, vel.xyz);
	float posDotVelSq = posDotVel * posDotVel;
	float velSq = dot(vel.xyz, vel.xyz);
	float R = 2. * blackHoleMass;
	float r2 = r * r;
	float r3 = r * r2;


	float _m_2r = blackHoleMass / (2. * r);
	float _1_plus_m_2r = 1. + _m_2r;
	float _1_minus_m_2r = 1. - _m_2r;
	
	float _1_plus_m_2r_sq = _1_plus_m_2r * _1_plus_m_2r;
	float _1_plus_m_2r_tothe3 = _1_plus_m_2r_sq * _1_plus_m_2r;
	float _1_plus_m_2r_tothe6 = _1_plus_m_2r_tothe3 * _1_plus_m_2r_tothe3;
	
	vec4 result;
	result.w = -blackHoleMass * 2. * vel.w * posDotVel / (_1_plus_m_2r * _1_minus_m_2r * r3);
	result.xyz = -blackHoleMass * (
		vel.w * vel.w * pos.xyz * _1_minus_m_2r / _1_plus_m_2r_tothe6 
		- 2. * vel.xyz * posDotVel
		+ velSq * pos.xyz 
	) / (_1_plus_m_2r * r3);
	return result;
}
*/}));
							} else if (objType == 'Alcubierre Warp Drive Bubble') {
								fsh.push(mlstr(function(){/*
uniform float warpBubbleThickness;
uniform float warpBubbleRadius;
uniform float warpBubbleVelocity;

vec4 accel(vec4 pos, vec4 vel) {
	float r = length(pos.xyz);
	float sigmaFront = warpBubbleThickness * (r + warpBubbleRadius);
	float sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	float sigmaBack = warpBubbleThickness * (r - warpBubbleRadius);
	float tanhSigmaCenter = tanh(sigmaCenter);
	float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanhSigmaCenter);
	float sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
	float dfScalar = sechDiff / (2. * r * tanhSigmaCenter);
	vec4 df;
	df.xyz = warpBubbleThickness * pos.xyz * dfScalar;
	df.w = -warpBubbleVelocity * warpBubbleThickness * pos.x * dfScalar;

	float u = f * warpBubbleVelocity;
	vec4 du = df * warpBubbleVelocity;
	
	vec4 result;
	result.w = 
		- du.x * vel.x * vel.x
		- du.y * vel.x * vel.y
		- du.z * vel.x * vel.z
		
		+ 2. * du.x * vel.w * vel.x * u
		- du.x * vel.w * vel.w * u * u
		+ du.y * vel.w * vel.y * u
		+ du.z * vel.w * vel.z * u
	;
	
	result.x = 
		du.w * vel.w * vel.w
		+ du.y * vel.w * vel.y
		+ du.y * vel.w * vel.y * u * u
		- du.y * vel.x * vel.y * u
		+ du.z * vel.w * vel.z
		+ du.z * vel.w * vel.z * u * u
		- du.z * vel.x * vel.z * u
		+ du.x * vel.w * vel.w * u
		- du.x * vel.w * vel.w * u * u * u
		+ 2. * du.x * vel.w * vel.x * u * u
		- du.x * vel.x * vel.x * u
	;
	result.y = -du.y * vel.w * (vel.x - vel.w * u);
	result.z = -du.z * vel.w * (vel.x - vel.w * u);
	return result;
}
*/}));
							}
						}
					
						fsh.push(mlstr(function(){/*
void main() {
	vec4 pos = texture2D(lightPosTex, uv);
	vec4 vel = texture2D(lightVelTex, uv);

*/}));
						if (pos_or_vel == 'pos') {
							fsh.push(mlstr(function(){/*
	gl_FragColor = pos + deltaLambda * vel;
}
*/}));
						} else if (pos_or_vel == 'vel') {
							fsh.push(mlstr(function(){/*
	gl_FragColor = vel + deltaLambda * accel(pos, vel);
}
*/}));
						}
					}

					var args = {
						vertexCode : vsh.join('\n'),
						fragmentCode : fsh.join('\n'),
					};

					args.context = thiz.glutil.context;
					args.uniforms = {
						lightPosTex : 0,
						lightVelTex : 1
					};
					args.vertexPrecision = 'best';
					args.fragmentPrecision = 'best';
					channel.shaders[objType][shaderType] = new thiz.glutil.ShaderProgram(args);
				});
			});
		});

		var cubeShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			vertexCode : mlstr(function(){/*
attribute vec2 vertex;
varying vec2 uv;
uniform mat4 projMat;
uniform vec4 angle;
const mat3 viewMatrix = mat3(
	0., 0., -1.,
	1., 0., 0., 
	0., -1., 0.);
vec3 quatRotate(vec4 q, vec3 v) { 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	uv = vertex;
	vec3 vtx3 = vec3(vertex * 2. - 1., 1.);
	vtx3 = quatRotate(vec4(angle.xyz, -angle.w), vtx3);
	vtx3 = viewMatrix * vtx3;
	vec4 vtx4 = vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
}
*/}),
			fragmentPrecision : 'best',
			fragmentCode : shaderCommonCode + mlstr(function(){/*
varying vec2 uv;
uniform samplerCube skyTex;
uniform sampler2D lightVelTex;
uniform vec4 viewAngle;
const mat3 viewMatrixInv = mat3(
	0., 1., 0.,
	0., 0., -1.,
	-1., 0., 0.);
vec3 quatRotate(vec4 q, vec3 v) { 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	vec3 dir = texture2D(lightVelTex, uv).xyz;
	dir = viewMatrixInv * viewMatrixInv * dir;
	dir = quatRotate(viewAngle, dir);
	gl_FragColor.xyz = textureCube(skyTex, dir).xyz;
	gl_FragColor.w = 1.; 
}*/}),
			uniforms : {
				skyTex : 0,
				lightVelTex : 1
			}
		});

		//scene graph, 6 quads oriented in a cube
		//I would use a cubemap but the Mali-400 doesn't seem to want to use them as 2D FBO targets ...
		this.cubeSides = [];
		for (var side = 0; side < 6; ++side) {
			this.cubeSides[side] = new this.glutil.SceneObject({
				//geometry : glutil.unitQuad.geometry,
				mode : gl.TRIANGLE_STRIP,
				attrs : {
					vertex : glutil.unitQuadVertexBuffer
				},
				uniforms : {
					viewAngle : this.glutil.view.angle,
					angle : angleForSide[side]
				},
				shader : cubeShader,
				texs : [skyTex, this.lightPosVelChannels[1].texs[0][side]],
			});
		}


	},

	resetField : function() {
		gl.viewport(0, 0, this.lightTexWidth, this.lightTexHeight);
		$.each(this.lightPosVelChannels, function(_,channel) {
			for (var side = 0; side < 6; ++side) {
				var shader = channel.shaders[objectType].reset;
				
				var uniforms = {};
				if (channel.uniforms !== undefined) {
					for (var i = 0; i < channel.uniforms.length; ++i) {
						var uniformName = channel.uniforms[i];
						uniforms[uniformName] = window[uniformName];
					}
				}
				var rotationMat3 = mat3.create();
				mat3.fromQuat(rotationMat3, angleForSide[side]);
				uniforms.rotation = rotationMat3;
				
				var fbo = channel.fbos[0][side];
				fbo.draw({
					callback:function(){
						glutil.unitQuad.draw({
							shader:shader,
							uniforms:uniforms
						});
					}
				});
			}
		});
		gl.viewport(0, 0, this.glutil.canvas.width, this.glutil.canvas.height);
	},
	
	updateLightPosTex : function() {	
		var thiz = this;
		gl.viewport(0, 0, this.lightTexWidth, this.lightTexHeight);
		$.each(this.lightPosVelChannels, function(_,channel) {
			for (var side = 0; side < 6; ++side) {
				var shader = channel.shaders[objectType].iterate;
				
				var uniforms = {};
				if (channel.uniforms !== undefined) {
					for (var i = 0; i < channel.uniforms.length; ++i) {
						var uniformName = channel.uniforms[i];
						uniforms[uniformName] = window[uniformName];
					}
				}
				var rotationMat3 = mat3.create();
				mat3.fromQuat(rotationMat3, angleForSide[side]);
				uniforms.rotation = rotationMat3;
				uniforms.lightPosTex = 0;
				uniforms.lightVelTex = 1;
				
				var fbo = channel.fbos[1][side];
				fbo.draw({
					callback:function(){
						glutil.unitQuad.draw({
							shader:shader,
							uniforms:uniforms,
							texs:[
								thiz.lightPosVelChannels[0].texs[0][side],
								thiz.lightPosVelChannels[1].texs[0][side]
							]
						});
					}
				});
			}
		});
		gl.viewport(0, 0, this.glutil.canvas.width, this.glutil.canvas.height);
		$.each(this.lightPosVelChannels, function(_,channel) {
			var tmp;
			tmp = channel.texs[0];
			channel.texs[0] = channel.texs[1];
			channel.texs[1] = tmp;
			tmp = channel.fbos[0];
			channel.fbos[0] = channel.fbos[1];
			channel.fbos[1] = tmp;
		});
	},

	update : function() {
		this.updateLightPosTex();
		
		//turn on magnification filter
		for (var side = 0; side < 6; ++side) {
			var tex = this.lightPosVelChannels[1].texs[0][side].obj;
			gl.bindTexture(gl.TEXTURE_2D, tex);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
			//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
			//gl.generateMipmap(gl.TEXTURE_2D);
			gl.bindTexture(gl.TEXTURE_2D, null);
		}
		
		for (var side = 0; side < 6; ++side) {
			//texs[0] is skyTex
			this.cubeSides[side].texs[1] = this.lightPosVelChannels[1].texs[0][side];
		}

		this.glutil.draw();	
		
		//turn off magnification filter
		for (var side = 0; side < 6; ++side) {
			var tex = this.lightPosVelChannels[1].texs[0][side].obj;
			gl.bindTexture(gl.TEXTURE_2D, tex);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
			//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
			gl.bindTexture(gl.TEXTURE_2D, null);
		}	
	}
});

