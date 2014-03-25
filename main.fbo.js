var canvas;
var gl;
var mouse;
var cubeSides;
var objectTypes = ['Black Hole', 'Alcubierre Warp Drive Bubble'];
var objectType = objectTypes[0];
var objectDist = 10;
var blackHoleMass = 1;
var warpBubbleThickness = 1;
var warpBubbleVelocity = .5;
var warpBubbleRadius = 2;
var deltaLambda = .1;	//ray forward iteration
//1024x1024 is just too big for my current card to handle
var lightTexWidth = 512;
var lightTexHeight = 512;

var ident4 = mat4.create();

/*
flags: flags used to find what shader to associate with this
	they are combined with each object type and put in 'shaders'
	from that, shaders[objectType] determines the shader to run
*/
var shaderTypes = ['reset', 'iterate'];
var lightPosVelChannels = [
	{flags:['pos']},
	{flags:['vel']}
];

var shaderCommonCode = mlstr(function(){/*
float tanh(float x) {
	float exp2x = exp(2. * x);
	return (exp2x - 1.) / (exp2x + 1.);
}

float sech(float x) {
	float expx = exp(x);
	return 2. * expx / (expx * expx + 1.);
}

float sechSq(float x) {
	float y = sech(x);
	return y * y;
}
*/});


function stupidPrint(s) {
	return;
	$.each(s.split('\n'), function(_,l) {
		console.log(l);
	});
}

var shaderScriptParts = [
	{
		flags : 'Black_Hole:Alcubierre_Warp_Drive_Bubble:reset:pos:vel:vsh',
		code : mlstr(function(){/*
attribute vec2 vertex;
varying vec3 vtxv;
uniform mat3 rotation;
void main() {
	vtxv = rotation * vec3(vertex * 2. - 1., 1.);
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/})
	},
	{
		flags : 'Black_Hole:reset:pos:fsh',
		code : mlstr(function(){/*
varying vec3 vtxv;
uniform float objectDist;
void main() {
	gl_FragColor = vec4(0.);
	gl_FragColor.x -= objectDist;
}
*/})
	},
	{
		flags : 'Alcubierre_Warp_Drive_Bubble:reset:pos:fsh',
		code : mlstr(function(){/*
varying vec3 vtxv;
void main() {
	gl_FragColor = vec4(0.);
}
*/})
	},
	{
		flags : 'Black_Hole:reset:vel:fsh',
		code : mlstr(function(){/*
uniform float objectDist;
uniform float blackHoleMass;
*/})
	},
	{
		flags : 'Alcubierre_Warp_Drive_Bubble:reset:vel:fsh',
		code : mlstr(function(){/*
uniform float objectDist;
uniform float warpBubbleThickness;
uniform float warpBubbleVelocity;
uniform float warpBubbleRadius;
*/})
	},
	{
		flags : 'Black_Hole:Alcubierre_Warp_Drive_Bubble:reset:vel:fsh',
		code : mlstr(function(){/*
varying vec3 vtxv;
void main() {
	vec3 vel = normalize(vtxv);
	gl_FragColor.xyz = vel;
*/})
	},
	{
		flags : 'Black_Hole:reset:vel:fsh',
		code : mlstr(function(){/*
//when initializing our metric:
//g_ab v^a v^b = 0 for our metric g
// (-1 + 2M/r) vt^2 + (vx^2 + vy^2 + vz^2) / (1 - 2M/r) = 0
// (1 - 2M/r) vt^2 = (vx^2 + vy^2 + vz^2) / (1 - 2M/r)
// vt^2 = (vx^2 + vy^2 + vz^2) / (1 - 2M/r)^2
// vt = ||vx,vy,vz|| / (1 - 2M/r)
	vec3 pos = vec3(-objectDist, 0., 0.);
	float r = length(pos);
	float oneMinus2MOverR = 1. - 2. * blackHoleMass / r;
	gl_FragColor.w = 1. / oneMinus2MOverR;
}
*/})
	},
	{
		flags : 'Alcubierre_Warp_Drive_Bubble:reset:vel:fsh',
		code : mlstr(function(){/*
	vec3 pos = vec3(0.);
	float rs = length(pos - vec3(objectDist, 0., 0.));
	float sigmaFront = warpBubbleThickness * (rs + warpBubbleRadius);
	float sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	float sigmaBack = warpBubbleThickness * (rs - warpBubbleRadius);
	float tanhSigmaCenter = tanh(sigmaCenter);
	float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanhSigmaCenter);
	float sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
	float dfScalar = sechDiff / (2. * rs * tanhSigmaCenter);				
	float vf = f * warpBubbleVelocity;
	float vf2 = vf * vf;
	gl_FragColor.w = 
		(warpBubbleVelocity * f + sqrt(
			vf2 * (1. + vel.x * vel.x) + 1.
		)) / (-1. + vf2);
}
*/})
	},
	{
		flags : 'Black_Hole:Alcubierre_Warp_Drive_Bubble:iterate:pos:vel:vsh',
		code : mlstr(function(){/*
attribute vec2 vertex;
varying vec2 uv;
varying vec3 vtxv;
uniform mat3 rotation;
void main() {
	uv = vertex;
	vtxv = rotation * vec3(vertex * 2. - 1., 1.);
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/})
	},
	{
		flags : 'Black_Hole:Alcubierre_Warp_Drive_Bubble:iterate:pos:vel:fsh',
		code : mlstr(function(){/*
varying vec2 uv;
varying vec3 vtxv;
uniform float deltaLambda;
uniform float objectDist;
uniform sampler2D lightPosTex;
uniform sampler2D lightVelTex;
*/})
	},
	{
		flags : 'Black_Hole:iterate:vel:fsh',
		code : mlstr(function(){/*
uniform float blackHoleMass;
*/})
	},
	{
		flags : 'Alcubierre_Warp_Drive_Bubble:iterate:vel:fsh',
		code : mlstr(function(){/*
uniform float warpBubbleThickness;
uniform float warpBubbleRadius;
uniform float warpBubbleVelocity;
*/})
	},
	{
		flags : 'Black_Hole:Alcubierre_Warp_Drive_Bubble:iterate:pos:vel:fsh',
		code : mlstr(function(){/*
void main() {
	vec4 pos = texture2D(lightPosTex, uv);
	vec4 vel = texture2D(lightVelTex, uv);

*/})
	},
	{
		flags : 'Black_Hole:Alcubierre_Warp_Drive_Bubble:iterate:pos:fsh',
		code : mlstr(function(){/*
	gl_FragColor = pos + deltaLambda * vel;
}
*/})
	},
	{
		flags : 'Black_Hole:iterate:vel:fsh',
		code : mlstr(function(){/*
	float r = length(pos.xyz);
	float oneMinus2MOverR = 1. - 2.*blackHoleMass/r;			
	float posDotVel = dot(pos.xyz, vel.xyz);
	float velDotVel = dot(vel.xyz, vel.xyz);
	float r2 = r * r;
	float invR2M = 1. / (r * oneMinus2MOverR);
	float rMinus2MOverR2 = oneMinus2MOverR / r;
	float MOverR2 = blackHoleMass / r2;
	gl_FragColor.x = vel.x - deltaLambda * MOverR2 * (rMinus2MOverR2 * pos.x * vel.w * vel.w + invR2M * (pos.x * velDotVel - 2. * vel.x * posDotVel));
	gl_FragColor.y = vel.y - deltaLambda * MOverR2 * (rMinus2MOverR2 * pos.y * vel.w * vel.w + invR2M * (pos.y * velDotVel - 2. * vel.y * posDotVel));
	gl_FragColor.z = vel.z - deltaLambda * MOverR2 * (rMinus2MOverR2 * pos.z * vel.w * vel.w + invR2M * (pos.z * velDotVel - 2. * vel.z * posDotVel));
	gl_FragColor.w = vel.w + deltaLambda * 2. * MOverR2 * invR2M * posDotVel * vel.w;
}
*/})
	},
	{
		flags : 'Alcubierre_Warp_Drive_Bubble:iterate:vel:fsh',
		code : mlstr(function(){/*
	float rs = sqrt(dot(pos,pos) + objectDist * (-2. * pos.x + objectDist) );
	//float rs = length(pos.xyz - vec3(objectDist, 0., 0.));
	float sigmaFront = warpBubbleThickness * (rs + warpBubbleRadius);
	float sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	float sigmaBack = warpBubbleThickness * (rs - warpBubbleRadius);
	float tanhSigmaCenter = tanh(sigmaCenter);
	float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanhSigmaCenter);
	float sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
	float dfScalar = sechDiff / (2. * rs * tanhSigmaCenter);
	vec4 fv;
	fv.xyz = warpBubbleThickness * pos.xyz * dfScalar;
	fv.w = -warpBubbleVelocity * warpBubbleThickness * pos.x * dfScalar;
	gl_FragColor.w = vel.w - deltaLambda * (f * f * fv.x * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * vel.w * vel.w
		- 2. * f * fv.x * warpBubbleVelocity * warpBubbleVelocity * vel.w * vel.x
		- 2. * f * fv.y * warpBubbleVelocity * warpBubbleVelocity / 2. * vel.w * vel.y
		- 2. * f * fv.z * warpBubbleVelocity * warpBubbleVelocity / 2. * vel.w * vel.z
		+ fv.x * warpBubbleVelocity * vel.x * vel.x
		+ 2. * fv.y * warpBubbleVelocity / 2. * vel.x * vel.y
		+ 2. * fv.z * warpBubbleVelocity / 2. * vel.x * vel.z
	);
	gl_FragColor.x = vel.x - deltaLambda * ((f * f * f * fv.x * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity - f * fv.x * warpBubbleVelocity * warpBubbleVelocity - fv.w * warpBubbleVelocity) * vel.w * vel.w
		- 2. * f * f * fv.x * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * vel.w * vel.x
		- 2. * (f * f * fv.y * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity + fv.y * warpBubbleVelocity) / 2. * vel.w * vel.y
		- 2. * (f * f * fv.z * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity + fv.z * warpBubbleVelocity) / 2. * vel.w * vel.z
		+ f * fv.x * warpBubbleVelocity * warpBubbleVelocity * vel.x * vel.x
		+ 2. * f * fv.y * warpBubbleVelocity * warpBubbleVelocity / 2. * vel.x * vel.y
		+ 2. * f * fv.z * warpBubbleVelocity * warpBubbleVelocity / 2. * vel.x * vel.z
	);
	gl_FragColor.y = vel.y + deltaLambda * (f * fv.y * warpBubbleVelocity * warpBubbleVelocity * vel.w * vel.w
		+ 2. * fv.y * warpBubbleVelocity / 2. * vel.w * vel.x
	);
	gl_FragColor.z = vel.z + deltaLambda * (f * fv.z * warpBubbleVelocity * warpBubbleVelocity * vel.w * vel.w
		+ 2. * fv.z * warpBubbleVelocity / 2. * vel.w * vel.x
	);
}
*/})
	}
];

function getScriptForFlags(flags) {
	flags = flags.clone();
	for (var i = 0; i < flags.length; ++i) {
		flags[i] = flags[i].replace(new RegExp(' ', 'g'), '_');
	}
	var code = '';
	$.each(shaderScriptParts, function(index, shaderScriptPart) {
		var parts = shaderScriptPart.flags.split(':');
		//if all flags are found in parts then use this shader
		//multiple concatenations? maybe later
		var failed = false;
		$.each(flags, function(_,flag) {
			if (parts.indexOf(flag) == -1) {
				failed = true;
				return false;	//break;
			}
		});
		if (failed) {
			return;	//continue;
		}
		var text = shaderScriptPart.code;
		stupidPrint('adding '+text);
		code += text; 
	});
	if (code == '') throw "couldn't find code for flags "+flags.join(':');
	stupidPrint('building '+flags.join(':')+' and getting '+code);
	return code;
}

function getShaderProgramArgsForFlags(flags) {
	var vshFlags = flags.clone();
	vshFlags.push('vsh');
	var vertexCode = getScriptForFlags(vshFlags);
	
	var fshFlags = flags.clone();
	fshFlags.push('fsh');
	var fragmentCode = getScriptForFlags(fshFlags);

	fragmentCode = shaderCommonCode + fragmentCode;

	return {
		vertexCode : vertexCode,
		fragmentCode : fragmentCode
	};
}

var SQRT_1_2 = Math.sqrt(.5);
var angleForSide = [
	[0, -SQRT_1_2, 0, -SQRT_1_2],
	[0, SQRT_1_2, 0, -SQRT_1_2],
	[SQRT_1_2, 0, 0, -SQRT_1_2],
	[SQRT_1_2, 0, 0, SQRT_1_2],
	[0, 0, 0, -1],
	[0, -1, 0, 0]
];

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	GL.resize();
}


var rotationMat3 = mat3.create();
function resetField() {
	gl.viewport(0, 0, lightTexWidth, lightTexHeight);
	$.each(lightPosVelChannels, function(_,channel) {
		for (var side = 0; side < 6; ++side) {
			var shader = channel.shaders[objectType].reset;
			
			var uniforms = {};
			if (channel.uniforms !== undefined) {
				for (var i = 0; i < channel.uniforms.length; ++i) {
					var uniformName = channel.uniforms[i];
					uniforms[uniformName] = window[uniformName];
				}
			}
			mat3.fromQuat(rotationMat3, angleForSide[side]);
			uniforms.rotation = rotationMat3;
			
			var fbo = channel.fbos[0][side];
			fbo.draw({
				callback:function(){
					GL.unitQuad.draw({
						shader:shader,
						uniforms:uniforms
					});
				}
			});
		}
	});
	gl.viewport(0, 0, GL.canvas.width, GL.canvas.height);
}


//update the uint8 array from the float array
//then upload the uint8 array to the gpu
function updateLightPosTex() {	
	gl.viewport(0, 0, lightTexWidth, lightTexHeight);
	$.each(lightPosVelChannels, function(_,channel) {
		for (var side = 0; side < 6; ++side) {
			var shader = channel.shaders[objectType].iterate;
			
			var uniforms = {};
			if (channel.uniforms !== undefined) {
				for (var i = 0; i < channel.uniforms.length; ++i) {
					var uniformName = channel.uniforms[i];
					uniforms[uniformName] = window[uniformName];
				}
			}
			mat3.fromQuat(rotationMat3, angleForSide[side]);
			uniforms.rotation = rotationMat3;
			uniforms.lightPosTex = 0;
			uniforms.lightVelTex = 1;
			
			var fbo = channel.fbos[1][side];
			fbo.draw({
				callback:function(){
					GL.unitQuad.draw({
						shader:shader,
						uniforms:uniforms,
						texs:[
							lightPosVelChannels[0].texs[0][side],
							lightPosVelChannels[1].texs[0][side]
						]
					});
				}
			});
		}
	});
	gl.viewport(0, 0, GL.canvas.width, GL.canvas.height);
	$.each(lightPosVelChannels, function(_,channel) {
		var tmp;
		tmp = channel.texs[0];
		channel.texs[0] = channel.texs[1];
		channel.texs[1] = tmp;
		tmp = channel.fbos[0];
		channel.fbos[0] = channel.fbos[1];
		channel.fbos[1] = tmp;
	});
}

// render loop
function update() {
	updateLightPosTex();
	
	//turn on magnification filter
	for (var side = 0; side < 6; ++side) {
		var tex = lightPosVelChannels[1].texs[0][side].obj;
		gl.bindTexture(gl.TEXTURE_2D, tex);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
		//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
		//gl.generateMipmap(gl.TEXTURE_2D);
		gl.bindTexture(gl.TEXTURE_2D, null);
	}
	
	for (var side = 0; side < 6; ++side) {
		//texs[0] is skyTex
		cubeSides[side].texs[1] = lightPosVelChannels[1].texs[0][side];
	}

	GL.draw();	
	
	//turn off magnification filter
	for (var side = 0; side < 6; ++side) {
		var tex = lightPosVelChannels[1].texs[0][side].obj;
		gl.bindTexture(gl.TEXTURE_2D, tex);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
		//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
		gl.bindTexture(gl.TEXTURE_2D, null);
	}	

	requestAnimFrame(update);
};

$(document).ready(main1);
	
var skyTexFilenames = [
	'skytex/sky-visible-cube-xp.png',
	'skytex/sky-visible-cube-xn.png',
	'skytex/sky-visible-cube-yp.png',
	'skytex/sky-visible-cube-yn.png',
	'skytex/sky-visible-cube-zp.png',
	'skytex/sky-visible-cube-zn.png'
];
function main1() {
	panel = $('#panel');
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()

	var objectTypeParamDivs = {};
	var refreshObjectTypeParamDivs = function() {
		$.each(objectTypeParamDivs, function(divObjectType,objectTypeParamDiv) {
			if (divObjectType == objectType) {
				objectTypeParamDiv.show();
			} else {
				objectTypeParamDiv.hide();
			}
		});
	};
	$.each(objectTypes, function(k,v) {
		objectTypeParamDivs[v] = $('#'+v.replace(new RegExp(' ', 'g'), '_')+'_params');
		var option = $('<option>', {text:v});
		option.appendTo($('#objectTypes'));
		if (v == objectType) {
			option.attr('selected', 'true');
		}
	});
	$('#objectTypes').change(function() {
		objectType = $('#objectTypes').val();
		refreshObjectTypeParamDivs();
		resetField();
	});
	refreshObjectTypeParamDivs();

	$.each([
		'deltaLambda',
		'objectDist',
		'blackHoleMass',
		'warpBubbleThickness',
		'warpBubbleVelocity',
		'warpBubbleRadius'
	], function(k,v) {
		var id = '#' + v;
		$(id).val(window[v]);
		$(id).change(function() {
			window[v] = $(id).val()*1;
			$(id).blur();
		});
	});
	
	$('#reset').click(function() {
		resetField();
	});

/* async for doesnt lkke pauses...
	$('#pause').click(function() {
		if (updateInterval === undefined) {
		} else {
		}
	});
*/

	try {
		gl = GL.init(canvas);
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	
	if (!gl.getExtension('OES_texture_float')) {
		throw 'This requires OES_texture_float';
	}

	gl.disable(gl.DITHER);

	$(skyTexFilenames).preload(main2);
}

var main2Initialized = false;
function main2() {
	if (main2Initialized) {
		console.log("main2 got called twice again.  check the preloader.");
		return;
	}
	main2Initialized = true; 
	
	GL.view.zNear = .1;
	GL.view.zFar = 100;
	GL.view.fovY = 90;
	quat.mul(GL.view.angle, /*90' x*/[SQRT_1_2,0,0,SQRT_1_2], /*90' -y*/[0,-SQRT_1_2,0,SQRT_1_2]);

	console.log('creating skyTex');
	var skyTex = new GL.TextureCube({
		flipY : true,
		generateMipmap : true,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR_MIPMAP_LINEAR,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		},
		urls : skyTexFilenames,
		done : function() {
			main3(this);
		}
	});
}

function main3(skyTex) {
	$.each(lightPosVelChannels, function(_,channel) {
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
				var tex = new GL.Texture2D({
					internalFormat : gl.RGBA,
					format : gl.RGBA,
					type : gl.FLOAT,
					width : lightTexWidth,
					height : lightTexHeight,
					magFilter : gl.NEAREST,
					minFilter : gl.NEAREST,
					wrap : {
						s : gl.CLAMP_TO_EDGE,
						t : gl.CLAMP_TO_EDGE
					}
				});
				texs[side] = tex;
				
				var fbo = new GL.Framebuffer();
				gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.obj);
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex.obj, 0);
				gl.bindFramebuffer(gl.FRAMEBUFFER, null);
				fbos[side] = fbo;
			}
		}
		
		channel.shaders = {};
		$.each(objectTypes, function(_,objType) {
			var flags = channel.flags.clone();
			flags.push(objType);
			channel.shaders[objType] = {};
			$.each(shaderTypes, function(_,shaderType) {
				var shflags = flags.clone();
				shflags.push(shaderType);
				var args = getShaderProgramArgsForFlags(shflags);
				args.uniforms = {
					lightPosTex : 0,
					lightVelTex : 1
				};
				args.vertexPrecision = 'best';
				args.fragmentPrecision = 'best';
				channel.shaders[objType][shaderType] = new GL.ShaderProgram(args);
			});
		});
	});

	var cubeShader = new GL.ShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;
varying vec2 uv;
uniform mat4 projMat;
uniform vec4 angle;
const mat3 viewMatrix = mat3(
	0., 0., -1.,
	1., 0., 0., 
	0., -1., 0.);
vec3 qtransform( vec4 q, vec3 v ){ 
	return v + 2.0*cross(cross(v, q.xyz ) + q.w*v, q.xyz);
}
void main() {
	uv = vertex;
	vec3 vtx3 = vec3(vertex * 2. - 1., 1.);
	vtx3 = qtransform(angle, vtx3);
	vtx3 = viewMatrix * vtx3;
	vec4 vtx4 = vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
}
*/}),
		fragmentCode : shaderCommonCode + mlstr(function(){/*
varying vec2 uv;
uniform samplerCube skyTex;
uniform sampler2D lightVelTex;
uniform vec4 viewAngle;
const mat3 viewMatrixInv = mat3(
	0., 1., 0.,
	0., 0., -1.,
	-1., 0., 0.);
vec3 qtransform( vec4 q, vec3 v ){ 
	return v + 2.0*cross(cross(v, q.xyz) + q.w*v, q.xyz);
}
void main() {
	vec3 dir = texture2D(lightVelTex, uv).xyz;
	dir = viewMatrixInv * viewMatrixInv * dir;
	dir = qtransform(vec4(viewAngle.xyz, -viewAngle.w), dir);
	gl_FragColor.xyz = textureCube(skyTex, dir).xyz;
	gl_FragColor.w = 1.; 
}*/}),
		uniforms : {
			skyTex : 0,
			lightVelTex : 1
		},
		vertexPrecision : 'best',
		fragmentPrecision : 'best'
	});

	//scene graph, 6 quads oriented in a cube
	//I would use a cubemap but the Mali-400 doesn't seem to want to use them as 2D FBO targets ...
	cubeSides = [];
	for (var side = 0; side < 6; ++side) {
		cubeSides[side] = new GL.SceneObject({
			//geometry : GL.unitQuad.geometry,
			mode : gl.TRIANGLE_STRIP,
			attrs : {
				vertex : GL.unitQuadVertexBuffer
			},
			uniforms : {
				viewAngle : GL.view.angle,
				angle : angleForSide[side]
			},
			shader : cubeShader,
			texs : [skyTex, lightPosVelChannels[1].texs[0][side]],
		});
	}

	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [dy, dx, 0], rotAngle);

			quat.mul(GL.view.angle, GL.view.angle, tmpQ);
			quat.normalize(GL.view.angle, GL.view.angle);
		},
		zoom : function(dz) {
			GL.view.fovY *= Math.exp(-.0003 * dz);
			GL.view.fovY = Math.clamp(GL.view.fovY, 1, 179);
			GL.updateProjection();
		}
	});

	resetField();

	$(window).resize(resize);
	resize();
	update();
}


