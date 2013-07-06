/*
The original rendered a viewport quad and subsequently iterated over it to raytrace through the scene. Worked great until you moved the camera.
I wanted to make a subsequent one that FBO'd six sides of a cubemap so you could at least turn your head (would still require redraw if you translated the camera).
But my current hardware doesn't even support float buffers, so I'll have to be creative about my current setup.
I could encode a cubemap as rgb => xyz, keep it normalized, and iterate through spacetime with that.  resolution would be low.  interpolation could help, or a 2nd tex for extra precision.
*/

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
var lightTexWidth = 256;
var lightTexHeight = 256;

var ident4 = mat4.create();

/*
flags: flags used to find what shader to associate with this
	they are combined with each object type and put in 'shaders'
	from that, shaders[objectType] determines the shader to run
*/
var shaderTypes = ['reset', 'iterate'];
var lightPosVelChannels = [
	{flags:['pos', 'x']},
	{flags:['pos', 'y']},
	{flags:['pos', 'z']},
	{flags:['pos', 't']},
	{flags:['vel', 'x']},
	{flags:['vel', 'y']},
	{flags:['vel', 'z']},
	{flags:['vel', 't']},
];
var unitQuadVertexBuffer;
var fboQuad;

function stupidPrint(s) {
	return;
	$.each(s.split('\n'), function(_,l) {
		console.log(l);
	});
}

function getScriptForFlags(flags) {
	flags = flags.clone();
	for (var i = 0; i < flags.length; ++i) {
		flags[i] = flags[i].replace(new RegExp(' ', 'g'), '_');
	}
	var code = '';
	$('script').each(function(index) {
		var id = $(this).attr('id');
		if (id === undefined) return;	//continue;
		var parts = id.split(':');
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
		var text = $(this).text();
		stupidPrint('adding '+text);
		code += text; 
		//return false;	//break;
	});
	if (code == '') throw "couldn't find code for flags "+flags.join(':');
	stupidPrint('building '+flags.join(':')+' and getting '+code);
	return code;
}

function buildShaderForFlags(flags) {
	var vshFlags = flags.clone();
	vshFlags.push('vsh');
	var vertexCode = getScriptForFlags(vshFlags);
	
	var fshFlags = flags.clone();
	fshFlags.push('fsh');
	var fragmentCode = getScriptForFlags(fshFlags);

	fragmentCode = 'precision mediump float;\n' + $('#shader-common').text() + fragmentCode;

	return new GL.ShaderProgram({
		vertexCode : vertexCode,
		fragmentCode : 
			(useFloatTextures ? '#define USE_TEXTURE_FLOAT\n' : '')
			+ fragmentCode
	});
}

var SQRT_1_2 = Math.sqrt(.5);
var angleForSide = [
	[SQRT_1_2,0,SQRT_1_2,0],
	[SQRT_1_2,0,-SQRT_1_2,0],
	[-SQRT_1_2,0,0,SQRT_1_2],
	[SQRT_1_2,0,0,SQRT_1_2],
	[1,0,0,0],
	[0,0,1,0]
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
					fboQuad.draw({
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
			uniforms.lightPosXTex = 0;
			uniforms.lightPosYTex = 1;
			uniforms.lightPosZTex = 2;
			uniforms.lightPosTTex = 3;
			uniforms.lightVelXTex = 4;
			uniforms.lightVelYTex = 5;
			uniforms.lightVelZTex = 6;
			uniforms.lightVelTTex = 7;
			
			var fbo = channel.fbos[1][side];
			fbo.draw({
				callback:function(){
					fboQuad.draw({
						shader:shader,
						uniforms:uniforms,
						texs:[
							lightPosVelChannels[0].texs[0][side],
							lightPosVelChannels[1].texs[0][side],
							lightPosVelChannels[2].texs[0][side],
							lightPosVelChannels[3].texs[0][side],
							lightPosVelChannels[4].texs[0][side],
							lightPosVelChannels[5].texs[0][side],
							lightPosVelChannels[6].texs[0][side],
							lightPosVelChannels[7].texs[0][side]
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
		for (var j = 4; j <= 6; ++j) {
			lightPosVelChannels[j].texs[0][side]
				.bind()
				.setArgs({magFilter:gl.LINEAR})
				.unbind();
		}
	}

	for (var side = 0; side < 6; ++side) {
		//texs[0] is skyTex
		cubeSides[side].texs[1] = lightPosVelChannels[4].texs[0][side];
		cubeSides[side].texs[2] = lightPosVelChannels[5].texs[0][side];
		cubeSides[side].texs[3] = lightPosVelChannels[6].texs[0][side];
	}

	//turn off magnification filter
	for (var side = 0; side < 6; ++side) {
		for (var j = 4; j <= 6; ++j) {
			lightPosVelChannels[j].texs[0][side]
				.bind()
				.setArgs({magFilter:gl.NEAREST})
				.unbind();
		}
	}	
	
	GL.draw();	
	requestAnimFrame(update);
};

var useFloatTextures = false; 
$(document).ready(function(){
	panel = $('#panel');	
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()

	$.each(objectTypes, function(k,v) {
		var option = $('<option>', {text:v});
		option.appendTo($('#objectTypes'));
		if (v == objectType) {
			option.attr('selected', 'true');
		}
	});
	$('#objectTypes').change(function() {
		objectType = $('#objectTypes').val();
		resetField();
	});

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
		gl = GL.init(canvas, {debug:true});
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	
	useFloatTextures = gl.getExtension('OES_texture_float');

	gl.disable(gl.DITHER);

	GL.view.zNear = .1;
	GL.view.zFar = 100;
	GL.view.fovY = 45;
	quat.mul(GL.view.angle, [SQRT_1_2,0,SQRT_1_2,0], [-SQRT_1_2,SQRT_1_2,0,0]);

	var skyTex = new GL.TextureCube({
		flipY : true,
		/*
		generateMipmap : true,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR_MIPMAP_LINEAR,
		*/
		magFilter : gl.LINEAR,
		minFilter : gl.NEAREST,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		},
		urls : [
			'skytex/sky-infrared-cube-xp.png',
			'skytex/sky-infrared-cube-xn.png',
			'skytex/sky-infrared-cube-yp.png',
			'skytex/sky-infrared-cube-yp.png',
			'skytex/sky-infrared-cube-zn.png',
			'skytex/sky-infrared-cube-zn.png'
		]
	});

	$.each(lightPosVelChannels, function(_,channel) {
		//working on how to organize this
		channel.uniforms = ['blackHoleMass', 'warpDriveThickness', 'warpDriveRadius', 'warpDriveVelocity', 'objectDist', 'deltaLambda'];
		
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
					type : useFloatTextures ? gl.FLOAT : gl.UNSIGNED_BYTE,
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
				channel.shaders[objType][shaderType] = buildShaderForFlags(shflags);
			});
		});
	});

	unitQuadVertexBuffer = new GL.ArrayBuffer({
		dim : 2,
		data : [0, 0, 1, 0, 0, 1, 1, 1]
	});

	//used for off-screen rendering and not part of the scene graph
	fboQuad = new GL.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		vertexBuffer : unitQuadVertexBuffer,
		parent : null
	});	
	
	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCode : 
			(useFloatTextures ? '#define USE_TEXTURE_FLOAT\n' : '')
			+ 'precision mediump float;\n'
			+ $('#shader-common').text()
			+ $('#cube-fsh').text(),
		uniforms : {
			skyTex : 0,
			lightVelXTex : 1,
			lightVelYTex : 2,
			lightVelZTex : 3
		}
	});


	//scene graph, 6 quads oriented in a cube
	//I would use a cubemap but the Mali-400 doesn't seem to want to use them as 2D FBO targets ...
	cubeSides = [];
	for (var side = 0; side < 6; ++side) {
		cubeSides[side] = new GL.SceneObject({
			mode : gl.TRIANGLE_STRIP,
			vertexBuffer : unitQuadVertexBuffer,
			shader : cubeShader,
			texs : [skyTex],
			angle : angleForSide[side]
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
});

