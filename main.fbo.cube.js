var canvas;
var gl;
var mouse;
var cubeObj;
var objectTypes = ['Black Hole', 'Alcubierre Warp Drive Bubble'];
var objectType = objectTypes[0];
var objectDist = 10;
var blackHoleMass = 1;
var warpBubbleThickness = 1;
var warpBubbleVelocity = .5;
var warpBubbleRadius = 2;
var deltaLambda = .1;	//ray forward iteration
var lightTexWidth = 512;
var lightTexHeight = 512;

var ident4 = mat4.create();

/*
flags: flags used to find what shader to associate with this
	they are combined with each object type and put in 'shaders'
	from that, shaders[objectType] determines the shader to run
*/
var lightPosTex, lightVelTex;
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

function getShaderProgramArgsForFlags(flags) {
	var vshFlags = flags.clone();
	vshFlags.push('vsh');
	var vertexCode = getScriptForFlags(vshFlags);
	
	var fshFlags = flags.clone();
	fshFlags.push('fsh');
	var fragmentCode = getScriptForFlags(fshFlags);

	fragmentCode = 'precision highp float;\n' + $('#shader-common').text() + fragmentCode;

	return {
		vertexCode : vertexCode,
		fragmentCode : fragmentCode
	};
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
	
	for (var side = 0; side < 6; ++side) {
		var shader = lightPosTex.shaders[objectType].reset;
		
		if (lightPosTex.uniforms !== undefined) {
			for (var i = 0; i < lightPosTex.uniforms.length; ++i) {
				var uniformName = lightPosTex.uniforms[i];
				lightPosTex.uniforms[uniformName] = window[uniformName];
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
					fboQuad.draw({
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
		lightPosVelChannels[1].texs[0][side]
			.bind()
			.setArgs({magFilter:gl.LINEAR})
			.unbind();
	}
	
	//texs[0] is skyTex
	cubeObj.texs[1] = lightVelTex.texs[0][side];

	//turn off magnification filter
	for (var side = 0; side < 6; ++side) {
		lightPosVelChannels[1].texs[0][side]
			.bind()
			.setArgs({magFilter:gl.NEAREST})
			.unbind();
	}	
	
	GL.draw();	
	requestAnimFrame(update);
};

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
			'skytex/sky-visible-cube-xp.png',
			'skytex/sky-visible-cube-xn.png',
			'skytex/sky-visible-cube-yp.png',
			'skytex/sky-visible-cube-yp.png',
			'skytex/sky-visible-cube-zn.png',
			'skytex/sky-visible-cube-zn.png'
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
				channel.shaders[objType][shaderType] = new GL.ShaderProgram(args);
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
		attrs : {
			vertex : unitQuadVertexBuffer
		},
		parent : null
	});	
	
	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCode : 
			'precision highp float;\n'
			+ $('#shader-common').text()
			+ $('#cube-fsh').text(),
		uniforms : {
			skyTex : 0,
			lightVelTex : 1,
		}
	});

	var cubeVtxArray = new Float32Array(3*8);
	for (var i = 0; i < 8; i++) {
		cubeVtxArray[0+3*i] = 2*(i&1)-1;
		cubeVtxArray[1+3*i] = 2*((i>>1)&1)-1;
		cubeVtxArray[2+3*i] = 2*((i>>2)&1)-1;
	}

	cubeVtxBuf = new GL.ArrayBuffer({
		data : cubeVtxArray 
	});

	var cubeIndexBuf = new GL.ElementArrayBuffer({
		data : [
			5,7,3,3,1,5,		// <- each value has the x,y,z in the 0,1,2 bits (off = 0, on = 1)
			6,4,0,0,2,6,
			2,3,7,7,6,2,
			4,5,1,1,0,4,
			6,7,5,5,4,6,
			0,1,3,3,2,0
		]
	});

	cubeObj = new GL.SceneObject({
		mode : gl.TRIANGLES,
		attrs : {
			vertex : cubeVtxBuf
		},
		indexes : cubeIndexBuf,
		shader : cubeShader,
		texs : [skyTex, lightVelTex],
		static : false
	});

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

