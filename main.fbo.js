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
var updateInterval = undefined; 
//I'm updating in software, so on my tablet 512 is a bit slow
//if I stop iteration after reaching a steady state, maybe I'll up resolution then
//...or do some sort of adaptive thing ...
var lightTexWidth = 256;
var lightTexHeight = 256;

var ident4 = mat4.create();
mat4.identity(ident4);

/*
flags: flags used to find what shader to associate with this
	they are combined with each object type and put in 'shaders'
	from that, shaders[objectType] determines the shader to run
*/
var shaderTypes = ['reset'];
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
		code += text; 
		return false;	//break;
	});
	if (code == '') throw "couldn't find code for flags "+flags.join(':');
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
		fragmentCode : fragmentCode
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


function resetField() {
	if (updateInterval !== undefined) {
		clearInterval(updateInterval); 
	}
	updateInterval = undefined;

	var rotationMat3 = mat3.create();
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

	//updateLightPosTex();
}


//update the uint8 array from the float array
//then upload the uint8 array to the gpu
function updateLightPosTex() {	
	//var progress = $('#update-progress');
	//progress.attr('value', 0);
	updateInterval = asyncfor({
		start : 0,
		end : 6,
		callback : function(side) {
			var srci = side * lightTexWidth * lightTexHeight * 8;
			var dsti = 0;
			for (var i = 0; i < lightTexWidth * lightTexHeight; ++i) {
				//read positions and velocities
				var oldPx = lightBuf[srci+0];
				var oldPy = lightBuf[srci+1];
				var oldPz = lightBuf[srci+2];
				var oldPt = lightBuf[srci+3];
				var oldVx = lightBuf[srci+4];
				var oldVy = lightBuf[srci+5];
				var oldVz = lightBuf[srci+6];
				var oldVt = lightBuf[srci+7];
				
				//cache change in positions by velocities
				var newPx = oldPx + oldVx * deltaLambda;
				var newPy = oldPy + oldVy * deltaLambda;
				var newPz = oldPz + oldVz * deltaLambda;
				var newPt = oldPt + oldVt * deltaLambda;
				var newVx, newVy, newVz, newVt;
				
				//update velocity by geodesic equation
				if (objectType == 'Black Hole') {
					// Schwarzschild Cartesian metric
					//aux variables:
					var r = Math.sqrt(oldPx * oldPx + oldPy * oldPy + oldPz * oldPz);
					var oneMinus2MOverR = 1 - 2*blackHoleMass/r;			
					var posDotVel = oldPx * oldVx + oldPy * oldVy + oldPz * oldVz;
					var velDotVel = oldVx * oldVx + oldVy * oldVy + oldVz * oldVz;
					var r2 = r * r;
					var invR2M = 1 / (r * oneMinus2MOverR);
					var rMinus2MOverR2 = oneMinus2MOverR / r;
					var MOverR2 = blackHoleMass / r2;
					newVx = oldVx - deltaLambda * MOverR2 * (rMinus2MOverR2 * oldPx * oldVt * oldVt + invR2M * (oldPx * velDotVel - 2 * oldVx * posDotVel));
					newVy = oldVy - deltaLambda * MOverR2 * (rMinus2MOverR2 * oldPy * oldVt * oldVt + invR2M * (oldPy * velDotVel - 2 * oldVy * posDotVel));
					newVz = oldVz - deltaLambda * MOverR2 * (rMinus2MOverR2 * oldPz * oldVt * oldVt + invR2M * (oldPz * velDotVel - 2 * oldVz * posDotVel));
					newVt = oldVt + deltaLambda * 2 * MOverR2 * invR2M * posDotVel * oldVt;
				} else if (objectType == 'Alcubierre Warp Drive Bubble') {
					var r = Math.sqrt(oldPx * oldPx + oldPy * oldPy + oldPz * oldPz);
					var sigmaFront = warpBubbleThickness * (r + warpBubbleRadius);
					var sigmaCenter = warpBubbleThickness * r;
					var sigmaBack = warpBubbleThickness * (r - warpBubbleRadius);
					var tanhSigmaCenter = tanh(sigmaCenter);
					var f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2 * tanhSigmaCenter);
					var sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
					var dfScalar = sechDiff / (2 * r * tanhSigmaCenter);
					var ft = -warpBubbleVelocity * warpBubbleThickness * oldPx * dfScalar;
					var fx = warpBubbleThickness * oldPx * dfScalar;
					var fy = warpBubbleThickness * oldPy * dfScalar;
					var fz = warpBubbleThickness * oldPz * dfScalar;
			
					//if I ever choose to keep track of v^t...
					newVt = oldVt - deltaLambda * (f * f * fx * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * oldVt * oldVt
						- 2. * f * fx * warpBubbleVelocity * warpBubbleVelocity * oldVt * oldVx
						- 2. * f * fy * warpBubbleVelocity * warpBubbleVelocity / 2. * oldVt * oldVy
						- 2. * f * fz * warpBubbleVelocity * warpBubbleVelocity / 2. * oldVt * oldVz
						+ fx * warpBubbleVelocity * oldVx * oldVx
						+ 2. * fy * warpBubbleVelocity / 2. * oldVx * oldVy
						+ 2. * fz * warpBubbleVelocity / 2. * oldVx * oldVz
					);
					newVx = oldVx - deltaLambda * ((f * f * f * fx * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity - f * fx * warpBubbleVelocity * warpBubbleVelocity - ft * warpBubbleVelocity) * oldVt * oldVt
						- 2. * f * f * fx * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * oldVt * oldVx
						- 2. * (f * f * fy * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity + fy * warpBubbleVelocity) / 2. * oldVt * oldVy
						- 2. * (f * f * fz * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity + fz * warpBubbleVelocity) / 2. * oldVt * oldVz
						+ f * fx * warpBubbleVelocity * warpBubbleVelocity * oldVx * oldVx
						+ 2. * f * fy * warpBubbleVelocity * warpBubbleVelocity / 2. * oldVx * oldVy
						+ 2. * f * fz * warpBubbleVelocity * warpBubbleVelocity / 2. * oldVx * oldVz
					);
					newVy = oldVy + deltaLambda * (f * fy * warpBubbleVelocity * warpBubbleVelocity * oldVt * oldVt
						+ 2. * fy * warpBubbleVelocity / 2. * oldVt * oldVx
					);
					newVz = oldVz + deltaLambda * (f * fz * warpBubbleVelocity * warpBubbleVelocity * oldVt * oldVt
						+ 2. * fz * warpBubbleVelocity / 2. * oldVt * oldVx
					);
				}
				// write back results
				lightBuf[srci++] = newPx;
				lightBuf[srci++] = newPy;
				lightBuf[srci++] = newPz;
				lightBuf[srci++] = newPt;
				lightBuf[srci++] = newVx;
				lightBuf[srci++] = newVy;
				lightBuf[srci++] = newVz;
				lightBuf[srci++] = newVt;
				//don't bother update vw, I don't store it and just reset it afterwards
				//copy floats to uint8 texture
				var s = Math.sqrt(newVx * newVx + newVy * newVy + newVz * newVz);
				lightPosTexData[side][dsti++]  = 255 * (newVx / s * .5 + .5);
				lightPosTexData[side][dsti++]  = 255 * (newVy / s * .5 + .5);
				lightPosTexData[side][dsti++]  = 255 * (newVz / s * .5 + .5);
			}
		},
		done : function() {
			lightPosTex.bind();
			for (var side = 0; side < 6; ++side) {
				assertEquals(lightTexWidth * lightTexHeight * 3, lightPosTexData[side].length); 
				//gl.texSubImage2D(gl.TEXTURE_CUBE_MAP_POSITIVE_X + side, 0, 0, 0, lightTexWidth, lightTexHeight, gl.RGB, gl.UNSIGNED_BYTE, lightPosTexData[side]);
				gl.texImage2D(gl.TEXTURE_CUBE_MAP_POSITIVE_X + side, 0, gl.RGB, lightTexWidth, lightTexHeight, 0, gl.RGB, gl.UNSIGNED_BYTE, lightPosTexData[side]);
			}
			lightPosTex.unbind();		
			
			updateLightPosTex();	
		}
	});
}

// render loop

var tmpRotMat = mat4.create();	
function update() {
	for (var side = 0; side < 6; ++side) {
		//texs[0] is skyTex
		cubeSides[side].texs[1] = lightPosVelChannels[4].texs[0][side];
		cubeSides[side].texs[2] = lightPosVelChannels[5].texs[0][side];
		cubeSides[side].texs[3] = lightPosVelChannels[6].texs[0][side];
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
		gl = GL.init(canvas, {debug:true});
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
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
			'skytex/sky-infrared-cube-xp.png',
			'skytex/sky-infrared-cube-xn.png',
			'skytex/sky-infrared-cube-yp.png',
			'skytex/sky-infrared-cube-yp.png',
			'skytex/sky-infrared-cube-zn.png',
			'skytex/sky-infrared-cube-zn.png'
		]
	});

	$.each(lightPosVelChannels, function(_,channel) {
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
					type : gl.UNSIGNED_BYTE,
					width : lightTexWidth,
					height : lightTexHeight,
					magFilter : gl.LINEAR,
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
		uniforms : {
			projMat : ident4,
			mvMat : ident4 
		},
		parent : null
	});	
	
	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCode : 
			'precision mediump float;\n'
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
		var angle = angleForSide[side];
		var m3 = mat3.create();
		mat3.fromQuat(m3, angle);
		//negative z is forward
		var v3 = vec3.create();
		v3[0] = -.5 * m3[6];
		v3[1] = -.5 * m3[7];
		v3[2] = -.5 * m3[8];

		cubeSides[side] = new GL.SceneObject({
			mode : gl.TRIANGLE_STRIP,
			vertexBuffer : unitQuadVertexBuffer,
			shader : cubeShader,
			uniforms : {
				rotMat3 : m3 
			},
			//angle : angle, 
			//pos : v3,
			//here's a breaking point of the current structure:
			//for static objects, the scene graph "optimizes" and uses the parent if no transformation information is provided
			//for non-static objects, the mvMat member is recomputed each draw
			//... what if we want a member mat that is modified in ways other than those recomputed ?  like scale?
			//  currently no dice.
			//texs will look like [skyTex, lightPosVelChannels[4, 5, and 6].texs[0][side]],
			// but those textures are ever-changing and will need to be updated each frame
			texs : [skyTex],
			angle : angleForSide[side],
		});
		
		//...fixing the above mentioned problem...
		//cubeSides[side].uniforms.projMat = GL.projMat;
		//var mvMat = mat4.create();
		//mat4.fromQuat(mvMat, angleForSide[side]);
		//cubeSides[side].uniforms.mvMat = mvMat;
		//for some reason i can't access the same uniform from within my vertex and fragment shader ...
		//cubeSides[side].uniforms.mvMat2 = mvMat;
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

