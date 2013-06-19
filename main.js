/*
The original rendered a viewport quad and subsequently iterated over it to raytrace through the scene. Worked great until you moved the camera.
I wanted to make a subsequent one that FBO'd six sides of a cubemap so you could at least turn your head (would still require redraw if you translated the camera).
But my current hardware doesn't even support float buffers, so I'll have to be creative about my current setup.
I could encode a cubemap as rgb => xyz, keep it normalized, and iterate through spacetime with that.  resolution would be low.  interpolation could help, or a 2nd tex for extra precision.
*/

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
var dLambda = .1;	//ray forward iteration
var updateInterval = undefined; 
//I'm updating in software, so on my tablet 512 is a bit slow
//if I stop iteration after reaching a steady state, maybe I'll up resolution then
//...or do some sort of adaptive thing ...
var lightTexWidth = 256;
var lightTexHeight = 256;
var lightBuf; 
var lightPosTexData = [];
for (var side = 0; side < 6; ++side) {
	lightPosTexData[side] = new Uint8Array(3 * lightTexWidth * lightTexHeight);
}
var lightPosTex; 

function tanh(x) {
	var exp2x = Math.exp(2 * x);
	return (exp2x - 1) / (exp2x + 1);
}

function sech(x) {
	var expx = Math.exp(x);
	return 2. * expx / (expx * expx + 1.);
}

function sechSq(x) {
	var y = sech(x);
	return y * y;
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

	var v3 = vec4.create();
	//lightBuf[side][u][v][x,y,z,w,vx,vy,vz,vw]
	var i = 0;
	for (var side = 0; side < 6; ++side) {
		for (var v = 0; v < lightTexHeight; ++v) {
			for (var u = 0; u < lightTexWidth; ++u) {
				v3[0] = (u + .5) / lightTexWidth * 2 - 1;
				v3[1] = 1 - (v + .5) / lightTexHeight * 2;
				v3[2] = 1;
				vec3.transformQuat(v3,v3,angleForSide[side]);

				//[3] will be the time (0'th) coordinate
				//light velocities must be unit in minkowski space
				//-vt^2 + vx^2 + vy^2 + vz^2 = -1
				//vt^2 = vx^2 + vy^2 + vz^2 + 1
				//v3[3] = Math.sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] + 1);
				//this still allows a degree of freedom in the spatial length ... (why is this?)
				//let's keep it normalized
				vec3.normalize(v3, v3);
				
				//position (relative to black hole)
				lightBuf[i++] = v3[0] - objectDist;
				lightBuf[i++] = v3[1];
				lightBuf[i++] = v3[2];
				//velocity
				lightBuf[i++] = v3[0];
				lightBuf[i++] = v3[1];
				lightBuf[i++] = v3[2];
			}
		}
	}
	updateLightPosTex();
}


//update the uint8 array from the float array
//then upload the uint8 array to the gpu
function updateLightPosTex() {	
	//var progress = $('#update-progress');
	//progress.attr('value', 0);
	var p = vec4.create();
	var v = vec4.create();
	updateInterval = asyncfor({
		start : 0,
		end : 6,
		callback : function(side) {
			var srci = side * lightTexWidth * lightTexHeight * 6;
			var dsti = 0;
			for (var i = 0; i < lightTexWidth * lightTexHeight; ++i) {
				//read positions and velocities
				p[0] = lightBuf[srci++];
				p[1] = lightBuf[srci++];
				p[2] = lightBuf[srci++];
				v[0] = lightBuf[srci++];
				v[1] = lightBuf[srci++];
				v[2] = lightBuf[srci++];
				
				//update positions by velocities
				p[0] += v[0] * dLambda;
				p[1] += v[1] * dLambda;
				p[2] += v[2] * dLambda;
			
				//update velocity by geodesic equation
				if (objectType == 'Black Hole') {
					// Schwarzschild Cartesian metric
					//aux variables:
					var r = vec3.length(p);
					var oneMinus2MOverR = 1 - 2*blackHoleMass/r;
					//g_ab v^a v^b = 0 for our metric g
					// (-1 + 2M/r) vt^2 + (vx^2 + vy^2 + vz^2) / (1 - 2M/r) = 0
					// (1 - 2M/r) vt^2 = (vx^2 + vy^2 + vz^2) / (1 - 2M/r)
					// vt^2 = (vx^2 + vy^2 + vz^2) / (1 - 2M/r)^2
					// vt = ||vx,vy,vz|| / (1 - 2M/r)
					var vtSq = vec3.dot(v, v) / (oneMinus2MOverR * oneMinus2MOverR);
			
					var posDotVel = p[0] * v[0] + p[1] * v[1] + p[2] * v[2];
					var posDotVelSq = p[0] * v[0] * v[0] + p[1] * v[1] * v[1] + p[2] * v[2] * v[2];
					var r2 = r * r;
					var invR2M = 1 / (r * oneMinus2MOverR);
					var rMinus2MOverR2 = oneMinus2MOverR / r;
					var MOverR2 = blackHoleMass / r2;
					v[0] -= dLambda * MOverR2 * (rMinus2MOverR2 * p[0] * vtSq + invR2M * (2 * v[0] * posDotVel - posDotVelSq));
					v[1] -= dLambda * MOverR2 * (rMinus2MOverR2 * p[1] * vtSq + invR2M * (2 * v[1] * posDotVel - posDotVelSq));
					v[2] -= dLambda * MOverR2 * (rMinus2MOverR2 * p[2] * vtSq + invR2M * (2 * v[2] * posDotVel - posDotVelSq));
				} else if (objectType == 'Alcubierre Warp Drive Bubble') {
					var vt = vec3.length(v);
					
					var r = vec3.length(p);
					var sigmaFront = warpBubbleThickness * (r + warpBubbleRadius);
					var sigmaCenter = warpBubbleThickness * r;
					var sigmaBack = warpBubbleThickness * (r - warpBubbleRadius);
					var tanhSigmaCenter = tanh(sigmaCenter);
					var f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2 * tanhSigmaCenter);
					var sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
					var dfScalar = sechDiff / (2 * r * tanhSigmaCenter);
					var ft = -warpBubbleVelocity * warpBubbleThickness * p[0] * dfScalar;
					var fx = warpBubbleThickness * p[0] * dfScalar;
					var fy = warpBubbleThickness * p[1] * dfScalar;
					var fz = warpBubbleThickness * p[2] * dfScalar;
			
					//if I ever choose to keep track of v^t...
					//negRelDiff2.w = f * f * fx * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * vt * vt
					//	- 2. * f * fx * warpBubbleVelocity * warpBubbleVelocity * vt * v[0]
					//	- 2. * f * fy * warpBubbleVelocity * warpBubbleVelocity / 2. * vt * v[1]
					//	- 2. * f * fz * warpBubbleVelocity * warpBubbleVelocity / 2. * vt * v[2]
					//	+ fx * warpBubbleVelocity * v[0] * v[0]
					//	+ 2. * fy * warpBubbleVelocity / 2. * v[0] * v[1]
					//	+ 2. * fz * warpBubbleVelocity / 2. * v[0] * v[2]
					//;
					v[0] -= dLambda * ((f * f * f * fx * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity - f * fx * warpBubbleVelocity * warpBubbleVelocity - ft * warpBubbleVelocity) * vt * vt
						- 2. * f * f * fx * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity * vt * v[0]
						- 2. * (f * f * fy * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity + fy * warpBubbleVelocity) / 2. * vt * v[1]
						- 2. * (f * f * fz * warpBubbleVelocity * warpBubbleVelocity * warpBubbleVelocity + fz * warpBubbleVelocity) / 2. * vt * v[2]
						+ f * fx * warpBubbleVelocity * warpBubbleVelocity * v[0] * v[0]
						+ 2. * f * fy * warpBubbleVelocity * warpBubbleVelocity / 2. * v[0] * v[1]
						+ 2. * f * fz * warpBubbleVelocity * warpBubbleVelocity / 2. * v[0] * v[2]
					);
					v[1] += dLambda * (f * fy * warpBubbleVelocity * warpBubbleVelocity * vt * vt
						+ 2. * fy * warpBubbleVelocity / 2. * vt * v[0]
					);
					v[2] += dLambda * (f * fz * warpBubbleVelocity * warpBubbleVelocity * vt * vt
						+ 2. * fz * warpBubbleVelocity / 2. * vt * v[0]
					);
				}
				// write back results
				srci -= 6;
				lightBuf[srci++] = p[0];
				lightBuf[srci++] = p[1];
				lightBuf[srci++] = p[2];
				lightBuf[srci++] = v[0];
				lightBuf[srci++] = v[1];
				lightBuf[srci++] = v[2];
				//don't bother update v[3], I don't store it and just reset it afterwards
				//copy floats to uint8 texture
				var s = vec3.length(v);
				lightPosTexData[side][dsti++]  = 255 * (v[0] / s * .5 + .5);
				lightPosTexData[side][dsti++]  = 255 * (v[1] / s * .5 + .5);
				lightPosTexData[side][dsti++]  = 255 * (v[2] / s * .5 + .5);
			}
			//progress.attr('value', 100 * (side + 1) / 6);
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
	
	lightPosTex = new GL.TextureCube({
		internalFormat : gl.RGB,
		format : gl.RGB,
		type : gl.UNSIGNED_BYTE,
		width : lightTexWidth,
		height : lightTexHeight,
		data : lightPosTexData[side],
		magFilter : gl.LINEAR,
		minFilter : gl.NEAREST,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		}
	});

	lightBuf = new Float32Array(6 * 3 * 2 * lightTexWidth * lightTexHeight);
	resetField();

	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCodeID : 'cube-fsh',
		uniforms : {
			skyTex : 0,
			lightPosTex : 1,
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
		vertexBuffer : cubeVtxBuf,
		indexBuffer : cubeIndexBuf,
		shader : cubeShader,
		texs : [skyTex, lightPosTex],
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

	$(window).resize(resize);
	resize();
	update();
});

