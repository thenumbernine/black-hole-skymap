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
var deltaLambda = .1;	//ray forward iteration
var updateInterval = undefined; 
//I'm updating in software, so on my tablet 512 is a bit slow
//if I stop iteration after reaching a steady state, maybe I'll up resolution then
//...or do some sort of adaptive thing ...
var lightTexWidth = 256;
var lightTexHeight = 256;
var lightBuf; 
var lightVelTexData = [];
for (var side = 0; side < 6; ++side) {
	lightVelTexData[side] = new Uint8Array(3 * lightTexWidth * lightTexHeight);
}
var lightVelTex; 

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

	var vel = vec4.create();
	var pos = vec4.create();
	//lightBuf[side][u][v][x,y,z,w,vx,vy,vz,vw]
	var i = 0;
	for (var side = 0; side < 6; ++side) {
		for (var v = 0; v < lightTexHeight; ++v) {
			for (var u = 0; u < lightTexWidth; ++u) {
				vel[0] = (u + .5) / lightTexWidth * 2 - 1;
				vel[1] = 1 - (v + .5) / lightTexHeight * 2;
				vel[2] = 1;
				vec3.transformQuat(vel,vel,angleForSide[side]);

				//[3] will be the time (0'th) coordinate
				//light velocities must be unit in minkowski space
				//-vt^2 + vx^2 + vy^2 + vz^2 = -1
				//vt^2 = vx^2 + vy^2 + vz^2 + 1
				//vel[3] = Math.sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2] + 1);
				//this still allows a degree of freedom in the spatial length ... (why is this?)
				//let's keep it normalized
				vec3.normalize(vel, vel);
				//vec4.copy(pos, vel);
				//pos[3] = 0;
				vec4.set(pos, 0,0,0,0);

				if (objectType == 'Black Hole') {
					//when initializing our metric:
					//g_ab v^a v^b = 0 for our metric g
					// (-1 + 2M/r) vt^2 + (vx^2 + vy^2 + vz^2) / (1 - 2M/r) = 0
					// (1 - 2M/r) vt^2 = (vx^2 + vy^2 + vz^2) / (1 - 2M/r)
					// vt^2 = (vx^2 + vy^2 + vz^2) / (1 - 2M/r)^2
					// vt = ||vx,vy,vz|| / (1 - 2M/r)
					pos[0] -= objectDist;	//center of metric is black hole origin
					var r = vec3.length(pos);
					var oneMinus2MOverR = 1 - 2*blackHoleMass/r;
					vel[3] = 1 / oneMinus2MOverR;
				} else if (objectType == 'Alcubierre Warp Drive Bubble') {
					//g_ab v^a v^b = 0
					//... later
					var rs = vec3.dist(pos, vec3.fromValues(objectDist,0,0));	//center is view pos, but rs is measured from bubble origin
					var sigmaFront = warpBubbleThickness * (rs + warpBubbleRadius);
					var sigmaCenter = warpBubbleThickness * rs;
					var sigmaBack = warpBubbleThickness * (rs - warpBubbleRadius);
					var tanhSigmaCenter = tanh(sigmaCenter);
					var f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2 * tanhSigmaCenter);
				
					var vf = f * warpBubbleVelocity;
					var vf2 = vf * vf;
					vel[3] = 
						(warpBubbleVelocity * f + Math.sqrt(
							vf2 * (1 + vel[0] * vel[0]) + 1
						)) / (-1 + vf2);
				}
				
				//position (relative to black hole)
				lightBuf[i++] = pos[0];
				lightBuf[i++] = pos[1];
				lightBuf[i++] = pos[2];
				lightBuf[i++] = pos[3];
				//velocity
				lightBuf[i++] = vel[0];
				lightBuf[i++] = vel[1];
				lightBuf[i++] = vel[2];
				lightBuf[i++] = vel[3]; 
			}
		}
	}
	updateLightVelTex();
}


//update the uint8 array from the float array
//then upload the uint8 array to the gpu
function updateLightVelTex() {	
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
				lightVelTexData[side][dsti++]  = 255 * (newVx / s * .5 + .5);
				lightVelTexData[side][dsti++]  = 255 * (newVy / s * .5 + .5);
				lightVelTexData[side][dsti++]  = 255 * (newVz / s * .5 + .5);
			}
		},
		done : function() {
			lightVelTex.bind();
			for (var side = 0; side < 6; ++side) {
				assertEquals(lightTexWidth * lightTexHeight * 3, lightVelTexData[side].length); 
				//gl.texSubImage2D(gl.TEXTURE_CUBE_MAP_POSITIVE_X + side, 0, 0, 0, lightTexWidth, lightTexHeight, gl.RGB, gl.UNSIGNED_BYTE, lightVelTexData[side]);
				gl.texImage2D(gl.TEXTURE_CUBE_MAP_POSITIVE_X + side, 0, gl.RGB, lightTexWidth, lightTexHeight, 0, gl.RGB, gl.UNSIGNED_BYTE, lightVelTexData[side]);
			}
			lightVelTex.unbind();		
			
			updateLightVelTex();	
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
	
	GL.view.zNear = .1;
	GL.view.zFar = 100;
	GL.view.fovY = 45;
	quat.mul(GL.view.angle, [SQRT_1_2,0,SQRT_1_2,0], [-SQRT_1_2,SQRT_1_2,0,0]);

	var skyTex = new GL.TextureCube({
		flipY : true,
		generateMipmap : true,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR_MIPMAP_LINEAR,
		/*
		magFilter : gl.LINEAR,
		minFilter : gl.NEAREST,
		*/
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
	
	lightVelTex = new GL.TextureCube({
		internalFormat : gl.RGB,
		format : gl.RGB,
		type : gl.UNSIGNED_BYTE,
		width : lightTexWidth,
		height : lightTexHeight,
		data : lightVelTexData[side],
		magFilter : gl.LINEAR,
		minFilter : gl.NEAREST,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		}
	});

	lightBuf = new Float32Array(6 * 4 * 2 * lightTexWidth * lightTexHeight);
	resetField();

	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCodeID : 'cube-fsh',
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
		uniforms : {
			viewAngle : GL.view.angle
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

	$(window).resize(resize);
	resize();
	update();
});

