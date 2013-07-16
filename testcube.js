/*
simple skybox viewer.  like the volumizer, I should really make this its own utility ...
*/

var canvas;
var gl;
var mouse;
var cubeObj;

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

function update() {
	GL.draw();
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

var glMaxCubeMapTextureSize;
function main1(){
	console.log('main1');
	panel = $('#panel');	
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()

	try {
		gl = GL.init(canvas);
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	
	glMaxCubeMapTextureSize = gl.getParameter(gl.MAX_CUBE_MAP_TEXTURE_SIZE);

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
	GL.view.fovY = 45;
	quat.mul(GL.view.angle, [SQRT_1_2,0,SQRT_1_2,0], [-SQRT_1_2,SQRT_1_2,0,0]);

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
		onload : function(side,url,image) {
			if (image.width > glMaxCubeMapTextureSize || image.height > glMaxCubeMapTextureSize) {
				throw "cube map size "+image.width+"x"+image.height+" cannot exceed "+glMaxCubeMapTextureSize;
			}
		},
		done : function() {
			main3(this);
		}
	});
}

function main3(skyTex) {
	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCodeID : 'cube-fsh',
		uniforms : {
			skyTex : 0
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
		texs : [skyTex],
		static : false
	});
	
	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [-dy, dx, 0], rotAngle);

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
}

