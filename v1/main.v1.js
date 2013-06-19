/*
The original rendered a viewport quad and subsequently iterated over it to raytrace through the scene. Worked great until you moved the camera.
I wanted to make a subsequent one that FBO'd six sides of a cubemap so you could at least turn your head (would still require redraw if you translated the camera).
But my current hardware doesn't even support float buffers, so I'll have to be creative about my current setup.
I could encode a cubemap as rgb => xyz, keep it normalized, and iterate through spacetime with that.  resolution would be low.  interpolation could help, or a 2nd tex for extra precision.
*/

var canvas;
var gl;
var mouse;
var skyTex;
var doIteration = 2;
var lightInitialized = false;
var drawLightShader;
var initLightPosShader, initLightVelShader;
var iterateLightPosShader, iterateLightVelShader;
var lightPosTexs = [];
var lightVelTexs = [];
var fbo;
var viewVtx, viewVtxBuf, viewQuad;
var unitQuad;
var unitOrthoMat = mat4.create();
mat4.ortho(unitOrthoMat,0,1,0,1,-1,1);
var identMat = mat4.create();

var uvs = [
	[0,0],
	[1,0],
	[0,1],
	[1,1]
];

function onresize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	GL.resize();
}

// render loop

var update = function() {
	
	if (!lightInitialized) {
		lightInitialized = true;
	
		gl.viewport(0, 0, fbo.width, fbo.height);
		var aspectRatio = canvas.width / canvas.height;
		
		$.each([{
			tex : lightPosTexs[0], 
			shader : initLightPosShader
		},{
			tex : lightVelTexs[0],
			shader : initLightVelShader
		}], function(_,args){
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, args.tex.target, args.tex.obj, 0); 
			fbo.check();

			$.each(uvs, function(i,uv){
				viewVtx[0+3*i] = (uv[0] - .5) * 2 * aspectRatio * GL.view.fovY;
				viewVtx[1+3*i] = (uv[1] - .5) * 2 * GL.view.fovY;
				viewVtx[2+3*i] = -1;
			});
			viewVtxBuf.updateData(viewVtx);
			viewQuad.draw({
				shader : args.shader
			});
			
			fbo.unbind();
		});
	}

	if (doIteration != 0) {
		if (doIteration == 1) doIteration = 0;
	
		$.each([{
			tex:lightPosTexs[1],
			shader:iterateLightPosShader
		},{
			tex:lightVelTexs[1],
			shader:iterateLightVelShader
		}], function(_,args){
			fbo.bind();
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, args.tex.target, args.tex.obj, 0); 
			fbo.check();
			unitQuad.draw({
				shader : args.shader,
				texs : [lightPosTexs[0], lightVelTexs[0]],
				uniforms : {
					projMat : unitOrthoMat,
					mvMat : identMat
				}
			});
			fbo.unbind();
		});
		
		var tmp;
		
		tmp = lightPosTexs[0];
		lightPosTexs[0] = lightPosTexs[1];
		lightPosTexs[1] = tmp;
		
		tmp = lightVelTexs[0];
		lightVelTexs[0] = lightVelTexs[1];
		lightVelTexs[1] = tmp;
	}

	gl.viewport(0, 0, canvas.width, canvas.height);
	gl.clear(gl.COLOR_BUFFER_BIT);
	
	unitQuad.draw({
		shader : drawLightShader,
		texs : [lightPosTexs[0], skyTex],
		uniforms : {
			projMat : unitOrthoMat,
			mvMat : identMat
		}
	});

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
	var sqrt1_2 = Math.sqrt(.5);
	quat.mul(GL.view.angle, [0,sqrt1_2,0,sqrt1_2], [sqrt1_2,0,0,-sqrt1_2]);

	skyTex = new GL.TextureCube({
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
			'sky-infrared-cube-xp.png',
			'sky-infrared-cube-xn.png',
			'sky-infrared-cube-yp.png',
			'sky-infrared-cube-yp.png',
			'sky-infrared-cube-zn.png',
			'sky-infrared-cube-zn.png'
		]
	});

/* the original code */
	var shaderDefs = $('#shader-defs').text();
	var initLightShaderVertexCode = shaderDefs + $('#init-light-vsh').text();
	var initLightShaderFragmentCode = 'precision mediump float;\n' + shaderDefs + $('#init-light-fsh').text(); 
	
	var initLightPosShader = new GL.ShaderProgram({
		vertexCode : initLightShaderVertexCode,
		fragmentCode : initLightShaderFragmentCode.replace('$assign', 'gl_FragColor = rel;')
	});
		
	var initLightVelShader = new GL.ShaderProgram({
		vertexCode : initLightShaderVertexCode,
		fragmentCode : initLightShaderFragmentCode.replace('$assign', 'gl_FragColor = relDiff;')
	})

	var lightRes = 1024;
	$.each([lightPosTexs, lightVelTexs], function(_,texs){
		for (var i = 0; i < 2; i++) {	
			texs[i] = new GL.Texture2D({
				width : lightRes,
				height : lightRes,
				type : gl.FLOAT,
				minFilter : gl.NEAREST,
				magFilter : gl.NEAREST
			});
		}
	});
	fbo = new GL.Framebuffer({
		width : lightRes, 
		height : lightRes
	});

	var iterateLightShaderVertexCode = $('#iterate-light-vsh').text();
	var iterateLightShaderFragmentCode = 'precision mediump float;\n' + shaderDefs + $('#iterate-light-fsh').text();
	
	var iterateLightShaderUniforms = {
		posTex : 0,
		velTex : 1
	};

	iterateLightPosShader = new GL.ShaderProgram({
		vertexCode : iterateLightShaderVertexCode,
		fragmentCode : iterateLightShaderFragmentCode.replace('$assign', 'gl_FragColor = rel;'),
		uniforms : iterateLightShaderUniforms
	});
		
	iterateLightVelShader = new GL.ShaderProgram({
		vertexCode : iterateLightShaderVertexCode,
		fragmentCode : iterateLightShaderFragmentCode.replace('$assign', 'gl_FragColor = relDiff;'),
		uniforms : iterateLightShaderUniforms
	});
	
/*
how it'll work ...
1) start with the normalized world vector
2) iterate geodesic (backwards) a few times
3) use the resulting coordinates

you could do the iteration in the shader loop, or you could use a fbo to store state information ...
*/
	drawLightShader = new GL.ShaderProgram({
		vertexCodeID : 'draw-light-vsh',
		fragmentCode : 'precision mediump float;\n' + shaderDefs + $('#draw-light-fsh').text(),
		uniforms : {
			posTex : 0,
			cubeTex : 1
		}
	});
	
	gl.clearColor(.3, .3, .3, 1);

	viewVtx = new Float32Array(12);
	$.each(uvs, function(i,uv) {
		viewVtx[0+3*i] = uv[0]*2-1;
		viewVtx[1+3*i] = uv[1]*2-1;
		viewVtx[2+3*i] = -1;
	});
	viewVtxBuf = new GL.ArrayBuffer({
		data : viewVtx,
		usage : gl.DYNAMIC_DRAW
	});
	viewQuad = new GL.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		vertexBuffer : viewVtxBuf,
		hidden : true	//not a part of the scene, used for manually drawing
	});

	var unitVtx = new Float32Array(8);
	$.each(uvs, function(i,uv) {
		unitVtx[0+2*i] = uv[0];
		unitVtx[1+2*i] = uv[1];
	});
	unitQuad = new GL.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		vertexBuffer : new GL.ArrayBuffer({data : unitVtx, dim : 2}),
		hidden : true	//not a part of the scene, used for manually drawing
	});


	onresize();
	update();

/* debugging * /

	var cubeShader = new GL.ShaderProgram({
		vertexCodeID : 'cube-vsh',
		fragmentCodeID : 'cube-fsh',
		uniforms : {
			tex : 0
		}
	});

	cubeArrayBuf = new Float32Array(3*8);
	for (var i = 0; i < 8; i++) {
		cubeArrayBuf[3*i+0] = 2*(i&1)-1;
		cubeArrayBuf[3*i+1] = 2*((i>>1)&1)-1;
		cubeArrayBuf[3*i+2] = 2*((i>>2)&1)-1;
	}

	cubeTriArrayBuf = new Float32Array(3*3*2*2*3);
	for (var d = 0; d < 3; d++) {	//dimension
		for (var s = 0; s < 2; s++) {	//side
			for (var f = 0; f < 2; f++) {	//face tri 0/1 per side
				for (var i = 0; i < 3; i++) {	//index in tri 0-2
					var j = i + f * 2;	//index in face quad 0-3
					j &= 3;
					if (s) j = 4 - j;	//spin plus faces
					j &= 3;
					j ^= (j & 2) >> 1; //switch from tri strip to quad indexes
					var e = 0;	//cube vtx index
					if (s) e |= 1 << d;	//if plus side, set the dim bit
					//j bits coincide with d+1 and d+2 (mod 3) bits
					var d1 = (d + 1) % 3;
					var d2 = (d + 2) % 3;
					e |= (j & 1) << d1;
					e |= ((j >> 1) & 1) << d2;
					for (var v = 0; v < 3; v++) {	//vtx dim (xyz)
						cubeTriArrayBuf[v+3*(i+3*(f+2*(s+2*(d))))] = cubeArrayBuf[v+3*e];
					}
				}
			}
		}
	}

	cubeVtxBuf = new GL.ArrayBuffer({
		data : cubeTriArrayBuf
	});


	cubeObj = new GL.SceneObject({
		mode : gl.TRIANGLES,
		vertexBuffer : cubeVtxBuf,
		shader : cubeShader,
		texs : [skyTex],
		static : false
	});

	onresize();
	var draw = function() {
		gl.clear(gl.COLOR_BUFFER_BIT);
		unitQuad.draw({
			shader : cubeShader,
			texs : [skyTex],
			uniforms : {
				projMat : unitOrthoMat,
				mvMat : identMat
			}
		});
		requestAnimFrame(draw);
	};
	draw();
	/**/

	var tmpRotMat = mat4.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			mat4.identity(tmpRotMat);
			mat4.rotate(tmpRotMat, tmpRotMat, 
				Math.PI / 180 * Math.sqrt(dx*dx + dy*dy),
				[dy, dx, 0]);
			//mat4.translate(mvMat, mvMat, [10*dx/canvas.width, -10*dy/canvas.height, 0]);
			mat4.mul(GL.mvMat, tmpRotMat, GL.mvMat);

		}
	});
});

