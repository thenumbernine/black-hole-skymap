GeodesicTestCubeRenderer = makeClass({
	testInit : function() {},
	initScene : function(skyTex) {
		var cubeShader = new GL.ShaderProgram({
			vertexPrecision : 'best',
			vertexCode : mlstr(function(){/*
attribute vec3 vertex;
varying vec3 vertexv;
uniform mat4 projMat;

void main() {
	vertexv = vertex;
	gl_Position = projMat * vec4(vertex, 1.);
}
*/}),
			fragmentPrecision : 'best',
			fragmentCode : mlstr(function(){/*
varying vec3 vertexv;
uniform samplerCube skyTex;

uniform vec4 viewAngle;
vec3 quatRotate(vec4 q, vec3 v) { 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	vec3 dir = vertexv;
	dir = quatRotate(viewAngle, dir);
	gl_FragColor = textureCube(skyTex, dir);
	gl_FragColor.w = 1.; 
}
*/}),
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

		var cubeVtxBuf = new GL.ArrayBuffer({
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

		var cubeObj = new GL.SceneObject({
			mode : gl.TRIANGLES,
			attrs : {
				vertex : cubeVtxBuf
			},
			uniforms : {
				viewAngle : GL.canvasRenderer.view.angle
			},
			indexes : cubeIndexBuf,
			shader : cubeShader,
			texs : [skyTex],
			static : false
		});
	},

	resetField : function() {},

	update : function() {
		GL.canvasRenderer.draw();
	}
});
