function makeGeodesicTestCubeRenderer(_G) {
const glutil = _G.glutil;
const gl = glutil.context;
class GeodesicTestCubeRenderer {
	initScene(skyTex) {
		const cubeShader = new glutil.Program({
			vertexCode : `
in vec3 vertex;
out vec3 vertexv;
uniform mat4 projMat;

void main() {
	vertexv = vertex;
	gl_Position = projMat * vec4(vertex, 1.);
}
`,
			fragmentCode : `
in vec3 vertexv;
uniform samplerCube skyTex;
uniform vec4 viewAngle;

vec4 quatConj(vec4 q) {
	return vec4(-q.xyz, q.w);
}
vec3 quatRotate(vec4 q, vec3 v) { 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
out vec4 fragColor;
void main() {
	vec3 dir = vertexv;
	dir = quatRotate(quatConj(viewAngle), dir);
	fragColor = texture(skyTex, dir);
	fragColor.w = 1.; 
}
`,
			uniforms : {
				skyTex : 0
			}
		});

		const cubeVtxArray = new Float32Array(3*8);
		for (let i = 0; i < 8; i++) {
			cubeVtxArray[0+3*i] = 2*(i&1)-1;
			cubeVtxArray[1+3*i] = 2*((i>>1)&1)-1;
			cubeVtxArray[2+3*i] = 2*((i>>2)&1)-1;
		}

		const cubeVtxBuf = new glutil.ArrayBuffer({
			data : cubeVtxArray 
		});

		const cubeIndexBuf = new glutil.ElementArrayBuffer({
			data : [
				5,7,3,3,1,5,		// <- each value has the x,y,z in the 0,1,2 bits (off = 0, on = 1)
				6,4,0,0,2,6,
				2,3,7,7,6,2,
				4,5,1,1,0,4,
				6,7,5,5,4,6,
				0,1,3,3,2,0
			]
		});

		const cubeObj = new glutil.SceneObject({
			mode : gl.TRIANGLES,
			attrs : {
				vertex : cubeVtxBuf
			},
			uniforms : {
				viewAngle : glutil.view.angle
			},
			indexes : cubeIndexBuf,
			shader : cubeShader,
			texs : [skyTex],
			static : false
		});
	}

	resetField() {}

	update() {
		glutil.draw();
	}
}

return GeodesicTestCubeRenderer;
}
export {makeGeodesicTestCubeRenderer};
