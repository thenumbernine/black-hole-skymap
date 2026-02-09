import {Canvas, Option} from '/js/dom.js';
import {mat4, quat} from '/js/gl-matrix-3.4.1/index.js';
import {assert, assertExists} from '/js/util.js';
import {clamp, getIDs, removeFromParent, preload, hide, show, hidden} from '/js/util.js';
import {GLUtil} from '/js/gl/gl.js';
import {Mouse3D} from '/js/mouse3d.js';
import {makeTextureCube} from '/js/gl/TextureCube.js';
import {makeGradient} from '/js/gl/Gradient.js';
import {makeGeodesicFBORenderer} from './renderer.fbo.js';
import {makeGeodesicSWRenderer} from './renderer.sw.js';
import {makeGeodesicTestCubeRenderer} from './renderer.testcube.js';

const ids = getIDs();
window.ids = ids;

const urlparams = new URLSearchParams(window.location.search);

// putting all these here cuz both the renderer and the UI callbacks want a named table
const _G = {};
window._G = _G;

let canvas;
let gl;
let glutil;

const objectTypes = [
	'Schwarzschild Black Hole',
	'Kerr Black Hole degeneracy',
	'Kerr Black Hole',
	'Alcubierre Warp Drive Bubble',
];
_G.objectTypes = objectTypes;
_G.objectType = objectTypes[0];
_G.objectDist = 10;
_G.blackHoleMass = .5;
_G.blackHoleCharge = 0.;
_G.blackHoleAngularVelocity = 0.;

_G.warpBubbleThickness = .01;
_G.warpBubbleVelocity = 1.5;
_G.warpBubbleRadius = 1;

_G.deltaLambda = .1;	//ray forward iteration
_G.simTime = 0;

const ident4 = mat4.create();

const tanh = x => {
	const exp2x = Math.exp(2 * x);
	return (exp2x - 1) / (exp2x + 1);
}

const sech = x => {
	const expx = Math.exp(x);
	return 2. * expx / (expx * expx + 1.);
}

const sechSq = x => {
	const y = sech(x);
	return y * y;
}

_G.shaderCommonCode = `
float sech(float x) {
	float expx = exp(x);
	return 2. * expx / (expx * expx + 1.);
}

float sechSq(float x) {
	float y = sech(x);
	return y * y;
}

vec3 quatRotate(vec4 q, vec3 v) {
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

vec4 quatConj(vec4 q) {
	return vec4(q.xyz, -q.w);
}
`;

const stupidPrint = s => {
/*
	s.split('\n').forEach(l => {
		console.log(l);
	});
*/
}

const SQRT_1_2 = Math.sqrt(.5);
//forward-transforming (object rotations)
_G.angleForSide = [
	[0, -SQRT_1_2, 0, -SQRT_1_2],
	[0, SQRT_1_2, 0, -SQRT_1_2],
	[SQRT_1_2, 0, 0, -SQRT_1_2],
	[SQRT_1_2, 0, 0, SQRT_1_2],
	[0, 0, 0, -1],
	[0, -1, 0, 0]
];




//names of all renderers
const skyboxRendererClassNames = [
	'GeodesicTestCubeRenderer',
	'GeodesicSWRenderer',
	'GeodesicFBORenderer',
];

let skyboxRenderer;
let skyboxRendererClassName;

//I would like to eventually instanciate all renderers and allow them to be toggled at runtime
//however courtesy of the scenegraph's globals (which I am not too happy about my current design), this will take a bit more work
//so in the mean time, this will take a page reset every time the glutil changes

const resize = () => {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	glutil.resize();

	const info = ids.info;
	const width = window.innerWidth
		- parseInt(info.style.paddingLeft)
		- parseInt(info.style.paddingRight);
	info.style.width = width+'px';
	const height = window.innerHeight
		- parseInt(info.style.paddingTop)
		- parseInt(info.style.paddingBottom);
	info.style.height = (height - 32)+'px';
}

// render loop

const update = () => {
	skyboxRenderer.update();
	requestAnimationFrame(update);
}

_G.mouseMethod = 'rotateCamera';
//_G.mouseMethod = 'rotateObject';

_G.drawMethod = 'background';

let mouse;

_G.objectAngle = quat.create();

const initAngle = [];
const initAngleInv = [];

const skyTexFilenames = [
	'skytex/sky-visible-cube-xp.png',
	'skytex/sky-visible-cube-xn.png',
	'skytex/sky-visible-cube-yp.png',
	'skytex/sky-visible-cube-yn.png',
	'skytex/sky-visible-cube-zp.png',
	'skytex/sky-visible-cube-zn.png',
];

ids.panelButton.addEventListener('click', e => {
	if (hidden(ids.panel)) {
		show(ids.panel);
		hide(ids.info);
	} else {
		hide(ids.panel);
	}
});
ids.infoButton.addEventListener('click', e => {
	if (hidden(ids.info)) {
		show(ids.info);
		hide(ids.panel);
	} else {
		hide(ids.info);
	}
});

canvas = Canvas({
	style : {
		left : 0,
		top : 0,
		position : 'absolute',
		userSelect : 'none',
	},
	prependTo : document.body,
});
window.canvas = canvas;

const objectTypeParamDivs = {};
const refreshObjectTypeParamDivs = () => {
	Object.entries(objectTypeParamDivs).forEach(entry => {
		const [divObjectType, objectTypeParamDiv] = entry;
		if (divObjectType == _G.objectType) {
			show(objectTypeParamDiv);
		} else {
			hide(objectTypeParamDiv);
		}
	});
};
objectTypes.forEach(v => {
	const id = v+' params';
	if (!(id in ids)) {
		console.log("couldn't find params for ", v, id);
		return;
	}
	objectTypeParamDivs[v] = ids[id];
	const option = Option({innerText:v, appendTo:ids.objectTypes});
	if (v == _G.objectType) {
		option.setAttribute('selected', 'true');
	}
});
ids.objectTypes.addEventListener('change', e => {
	_G.objectType = ids.objectTypes.value;
	refreshObjectTypeParamDivs();
	skyboxRenderer.resetField();
});
refreshObjectTypeParamDivs();

[
	'deltaLambda',
	'objectDist',
	'blackHoleMass',
	'blackHoleCharge',
	'blackHoleAngularVelocity',
	'warpBubbleThickness',
	'warpBubbleVelocity',
	'warpBubbleRadius'
].forEach(v => {
	const o = ids[v];
	o.value = _G[v];
	o.addEventListener('change', e => {
		_G[v] = o.value*1;
		o.blur();
	});
});

{
	skyboxRendererClassName = 'GeodesicFBORenderer';
	const classname = urlparams.get('renderer');
	if (classname) {
		skyboxRendererClassName = classname;
	}
	if (skyboxRendererClassNames.indexOf(skyboxRendererClassName) == -1) throw "unable to find skybox renderer named "+skyboxRendererClassName;
}

skyboxRendererClassNames.forEach(name => {
	if (!(name in ids)) {
		console.log("couldn't find radio for ", name);
		return;
	}
	const radio = ids[name];
	// TODO recmbine but with this param exchanged, like in 4d-renderer
	radio.addEventListener('click', e => {
		location.href = 'index.html?renderer=' + name;
	});
	if (name == skyboxRendererClassName) radio.checked = false;
});

try {
	glutil = new GLUtil({canvas:canvas});
	gl = glutil.context;
} catch (e) {
	removeFromParent(canvas);
	show(ids.webglfail);
	throw e;
}
show(ids.menu);
show(ids.panel);
glutil.import('Gradient', makeGradient);
const glMaxCubeMapTextureSize = gl.getParameter(gl.MAX_CUBE_MAP_TEXTURE_SIZE);
_G.glMaxCubeMapTextureSize = glMaxCubeMapTextureSize;

const hsvTex = new glutil.Gradient.HSVTexture(256);
hsvTex.bind();
gl.texParameteri(hsvTex.target, gl.TEXTURE_WRAP_S, gl.REPEAT);
hsvTex.unbind();
_G.hsvTex = hsvTex;

ids.reset.addEventListener('click', e => {
	skyboxRenderer.resetField();
});

ids.reset_view.addEventListener('click', e => {
	quat.copy(_G.objectAngle, quat.create());
	quat.copy(glutil.view.angle, initAngle);
	skyboxRenderer.resetField();
});

//classes themselves
const rendererClassGen = {
	GeodesicTestCubeRenderer : makeGeodesicTestCubeRenderer,
	GeodesicSWRenderer : makeGeodesicSWRenderer,
	GeodesicFBORenderer : makeGeodesicFBORenderer,
};
const classgen = assertExists(rendererClassGen, skyboxRendererClassName);

_G.glutil = glutil;
const cl = classgen(_G);
skyboxRenderer = new cl();

gl.disable(gl.DITHER);

preload(skyTexFilenames, () => {
	glutil.view.zNear = .1;
	glutil.view.zFar = 100;
	glutil.view.fovY = 90;
	quat.mul(glutil.view.angle, /*90' y*/[0,SQRT_1_2,0,SQRT_1_2], /*90' -x*/[-SQRT_1_2,0,0,SQRT_1_2]);
	quat.copy(initAngle, glutil.view.angle);
	quat.conjugate(initAngleInv, initAngle);

	const skyTex = new glutil.TextureCube({
		flipY : false,
		generateMipmap : true,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR_MIPMAP_LINEAR,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		},
		urls : skyTexFilenames,
		onload : (side,url,image) => {
			if (image.width > glMaxCubeMapTextureSize || image.height > glMaxCubeMapTextureSize) {
				throw "cube map size "+image.width+"x"+image.height+" cannot exceed "+glMaxCubeMapTextureSize;
			}
		},
		done : () => {
			skyboxRenderer.initScene(skyTex);
			document.querySelectorAll('input[name="mouseMethod"]').forEach(o => {
				o.addEventListener('click', e => {
					_G[o.name] = o.value;
				});
			});
			document.querySelectorAll('input[name="drawMethod"]').forEach(o => {
				o.addEventListener('click', e => {
					_G[o.name] = o.value;
				});
			});

			const tmpQ = quat.create();
			mouse = new Mouse3D({
				pressObj : canvas,
				move : (dx,dy) => {
					const rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
					quat.setAxisAngle(tmpQ, [dy, dx, 0], rotAngle);

					if (_G.mouseMethod == 'rotateCamera') {
						quat.conjugate(tmpQ, tmpQ);
						quat.mul(glutil.view.angle, tmpQ, glutil.view.angle);
						quat.normalize(glutil.view.angle, glutil.view.angle);
					} else if (_G.mouseMethod == 'rotateObject') {
						//rotate into view space
						quat.conjugate(tmpQ, tmpQ);
						quat.mul(tmpQ, tmpQ, initAngle);
						quat.mul(tmpQ, initAngleInv, tmpQ);

						quat.mul(_G.objectAngle, _G.objectAngle, tmpQ);
						quat.normalize(_G.objectAngle, _G.objectAngle);
						skyboxRenderer.resetField();
					}
				},
				zoom : dz => {
					glutil.view.fovY *= Math.exp(-.0003 * dz);
					glutil.view.fovY = clamp(glutil.view.fovY, 1, 179);
					glutil.updateProjection();
				}
			});

			const setRunning = e => { skyboxRenderer.runSimulation = ids.runSimulation.checked; };
			ids.runSimulation.addEventListener('click', setRunning);
			setRunning();

			skyboxRenderer.resetField();

			window.addEventListener('resize', resize);
			resize();
			update();
		},
	});
});
