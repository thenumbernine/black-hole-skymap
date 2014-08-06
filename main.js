var glMaxCubeMapTextureSize;
var canvas;
var gl;
var glutil;

var objectTypes = ['Black Hole', 'Alcubierre Warp Drive Bubble'];
var objectType = objectTypes[0];
var objectDist = 10;
var blackHoleMass = 1;
var warpBubbleThickness = 1;
var warpBubbleVelocity = .5;
var warpBubbleRadius = 2;
var deltaLambda = .1;	//ray forward iteration

var ident4 = mat4.create();

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

var shaderCommonCode = mlstr(function(){/*
float tanh(float x) {
	float exp2x = exp(2. * x);
	return (exp2x - 1.) / (exp2x + 1.);
}

float sech(float x) {
	float expx = exp(x);
	return 2. * expx / (expx * expx + 1.);
}

float sechSq(float x) {
	float y = sech(x);
	return y * y;
}
*/});

function stupidPrint(s) {
	return;
	$.each(s.split('\n'), function(_,l) {
		console.log(l);
	});
}

var SQRT_1_2 = Math.sqrt(.5);
//forward-transforming (object rotations)
var angleForSide = [
	[0, -SQRT_1_2, 0, -SQRT_1_2],
	[0, SQRT_1_2, 0, -SQRT_1_2],
	[SQRT_1_2, 0, 0, -SQRT_1_2],
	[SQRT_1_2, 0, 0, SQRT_1_2],
	[0, 0, 0, -1],
	[0, -1, 0, 0]
];




//names of all renderers
var skyboxRendererClassNames = [
	'GeodesicTestCubeRenderer',
	'GeodesicSWRenderer',
	'GeodesicFBORenderer'
];

var skyboxRenderer;
var skyboxRendererClassName;

//I would like to eventually instanciate all renderers and allow them to be toggled at runtime
//however courtesy of the scenegraph's globals (which I am not too happy about my current design), this will take a bit more work
//so in the mean time, this will take a page reset every time the glutil changes

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	glutil.resize();

	var info = $('#info');
	var width = window.innerWidth 
		- parseInt(info.css('padding-left'))
		- parseInt(info.css('padding-right'));
	info.width(width);
	var height = window.innerHeight
		- parseInt(info.css('padding-top'))
		- parseInt(info.css('padding-bottom'));
	info.height(height - 32);
}

// render loop

function update() {
	skyboxRenderer.update();
	requestAnimFrame(update);
};

var mouse;
function main3(skyTex) {
	skyboxRenderer.initScene(skyTex);

	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [dy, dx, 0], rotAngle);

			quat.mul(glutil.view.angle, glutil.view.angle, tmpQ);
			quat.normalize(glutil.view.angle, glutil.view.angle);
		},
		zoom : function(dz) {
			glutil.view.fovY *= Math.exp(-.0003 * dz);
			glutil.view.fovY = Math.clamp(glutil.view.fovY, 1, 179);
			glutil.updateProjection();
		}
	});

	skyboxRenderer.resetField();

	$(window).resize(resize);
	resize();
	update();
}


var main2Initialized = false;
function main2() {
	if (main2Initialized) {
		console.log("main2 got called twice again.  check the preloader.");
		return;
	}
	main2Initialized = true; 
	
	glutil.view.zNear = .1;
	glutil.view.zFar = 100;
	glutil.view.fovY = 90;
	quat.mul(glutil.view.angle, /*90' x*/[SQRT_1_2,0,0,SQRT_1_2], /*90' -y*/[0,-SQRT_1_2,0,SQRT_1_2]);

	console.log('creating skyTex');
	var skyTex = new glutil.TextureCube({
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


var skyTexFilenames = [
	'skytex/sky-visible-cube-xp.png',
	'skytex/sky-visible-cube-xn.png',
	'skytex/sky-visible-cube-yp.png',
	'skytex/sky-visible-cube-yn.png',
	'skytex/sky-visible-cube-zp.png',
	'skytex/sky-visible-cube-zn.png'
];

function main1() {
	$('#panelButton').click(function() {
		var panel = $('#panel');	
		if (panel.css('display') == 'none') {
			panel.show();
			$('#info').hide();
		} else {
			panel.hide();
		}
	});
	$('#infoButton').click(function() {
		var info = $('#info');
		if (info.css('display') == 'none') {
			info.show();
			$('#panel').hide();
		} else {
			info.hide();
		}
	});
	
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(canvas).disableSelection()

	var objectTypeParamDivs = {};
	var refreshObjectTypeParamDivs = function() {
		$.each(objectTypeParamDivs, function(divObjectType,objectTypeParamDiv) {
			if (divObjectType == objectType) {
				objectTypeParamDiv.show();
			} else {
				objectTypeParamDiv.hide();
			}
		});
	};
	$.each(objectTypes, function(k,v) {
		objectTypeParamDivs[v] = $('#'+v.replace(new RegExp(' ', 'g'), '_')+'_params');
		var option = $('<option>', {text:v});
		option.appendTo($('#objectTypes'));
		if (v == objectType) {
			option.attr('selected', 'true');
		}
	});
	$('#objectTypes').change(function() {
		objectType = $('#objectTypes').val();
		refreshObjectTypeParamDivs();
		skyboxRenderer.resetField();
	});
	refreshObjectTypeParamDivs();

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

	skyboxRendererClassName = 'GeodesicFBORenderer';
	var classname = $.url().param('renderer');
	if (classname) {
		skyboxRendererClassName = classname;
	}
	if (skyboxRendererClassNames.indexOf(skyboxRendererClassName) == -1) throw "unable to find skybox renderer named "+skyboxRendererClassName;

	$.each(skyboxRendererClassNames, function(i,name) {
		var radio = $('#' + name);
		radio.click(function() {
			location.href = 'index.html?renderer=' + name;
		});
		if (name == skyboxRendererClassName) radio.attr('checked', 'checked');
	});

	try {
		glutil = new GLUtil({canvas:canvas});
		gl = glutil.context;
	} catch (e) {
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	$('#menu').show();
	
	glMaxCubeMapTextureSize = gl.getParameter(gl.MAX_CUBE_MAP_TEXTURE_SIZE);

	$('#reset').click(function() {
		skyboxRenderer.resetField();
	});

	skyboxRenderer = new (window[skyboxRendererClassName])(glutil);

	gl.disable(gl.DITHER);

	$(skyTexFilenames).preload(main2);
}

$(document).ready(main1);

