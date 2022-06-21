/*
used by gl-util-font and gl-util-kernel
pretty simple
*/

if (!GLUtil) throw "require gl-util.js before gl-util-unitquad.js";

GLUtil.prototype.oninit.push(function() {
	var glutil = this;

	//2D tri strip front facing unit quad
	this.unitQuadVertexes = new Float32Array([
		0,0,
		1,0,
		0,1,
		1,1
	]);

	this.unitQuadVertexBuffer = new this.ArrayBuffer({
		dim : 2,
		data : this.unitQuadVertexes,
		keep : true
	});

	this.unitQuadGeom = new this.Geometry({
		mode : this.context.TRIANGLE_STRIP,
		vertexes : this.unitQuadVertexBuffer
	});

	this.unitQuad = new this.SceneObject({
		mode : this.context.TRIANGLE_STRIP,
		attrs : {
			vertex : this.unitQuadVertexBuffer
		},
		parent : null,
		static : true
	});
});
