if (!GLUtil) throw "require gl-util.js before gl-util-font.js";

GLUtil.prototype.oninit.push(function() {
	var glutil = this;
	
	/*
	args:
		context (optional)
		width
		colors
		dontRepeat
	*/
	var GradientTexture = makeClass({
		super : glutil.Texture2D,
		init : function(args) {
			this.context = args.context || glutil.context;
			var width = args.width;
			var colors = args.colors;
			var data = new Uint8Array(width * 3);
			for (var i = 0; i < width; i++) {
				var f = (i+.5)/width;
				if (args.dontRepeat) {
					f *= colors.length - 1;
				} else {
					f *= colors.length;
				}
				var ip = parseInt(f);
				f -= ip;
				var iq = (ip + 1) % colors.length;
				var g = 1. - f;	
				for (var k = 0; k < 3; k++) {
					data[k+3*i] = 255*(colors[ip][k] * g + colors[iq][k] * f);
				}
			}
			glutil.Texture2D.call(this, {
				context : this.context,
				width : width,
				height : 1,
				format : this.context.RGB,
				internalFormat : this.context.RGB,
				data : data,
				minFilter : this.context.NEAREST,
				magFilter : this.context.LINEAR,
				wrap : { s : this.context.CLAMP_TO_EDGE, t : this.context.CLAMP_TO_EDGE }
			});
		}
	});
	this.GradientTexture = GradientTexture;

	var HSVTexture = makeClass({
		super : glutil.Texture2D,
		init : function(width) {
			GradientTexture.call(this, {
				width : width,
				colors : [
					[1,0,0],
					[1,1,0],
					[0,1,0],
					[0,1,1],
					[0,0,1],
					[1,0,1]
				]
			});
		}
	});
	this.HSVTexture = HSVTexture;
});
