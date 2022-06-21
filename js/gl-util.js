/*
WebGL helper classes

current setup
GLUtil - main class holding it all 
	Context - holds a webgl-capable canvas
		also holds the View associated with this canvasRenderer
		also holds the Scene associated with this canvasRenderer
	View - holds the camera info / how to render the canvasRenderer
	Scene - holds all SceneGraph objects, including root
		currently also holds the mvMat associated with the View and projMat associated with the Context
	SceneObject - holds information pertaining to an object in the scene
		Geometry
		Texture is one of the following:
			Texture1D
			Texture2D
			TextureCube
		Shader

	Framebuffer
*/
GLUtil = makeClass(new function(){
	//in this block 'this' is the prototype / the class object
	var GLUtilPrototype = this;

	//what webgl names to search through
	this.webGLNames = ['webgl', 'experimental-webgl'];

	/*
	on-init callbacks to be called by each created Context.
	register these by plugins as follows: GL.oninit.push(function() { ... });
	*/
	this.oninit = [];

	/*
	this combines the GL context, render target, and viewport responsabilities
	so maybe it should be split up later
	
	create a new scene with associated canvas and view
	args:
		one of these:
			canvas = which canvas to use
			fullscreen = whether to create our own fullscreen canvas 
		canvasArgs = canvas.getContext arguments, including:
			premultipliedAlpha (default false)
			alpha (default false)
	*/
	this.init = function(args) {
		var glutil = this;
	
		this.canvas = args.canvas;
		if (args.fullscreen) {
			window.scrollTo(0,1);
			if (this.canvas === undefined) {
				this.canvas = $('<canvas>').prependTo(document.body).get(0);
			}
		
			$(this.canvas).css({
				left : 0,
				top : 0,
				position : 'absolute'
			});
			var resize = function() {
				glutil.canvas.width = window.innerWidth;
				glutil.canvas.height = window.innerHeight;
				glutil.resize();
			};
			$(window).bind('resize', resize);
			//also call resize after init is done
			setTimeout(resize, 0);
		}
		if (this.canvas === undefined) throw 'expected canvas or fullscreen';

		var canvasArgs = args.canvasArgs;
		if (canvasArgs === undefined) canvasArgs = {};
		
		/*
		this is supposed to save me from having to write 1's in the dest alpha channel
		to keep the dest image from being invisible
		but it is buggy in firefox and safari
		*/
		if (canvasArgs.premultipliedAlpha === undefined) canvasArgs.premultipliedAlpha = false;
		
		/*
		this is supposed to slow things down
		but it is also supposed to allow folks to take screenshots ...
		*/
		//args.preserveDrawingBuffer
		
		if (canvasArgs.alpha === undefined) canvasArgs.alpha = false;

		var gl = undefined;
		for (var i = 0; i < this.webGLNames.length; i++) {
			try {
				//console.log('trying to init gl context of type', this.webGLNames[i]);
				gl = this.canvas.getContext(this.webGLNames[i], canvasArgs);
			} catch (e) {
				//console.log('failed with exception', e);
			}
			if (gl) break;
		}
		if (gl === undefined) {
			throw "Couldn't initialize WebGL =(";
		}
	
		if (args !== undefined && args.debug) {
			gl = WebGLDebugUtils.makeDebugContext(gl);	
		}

		this.context = gl;

		//gather extensions
		gl.getExtension('OES_element_index_uint');
		gl.getExtension('OES_standard_derivatives');
		gl.getExtension('OES_texture_float');
		gl.getExtension('OES_texture_float_linear');

		//initialize variables based on the gl context object constants:

		this.wrapMap = {
			s : gl.TEXTURE_WRAP_S,
			t : gl.TEXTURE_WRAP_T
		};

		//detect precision used

		this.fragmentPrecision = 'precision mediump float;\n';
		this.vertexPrecision = '';

		var vtxhigh = gl.getShaderPrecisionFormat(gl.VERTEX_SHADER, gl.HIGH_FLOAT)
		if (vtxhigh.rangeMin !== 0 && vtxhigh.rangeMax !== 0 && vtxhigh.precision !== 0) {
			this.vertexPrecision = 'precision highp float;\n';
		}
		var fraghigh = gl.getShaderPrecisionFormat(gl.FRAGMENT_SHADER, gl.HIGH_FLOAT)
		if (fraghigh.rangeMin !== 0 && fraghigh.rangeMax !== 0 && fraghigh.precision !== 0) {
			this.fragmentPrecision = 'precision highp float;\n';
		}


		//define classes here, so they have access to the glutil





		/*
		view object, for matrix deduction
		*/
		var View = makeClass({
			/*
			args:
				zNear (optional)
				zFar
				fovY
				ortho
				pos
				angle
			*/
			init : function(args) {
				this.zNear = 1;
				this.zFar = 2000;
				this.fovY = 90;	// corresponding with 1:1 x:z
				this.ortho = false;
				this.pos = vec3.create();
				this.angle = quat.create();
				if (args !== undefined) {
					if (args.zNear !== undefined) this.zNear = args.zNear;
					if (args.zFar !== undefined) this.zFar = args.zFar;
					if (args.fovY !== undefined) this.fovY = args.fovY;
					if (args.ortho !== undefined) this.ortho == args.ortho;
					if (args.pos !== undefined) vec3.copy(this.pos, args.pos);
					if (args.angle !== undefined) quat.copy(this.angle, args.angle);
				}
			}
		});
		this.View = View;

		var SceneObject;
		var Scene = makeClass({
			init : function() {
				//traditional gl matrices
				this.projMat = mat4.create();
				this.mvMat = mat4.create();
			
				this.root = new SceneObject({
					scene : this,
					parent : undefined,
					geometry : undefined
				});
			},

			setupMatrices : (function(){
				var viewAngleInv = quat.create();
				var viewPosInv = vec3.create();
				return function() {
					quat.conjugate(viewAngleInv, glutil.view.angle);
					mat4.fromQuat(this.mvMat, viewAngleInv);
					vec3.negate(viewPosInv, glutil.view.pos);
					mat4.translate(this.mvMat, this.mvMat, viewPosInv);
				};
			})()
		});
		this.Scene = Scene;

		var Shader = makeClass({
			/*
			args:
				(appended in this order)
				code = shader code,
				id = the id of the DOM element containing the shader code
			*/
			init : function(args) {
				var code = '';
				if (args !== undefined && args.code !== undefined) {
					if (code === undefined) code = '';
					code += args.code;
				}
				if (args !== undefined && args.id !== undefined) {
					if (code === undefined) code = '';
					var src = $('#'+args.id);
					code += src.text();
				}
				if (code === undefined) throw "expected code or id";

				this.obj = gl.createShader(this.shaderType);
				gl.shaderSource(this.obj, code);
				gl.compileShader(this.obj);
				if (!gl.getShaderParameter(this.obj, gl.COMPILE_STATUS)) {
					//stupid grep for tablet aLogCat
					$.each(code.split('\n'), function(i,line) {
						console.log((i+1)+': '+line);
					});
					throw gl.getShaderInfoLog(this.obj);
				}
			}
		});
		this.Shader = Shader;

		var VertexShader = makeClass({
			super : Shader,
			shaderType : gl.VERTEX_SHADER
		});
		this.VertexShader = VertexShader;

		var FragmentShader = makeClass({
			super : Shader,
			shaderType : gl.FRAGMENT_SHADER
		});
		this.FragmentShader = FragmentShader;

		//returns an array of gl.uniform* functions to use with this uniform: 
		//0: used when multiple primitive values are passed
		//1: used when an array is passed
		//2: used when an array is passed as a matrix
		var getUniformSettersForGLType = function(gl, gltype) {
			switch (gltype) {
			case gl.FLOAT: 
				return {arg:gl.uniform1f, count:1};
			case gl.INT:
			case gl.BOOL:
			case gl.SAMPLER_2D:
			case gl.SAMPLER_CUBE: 
				return {arg:gl.uniform1i, count:1};
			case gl.FLOAT_VEC2: 
				return {arg:gl.uniform2f, count:2, vec:gl.uniform2fv};
			case gl.INT_VEC2:
			case gl.BOOL_VEC2:
				return {arg:gl.uniform2i, count:2, vec:gl.uniform2iv};
			case gl.FLOAT_VEC3: 
				return {arg:gl.uniform3f, count:3, vec:gl.uniform3fv};
			case gl.INT_VEC3:
			case gl.BOOL_VEC3:
				return {arg:gl.uniform3i, count:3, vec:gl.uniform3iv};
			case gl.FLOAT_VEC4: 
				return {arg:gl.uniform4f, count:4, vec:gl.uniform4fv};
			case gl.INT_VEC4:
			case gl.BOOL_VEC4:
				return {arg:gl.uniform4i, count:4, vec:gl.uniform4iv};
			case gl.FLOAT_MAT2:
				return {mat:gl.uniformMatrix2fv};
			case gl.FLOAT_MAT3:
				return {mat:gl.uniformMatrix3fv};
			case gl.FLOAT_MAT4:
				return {mat:gl.uniformMatrix4fv};
			}
		};

		var ArrayBuffer;
		var ShaderProgram = makeClass({
			/*
			args:
					one of the following:
				vertexShader = the VertexShader object to link with
				vertexCode = the vertex shader code
				vertexCodeID = the id of the DOM element containing the vertex shader code
					one of the following:
				fragmentShader = the FragmentShader object to link with
				fragmentCode = the fragment shader code
				fragmentCodeID = the id of the DOM element containing the fragment shader code
					and any of the following:
				vertexPrecision = set to 'best' for generating the best possible precision
				fragmentPrecision = set to 'best' for generating the best possible precision
				uniforms = a key-value map containing initial values of any uniforms
			*/
			init : function(args) {
				var thiz = this;
				this.vertexShader = args.vertexShader;
				if (this.vertexShader === undefined) {
					var vertexCode = args.vertexCode;
					if (args.vertexPrecision === 'best') {
						if (vertexCode === undefined) vertexCode = '';
						vertexCode = glutil.vertexPrecision + vertexCode;
					}
					this.vertexShader = new glutil.VertexShader({
						code : vertexCode,
						id : args.vertexCodeID
					});
				}
				if (!this.vertexShader) throw "expected vertexShader or vertexCode or vertexCodeID";

				this.fragmentShader = args.fragmentShader;
				if (this.fragmentShader === undefined) {
					var fragmentCode = args.fragmentCode;
					if (args.fragmentPrecision === 'best') {
						if (fragmentCode === undefined) fragmentCode = '';
						fragmentCode = glutil.fragmentPrecision + fragmentCode;
					}
					this.fragmentShader = new glutil.FragmentShader({
						code : fragmentCode,
						id : args.fragmentCodeID
					});
				}
				if (!this.fragmentShader) throw "expected fragmentShader or fragmentCode or fragmentCodeID";

				this.obj = gl.createProgram();
				gl.attachShader(this.obj, this.vertexShader.obj);
				gl.attachShader(this.obj, this.fragmentShader.obj);
				
				gl.linkProgram(this.obj);
				if (!gl.getProgramParameter(this.obj, gl.LINK_STATUS)) {
					//throw 'Link Error: '+gl.getShaderInfoLog(this.obj);	
					console.log('vertex code:');
					$.each((args.vertexCode || $('#'+args.vertexCodeID).text()).split('\n'), function(i,line) {
						console.log((i+1)+': '+line);
					});
					console.log('fragment code:');
					$.each((args.fragmentCode || $('#'+args.fragmentCodeID).text()).split('\n'), function(i,line) {
						console.log((i+1)+': '+line);
					});
					throw "Could not initialize shaders";
				}
				
				gl.useProgram(this.obj);
				
				this.uniforms = {};
				var maxUniforms = gl.getProgramParameter(this.obj, gl.ACTIVE_UNIFORMS);
				for (var i = 0; i < maxUniforms; i++) {
					var info = gl.getActiveUniform(this.obj, i);
					info.loc = gl.getUniformLocation(this.obj, info.name);
					info.setters = getUniformSettersForGLType(gl, info.type);
					this.uniforms[i] = info;
					this.uniforms[info.name] = info;
				}

				this.attrs = {};
				var maxAttrs = gl.getProgramParameter(this.obj, gl.ACTIVE_ATTRIBUTES);
				for (var i = 0; i < maxAttrs; i++) {
					var info = gl.getActiveAttrib(this.obj, i);
					info.loc = gl.getAttribLocation(this.obj, info.name);
					this.attrs[info.name] = info;
				}

				if (args.uniforms !== undefined) {
					this.setUniforms(args.uniforms);
				}

				gl.useProgram(null);
			},
			use : function() {
				gl.useProgram(this.obj);
				return this;
			},
			useNone : function() {
				gl.useProgram(null);
				return this;
			},
			setUniforms : function(uniforms) {
				for (k in uniforms) {
					this.setUniform(k, uniforms[k]);
				}
				return this;
			},
			/*
			type-detecting uniform setting
			currenly only handles unpacked arguments
			and currently calls everything through uniformf
			*/
			setUniform : function() {
				var name = arguments[0];
				var info = this.uniforms[name];
				if (info === undefined) return;	//throw?  but if a uniform isn't used it'll be removed, and its info will return null ... is this an error?
				var value = arguments[1];
				var isArray = typeof(value) == 'object';	//$.isArray(value);
				var setters = info.setters;
				var loc = info.loc;
				if (!isArray) {
					var setter = setters.arg;
					if (!setter) throw "failed to find non-array setter for uniform "+name;
					Array.prototype.splice.call(arguments, 0, 1, loc);
					if (arguments.length < setters.count) {
						throw 'setUniform('+name+') needed '+setters.count+' arguments';
					}
					setter.apply(gl, arguments);
				} else {
					if (setters.vec) {
						setters.vec.call(gl, loc, value);
					} else if (setters.mat) {
						setters.mat.call(gl, loc, false, value);
					} else {
						throw "failed to find array setter for uniform "+name;
					}
				}
			},
			setAttrs : function(attrs) {
				for (k in attrs) {
					this.setAttr(k, attrs[k]);
				}
			},
			setAttr : function(name, buffer) {
				var info = this.attrs[name];
				if (info === undefined) return;
				gl.enableVertexAttribArray(info.loc);
				// array buffer object, assume packed
				if (buffer.__proto__ === ArrayBuffer.prototype) {
					gl.bindBuffer(gl.ARRAY_BUFFER, buffer.obj);
					gl.vertexAttribPointer(info.loc, buffer.dim, gl.FLOAT, false, 0, 0);
				//table object, try to derive values
				} else {
					var attrInfo = buffer;
					buffer = assertExists(attrInfo, 'buffer');
					var size = attrInfo.size !== undefined ? attrInfo.size : buffer.dim;
					//TODO make underlying type modular, and store as a parameter of the buffer
					var type = gl.FLOAT;
					var normalized = attrInfo.normalized !== undefined ? attrInfo.normalized : false;
					var offset = attrInfo.offset !== undefined ? attrInfo.offset : 0;
					var stride = attrInfo.stride !== undefined ? attrInfo.stride : 0;
					gl.bindBuffer(gl.ARRAY_BUFFER, buffer.obj);
					gl.vertexAttribPointer(info.loc, size, type, normalized, stride, offset);
				}
			},
			removeAttrs : function(attrs) {
				for (k in attrs) {
					this.removeAttr(k);
				}
			},
			removeAttr : function(name) {
				var info = this.attrs[name];
				if (info === undefined) return;
				gl.disableVertexAttribArray(info.loc);
			}
		});
		this.ShaderProgram = ShaderProgram;

		var Texture = makeClass({
			/*
			args:
				glutil (optional)
				everything else handled by setArgs
			*/
			init : function(args) {
				this.obj = gl.createTexture();
				gl.bindTexture(this.target, this.obj);
				if (args !== undefined) this.setArgs(args);
				gl.bindTexture(this.target, null);
			},
			//target provided upon init 
			bind : function(unit) { 
				if (unit !== undefined) gl.activeTexture(gl.TEXTURE0 + unit);
				gl.bindTexture(this.target, this.obj);
				return this;
			},
			unbind : function(unit) { 
				if (unit !== undefined) gl.activeTexture(gl.TEXTURE0 + unit);
				gl.bindTexture(this.target, null); 
				return this;
			},
			/*
			args:
				flipY: use UNPACK_FLIP_Y_WEBGL
				dontPremultiplyAlpha: don't use UNPACK_PREMULTIPLY_ALPHA_WEBGL
				magFilter: texParameter TEXTURE_MAG_FILTER
				minFilter: texParameter TEXTURE_MIN_FILTER
				generateMipmap: specifies to call generateMipmap after loading the texture data
				url: url of image to load
				wrap : wrap info
				onload: any onload callback to be used with url
			*/
			setArgs : function(args) {
				var target = this.target;
				if (args.alignment) gl.pixelStorei(gl.UNPACK_ALIGNMENT, args.alignment);
				if (args.flipY === true) gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
				else if (args.flipY === false) gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
				if (!args.dontPremultiplyAlpha) gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
				if (args.magFilter) gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, args.magFilter);
				if (args.minFilter) gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, args.minFilter);
				if (args.wrap) {
					this.setWrap(args.wrap);
				}
				this.setData(args);
				return this;
			},
			setWrap : function(args) {
				for (var k in args) {
					gl.texParameteri(this.target, glutil.wrapMap[k] || k, args[k]);
				}
				return this;
			},
			//typically overwritten. default calls setImage if args.data is provided
			setData : function(args) {
				if (args.data) {
					this.setImage(args);
				}
			}
		});
		this.Texture = Texture;

		/*
		args match Texture2D.setArgs
		with the exception of:
			onload : called when the image loads, after the data is set to the texture
			onerror : called if the image errors
		*/
		var Texture2D = makeClass({
			super : Texture,
			target : gl.TEXTURE_2D,
			setData : function(args) {
				if (args.url) {
					var image = new Image();
					var thiz = this;
					image.onload = function() {
						args.data = image;
						gl.bindTexture(thiz.target, thiz.obj);
						thiz.setImage(args);
						gl.bindTexture(thiz.target, null);
						
						if (args.onload) args.onload.call(thiz, args.url, image);
					};
					image.onerror = args.onerror;
					image.src = args.url;
				} else {
					if (args.data === undefined) args.data = null;
					this.setImage(args);
				}
			},
			/*
			args:
				target (default this.target)
				level (default 0)
				internalFormat (default gl.RGBA)
				width (optional)
				height (optional)
				border (optional, default 0 if needed)
				format (default gl.RGBA)
				type (default gl.UNSIGNED_BYTE)
				data (required - ArrayBufferView, ImageData, HTMLImageElement, HTMLCanvasElement, HTMLVideoElement)
					- though only ArrayBufferView requires width, height, and optionally border
			*/
			setImage : function(args) {
				var target = args.target !== undefined ? args.target : this.target;
				var level = args.level !== undefined ? args.level : 0;
				var internalFormat = args.internalFormat !== undefined ? args.internalFormat : gl.RGBA;
				var format = args.format !== undefined ? args.format : gl.RGBA;
				var type = args.type !== undefined ? args.type : gl.UNSIGNED_BYTE;
				var width = args.width;
				var height = args.height;
				var border = args.border !== undefined ? args.border : 0;
				
				//store?  WebGL has no getTexParameteri(gl.TEXTURE_WIDTH) ...
				this.internalFormat = internalFormat;
				this.format = format;
				this.type = type;
				this.width = width;
				this.height = height;
				
				//NOTICE this method only works for ArrayBufferView.  maybe that should be my test
				//console.log('setting image target',target,'level',level,'internalFormat',internalFormat,'width',width,'height',height,'border',border,'format',format,'type',type,'data',args.data);
				if (width === undefined && height === undefined) {
					//assume it's an image
					gl.texImage2D(target, level, internalFormat, format, type, args.data);
				} else {
					if (typeof(args.data) != 'function') {
						//assume it's a buffer
						gl.texImage2D(target, level, internalFormat, width, height, border, format, type, args.data);
					} else {
						//procedural generation
						var i = 0;
						
						//TODO get number of channels for format, rather than overriding it...
						format = gl.RGBA;
						var channels = 4;

						var scale = undefined;
						var data = undefined;
						if (type == gl.UNSIGNED_BYTE) {
							data = new Uint8Array(width * height * channels);
							scale = 255;
						} else if (type == gl.FLOAT) {
							data = new Float32Array(width * height * channels);
							scale = 1;
						}
				
						for (var y = 0; y < height; ++y) {
							for (var x = 0; x < width; ++x) {
								var c = args.data(x,y);
								for (var ch = 0; ch < 4; ++ch) {
									var d = c[ch];
									if (d === undefined) d = 0;
									d *= scale;
									data[i] = d;
									++i;
								}
							}
						}
						gl.texImage2D(target, level, internalFormat, width, height, border, format, type, data);
					}
				}
				if (args.generateMipmap) {
					gl.generateMipmap(this.target);
				}
			}
		});
		this.Texture2D = Texture2D;
		
		var TextureCube = makeClass({
			super : Texture,
			target : gl.TEXTURE_CUBE_MAP,
			getTargetForSide : function(side) {	//static
				return gl.TEXTURE_CUBE_MAP_POSITIVE_X + side;
			},
			setArgs : function(args) {
				Texture.prototype.setArgs.call(this, args);
				if (args.urls) {
					var loadedCount = 0;
					//store 'generateMipmap' up front.  we can't set it per-loaded-face, we have to as a whole.
					var generateMipmap = args.generateMipmap;
					args.generateMipmap = undefined;
					var thiz = this;
					$.each(args.urls, function(side,url) {
						var image = new Image();
						image.onload = function() {
							args.data = image;
							args.target = thiz.getTargetForSide(side);
							gl.bindTexture(thiz.target, thiz.obj);
							Texture2D.prototype.setImage.call(thiz, args);
							gl.bindTexture(thiz.target, null);
						
							if (args.onload) args.onload.call(thiz,side,url,image);
						
							//provide an overall all-sides-loaded callback
							//TODO make it generic?  add 'done' and 'onload' to Texture2D too?
							loadedCount++;
							if (loadedCount == 6) {
								if (generateMipmap) {
									gl.bindTexture(gl.TEXTURE_CUBE_MAP, thiz.obj);
									gl.generateMipmap(gl.TEXTURE_CUBE_MAP);
									gl.bindTexture(gl.TEXTURE_CUBE_MAP, null);
								}
								if (args.done) args.done.call(thiz);
							}
						};
						image.src = url;
					});
				} 
			},
			setData : function(args) {
				if (args.data === undefined) return;
				
				//store 'generateMipmap' up front.  we can't set it per-loaded-face, we have to as a whole.
				var generateMipmap = args.generateMipmap;
				args.generateMipmap = undefined;
		
				gl.bindTexture(this.target, this.obj);
				var isArray = typeof(args.data) == 'object';	//$.isArray(value);
				//console.log('isArray?',isArray);
				if (isArray && args.data.length >= 6) {
					var srcdata = args.data;
					for (var side = 0; side < 6; ++side) {
						args.data = srcdata[side];
						args.target = this.getTargetForSide(side);
						//console.log('setting target',args.target,' to data ',args.data);
						Texture2D.prototype.setImage.call(this, args);
					}
				} else if (typeof(args.data) == 'function') {
					var srcdata = args.data;
					for (var side = 0; side < 6; ++side) {
						args.data = function(x,y) {
							return srcdata(x,y,side);
						};
						args.target = this.getTargetForSide(side);
						Texture2D.prototype.setImage.call(this, args);
					}
				}
				
				if (generateMipmap) {
					//console.log('generating mipmaps of data-driven cubemap');
					gl.generateMipmap(gl.TEXTURE_CUBE_MAP);
				}
				
				gl.bindTexture(this.target, null);
		
			}
		});
		this.TextureCube = TextureCube;

		ArrayBuffer = makeClass({
			/*
			args:
				one of the two:
					data = either a Float32Array object, or a constructor for a Float32Array object
					count = how many vertexes to create
				usage = gl.bufferData usage
				dim = dimension / # elements per vector in data. only used for attrs and calculating length. default 3
				keep = optional, default true, set to false to not retain data in .data
			*/
			init : function(args) {
				if (args.keep === undefined) args.keep = true;
				this.obj = gl.createBuffer();
				this.dim = args.dim !== undefined ? args.dim : 3;
				var data = args.data;
				if (data === undefined) {
					if (args.count !== undefined) {
						data = new Float32Array(this.dim * args.count);
					} else {
						throw "expected 'data' or 'count'";
					}
				}
				this.setData(data, args.usage || gl.STATIC_DRAW, args.keep);
			},
			setData : function(data, usage, keep) {
				if (data.constructor != Float32Array) {
					data = new Float32Array(data);
				}
				gl.bindBuffer(gl.ARRAY_BUFFER, this.obj);
				gl.bufferData(gl.ARRAY_BUFFER, data, usage);
				gl.bindBuffer(gl.ARRAY_BUFFER, null);
				this.count = data.length / this.dim;
				if (keep) this.data = data;
			},
			updateData : function(data, offset) {
				if (offset === undefined) offset = 0;
				if (data === undefined) data = this.data;
				if (data.constructor != Float32Array) {
					data = new Float32Array(data);
				}
				gl.bindBuffer(gl.ARRAY_BUFFER, this.obj);
				gl.bufferSubData(gl.ARRAY_BUFFER, offset, data);
				gl.bindBuffer(gl.ARRAY_BUFFER, null);
			}
		});
		this.ArrayBuffer = ArrayBuffer;
		
		var ElementArrayBuffer = makeClass({
			/*
			args:
				data = either a Uint16Array object, or a constructor for a Uint16Array object
			*/
			init : function(args) {
				this.obj = gl.createBuffer();
				this.setData(args.data, args.usage || gl.STATIC_DRAW, args.keep);
			},
			setData : function(data, usage, keep) {
				if (data.constructor != Uint16Array && 
					data.constructor != Uint32Array) 
				{
					//in case of uint, default to uint
					// otherwise default to ushort
					var type = Uint16Array;
					if (gl.getExtension('OES_element_index_uint')) type = Uint32Array;
					data = new type(data);
				}
				this.type = data.constructor == Uint32Array ? gl.UNSIGNED_INT : gl.UNSIGNED_SHORT;
				
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.obj);
				gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, data, usage);
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
				this.count = data.length;

				if (keep) {
					this.data = data;
				}
			},
			updateData : function(data, offset) {
				if (offset === undefined) offset = 0;
				if (data.constructor != Uint16Array) {
					data = new Uint16Array(data);
				}
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.obj);
				gl.bufferSubData(gl.ELEMENT_ARRAY_BUFFER, offset, data);
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
			}
		});
		this.ElementArrayBuffer = ElementArrayBuffer;

		var Framebuffer = makeClass({
			/*
			args:
				width : framebuffer width.  required with depth.
				height : framebuffer height.  required with depth.
				useDepth : set to create a depth renderbuffer for this framebuffer.
			*/
			init : function(args) {
				this.width = args && args.width;
				this.height = args && args.height;
				this.obj = gl.createFramebuffer();
				gl.bindFramebuffer(gl.FRAMEBUFFER, this.obj);
				if (args !== undefined && args.useDepth) {
					this.depthObj = gl.createRenderbuffer();
					gl.bindRenderbuffer(gl.RENDERBUFFER, this.depthObj);
					gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, this.width, this.height);
					gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, this.depthObj);
					gl.bindRenderbuffer(gl.RENDERBUFFER, null);
				}
				gl.bindFramebuffer(gl.FRAMEBUFFER, null);
			},
			bind : function() {
				gl.bindFramebuffer(gl.FRAMEBUFFER, this.obj);
			},
			unbind : function() {
				gl.bindFramebuffer(gl.FRAMEBUFFER, null);
			},
			fboErrors : [
				'FRAMEBUFFER_INCOMPLETE_ATTACHMENT',
				'FRAMEBUFFER_INCOMPLETE_DIMENSIONS',
				'FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER',
				'FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT'
			],
			check : function() {
				var status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
				if (status != gl.FRAMEBUFFER_COMPLETE) {
					var errstr = 'glCheckFramebufferStatus status=' + status;
					$.each(this.fboErrors, function(i,fboError) {
						if (gl[fboError] == status) {
							errstr += ' error=' + fboError;
							return true;	//break;
						}
					});
					throw errstr;
				}
			},
			setColorAttachmentTex2D : function(index, tex, target, level) {
				if (index === undefined) index = 0;
				if (target === undefined) target = gl.TEXTURE_2D;
				if (level === undefined) level = 0;
				this.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + index, target, tex.obj, level);
				this.unbind();
			},
			setColorAttachmentTexCubeMapSide : function(index, tex, side, level) {
				if (side === undefined) side = index;
				if (level === undefined) level = 0;
				this.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + index, this.TextureCube.prototype.getTargetForSide(side), tex, level);
				this.unbind();
			},
	/* WebGL only supports one color attachment at a time ...
			setColorAttachmentTexCubeMap : function(tex, level) {
				this.bind();
				for (var i = 0; i < 6; i++) {
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + i, this.TextureCube.prototype.getTargetForSide(i), tex, level || 0);
				}
				this.unbind();
			},
	*/
	/* only available by extension ...
	function FBO:setColorAttachmentTex3D(index, tex, slice, target, level)
		if not tonumber(slice) then error("unable to convert slice to number: " ..tostring(slice)) end
		slice = tonumber(slice)
		self:bind()
		gl.glFramebufferTexture3D(gl.GL_FRAMEBUFFER, gl.GL_COLOR_ATTACHMENT0 + index, target or gl.GL_TEXTURE_3D, tex, level or 0, slice)
		self:unbind()
	end
	*/
			//general, object-based type-deducing
			setColorAttachment : function(index, tex) {
				if (typeof(tex) == 'object') {
					if (tex.constructor === Texture2D) {
						//javascript splice won't work, so array-clone it first or whatever needs to be done
						this.setColorAttachmentTex2D(index, tex.obj)	//, arguments.splice(2));
					// cube map? side or all at once?
					//elseif mt == Tex3D then
					//	self:setColorAttachmentTex3D(index, tex.id, ...)
					} else if (tex.constructor == WebGLTexture) {
						this.setColorAttachmentTex2D(index, tex);	// though this could be a 3d slice or a cube side...
					} else {
						throw "Can't deduce how to attach the object.  Try using an explicit attachment method";
					}
				} else {
					throw "Can't deduce how to attach the object.  Try using an explicit attachment method";
				}
			},
			
			/*
			if index is a number then it binds the associated color attachment at 'GL_COLOR_ATTACHMENT0+index' and runs the callback
			if index is a table then it runs through the ipairs,
				binding the associated color attachment at 'GL_COLOR_ATTACHMENT0+index[side]'
				and running the callback for each, passing the side as a parameter
			*/
			drawToCallback : function(callback/*, index*/) {
				this.bind();
				this.check();
				//no need to preserve the previous draw buffer in webgl
				//simply binding a framebuffer changes the render target to it
				callback();	
				this.unbind();
			},

			/*
			args:
				viewport
				shader
				uniforms
				texs
				callback
			*/
			draw : function(args) {
				var oldvp;
				if (args.viewport) {
					var vp = args.viewport;
					var oldvp = gl.getParameter(gl.VIEWPORT);
					gl.viewport.apply(gl, args.viewport);
				}
				//if (args.resetProjection) throw 'not supported in webgl';
				
				//something to consider:
				//the draw callback will most likely need a shader to bind its vertex attribute to
				//the fbo itself doesn't necessarily need one, nor does it store uniforms
				//so args.shader, args.uniforms, and args.texs might be moot here
				if (args.shader) {
					gl.useProgram(args.shader.obj);
					if (args.uniforms) {
						if (args.uniforms) {
							args.shader.setUniforms(args.uniforms);
						}
					}
				}
				if (args.texs) bindTextureSet(gl, args.texs);

				//if (args.color) throw 'color not supported in webgl';
				//if (args.dest) throw 'multiple color attachments not supported in webgl';
				
				// no one seems to use fbo:draw... at all...
				// so why preserve a function that no one uses?
				// why not just merge it in here?
				this.drawToCallback(args.callback/* || drawScreenQuad, args.colorAttachment || 0*/);
				
				if (args.texs) unbindTextureSet(gl, args.texs);
				if (args.shader) {
					gl.useProgram(null);
				}

				if (args.viewport) {
					gl.viewport.apply(gl, oldvp);
				}
			}
		});
		this.Framebuffer = Framebuffer;

		var bindTextureSet = function(gl, texs) {
			for (var k in texs) {
				if (texs.hasOwnProperty(k)) {
					gl.activeTexture(gl.TEXTURE0 + parseInt(k));
					var tex = texs[k];
					if (tex) {	//because java can enumerate through undefined values
						gl.bindTexture(tex.target, tex.obj);
					}
				}
			}
			gl.activeTexture(gl.TEXTURE0);
		};
		var unbindTextureSet = function(gl, texs) {
			for (var k in texs) {
				if (texs.hasOwnProperty(k)) {
					gl.activeTexture(gl.TEXTURE0 + parseInt(k));
					var tex = texs[k];
					if (tex) {	//because java can enumerate through undefined values
						gl.bindTexture(tex.target, null);
					}
				}
			}
			gl.activeTexture(gl.TEXTURE0);
		};

		var Attribute;
		var Geometry = makeClass({
			/*
			args:
				mode
				count (optional).  required unless 'indexes' or 'vertexes' is provided.
				indexes (optional).  specifies to use drawElements instead of drawArrays
				vertexes (optional).  either Attribute (holding ArrayBuffer) or ArrayBuffer.  solely used for providing 'count' when 'indexes' and 'count' is not used.
				offset (optional).  default 0
			*/
			init : function(args) {
				this.mode = args.mode;
				this.count = args.count;
				this.indexes = args.indexes;
				this.vertexes = args.vertexes;
				this.offset = args.offset !== undefined ? args.offset : 0;
			},
			/*
			args:
				mode : overrides mode
				count : overrides count
				offset : overrides offset
			*/
			draw : function(args) {
				var mode = this.mode;
				var count = this.count;
				var offset = this.offset;
				//allow overrides?  for which variables?
				if (args !== undefined) {
					if (args.mode !== undefined) mode = args.mode;
					if (args.count !== undefined) count = args.count;
					if (args.offset !== undefined) offset = args.offset;
				}
				if (this.indexes !== undefined) {
					gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.indexes.obj);
					if (count === undefined) {
						count = this.indexes.count;
					}
					gl.drawElements(mode, count, this.indexes.type, offset);
					gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
				} else {
					if (count === undefined && this.vertexes !== undefined) {
						if (this.vertexes.__proto__ !== Attribute.prototype) {
							count = this.vertexes.count;
						} else {
							count = this.vertexes.buffer.count;
						}
					}
					if (count > 0) {
						gl.drawArrays(mode, offset, count);
					}
				}
			}
		});
		this.Geometry = Geometry;

		Attribute = makeClass({
			/*
			args:
				buffer = ArrayBuffer object
				size = dimension of the buffer, default buffer.dim
				type = type of the buffer, default gl.FLOAT (soon to be buffer.type)
				normalize = whether to normalize the buffer 
				stride = stride of the buffer, default 0
				offset = offset of the buffer, default 0
			if args is a ArrayBuffer then it is treated as the buffer argument
			*/
			init : function(args) {
				if (args.__proto__ == ArrayBuffer.prototype) {
					this.buffer = args;
					this.size = this.buffer.dim;
					this.type = gl.FLOAT;
					this.normalize = false;
					this.stride = 0;
					this.offset = 0;
				} else {
					this.buffer = assertExists(args, 'buffer');
					this.size = args.size !== undefined ? args.size : this.buffer.dim;
					this.type = args.type !== undefined ? args.type : gl.FLOAT;	//TODO this.buffer.type
					this.normalize = args.normalize !== undefined ? args.normalize : false;
					this.stride = args.stride !== undefined ? args.stride : 0;
					this.offset = args.offset !== undefined ? args.offset : 0;
				}
			}
		});
		this.Attribute = Attribute;

		SceneObject = makeClass({
			/*
			args:
				scene

				geometry
					-or-
				mode
				count 	used to specify the number of elements to render.  not necessary if attrs.vertex is provided.
				index	(optional) used to specify drawElements instead of drawArrays	
				offset	offset into arrays to draw.  default: 0
				
				shader
				uniforms
				attrs:
					[attributeName] = glutil.Attribute object, or semblance of an Attribute object, or an ArrayBuffer.
									I was considering migrating this all to Attribute objects before I considered how many dereferences into the attrs table there are.  (Any time dynamic data is updated.)
									So until I switch everything over, I'll just code this, Shader.setAttr, and Geometry.draw to handle both.  Old code will run the same.  New code will run the new way.
					vertex = vertex attribute buffer.  used to override 'count'
				texs

				scenegraph / questionable vars:
					blend
					useDepth
					static
					parent
					pos
					angle
			*/
			init : function(args) {
				if (!args) args = {};
				this.scene = args.scene || glutil.scene;
				
				this.shader = args.shader;
				this.uniforms = args.uniforms || {};
				if (args.attrs) {
					this.attrs = {};
					for (var k in args.attrs) {
						var v = args.attrs[k];
						if (v.buffer !== undefined) {
							v = new Attribute(v);
						}
						this.attrs[k] = v;
					}
				}
				this.texs = args.texs;
				this.blend = args.blend;
				this.useDepth = args.useDepth;

				//TODO this should be bool-cast, and should probably be after the implicit assignment
				if ('static' in args) this.static = args.static;
				if (args.pos) {
					this.pos = vec3.clone(args.pos);
					this.static = false;
				}
				if (args.angle) {
					this.angle = quat.clone(args.angle);
					this.static = false;
				}

				if ('geometry' in args) {
					this.geometry = args.geometry;
				} else {
					this.geometry = new glutil.Geometry({
						mode : args.mode,
						count : args.count,
						offset : args.offset,
						indexes : args.indexes,
						vertexes : this.attrs !== undefined ? this.attrs.vertex : undefined
					});
				}

				if (!this.static) {
					if (this.pos === undefined) {
						this.pos = vec3.create();
					}

					if (this.angle === undefined) {
						this.angle = quat.create();
					}
				}
			
				this.children = [];
				
				if (args && 'parent' in args) {
					this.parent = args.parent;
				} else {
					this.parent = this.scene.root;
				}
				if (this.parent) {
					this.parent.children.push(this);
				}

				if (this.static) {
					if (this.parent) {
						this.targetMat = this.parent.targetMat;
					} else {
						this.targetMat = this.scene.mvMat;
					}
				} else {
					this.localMat = mat4.create();
					this.mvMat = mat4.create();
					this.targetMat = this.mvMat;
				}
			
				if (this.uniforms.projMat === undefined) this.uniforms.projMat = this.scene.projMat;
				if (this.uniforms.mvMat === undefined) this.uniforms.mvMat = this.targetMat;
				//hack to undo the model part of modelview.  TODO instead separate model and view ...
				if (this.localMat && this.uniforms.localMat === undefined) this.uniforms.localMat = this.localMat;
				if (this.uniforms.viewMatInv === undefined) this.uniforms.viewMatInv = this.scene.mvMat;
			},
		
			static : true, //default
			
			setupMatrices : function() {
				//TODO make matrix stuff optional?
				if (!this.static) {
					mat4.fromRotationTranslation(this.localMat, this.angle, this.pos);
					if (this.parent) {
						mat4.multiply(this.mvMat, this.parent.targetMat, this.localMat);
					} else {
						mat4.multiply(this.mvMat, this.scene.mvMat, this.localMat);
					}
				}
			},
			
			/*
			args: all optional and all overrides for args of constructor and shader constructor
				shader
				uniforms
				attrs
				mode
				count
				offset
			*/
			draw : function(args) {
				this.setupMatrices();

				//TODO push attrib anyone?
			
				var blend = this.blend || (args && args.blend);
				if (blend) {
					gl.blendFunc.apply(gl, blend);
					gl.enable(gl.BLEND);
				}

				if (this.useDepth === true) {
					gl.enable(gl.DEPTH_TEST);
				} else if (this.useDepth === false) {
					gl.disable(gl.DEPTH_TEST);
				}

				if (this.texs) bindTextureSet(gl, this.texs);
				if (args && args.texs) bindTextureSet(gl, args.texs);

				var shader = this.shader;
				if (args && args.shader) shader = args.shader;
				
				if (shader) {
					gl.useProgram(shader.obj);

					if (this.uniforms) shader.setUniforms(this.uniforms);
					if (args && args.uniforms) shader.setUniforms(args.uniforms);
					
					if (this.attrs) shader.setAttrs(this.attrs);
					if (args && args.attrs) shader.setAttrs(args.attrs);
				}
				
				if (this.geometry) {
					this.geometry.draw(args);
				}
				
				//nest within state & shader binding so children can inherit
				// they can also screw up state, mind you
				for (var i = 0; i < this.children.length; i++) {
					var child = this.children[i];
					if (!child.hidden) {
						child.draw();
					}
				}

				if (shader) {
					if (this.attrs) shader.removeAttrs(this.attrs);
					if (args && args.attrs) shader.removeAttrs(args.attrs);
					
					gl.useProgram(null);
				}
				
				if (args && args.texs) unbindTextureSet(gl, args.texs);
				if (this.texs) unbindTextureSet(gl, this.texs);
		
				if (blend) {
					gl.disable(gl.BLEND);
				}
			},
			
			remove : function() {
				if (!this.parent) return;
				this.parent.children.remove(this);
				this.parent = undefined;
			},
			
			appendTo : function(parent) {
				this.remove();
				this.parent = parent;
				this.parent.children.push(this);
			},
			
			prependTo : function(parent) {
				this.remove();
				this.parent = parent;
				this.parent.children.splice(0, 0, this);
			}
		});
		this.SceneObject = SceneObject;


		//create objects:

		//camera
		this.view = new this.View();

		//scenegraph
		this.scene = new this.Scene();

		//call on-init callbacks

		$.each(this.oninit, function(k,v) {
			v.call(glutil);
		});
	};

	var frames = 0;
	var lastTime = Date.now();
	this.draw = function() {
		if (this.onfps) {
			frames++;
			var thisTime = Date.now();
			if (thisTime - lastTime > 1000) {
				var fps = frames * 1000 / (thisTime - lastTime);
				if (this.onfps) this.onfps(fps);
				frames = 0;
				lastTime = thisTime;	
			}
		}

		this.scene.setupMatrices();

		//TODO modular?
		this.context.clear(this.context.COLOR_BUFFER_BIT | this.context.DEPTH_BUFFER_BIT);
	
		if (!this.scene.root.hidden) this.scene.root.draw();

		if (this.ondraw) this.ondraw();

		this.clearAlpha();
	};

	this.clearAlpha = function() {
		//work around canvas alpha crap
		this.context.colorMask(false,false,false,true);
		this.context.clear(this.context.COLOR_BUFFER_BIT);
		this.context.colorMask(true,true,true,true);
	};

	/*
	must be called manually 
	 (because it makes no assumptions of what the canvas should be resized to
	  or of whether the canvas resize callback fired before or after it did)
	*/
	this.resize = function() {
		this.context.viewport(0, 0, this.canvas.width, this.canvas.height);
		this.updateProjection();
		
		//auto draw on resize?
		//flag?
		//or leave it up to the caller?
		if (this.dontDrawOnResize) return;
		this.draw();
	};

	/*
	dir = unit direction result
	xf = fraction x coordinte
	fy = fraction y coordinate
	*/
	this.mouseDir = function(dir,xf,yf) {
		//basis: [0] = right, [1] = up, [2] = backwards
		var x = vec3.create(); vec3.quatXAxis(x, this.view.angle);
		var y = vec3.create(); vec3.quatYAxis(y, this.view.angle);
		var z = vec3.create(); vec3.quatZAxis(z, this.view.angle);
		var aspectRatio = this.canvas.width / this.canvas.height;
		var mxf = xf * 2 - 1;
		var myf = 1 - yf * 2;
		var tanFovY = Math.tan(this.view.fovY * Math.PI / 360);
		var px = this.view.pos[0];
		var py = this.view.pos[1];
		var pz = this.view.pos[2];
		dir[0] = -z[0] + tanFovY * (aspectRatio * mxf * x[0] + myf * y[0]);
		dir[1] = -z[1] + tanFovY * (aspectRatio * mxf * x[1] + myf * y[1]);
		dir[2] = -z[2] + tanFovY * (aspectRatio * mxf * x[2] + myf * y[2]);
	};

	/*
	must be manually called when any view projection matrix values change:
		aspectRatio, fovY, zNear, zFar
	*/
	this.updateProjection = function() {
		var projMat = this.scene.projMat;
		var aspectRatio = this.canvas.width / this.canvas.height;
		if (this.view.ortho) {
			var fovY = this.view.fovY;
			mat4.ortho(projMat,
				-aspectRatio * fovY,
				aspectRatio * fovY,
				-fovY,
				fovY,
				this.view.zNear,
				this.view.zFar);
		} else {
			var tanFovY = Math.tan(this.view.fovY * Math.PI / 360);
			mat4.frustum(projMat, 
				-aspectRatio * tanFovY * this.view.zNear, 
				aspectRatio * tanFovY * this.view.zNear, 
				-tanFovY * this.view.zNear, 
				tanFovY * this.view.zNear, 
				this.view.zNear, 
				this.view.zFar);
		}
	};

	//might require preserveDrawingBuffer ...
	this.screenshot = function() {
		/* download ... as a fixed-filename that can't be given an extension ... */
		var data = this.canvas.toDataURL('image/png');
		document.location.href = data.replace('image/png', 'image/octet');
		/**/
	
		/* download as a specified filename (by encoding in anchor element and simulating click) * /
		var mimeType = 'image/octet';
		var filename = 'download.png';
		var data = this.canvas.toDataURL(mimeType);
window.downloadData = data;
		var blob = new Blob([data], {type: mimeType});

		var downloadAnchor = document.createElement('a');
window.downloadAnchor = downloadAnchor;
		downloadAnchor.download = filename; 
		downloadAnchor.href = window.URL.createObjectURL(blob);
		downloadAnchor.textContent = 'Download Ready';

		downloadAnchor.dataset.downloadurl = [
			mimeType, 
			downloadAnchor.download, 
			downloadAnchor.href].join(':');
		downloadAnchor.dataset.disabled = false;

		document.body.appendChild(downloadAnchor);

		downloadAnchor.onclick = function(e) {
			if ('disabled' in this.dataset) {
				return false;
			}

			downloadAnchor.textContent = '';
			downloadAnchor.dataset.disabled = true;

			// Need a small delay for the revokeObjectURL to work properly.
			setTimeout(function() {
				window.URL.revokeObjectURL(downloadAnchor.href);
				//document.body.removeChild(downloadAnchor);
			}, 1500);
		};

		setTimeout(function() {
			$(downloadAnchor).trigger('click');
		}, 1000);
	*/
	};
});

vec3.quatXAxis = function(res, q) {
	var x = q[0], y = q[1], z = q[2], w = q[3];
	res[0] = 1 - 2 * (y * y + z * z); 
	res[1] = 2 * (x * y + z * w); 
	res[2] = 2 * (x * z - w * y); 
};

vec3.quatYAxis = function(res, q) {
	var x = q[0], y = q[1], z = q[2], w = q[3];
	res[0] = 2 * (x * y - w * z);
	res[1] = 1 - 2 * (x * x + z * z);
	res[2] = 2 * (y * z + w * x);
};

vec3.quatZAxis = function(res, q) {
	var x = q[0], y = q[1], z = q[2], w = q[3];
	res[0] = 2 * (x * z + w * y); 
	res[1] = 2 * (y * z - w * x); 
	res[2] = 1 - 2 * (x * x + y * y); 
};

