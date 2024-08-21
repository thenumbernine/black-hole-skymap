#!/usr/bin/env luajit

--[[
cube map render
shader that iterates along geodesic from the view plane backwards through time to infinity
then reprojects onto the cubemap

schwarzschild metric: ds^2 = -(1-R/r)dt^2 + 1/(1-R/r)dr^2 + r^2 (dθ^2 + sin(θ)^2 dφ^2)
where R = 2M is the schwarzschild radius

christoffel symbols:

conn^r_tt = R(r-R)/(2r^3)
conn^t_tr = -conn^r_rr = R/(2r(r-R))
conn^φ_rφ = conn^θ_rθ = 1/r
conn^r_φφ = -(r - R)
conn^θ_φθ = cos(φ)/sin(φ)
conn^r_θθ = -(r - R) sin(φ)^2
conn^φ_θθ = -sin(φ) cos(φ)

-t'' = R/(r(r-R)) t' r'
-r'' = R(r-R)/(2r^3) t'^2 - R/(2r(r-R)) r'^2 - (r - R) φ'^2 - (r - R) sin(φ)^2 θ'^2
-θ'' = 2/r r' θ' + 2 cos(φ)/sin(φ) φ' θ'
-φ'' = 2/r r' φ' - sin(φ) cos(φ) θ'^2
--]]

local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local glu = require 'ffi.req' 'glu'
local sdl = require 'ffi.req' 'sdl'
local ig = require 'imgui'
local vec2d = require 'vec-ffi.vec2d'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'
local matrix_ffi = require 'matrix.ffi'
local GLTex2D = require 'gl.tex2d'
local GLTexCube = require 'gl.texcube'
local GLProgram = require 'gl.program'
local FBO = require 'gl.fbo'
local glreport = require 'gl.report'

local skyTex
--local viewRot = quatd(0,math.sqrt(.5),0,math.sqrt(.5))*quatd(math.sqrt(.5),0,0,-math.sqrt(.5))
--local viewAngleAxis = quatd(0,0,0,1)
--local zNear = .1
--local zFar = 100
--local tanFov = 1
--local zoomFactor = .9
local leftButtonDown

local doIteration = 2
local lightInitialized = false
local drawLightShader
local initLightPosShader, initLightVelShader
local iterateLightPosShader, iterateLightVelShader
local lightPosTexs = {}
local lightVelTexs = {}
local fbo

-- notice that the order matches global 'sides'
cubeFaces = {
	{5,7,3,1};		-- <- each value has the x,y,z in the 0,1,2 bits (off = 0, on = 1)
	{6,4,0,2};
	{2,3,7,6};
	{4,5,1,0};
	{6,7,5,4};
	{0,1,3,2};
}

local uvs = {
	vec2d(0,0),
	vec2d(1,0),
	vec2d(1,1),
	vec2d(0,1),
}

-- _G for input
deltaLambdaPtr = 1
iterationsPtr = 1

local App = require 'imguiapp.withorbit'()
App.viewUseBuiltinMatrixMath = true
App.title = 'black hole raytracer'
App.viewDist = 20

function App:initGL()
	App.super.initGL(self)

	self.fboProjMat = matrix_ffi({4,4}, 'float'):zeros():setOrtho(0, 1, 0, 1, -1, 1)
	self.identMat = matrix_ffi({4,4}, 'float'):zeros():setIdent()

	self.view.angle =
		quatd(-math.sqrt(.5), 0, 0, -math.sqrt(.5))
		* quatd(0, -math.sqrt(.5), 0, math.sqrt(.5))

	skyTex = GLTexCube{
		filenames = {
			'../skytex/sky-infrared-cube-xp.png',
			'../skytex/sky-infrared-cube-xn.png',
			'../skytex/sky-infrared-cube-yp.png',
			'../skytex/sky-infrared-cube-yn.png',
			'../skytex/sky-infrared-cube-zp.png',
			'../skytex/sky-infrared-cube-zn.png',
		},
		wrap={
			s=gl.GL_CLAMP_TO_EDGE,
			t=gl.GL_CLAMP_TO_EDGE,
			r=gl.GL_CLAMP_TO_EDGE,
		},
		magFilter = gl.GL_LINEAR,
		minFilter = gl.GL_LINEAR,
	}

	local shaderDefs = [[
#define DLAMBDA .1
#define M_PI 3.14159265358979311599796346854418516159057617187500
#define ITERATIONS 1

#define SCHWARZSCHILD_CARTESIAN
//#define SCHWARZSCHILD_SPHERIC
//#define ALCUBIERRE

#ifdef SCHWARZSCHILD_CARTESIAN
#define EVENT_HORIZON_RADIUS	1.
#define COORDINATE_CENTER	vec3(20., 0., 0.);	//working but not without some numerical instabilities
#define TO_COORDINATES(v)	(v)
#define FROM_COORDINATES(v)	(v)

#endif

#ifdef SCHWARZSCHILD_SPHERIC
#define EVENT_HORIZON_RADIUS	1.
#define COORDINATE_CENTER	vec3(20., 0., 0.);	//working but not without some numerical instabilities
#define TO_COORDINATES(v)	euclidianToSpheric(v)
#define FROM_COORDINATES(v)	sphericToEuclidian(v)
#endif

#ifdef ALCUBIERRE
#define COORDINATE_CENTER		vec3(0., 3., 0.);
#define TO_COORDINATES(v)		(v)
#define FROM_COORDINATES(v)		(v)
#define WARP_BUBBLE_RADIUS		2.
#define WARP_BUBBLE_DISTANCE	2.
#define WARP_BUBBLE_VELOCITY	.9
#define	WARP_THICKNESS			.1
#endif

/*
x = r
y = polar angle
z = inclination angle
*/

vec4 euclidianToSpheric(vec4 xyzt) {
	vec4 rhpt;
	float r2 = length(xyzt.xy);	//length in xy
	rhpt.x = length(xyzt.xyz);	//radial coordinate / length in xyz
	rhpt.y = atan(xyzt.y, xyzt.x);			//polar coordinate / angle about z axis
	rhpt.z = atan(r2, xyzt.z);				//inclination coordinate / angle from positive z axis to negative z axis
	rhpt.w = xyzt.w;

	rhpt.y = mod(rhpt.y, 2. * M_PI);
	rhpt.z = mod(rhpt.z, M_PI);

	return rhpt;
}

vec4 sphericToEuclidian(vec4 rhpt) {
	vec4 xyzt;
	xyzt.x = rhpt.x * cos(rhpt.y) * sin(rhpt.z);
	xyzt.y = rhpt.x * sin(rhpt.y) * sin(rhpt.z);
	xyzt.z = rhpt.x * cos(rhpt.z);
	xyzt.w = rhpt.w;
	return xyzt;
}

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

]]

--[[
inverse camera parameters:
[R^-1 0] [1 -T]   [R^-1, -R^-1 T]
[0    1] [0  1] = [   0,    1   ]
--]]
	local initLightShaderVertexCode = [[
in vec4 vertex;
out vec3 eyePos;
out vec3 vtxPos;
uniform mat4 projMat, mvMat;
void main() {
	mat3 rot = transpose(mat3(
		mvMat[0].xyz,
		mvMat[1].xyz,
		mvMat[2].xyz));
	eyePos = -rot * mvMat[3].xyz;
	vtxPos = rot * vertex.xyz + eyePos;

	gl_Position = projMat * vertex;
}
]]

	local initLightShaderFragmentCode = shaderDefs .. [[
in vec3 eyePos;
in vec3 vtxPos;
out vec4 fragColor;

void main() {
	vec4 rel = vec4(eyePos, 1.);
	vec4 relDirEnd = vec4(vtxPos - eyePos, 1.);

	//convert to spherical coordinates for iteration
	rel = TO_COORDINATES(rel);

	//find the difference of a vector pointing in the initial direction, in the black hole's spherical coordinates
	relDirEnd = TO_COORDINATES(relDirEnd);

	//get the direction difference in spheric coordinates
	vec4 relDiff = relDirEnd - rel;

	//...and start us off in the view plane
	rel = relDirEnd;

	$assign
}
]]

	initLightPosShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = shaderDefs .. initLightShaderVertexCode,
		fragmentCode = initLightShaderFragmentCode:gsub('$assign', 'fragColor = rel;'),
	}:useNone()

	initLightVelShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = shaderDefs .. initLightShaderVertexCode,
		fragmentCode = initLightShaderFragmentCode:gsub('$assign', 'fragColor = relDiff;'),
	}:useNone()

	local lightRes = 1024
	for _,texs in ipairs{lightPosTexs, lightVelTexs} do
		for i=1,2 do
			texs[i] = GLTex2D{
				width=lightRes,
				height=lightRes,
				format=gl.GL_RGBA,
				type=gl.GL_FLOAT,
				internalFormat=gl.GL_RGBA32F,
				minFilter=gl.GL_NEAREST,
				magFilter=gl.GL_NEAREST,
			}
		end
	end

	fbo = FBO{width=lightRes, height=lightRes}:unbind()

	local iterateLightShaderVertexCode = [[
in vec4 vertex;
out vec2 tc;
uniform mat4 projMat, mvMat;
void main() {
	tc = vertex.xy;
	gl_Position = projMat * mvMat * vertex;
}
]]

	local iterateLightShaderFragmentCode = [[
in vec2 tc;
out vec4 fragColor;

uniform sampler2D posTex;
uniform sampler2D velTex;
void main() {

	vec4 rel = texture(posTex, tc);
	vec4 relDiff = texture(velTex, tc);

	//integrate along geodesic
	for (int i = 0; i < ITERATIONS; i++) {
		//-x''^a = conn^a_bc x'^b x'^c
		vec4 negRelDiff2;

#ifdef SCHWARZSCHILD_CARTESIAN
		float r = length(rel.xyz);
		float posDotVel = dot(rel.xyz, relDiff.xyz);
		float posDotVelSq = posDotVel * posDotVel;
		float relDiffSq = dot(relDiff.xyz, relDiff.xyz);
		float R = EVENT_HORIZON_RADIUS;
		float rSq = r * r;
		float relDiffTSq = relDiff.w * relDiff.w;
		float scale = (R / (r * rSq)) * (.5 * relDiffTSq * (1. - R / r) + relDiffSq - .5 * posDotVelSq * (3. * r - 2. * R) / (rSq * (r - R)));
		negRelDiff2.xyz = rel.xyz * scale;
		negRelDiff2.w = R / (rSq * (r - R)) * posDotVel * relDiff.w;
#endif

#ifdef SCHWARZSCHILD_SPHERIC
		rel.y = mod(rel.y, 2. * M_PI);
		rel.z = mod(rel.z, 2. * M_PI);

		float radius = rel.x;
		float radiusFromEventHorizon = radius - EVENT_HORIZON_RADIUS;
		float cosPhi = cos(rel.z);
		float sinPhi = sin(rel.z);
		float cotPhi = cosPhi / sinPhi;

		negRelDiff2.w = EVENT_HORIZON_RADIUS / (radius * radiusFromEventHorizon) * relDiff.w * relDiff.x;
		negRelDiff2.x = EVENT_HORIZON_RADIUS * radiusFromEventHorizon / (2. * radius * radius * radius) * relDiff.w * relDiff.w
			- EVENT_HORIZON_RADIUS / (2. * radius * radiusFromEventHorizon) * relDiff.x * relDiff.x
			- radiusFromEventHorizon * sinPhi * sinPhi * relDiff.y * relDiff.y
			- radiusFromEventHorizon * relDiff.z * relDiff.z
		;
		negRelDiff2.y = 2. / radius * relDiff.x * relDiff.y + 2. * cotPhi * relDiff.y * relDiff.z;
		negRelDiff2.z = 2. / radius * relDiff.x * relDiff.z - sinPhi * cosPhi * relDiff.y * relDiff.y;
#endif

#ifdef ALCUBIERRE
		float rs = length(rel.xyz - vec3(WARP_BUBBLE_DISTANCE, 0., 0.));
		float sigmaFront = WARP_THICKNESS * (rs + WARP_BUBBLE_RADIUS);
		float sigmaCenter = WARP_THICKNESS * rs;
		float sigmaBack = WARP_THICKNESS * (rs - WARP_BUBBLE_RADIUS);
		float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanh(sigmaCenter));
		float sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
		float dfScalar = sechDiff / (2. * rs * tanh(sigmaCenter));
		float ft = -WARP_BUBBLE_VELOCITY * WARP_THICKNESS * (rel.x - WARP_BUBBLE_DISTANCE) * dfScalar;
		float fx = WARP_THICKNESS * (rel.x - WARP_BUBBLE_DISTANCE) * dfScalar;
		float fy = WARP_THICKNESS * (rel.y) * dfScalar;
		float fz = WARP_THICKNESS * (rel.z) * dfScalar;

		negRelDiff2.w = f * f * fx * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * relDiff.w * relDiff.w
			- 2. * f * fx * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * relDiff.w * relDiff.x
			- 2. * f * fy * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY / 2. * relDiff.w * relDiff.y
			- 2. * f * fz * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY / 2. * relDiff.w * relDiff.z
			+ fx * WARP_BUBBLE_VELOCITY * relDiff.x * relDiff.x
			+ 2. * fy * WARP_BUBBLE_VELOCITY / 2. * relDiff.x * relDiff.y
			+ 2. * fz * WARP_BUBBLE_VELOCITY / 2. * relDiff.x * relDiff.z
		;
		negRelDiff2.x = (f * f * f * fx * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY - f * fx * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY - ft * WARP_BUBBLE_VELOCITY) * relDiff.w * relDiff.w
			- 2. * f * f * fx * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * relDiff.w * relDiff.x
			- 2. * (f * f * fy * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY + fy * WARP_BUBBLE_VELOCITY) / 2. * relDiff.w * relDiff.y
			- 2. * (f * f * fz * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY + fz * WARP_BUBBLE_VELOCITY) / 2. * relDiff.w * relDiff.z
			+ f * fx * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * relDiff.x * relDiff.x
			+ 2. * f * fy * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY / 2. * relDiff.x * relDiff.y
			+ 2. * f * fz * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY / 2. * relDiff.x * relDiff.z
		;
		negRelDiff2.y = -f * fy * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * relDiff.w * relDiff.w
			+ 2. * fy * WARP_BUBBLE_VELOCITY / 2. * relDiff.w * relDiff.x
		;
		negRelDiff2.z = -f * fz * WARP_BUBBLE_VELOCITY * WARP_BUBBLE_VELOCITY * relDiff.w * relDiff.w
			+ 2. * fz * WARP_BUBBLE_VELOCITY / 2. * relDiff.w * relDiff.x
		;
#endif
		rel += DLAMBDA * relDiff;
		relDiff -= DLAMBDA * negRelDiff2;
	}

	$assign
}
]]
	local iterateLightShaderUniforms = {
		posTex = 0,
		velTex = 1,
	}

	iterateLightPosShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = iterateLightShaderVertexCode,
		fragmentCode = shaderDefs .. iterateLightShaderFragmentCode:gsub('$assign', 'fragColor = rel;'),
		uniforms = iterateLightShaderUniforms,
	}:useNone()

	iterateLightVelShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = iterateLightShaderVertexCode,
		fragmentCode = shaderDefs .. iterateLightShaderFragmentCode:gsub('$assign', 'fragColor = relDiff;'),
		uniforms = iterateLightShaderUniforms,
	}:useNone()


--[[
how it'll work ...
1) start with the normalized world vector
2) iterate geodesic (backwards) a few times
3) use the resulting coordinates

you could do the iteration in the shader loop, or you could use a fbo to store state information ...
--]]
	drawLightShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = [[
in vec4 vertex;
out vec2 tc;
uniform mat4 projMat;
void main() {
	tc = vertex.xy;
	gl_Position = projMat * vertex;
}
]],
		fragmentCode = shaderDefs .. [[
in vec2 tc;
out vec4 fragColor;
uniform sampler2D posTex;
uniform samplerCube cubeTex;

void main() {
	vec4 rel = texture(posTex, tc);

#if 0
	//warn if the light ended up in the event horizon
	if (rel.x <= EVENT_HORIZON_RADIUS) {
		fragColor = vec4(tc.x, tc.y, 0., 1.);
	} else
#endif
	{
		//convert back to euclidian space
		vec3 result = FROM_COORDINATES(rel).xyz;

		//convert from black-hole-centered to view-centered
		//result += mvMat[3].xyz;
		result += COORDINATE_CENTER;
		//result += vec3(
		//	dot(mvMat[0].xyz, mvMat[3].xyz),
		//	dot(mvMat[1].xyz, mvMat[3].xyz),
		//	dot(mvMat[2].xyz, mvMat[3].xyz));

		//with distortion
		fragColor = texture(cubeTex, result);

		//no distortion
	//	fragColor = texture(cubeTex, pos);

		//difference
	//	vec4 color = texture(cubeTex, result) - texture(cubeTex, pos);
	//	fragColor = vec4(abs(color.x), abs(color.y), abs(color.z), abs(color.w));
	}
}
]],
		uniforms = {
			posTex = 0,
			cubeTex = 1,
		},
	}:useNone()

	self.wireframeShader = GLProgram{
		version = 'latest',
		precision = 'best',
		vertexCode = [[
in vec4 vertex;
uniform mat4 projMat, mvMat;
void main() {
	gl_Position = projMat * mvMat * vertex;
}
]],
		fragmentCode = [[
out vec4 fragColor;
void main() {
	fragColor = vec4(1., 1., 1., 1.);
}
]],
	}:useNone()

	gl.glClearColor(.3, .3, .3, 1)
end

function App:event(e)
	App.super.event(self, e)
	if e[0].type == sdl.SDL_MOUSEMOTION then
		if leftButtonDown then
			lightInitialized = false
		end
	elseif e[0].type == sdl.SDL_MOUSEBUTTONDOWN then
		if e[0].button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = true
		elseif e[0].button.button == sdl.SDL_BUTTON_WHEELUP then
			lightInitialized = false
		elseif e[0].button.button == sdl.SDL_BUTTON_WHEELDOWN then
			lightInitialized = false
		end
	elseif e[0].type == sdl.SDL_MOUSEBUTTONUP then
		if e[0].button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = false
		end
	elseif e[0].type == sdl.SDL_KEYDOWN then
		if e[0].key.keysym.sym == sdl.SDLK_i then
			doIteration = 2
		elseif e[0].key.keysym.sym == sdl.SDLK_j then
			doIteration = 1
		end
	end
end

function App:updateGUI()
	ig.igText'testing testing'
	ig.luatableInputFloat('delta lambda', _G, 'deltaLambdaPtr')
	ig.luatableInputInt('iterations', _G, 'iterationsPtr')
end

function App:update()
	local view = self.view
	local viewWidth, viewHeight = self.width, self.height
	local aspectRatio = viewWidth / viewHeight

	-- init light if necessary
	if not lightInitialized then
		lightInitialized = true

		gl.glViewport(0, 0, fbo.width, fbo.height)

		local tanFovY = math.tan(view.fovY / 2)
		local tanFovX = aspectRatio * tanFovY
local oldpos = view.pos
view.pos = vec3d(0,0,0)
view:setup(aspectRatio)
view.pos = oldpos

		for _,args in ipairs{
			{tex=lightPosTexs[1], shader=initLightPosShader},
			{tex=lightVelTexs[1], shader=initLightVelShader},
		} do
			fbo:draw{
				dest = args.tex,
				shader = args.shader,
				uniforms = {
					mvMat = view.mvMat.ptr,
					projMat = view.projMat.ptr,
				},
				callback = function()
					gl.glBegin(gl.GL_QUADS)
					for _,uv in ipairs(uvs) do
						gl.glVertex3f((uv.x - .5) * 2 * tanFovX,
							(uv.y - .5) * 2 * tanFovY,
							-1)
					end
					gl.glEnd()
				end,
			}
		end
	end

	if doIteration ~= 0 then
		if doIteration == 1 then doIteration = 0 end
		--print('iterating...')
		gl.glViewport(0, 0, fbo.width, fbo.height)
		gl.glClear(gl.GL_COLOR_BUFFER_BIT)

		for _,args in ipairs{
			{tex=lightPosTexs[2].id, shader=iterateLightPosShader},
			{tex=lightVelTexs[2].id, shader=iterateLightVelShader},
		} do
			fbo:draw{
				dest = args.tex,
				shader = args.shader,
				uniforms = {
					projMat = self.fboProjMat.ptr,
					mvMat = self.identMat.ptr,
				},
				texs = {lightPosTexs[1], lightVelTexs[1]},
				callback = function()
					gl.glBegin(gl.GL_QUADS)
					for _,uv in ipairs(uvs) do
						gl.glVertex2f(uv:unpack())
					end
					gl.glEnd()
				end,
			}
		end

		-- swap
		lightPosTexs[1], lightPosTexs[2] = lightPosTexs[2], lightPosTexs[1]
		lightVelTexs[1], lightVelTexs[2] = lightVelTexs[2], lightVelTexs[1]
	end

	gl.glViewport(0,0,viewWidth, viewHeight)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)

	view:setupModelView()
	--gl.glLoadIdentity()
	--gl.glTranslatef(20, 0, 0)

	drawLightShader:use()
	drawLightShader:setUniforms{
		projMat = self.fboProjMat.ptr,
		mvMat = view.mvMat.ptr,
	}
	lightPosTexs[1]:bind(0)
	skyTex:bind(1)
	gl.glBegin(gl.GL_QUADS)
	for _,uv in ipairs(uvs) do
		gl.glVertex2f(uv:unpack())
	end
	gl.glEnd()
	drawLightShader:useNone()

	gl.glLoadIdentity()
	gl.glActiveTexture(gl.GL_TEXTURE0)

	-- [[
	local oldpos = view.pos
	view:setup(aspectRatio)

	self.wireframeShader:use()
	self.wireframeShader:setUniforms{
		projMat = view.projMat.ptr,
		mvMat = view.mvMat.ptr,
	}

	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_LINE)
	gl.glBegin(gl.GL_QUADS)
	local s = 1
	for _,face in ipairs(cubeFaces) do
		for _,i in ipairs(face) do
			local x = bit.band(i, 1)
			local y = bit.band(bit.rshift(i, 1), 1)
			local z = bit.band(bit.rshift(i, 2), 1)
			gl.glVertex3d(s*(x*2-1),s*(y*2-1),s*(z*2-1))
		end
	end
	gl.glEnd()
	gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
	self.wireframeShader:useNone()
	--]]

	glreport('update done')
	App.super.update(self)
end

return App():run()
