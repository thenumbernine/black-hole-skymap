#!/usr/bin/env luajit

--[[
cube map render
shader that iterates along geodesic from the view plane backwards through time to infinity
then reprojects onto the cubemap

schwarzschild metric: ds^2 = -(1-R/r)dt^2 + 1/(1-R/r)dr^2 + r^2 (dtheta^2 + sin(theta)^2 dphi^2)
where R = 2M is the schwarzschild radius

christoffel symbols:

conn^r_t_t = R(r-R)/(2r^3)
conn^t_t_r = -conn^r_r_r = R/(2r(r-R))
conn^phi_r_phi = conn^theta_r_theta = 1/r
conn^r_phi_phi = -(r - R)
conn^theta_phi_theta = cos(phi)/sin(phi)
conn^r_theta_theta = -(r - R) sin(phi)^2
conn^phi_theta_theta = -sin(phi) cos(phi)

-t'' = R/(r(r-R)) t' r'
-r'' = R(r-R)/(2r^3) t'^2 - R/(2r(r-R)) r'^2 - (r - R) phi'^2 - (r - R) sin(phi)^2 theta'^2
-theta'' = 2/r r' theta' + 2 cos(phi)/sin(phi) phi' theta'
-phi'' = 2/r r' phi' - sin(phi) cos(phi) theta'^2
--]]

require 'ext'
local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local glu = require 'ffi.glu'
local sdl = require 'ffi.sdl'
local ig = require 'imgui'
local vec2d = require 'vec-ffi.vec2d'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'
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

App.title = 'black hole raytracer'
App.viewDist = 20

function App:initGL()
	App.super.initGL(self)

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
#version 120

#define DLAMBDA 0.1
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
	local initLightShaderVertexCode = shaderDefs .. [[
varying vec3 eyePos;
varying vec3 vtxPos;
void main() {
	mat3 rot = transpose(mat3(
		gl_ModelViewMatrix[0].xyz,
		gl_ModelViewMatrix[1].xyz,
		gl_ModelViewMatrix[2].xyz));
	eyePos = -rot * gl_ModelViewMatrix[3].xyz;
	vtxPos = rot * gl_Vertex.xyz + eyePos;

	gl_Position = gl_ProjectionMatrix * gl_Vertex;
}
]]

	local initLightShaderFragmentCode = shaderDefs .. [[
varying vec3 eyePos;
varying vec3 vtxPos;

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
		vertexCode = initLightShaderVertexCode,
		fragmentCode = initLightShaderFragmentCode:gsub('$assign', 'gl_FragColor = rel;'),
	}
	
	initLightVelShader = GLProgram{
		vertexCode = initLightShaderVertexCode,
		fragmentCode = initLightShaderFragmentCode:gsub('$assign', 'gl_FragColor = relDiff;'),
	}
	
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

	fbo = FBO{width=lightRes, height=lightRes}

	local iterateLightShaderVertexCode = [[
varying vec2 tc;
void main() {
	tc = gl_Vertex.xy;
	gl_Position = ftransform();
}
]]

	local iterateLightShaderFragmentCode = shaderDefs .. [[

uniform sampler2D posTex;
uniform sampler2D velTex;
varying vec2 tc;
void main() {

	vec4 rel = texture2D(posTex, tc);
	vec4 relDiff = texture2D(velTex, tc);
	
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
		vertexCode = iterateLightShaderVertexCode,
		fragmentCode = iterateLightShaderFragmentCode:gsub('$assign', 'gl_FragColor = rel;'),
		uniforms = iterateLightShaderUniforms,
	}
	
	iterateLightVelShader = GLProgram{
		vertexCode = iterateLightShaderVertexCode,
		fragmentCode = iterateLightShaderFragmentCode:gsub('$assign', 'gl_FragColor = relDiff;'),
		uniforms = iterateLightShaderUniforms,
	}
	
		
--[[
how it'll work ...
1) start with the normalized world vector
2) iterate geodesic (backwards) a few times
3) use the resulting coordinates

you could do the iteration in the shader loop, or you could use a fbo to store state information ...
--]]
	drawLightShader = GLProgram{
		vertexCode=[[
varying vec2 tc;
void main() {
	tc = gl_Vertex.xy;
	gl_Position = gl_ProjectionMatrix * gl_Vertex;
}
]],
		fragmentCode = shaderDefs .. [[
uniform sampler2D posTex;
uniform samplerCube cubeTex;
varying vec2 tc;

void main() {
	vec4 rel = texture2D(posTex, tc);
	
#if 0
	//warn if the light ended up in the event horizon
	if (rel.x <= EVENT_HORIZON_RADIUS) {
		gl_FragColor = vec4(tc.x, tc.y, 0., 1.);
	} else 
#endif
	{
		//convert back to euclidian space
		vec3 result = FROM_COORDINATES(rel).xyz;
		
		//convert from black-hole-centered to view-centered
		//result += gl_ModelViewMatrix[3].xyz;
		result += COORDINATE_CENTER; 
		//result += vec3(
		//	dot(gl_ModelViewMatrix[0].xyz, gl_ModelViewMatrix[3].xyz),
		//	dot(gl_ModelViewMatrix[1].xyz, gl_ModelViewMatrix[3].xyz),
		//	dot(gl_ModelViewMatrix[2].xyz, gl_ModelViewMatrix[3].xyz));
		
		//with distortion
		gl_FragColor = textureCube(cubeTex, result);
		
		//no distortion
	//	gl_FragColor = textureCube(cubeTex, pos);

		//difference
	//	vec4 color = textureCube(cubeTex, result) - textureCube(cubeTex, pos);
	//	gl_FragColor = vec4(abs(color.x), abs(color.y), abs(color.z), abs(color.w));
	}
}
]],
		uniforms = {
			posTex = 0,
			cubeTex = 1,
		},
	}
	
	gl.glClearColor(.3, .3, .3, 1)		
end

function App:event(event, eventPtr)
	App.super.event(self, event, eventPtr)
	if event.type == sdl.SDL_MOUSEMOTION then
		if leftButtonDown then
			lightInitialized = false
		end
	elseif event.type == sdl.SDL_MOUSEBUTTONDOWN then
		if event.button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = true
		elseif event.button.button == sdl.SDL_BUTTON_WHEELUP then
			lightInitialized = false
		elseif event.button.button == sdl.SDL_BUTTON_WHEELDOWN then
			lightInitialized = false
		end
	elseif event.type == sdl.SDL_MOUSEBUTTONUP then
		if event.button.button == sdl.SDL_BUTTON_LEFT then
			leftButtonDown = false
		end
	elseif event.type == sdl.SDL_KEYDOWN then
		if event.key.keysym.sym == sdl.SDLK_i then
			doIteration = 2
		elseif event.key.keysym.sym == sdl.SDLK_j then
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
	App.super.update(self)
	
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

		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glOrtho(0,1,0,1,-1,1)
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()
		
		for _,args in ipairs{
			{tex=lightPosTexs[2].id, shader=iterateLightPosShader},
			{tex=lightVelTexs[2].id, shader=iterateLightVelShader},
		} do
			fbo:draw{
				dest = args.tex,
				shader = args.shader,
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

	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(0,1,0,1,-1,1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	view:setupModelView()
	--gl.glLoadIdentity()
	--gl.glTranslatef(20, 0, 0)

	drawLightShader:use()
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
	--]]

	glreport('update done')
end

App():run()
