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

local openglapp = require 'openglapp'
local bit = require 'bit'
local ffi = require 'ffi'
local gl = require 'gl'
local glu = require 'ffi.req' 'glu'
local sdl = require 'sdl'

--local tw = require 'ffi.req' 'anttweakbar'

local Quat = require 'vec.quat'
require 'glutil.tex'
require 'glutil.shader'
require 'glutil.pingpong'

local skyTex
local viewRot = Quat(0,math.sqrt(.5),0,math.sqrt(.5))*Quat(math.sqrt(.5),0,0,-math.sqrt(.5))
local zNear = .1
local zFar = 100
local tanFov = 1
local zoomFactor = .9
local leftButtonDown

local bendLightShader
local initLightShader
local iterateLightShader

local cosTex	-- map from 0,2pi to value
local sinTex	-- map form 0,2pi to value
local atanTex	-- map from -1,1 to value
local sqrtTex	-- map from 0,100 to value

-- notice that the order matches global 'sides'
cubeFaces = {
	{5,7,3,1};		-- <- each value has the x,y,z in the 0,1,2 bits (off = 0, on = 1)
	{6,4,0,2};
	{2,3,7,6};
	{4,5,1,0};
	{6,7,5,4};
	{0,1,3,2};
}

local deltaLambdaPtr = ffi.new('double[1]', 1.)
local iterationsPtr = ffi.new('int[1]', 1)

local sdlVersion = ffi.new('SDL_version[1]')
local sdlEventCopy = ffi.new('SDL_Event[1]')	-- for anttweakbar

openglapp:run{
	init = function()
		sdl.SDL_GetVersion(sdlVersion)

		if tw then
			tw.TwInit(tw.TW_OPENGL, nil)

			local viewWidth, viewHeight = openglapp:size()
			tw.TwWindowSize(viewWidth, viewHeight)
		
			tw.TwDefine("GLOBAL help='testing testing'");
			local bar = tw.TwNewBar('Controls')
			tw.TwAddVarRW(bar, 'deltalambda', tw.TW_TYPE_DOUBLE, deltaLambdaPtr, " label='delta lambda' min=0 max=10 step=0.01 keyIncr=l keyDecr=L help='delta lambda'")
			tw.TwAddVarRW(bar, 'iterations', tw.TW_TYPE_INT32, iterationsPtr, " label='iterations' min=0 max=1000 step=1 keyIncr=i keyDecr=I help='iterations'")
		end
		
		skyTex = TexCube{
			filenames = {
				'skytex/sky-infrared-cube-xp.png',
				'skytex/sky-infrared-cube-xn.png',
				'skytex/sky-infrared-cube-yp.png',
				'skytex/sky-infrared-cube-yn.png',
				'skytex/sky-infrared-cube-zp.png',
				'skytex/sky-infrared-cube-zn.png',
			},
			wrap={
				s=gl.GL_CLAMP_TO_EDGE,
				t=gl.GL_CLAMP_TO_EDGE,
				r=gl.GL_CLAMP_TO_EDGE,
			},
			magFilter = gl.GL_LINEAR,
			minFilter = gl.GL_LINEAR,
		}
		
		initLightShader = ShaderProgram{
			vertexCode=[[
varying vec3 dir;			
void main() {
	vec4 mvpos = gl_ModelViewMatrix * gl_Vertex;
	vec4 center = gl_ModelViewMatrix * vec4(0., 0., 0., 1.);
	dir = mvpos.xyz - center.xyz;
	gl_Position = gl_ProjectionMatrix * mvpos;
}
]],
			fragmentCode=[[
varying vec3 dir;
void main() {
	gl_FragColor.xyz = normalize(dir);
	gl_FragColor.w = 1.;
}
]],
		}
		iterateLightShader = ShaderProgram{
		}
		
--[[
how it'll work ...
1) start with the normalized world vector
2) iterate geodesic (backwards) a few times
3) use the resulting coordinates

you could do the iteration in the shader loop, or you could use a fbo to store state information ...
--]]
		bendLightShader = ShaderProgram{
			vertexCode=[[
varying vec3 pos;
void main() {
	pos = normalize(gl_Vertex.xyz);
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
]],
			fragmentCode=[[
uniform samplerCube cubeTex;
varying vec3 pos;

#define M_PI 3.14159265358979311599796346854418516159057617187500
#define ITERATIONS 30
#define DLAMBDA 1.
#define EVENT_HORIZON_RADIUS 2.

vec4 euclidianToSpheric(vec4 xyzt) {
	vec4 rhpt;
	float r2 = sqrt(dot(xyzt.xy,xyzt.xy));
	rhpt.x = sqrt(r2*r2 + xyzt.z*xyzt.z);
	rhpt.y = atan(xyzt.y, xyzt.x);
	rhpt.z = atan(r2, xyzt.z);
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

void main() {
	const vec3 blackHoleCenter = vec3(20., 0., 0.);

	vec4 relDirEnd = vec4(normalize(pos), 1.);
	vec4 rel = vec4(0., 0., 0., 0.);

	//center coordinates at the stellar body center
	rel.xyz -= blackHoleCenter;
	relDirEnd.xyz -= blackHoleCenter;
	
	//convert to spherical coordinates for iteration
	rel = euclidianToSpheric(rel);
	
	//find the difference of a vector pointing in the initial direction, in the black hole's spherical coordinates
	relDirEnd = euclidianToSpheric(relDirEnd);

	//get the direction difference in spheric coordinates
	vec4 relDiff = relDirEnd - rel;
	
	//...and start us off in the view plane
	rel = relDirEnd;
	
	//TODO in the y- direction, iteration is going screwy
	
	//integrate along geodesic
	for (int i = 0; i < ITERATIONS; i++) {
		//-x''^a = conn^a_bc x'^b x'^c

		float radius = rel.x;
		float radiusFromEventHorizon = radius - EVENT_HORIZON_RADIUS;
		float cosPhi = cos(rel.z);
		float sinPhi = sin(rel.z);
		
		vec4 negRelDiff2;
		negRelDiff2.w = EVENT_HORIZON_RADIUS / (radius * radiusFromEventHorizon) * relDiff.w * relDiff.x;
		negRelDiff2.x = EVENT_HORIZON_RADIUS * radiusFromEventHorizon / (2. * radius * radius * radius) * relDiff.w * relDiff.w
			- EVENT_HORIZON_RADIUS / (2. * radius * radiusFromEventHorizon) * relDiff.x * relDiff.x
			- radiusFromEventHorizon * sinPhi * sinPhi * relDiff.y * relDiff.y
			- radiusFromEventHorizon * relDiff.z * relDiff.z;
		negRelDiff2.y = 2. / radius * relDiff.x * relDiff.y + 2. * cosPhi / sinPhi * relDiff.y * relDiff.z;
		negRelDiff2.z = 2. / radius * relDiff.x * relDiff.z - sinPhi * cosPhi * relDiff.y * relDiff.y;
	
		rel += DLAMBDA * relDiff;
		relDiff -= DLAMBDA * negRelDiff2;
		
		rel.y = mod(rel.y, 2. * M_PI);
		rel.z = mod(rel.z, M_PI);
		
//		rel.w = 0.;
//		relDiff.w = 1.;
	}
	
	//warn if the light ended up in the event horizon
	if (rel.x <= EVENT_HORIZON_RADIUS) {
		gl_FragColor = vec4(1., 0., 0., 1.);
	} else 
	{
		
		//convert back to euclidian space
		vec3 result = sphericToEuclidian(rel).xyz;
		
		//convert from black-hole-centered to view-centered
		result += blackHoleCenter;
		
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
			uniforms={
				cubeTex = 0,
			},
		}
		
		gl.glEnable(gl.GL_DEPTH_TEST)
		gl.glClearColor(.3, .3, .3, 1)		
	end,
	event = function(event)
		sdlEventCopy[0] = event	-- memcpy I hope
		if tw and 0 ~= tw.TwEventSDL(sdlEventCopy, sdlVersion[0].major, sdlVersion[0].minor) then return end
		if event.type == sdl.SDL_MOUSEMOTION then
			if leftButtonDown then
				local idx = event.motion.xrel
				local idy = event.motion.yrel
				local magn = math.sqrt(idx * idx + idy * idy)
				local dx = idx / magn
				local dy = idy / magn
				local r = Quat():fromAngleAxis(dy, dx, 0, -magn)
				viewRot = (r * viewRot):normalize()
			end
		elseif event.type == sdl.SDL_MOUSEBUTTONDOWN then
			if event.button.button == sdl.SDL_BUTTON_LEFT then
				leftButtonDown = true
			elseif event.button.button == sdl.SDL_BUTTON_WHEELUP then
				tanFov = tanFov * zoomFactor
			elseif event.button.button == sdl.SDL_BUTTON_WHEELDOWN then
				tanFov = tanFov / zoomFactor
			end
		elseif event.type == sdl.SDL_MOUSEBUTTONUP then
			if event.button.button == sdl.SDL_BUTTON_LEFT then
				leftButtonDown = false
			end
		end
	end,
	update = function()
		local viewWidth, viewHeight = openglapp:size()
		if tw then tw.TwWindowSize(viewWidth, viewHeight) end
		local aspectRatio = viewWidth / viewHeight

		gl.glViewport(0,0,viewWidth, viewHeight)
		gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
				
		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glFrustum(-zNear * aspectRatio * tanFov, zNear * aspectRatio * tanFov, -zNear * tanFov, zNear * tanFov, zNear, zFar);
		
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()
		
		local aa = viewRot:toAngleAxis()
		gl.glRotatef(aa[4], aa[1], aa[2], aa[3])
		
		bendLightShader:use()
		skyTex:bind(0)

		gl.glActiveTexture(gl.GL_TEXTURE0)
		gl.glBegin(gl.GL_QUADS)
		for _,face in ipairs(cubeFaces) do
			for _,vtx in ipairs(face) do
				local x = bit.band(vtx,1)*2-1
				local y = bit.band(bit.rshift(vtx,1),1)*2-1
				local z = bit.band(bit.rshift(vtx,2),1)*2-1
				gl.glTexCoord3f(x,y,z)
				gl.glVertex3f(x,y,z)
			end
		end
		gl.glEnd()
		bendLightShader:useNone()
		
		if tw then tw.TwDraw() end
		
		glreport('update done')
	end,
	exit = function()
		if tw then tw.TwTerminate() end
	end,
}
