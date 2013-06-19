local openglapp = require 'openglapp'
local bit = require 'bit'
local gl = require 'ffi.gl'
local glu = require 'ffi.glu'
local sdl = require 'ffi.sdl'
local Quat = require 'vec.quat'
local vec2 = require 'vec.vec2'
require 'glutil.tex'
require 'glutil.glutil'

local quadOffsets = {
	vec2(0,0),
	vec2(1,0),
	vec2(1,1),
	vec2(0,1),
}

local leftButtonDown
local projectionDisplayList = {}
local cubeMapDisplayList = {}
local viewRot = Quat()
local zoomScale = 1
local zoomFactor = .9
local zNear = .1
local zFar = 100
local tanFovX = 1
local tanFovY = 1

local tex
local cubeTex

local saved
		
local thetaStep = 1
local phiStep = 1

function drawSkySphere()
	tex:enable()
	tex:bind()
	gl.glBegin(gl.GL_QUADS)
	for baseTheta=-180,180-thetaStep,thetaStep do
		for basePhi=-90,90-phiStep,phiStep do
			local r = 1
			for _,quadOffset in ipairs(quadOffsets) do
				local theta = math.rad(baseTheta + thetaStep * quadOffset[1])
				local phi = math.rad(basePhi + phiStep * quadOffset[2])
				
				local spherePhi = phi + math.pi/2
				local x = r * math.cos(theta) * math.sin(spherePhi)
				local y = r * math.sin(theta) * math.sin(spherePhi)
				local z = r * math.cos(spherePhi)
				--[[
				local u = theta / 360
				local v = phi / 180
				--]]
				--[[
				we have theta and phi
				we need fwd-project to find x and y
				then remap that to the texture
				--]]
				local mollPhi = phi
				local mollTheta = 2 * math.asin(math.clamp(-1, 2 * mollPhi / math.pi, 1))
				for i=1,100 do
					local denom = 1 + math.cos(mollTheta)
					if math.abs(denom) < 1e-20 then break end
					local deltaMollTheta = -(mollTheta + math.sin(mollTheta) - math.pi * math.sin(mollPhi)) / denom
					mollTheta = mollTheta + deltaMollTheta
				end
				mollTheta = mollTheta * .5
				local u = .5 + .495 * theta * math.cos(mollTheta) / math.pi
				local v = .5 + .495 * math.sin(mollTheta)
				
				if not uMin or u < uMin then uMin = u end
				if not uMax or u > uMax then uMax = u end
				if not vMin or v < vMin then vMin = v end
				if not vMax or v > vMax then vMax = v end
				
				gl.glTexCoord2f(u,v)
				--gl.glColor3f(u,v,1)
				gl.glVertex3f(x,y,z)
				--gl.glVertex2f(u,v)
			end
		end
	end
	gl.glEnd()
	tex:disable()
end

-- only do this once
function drawSkyCube()
	cubeTex:enable()
	cubeTex:bind()
	gl.glBegin(gl.GL_QUADS)
	for baseTheta=-180,180-thetaStep,thetaStep do
		for basePhi=-90,90-phiStep,phiStep do
			local r = 1
			for _,quadOffset in ipairs(quadOffsets) do
				local theta = math.rad(baseTheta + thetaStep * quadOffset[1])
				local phi = math.rad(basePhi + phiStep * quadOffset[2])
				
				local spherePhi = phi + math.pi/2
				local x = r * math.cos(theta) * math.sin(spherePhi)
				local y = r * math.sin(theta) * math.sin(spherePhi)
				local z = r * math.cos(spherePhi)

				gl.glTexCoord3f(x,y,z)
				gl.glVertex3f(x,y,z)
			end
		end
	end
	gl.glEnd()
	cubeTex:disable()
end

function saveCubeSides()	
	local viewWidth, viewHeight = openglapp:size()
	local cubeSide = math.min(viewWidth, viewHeight)
	cubeSide = math.floor(2^math.floor(math.log(cubeSide)/math.log(2)) + .1)
	print('cubeSide',cubeSide)
	gl.glViewport(0,0,cubeSide,cubeSide)
	
	local im = img.new{width=cubeSide, height=cubeSide, channels=3}
	--[[
	now the first time we should render cubemaps and export them
	--]]
	for i=1,6 do
		gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
		
		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glFrustum(-zNear, zNear, -zNear, zNear, zNear, zFar);
		
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()
		
		local suffix
		local side = gl.GL_TEXTURE_CUBE_MAP_POSITIVE_X + i-1
		if side == gl.GL_TEXTURE_CUBE_MAP_POSITIVE_X then
			suffix = 'xp'
			gl.glRotatef(90,0,1,0)
			gl.glRotatef(180,1,0,0)
		elseif side == gl.GL_TEXTURE_CUBE_MAP_NEGATIVE_X then
			suffix = 'xn'
			gl.glRotatef(-90,0,1,0)
			gl.glRotatef(180,1,0,0)
		elseif side == gl.GL_TEXTURE_CUBE_MAP_POSITIVE_Y then
			suffix = 'yp'
			gl.glRotatef(-90,1,0,0)
		elseif side == gl.GL_TEXTURE_CUBE_MAP_NEGATIVE_Y then
			suffix = 'yn'
			gl.glRotatef(90,1,0,0)
		elseif side == gl.GL_TEXTURE_CUBE_MAP_POSITIVE_Z then
			suffix = 'zp'
			gl.glRotatef(180,1,0,0)
		elseif side == gl.GL_TEXTURE_CUBE_MAP_NEGATIVE_Z then
			suffix = 'zn'
			gl.glRotatef(180,0,0,1)
		end
		
		glCallOrRun(projectionDisplayList, drawSkySphere)
		
		gl.glReadPixels(0,0,cubeSide,cubeSide,gl.GL_RGB,gl.GL_UNSIGNED_BYTE, im:dataptr())
		im:save('sky-infrared-cube-'..suffix..'.png')
	end
end

openglapp:run{
	init = function()	
		gl.glEnable(gl.GL_DEPTH_TEST)
		gl.glClearColor(.3, .3, .3, 1)
	end,
	event = function(event)
		if event.type == sdl.SDL_MOUSEMOTION then
			if leftButtonDown then
				local idx = event.motion.xrel
				local idy = event.motion.yrel
				local magn = math.sqrt(idx * idx + idy * idy)
				local dx = idx / magn
				local dy = idy / magn
				local r = Quat():fromAngleAxis(dy, dx, 0, magn)
				viewRot = (r * viewRot):normalize()
			end
		elseif event.type == sdl.SDL_MOUSEBUTTONDOWN then
			if event.button.button == sdl.SDL_BUTTON_LEFT then
				leftButtonDown = true
			elseif event.button.button == sdl.SDL_BUTTON_WHEELUP then
				zoomScale = zoomScale / zoomFactor
			elseif event.button.button == sdl.SDL_BUTTON_WHEELDOWN then
				zoomScale = zoomScale * zoomFactor
			end
		elseif event.type == sdl.SDL_MOUSEBUTTONUP then
			if event.button.button == sdl.SDL_BUTTON_LEFT then
				leftButtonDown = false
			end
		end
	end,
	update = function()
		local viewWidth, viewHeight = openglapp:size()
		local aspectRatio = viewWidth / viewHeight

		gl.glViewport(0,0,viewWidth, viewHeight)
		gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
				
		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glFrustum(-zNear * aspectRatio * tanFovX, zNear * aspectRatio * tanFovX, -zNear * tanFovY, zNear * tanFovY, zNear, zFar);
		
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()
		gl.glTranslatef(0,0,-3)
		
		local aa = viewRot:toAngleAxis()
		gl.glRotatef(aa[4], aa[1], aa[2], aa[3])
		gl.glScalef(zoomScale,zoomScale,zoomScale)
		
		gl.glTranslatef(-1.5,0,0)
		if not projectionDisplayList.id then
			tex = Tex2D.load{
				filename='sky-infrared.png',
				magFilter=gl.GL_LINEAR,
				minFilter=gl.GL_LINEAR_MIPMAP_LINEAR,
				wrap={
					s=gl.GL_CLAMP,
					t=gl.GL_CLAMP,
				}
			}
		end
		glCallOrRun(projectionDisplayList, drawSkySphere)
		
		-- [[
		if not saved then
			saveCubeSides()
			saved = true
		end
		--]]
		
		gl.glTranslatef(3,0,0)
		if not cubeMapDisplayList.id then
			cubeTex = TexCube{
				filenames = {
					'sky-infrared-cube-xp.png',
					'sky-infrared-cube-xn.png',
					'sky-infrared-cube-yp.png',
					'sky-infrared-cube-yn.png',
					'sky-infrared-cube-zp.png',
					'sky-infrared-cube-zn.png',
				},
				wrap={
					s=gl.GL_CLAMP_TO_EDGE,
					t=gl.GL_CLAMP_TO_EDGE,
					r=gl.GL_CLAMP_TO_EDGE,
				},
				magFilter = gl.GL_LINEAR,
				minFilter = gl.GL_LINEAR,
			}
		end
		glCallOrRun(cubeMapDisplayList, drawSkyCube)
	end,
}

print('u range',uMin,uMax)
print('v range',vMin,vMax)
