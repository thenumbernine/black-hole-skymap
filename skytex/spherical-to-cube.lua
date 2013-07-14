local vec3 = require 'vec.vec3'
local Quat = require 'vec.quat'
local Image = require 'image'

local dstImgWidth = 1024
local dstImgHeight = 1024

local SQRT_1_2 = math.sqrt(.5)
local sides = {
	{name='xp', angle=Quat(SQRT_1_2,0,SQRT_1_2,0)},
	{name='xn', angle=Quat(SQRT_1_2,0,-SQRT_1_2,0)},
	{name='yp', angle=Quat(-SQRT_1_2,0,0,SQRT_1_2)},
	{name='yn', angle=Quat(SQRT_1_2,0,0,SQRT_1_2)},
	{name='zp', angle=Quat(1,0,0,0)},
	{name='zn', angle=Quat(0,0,1,0)},
}

local srcImg = Image('eso0932a.png')
local srcImgWidth, srcImgHeight
srcImgWidth, srcImgHeight = srcImg:size()

local startTime = os.time()
io.write('    ')
for sideIndex,side in ipairs(sides) do
	local dstImg = Image(dstImgWidth, dstImgHeight)
	for j=0,dstImgHeight-1 do
		local v = (j+.5)/dstImgHeight
		for i=0,dstImgWidth-1 do
			local u = (i+.5)/dstImgWidth
			-- start in screen space, 90 degree fov
			local vec = vec3(2*u-1,2*v-1,1)
			-- rotate to cube side in world space
			vec = side.angle:rotate(vec)
			-- pick apart spherical coordinates
			local r = vec:length()
			local x,y,z = unpack(vec)
			local phi = math.atan2(y,x)
			local theta = math.acos(z,r)
			local srcu = math.clamp(phi / math.pi *.5 + .5, 0, 1)
			local srcv = math.clamp(theta / math.pi, 0, 1)
			local srci = math.floor(srcu * (srcImgWidth - 1))
			local srcj = math.floor(srcv * (srcImgHeight - 1))
			dstImg(i,j,srcImg(srci,srcj))
			local thisTime = os.time()
			if thisTime ~= startTime then
				startTime = thisTime
				local percent = 100 * (i + dstImgWidth * (j + dstImgHeight)) / (dstImgWidth * dstImgHeight)
				io.write(("\008\008\008\008%3.f%%"):format(percent))
				io.stdout:flush()
			end
		end
	end
	-- currently using sdl_image, which I have only writing to bmp
	local dstImgFn = 'sky-visible-cube-'..side.name..'.png'
	print(' saving',dstImgFn)
	dstImg:save(dstImgFn)
end
