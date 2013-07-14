--[[
no gl required, but since the dest coords are mapped to the source coords, parts in the source are skipped, esp in regions of high intrinsic curvature
--]]
local args = {...}
local vec3 = require 'vec.vec3'
local Quat = require 'vec.quat'
local Image = require 'image'
require 'ext'

local destFilenamePrefix = 'sky-visible-cube-'
local destImageWidth = 1024
local destImageHeight = 1024
--[[ for testing algorithms
local destFilenamePrefix = 'sky-visible-cube-testsmall-'
local destImageWidth = 128
local destImageHeight = 128
--]]

local SQRT_1_2 = math.sqrt(.5)
local sides = table{
	{name='xp', angle=Quat(SQRT_1_2,0,SQRT_1_2,0)},
	{name='xn', angle=Quat(SQRT_1_2,0,-SQRT_1_2,0)},
	{name='yp', angle=Quat(-SQRT_1_2,0,0,SQRT_1_2)},
	{name='yn', angle=Quat(SQRT_1_2,0,0,SQRT_1_2)},
	{name='zp', angle=Quat(1,0,0,0)},
	{name='zn', angle=Quat(0,0,1,0)},
}

local useSides = table(args)
for _,useSide in ipairs(useSides) do
	if not sides:find(nil, function(side) return side.name == useSide end) then
		error('specified to build unknown side '..useSide)
	end
end
if #useSides == 0 then
	useSides = sides:map(function(side) return side.name end)
end

local srcImg = Image('eso0932a.png')
local srcImgWidth, srcImgHeight
srcImgWidth, srcImgHeight = srcImg:size()

-- build our filter
local function sinc(x) return x == 0 and 1 or math.sin(x) / x end
local filterSize = 3
local filter = {}
for i=1,2*filterSize+1 do
	local x = i - filterSize - 1
	if x == filterSize+1 then
		filter[i] = 1
	else
		filter[i] = sinc(x) * sinc(x/filterSize)
	end
end

local M_LOG_2 = math.log(2)
local function log2(x)
	return math.log(x) / M_LOG_2
end

local startTime = os.time()
local srcMipMaps = {srcImg}
do
	local mipMapWidth = math.floor(srcImgWidth / 2)
	local mipMapHeight = math.floor(srcImgHeight / 2)
	local lastMipMap = srcMipMaps[#srcMipMaps]
	local lastMipMapWidth, lastMipMapHeight = lastMipMap:size()
	for level=2,log2(math.max(srcImgWidth,srcImgHeight)) do
		local mipMapFilename = 'mipmap-'..level..'.png'
		if mipMapWidth == 0 or mipMapHeight == 0 then break end
		local newMipMap
		if io.fileexists(mipMapFilename) then
			newMipMap = Image(mipMapFilename)
		else
			io.write('saving '..mipMapFilename..'     ') io.flush()
			newMipMap = Image(mipMapWidth, mipMapHeight)
			-- TODO filter	
			for dstj=0,mipMapHeight-1 do
				for dsti=0,mipMapWidth-1 do
					local srci = math.floor(dsti*(lastMipMapWidth-1)/(mipMapWidth-1))
					local srcj = math.floor(dstj*(lastMipMapHeight-1)/(mipMapHeight-1))
					local tr, tg, tb, tw = 0, 0, 0, 0
					for fsj = math.max(srcj - filterSize, 0), math.min(srcj + filterSize, lastMipMapHeight-1) do
						for fsi = math.max(srci - filterSize, 0), math.min(srci + filterSize, lastMipMapWidth-1) do
							local r,g,b = lastMipMap(fsi,fsj)
							local w = filter[fsi - srci + filterSize + 1] * filter[fsj - srcj + filterSize + 1]
							tr = tr + r * w
							tg = tg + g * w
							tb = tb + b * w
							tw = tw + w
						end
					end
					if math.abs(tw) > 1e-9 then
						tr = tr / tw
						tg = tg / tw
						tb = tb / tw
					end
					newMipMap(dsti,dstj,tr,tg,tb)
					local thisTime = os.time()
					if thisTime ~= startTime then
						startTime = thisTime
						local percent = 100 * (dsti + mipMapWidth * dstj) / (mipMapWidth * mipMapHeight)
						io.write(("\008\008\008\008%3.f%%"):format(percent)) io.flush()
					end
				end
			end
			print()
			table.insert(srcMipMaps,newMipMap)
			mipMapWidth = math.floor(mipMapWidth / 2)
			mipMapHeight = math.floor(mipMapHeight / 2)
			newMipMap:save(mipMapFilename)
		end
	end
end

local automagically  = false	-- differentiation magnitude <=> miplevel is coming out constant across the whole image.  no good.
if automagically then
	require 'symmath'
	symmath.simplifyConstantPowers = true
end

--[[
mipmapping:
d/d[destx] ( srcu, srcv )
d/d[desty] ( srcu, srcv )
level = .5 log2(max( |d/d[destx]|^2, |d/d[desty]|^2 ))
--]]

local minReadMipLevel, maxReadMipLevel
local startTime = os.time()
for sideIndex,side in ipairs(sides) do
	if useSides:find(side.name) then
		print('building side '..side.name)
		local srciFunc
		local srcjFunc
		local dSrciDiFunc
		local dSrciDjFunc
		local dSrcjDiFunc
		local dSrcjDjFunc
		if automagically then
			local i = symmath.Variable('i')
			local j = symmath.Variable('j')
			local coords = {i,j}
			local u = (i + .5) / destImageWidth
			local v = (j + .5) / destImageHeight
			local vec = {2 * u - 1, 2 * v - 1, symmath.Constant(1)}
			-- I don't trust the quaternion rotate function to work on symmath so straightforward ...  maybe it would?
			local mat = side.angle:toMatrix()
			local rotvec = {}
			rotvec[1] = mat[1][1] * vec[1] + mat[2][1] * vec[2] + mat[3][1] * vec[3]
			rotvec[2] = mat[1][2] * vec[1] + mat[2][2] * vec[2] + mat[3][2] * vec[3]
			rotvec[3] = mat[1][3] * vec[1] + mat[2][3] * vec[2] + mat[3][3] * vec[3]
			local x, y, z = unpack(rotvec)
			local r = symmath.sqrt(x^2 + y^2 + z^2)
			local phi = symmath.atan2(y,x)
			local theta = symmath.acos(z/r)
			local srcu = phi / math.pi * .5 + .5
			local srcv = theta / math.pi
			-- right now, for accuracy/appearance's sake, division and constants are not simplified
			-- how about a flag to circumvent this?
			local srci = symmath.simplify(srcu * (srcImgWidth - 1))
			local srcj = symmath.simplify(srcv * (srcImgHeight - 1))
			local dSrciDi = symmath.simplify(symmath.diff(srci, i))
			local dSrciDj = symmath.simplify(symmath.diff(srci, j))
			local dSrcjDi = symmath.simplify(symmath.diff(srcj, i))
			local dSrcjDj = symmath.simplify(symmath.diff(srcj, j))
			local cmd
			srciFunc, cmd = symmath.compile(srci, coords)	-- still has to be clamped ...
			print('compiling srci as '..cmd)
			srcjFunc, cmd = symmath.compile(srcj, coords)
			print('compiling srcj as '..cmd)
			dSrciDiFunc, cmd = symmath.compile(dSrciDi, coords)
			print('compiling dSrciDi as '..cmd)
			dSrciDjFunc, cmd = symmath.compile(dSrciDj, coords)
			print('compiling dSrciDj srcj as '..cmd)
			dSrcjDiFunc, cmd = symmath.compile(dSrcjDi, coords)
			print('compiling dSrcjDi as '..cmd)
			dSrcjDjFunc, cmd = symmath.compile(dSrcjDj, coords)
			print('compiling dSrcjDj as '..cmd)
		end
		
		local dstImgFn = destFilenamePrefix..side.name..'.png'
		io.write('saving '..dstImgFn..'     ') io.flush()
		local dstImg = Image(destImageWidth, destImageHeight)
		for j=0,destImageHeight-1 do
			local v = (j+.5)/destImageHeight
			for i=0,destImageWidth-1 do
				local srci, srcj, dSrciDi, dSrciDj, dSrcjDi, dSrcjDj 	
				if not automagically then -- do the math by hand
					local u = (i+.5)/destImageWidth
					-- start in screen space, 90 degree fov
					local vec = vec3(2*u-1,2*v-1,1)
					-- rotate to cube side in world space
					vec = side.angle:rotate(vec)
					-- pick apart spherical coordinates
					local r = vec:length()
					local x,y,z = unpack(vec)
					local phi = math.atan2(y,x)
					local theta = math.acos(z/r)
					local srcu = math.clamp(phi / math.pi *.5 + .5, 0, 1)
					local srcv = math.clamp(theta / math.pi, 0, 1)
					--[[
					local du_dx = d/dx ( atan(y/x) / math.pi * .5 )
					local dv_dx = d/dx ( acos(z/r) / math.pi )
					local d_dx_sq = du_dx * du_dx + dv_dx * dv_dx
					local du_dy = d/dy ( atan(y/x) / math.pi * .5 )
					local dv_dy = d/dy ( acos(z/r) / math.pi )
					local d_dy_sq = du_dy * du_dy + dv_dy * dv_dy
					--]]
					srci = math.floor(srcu * (srcImgWidth - 1))
					srcj = math.floor(srcv * (srcImgHeight - 1))
					--[[
					dSrciDi = d(math.atan2(y,x))/di .5 / math.pi * (srcImgWidth - 1)
					dSrciDj = d(math.atan2(y,x))/di .5 / math.pi * dsrcu/dj * (srcImgWidth - 1)
					dSrcjDi = d(math.acos(z/r))/di / math.pi * (srcImgHeight - 1)
					dSrcjDj = d(math.acos(z/r))/dj / math.pi * (srcImgHeight - 1)
					--]]
				else -- automagically
					srci = math.clamp(0, math.floor(srciFunc(i,j)), srcImgWidth-1)
					srcj = math.clamp(0, math.floor(srcjFunc(i,j)), srcImgHeight-1)
					dSrciDi = dSrciDiFunc(i,j)
					dSrciDj = dSrciDjFunc(i,j)
					dSrcjDi = dSrcjDiFunc(i,j)
					dSrcjDj = dSrcjDjFunc(i,j)
				end
				--[[
				local d_di_lenSq = dSrciDi * dSrciDi + dSrcjDi * dSrcjDi
				local d_dj_lenSq = dSrciDj * dSrciDj + dSrcjDj * dSrcjDj
				local readMipLevel = .5 * log2(math.max(d_di_lenSq, d_dj_lenSq))
				--]]
				
				minReadMipLevel = minReadMinLevel and math.min(minReadMipLevel, readMipLevel) or readMipLevel
				maxReadMipLevel = maxReadMinLevel and math.max(maxReadMipLevel, readMipLevel) or readMipLevel
				dstImg(i,j,srcImg(srci,srcj))
				local thisTime = os.time()
				if thisTime ~= startTime then
					startTime = thisTime
					local percent = 100 * (i + destImageWidth * j) / (destImageWidth * destImageHeight)
					io.write(("\008\008\008\008%3.f%%"):format(percent)) io.stdout:flush()
				end
			end
		end
		print()
		-- currently using sdl_image, which I have only writing to bmp
		dstImg:save(dstImgFn)
	end
end

print('readMipLevel min',minReadMipLevel,'max',maxReadMipLevel)

