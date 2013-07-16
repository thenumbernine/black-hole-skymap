--[[
final winrar i think:
[0, -SQRT_1_2, 0, -SQRT_1_2],
[0, SQRT_1_2, 0, -SQRT_1_2],
[SQRT_1_2, 0, 0, -SQRT_1_2],
[SQRT_1_2, 0, 0, SQRT_1_2],
[0, 0, 0, -1],
[0, -1, 0, 0]
--]]

local Quat = require 'vec.quat'

local checks = {['0']=0,SQRT_1_2=math.sqrt(.5),['1']=1}
local function qprint(q)
	for i=1,4 do
		for k,v in pairs(checks) do
			if math.abs(q[i] - v) < 1e-5 then q[i] = k break end
			if math.abs(-q[i] - v) < 1e-5 then q[i] = '-'..k break end
		end
	end
	print('['..table.concat(q,',')..'],')
end


-- I should add these to gl-util.js GL.TextureCube somewhere ...

--xp
qprint(Quat:fromAngleAxis(1,0,0,180) * Quat(math.sqrt(.5),0,math.sqrt(.5),0))
--xn
qprint(Quat:fromAngleAxis(1,0,0,180) * Quat(math.sqrt(.5),0,-math.sqrt(.5),0))
--yp
qprint(Quat:fromAngleAxis(1,0,0,180) * Quat(math.sqrt(.5),0,0,math.sqrt(.5)))
--yn
qprint(Quat:fromAngleAxis(1,0,0,180) * Quat(-math.sqrt(.5),0,0,math.sqrt(.5)))
--zp
qprint(Quat:fromAngleAxis(1,0,0,180) * Quat(1,0,0,0))
--zn
qprint(Quat:fromAngleAxis(1,0,0,180) * Quat(0,0,1,0))

