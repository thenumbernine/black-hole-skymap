local Quat = require 'vec.quat'

local checks = {['0']=0,sqrt1_2=math.sqrt(.5),['1']=1}
local function qprint(q)
	for i=1,4 do
		for k,v in pairs(checks) do
			if math.abs(q[i] - v) < 1e-5 then q[i] = k break end
			if math.abs(-q[i] - v) < 1e-5 then q[i] = '-'..k break end
		end
	end
	print('['..table.concat(q,',')..'],')
end

--xp
qprint((Quat():fromAngleAxis(1,0,0,0) * Quat():fromAngleAxis(0,1,0,90)))
-- SQRT_1_2,0,SQRT_1_2,0
--xn
qprint((Quat():fromAngleAxis(1,0,0,0) * Quat():fromAngleAxis(0,1,0,-90)))
-- SQRT_1_2,0,-SQRT_1_2,0
--yp
qprint(Quat():fromAngleAxis(1,0,0,-90))
-- SQRT_1_2,-0,-0,SQRT_1_2
--yn
qprint(Quat():fromAngleAxis(1,0,0,90))
-- SQRT_1_2,0,0,SQRT_1_2
--zp
qprint(Quat():fromAngleAxis(1,0,0,180))
-- 1,0,0,0
--zn
qprint(Quat():fromAngleAxis(0,0,1,180))
-- 0,0,1,0

