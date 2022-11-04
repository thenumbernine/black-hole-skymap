#!/usr/bin/env lua
local range = require 'ext.range'
local tolua = require 'ext.tolua'
local gnuplot = require 'gnuplot'

local M = 1
local xmin = 0
local xmax = 10
local ymin = 0
local ymax = 10
local n = 1000
local rhos = range(n):mapi(function(i) return (i-1)/n*(xmax-xmin) + xmin end)
local rs = rhos:mapi(function(rho)
	return rho * (1 + M / (2 * rho))^2
end)

-- calculate an inverse using Newtons method
local function calcInv(y, f, df_dx)
	local x = y
	local epsilon = 1e-7
	local maxiter = 1000
	for iter=1,maxiter do
		local dx = f(x) / df_dx(x)
		x = x - dx
		if math.abs(dx) < epsilon then return x, iter end
	end
	return x, maxiter
end

-- now to solve rho as a function of r
-- r(rho) = rho * (1 + M / (2 * rho))^2
-- r(rho) - rho * (1 + M / (2 * rho))^2 = 0
-- y = f(x)
-- y - f(x) = 0
-- finv(y) = x
local rhosFromRs = rs:mapi(function(r)
	local function r_of_rho(rho)
		return rho * (1 + M / (2 * rho))^2 - r
	end
	local function dr_drho(rho)
		return 1 - M * M / (4 * rho * rho)
	end
	--return (select(2, calcInv(r, r_of_rho, dr_drho)))	-- how many iterations to converge?  could be up to 30...
	return (calcInv(r, r_of_rho, dr_drho))
end)
-- y = x^2 = const
-- x^2 - y = 0 for x = sqrt(y)

gnuplot{
	persist = true,
	xrange = {xmin, xmax},
	--yrange = {ymin, ymax},
	data = {rhos, rs, rhosFromRs},
	xtitle = 'rho',
	ytitle = 'r',
	style = 'data lines',
	{using='1:2', notitle=true},
	{using='3:2', notitle=true},
}
