#!/usr/bin/env luajit
--[[
horizons is in ecliptic frame
ephemeris is in equatorial frame
the transform from ecliptic to equatorial is:

	[ 1 0 0 ]
	[ 0 cos(h) -sin(h) ]
	[ 0 sin(h)  cos(h) ]
	
	for h = 23.4365472133 degrees
--]]

local path = require 'ext.path'
local table = require 'ext.table'
local range = require 'ext.range'
local string = require 'ext.string'
local json = require 'dkjson'

local matrix_ffi = require 'matrix.ffi'
matrix_ffi.real = 'double' 

local function jsonstrtovec3(t)
	return matrix_ffi(table.mapi(t, function(x)
		return tonumber(string.trim(x))
	end))
end

local function vec3tojsonstr(...)
	return table{...}:mapi(function(x) return tostring(x) end):setmetatable(nil)
end

local dynamicVarsFileName = 'dynamic-vars-orig.json'
local dynamicVarsJSON = assert(path(dynamicVarsFileName):read())	-- really is js, has a var = at the front
local dynamicVars = assert(json.decode((dynamicVarsJSON:match'^.-=(.*)$')))

local _, hsun = table.find(dynamicVars.coords, nil, function(o) return o.name:lower() == 'sun' end)
assert(hsun, "couldn't find sun in horizons dynamicVars")

local earthTiltInDeg = -23.4365472133
local cosTilt = math.cos(math.rad(earthTiltInDeg))
local sinTilt = math.sin(math.rad(earthTiltInDeg))
local EclipToEquat = matrix_ffi{
	{1,0,0},
	{0,cosTilt,-sinTilt},
	{0,sinTilt,cosTilt},
}

for _,o in ipairs(dynamicVars.coords)  do
	local hpos = jsonstrtovec3(o.pos)
	local hvel = jsonstrtovec3(o.vel)
--	hpos = hpos + jsonstrtovec3(hsun.pos)
--	hvel = hvel + jsonstrtovec3(hsun.vel)
	local hposInEquat = EclipToEquat * hpos
	local hvelInEquat = EclipToEquat * hvel

	print(
		(hposInEquat - hpos):norm() / hpos:norm(), 
		(hvelInEquat - hvel):norm() / hvel:norm()
	) 

	o.pos = vec3tojsonstr(hposInEquat:unpack())
	o.vel = vec3tojsonstr(hvelInEquat:unpack())
end
local jsonData = 'horizonsDynamicData = '..json.encode(dynamicVars, {indent=true})..';'

path'dynamic-vars.json':write(jsonData)
