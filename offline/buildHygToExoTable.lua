#!/usr/bin/env luajit

-- has to go before any line comments
print'{'

require 'ext'

local CSV = require 'csv'
print'-- reading hyg'
local hyg = CSV.file'../hyg/hygdata_v3.csv'
print'-- done reading hyg'
hyg:setColumnNames(hyg.rows:remove(1))

local rowForHip = {}
local rowForProper = {}
local rowForFlamsteed = {}
for _,row in ipairs(hyg.rows) do
	if row.hip ~= '' then
		rowForHip[row.hip] = row
	end
	if row.proper ~= '' then
		rowForProper[row.proper] = row
	end
	if row.flam ~= ''
	and row.con ~= ''
	then
		local flamsteed = row.flam..' '..row.con
		rowForFlamsteed[flamsteed] = row 
	end
end

local exoplanets = fromlua(file'../exoplanet/openExoplanetCatalog.lua':read())

-- [[ make sure there's no overlap of names
--[=[ duplicate names I found:

"SWEEPS-11" "SWEEPS-4", "BD+20 2184", "2MASS J08414382+2013368", "TrES-1", "TrES-3", "TrES-5", "OGLE-2016-BLG-1928", "TOI-1640", "TOI-130"
are names for both the star and its orbiting planet.  The planet always has an alternative name suffix 'b' ... but the star never has a suffix 'A'
does suffix 'A' imply a multiple-star system?

"TrES-2", "TrES-4", "WASP-2", "GJ 229"
are a planets and a *barycenter* with multiple stars, and the two stars around it are labeled A and B
"WASP-2"'s main star has no 'A' suffix, but still orbits the barycenter, so both the star and the barycenter have a matching name.

"Kepler-138 c", "KIC 7603200 b"
seems like a big mixup of 'b's and 'c's between different exoplanet classifications. 
--]=]
local allSystemNames = {}
local allBodyNames = {}
for _,exoplanet in ipairs(exoplanets) do
	local function addNames(names, to, what)
		for _,name in ipairs(names) do
			if to[name] then
				print("-- found duplicate of "..what..' name '..name)
			end
			to[name] = true
		end
	end
	addNames(exoplanet.names, allSystemNames, 'system')
	-- so system names and body names can overlap
	if exoplanet.bodies then
		for _,body in ipairs(exoplanet.bodies) do
			if body.names then
				addNames(body.names, allBodyNames, 'body')
			end
		end
	end
end
--]]


local function findHygForName(name, obj)
	local row = rowForProper[name] 
	if row then return row end

	row = rowForFlamsteed[name]
	if row then return row end

	local parts = name:split' '
	if parts[1] == 'HIP' then
		-- parts[2] should be a number
		local hip = parts[2]
		
		local base, sub = hip:match'^(%d+)/(%d+)$'
		if base then
			return base
		end

		if not hip:match'^%d+$' then
			error("got a HIP with a non-number "..name)
		end
		-- for a single-star system optional parts[3] should be a-z 
		-- for a multi-star system optional parts[3] should be A-Z and optional parts[4] should be a-z
		assert(#parts >= 2 and #parts <= 4)
		if #parts >= 3 then 
			-- can be multi-letter i.e. multi-star barycenter
			if not parts[3]:match'^%a+$' then error("got a HIP star non-letter "..name) end
			if #parts[3] > 1 then assert(obj.type == 'barycenter') end
		end
		if #parts >= 4 then 
			if not parts[3]:match'^%a+$' then error("got a HIP planet non-letter "..name) end
			if #parts[4] > 1 then assert(obj.type == 'barycenter') end
		end
	-- else
	-- sometimes a system's name will match a HYG name,
	-- and its planets will match HIP rows, 
	-- but the system itself won't 
	-- so TODO look for the name among other HYG info
		return rowForHip[hip]
	end
end

local function findHygForObj(obj)
	if obj.names then
		for _,name in ipairs(obj.names) do
			local row = findHygForName(name, obj)
			if row then return row end
		end
	end

	if obj.bodies then
		for _,body in ipairs(obj.bodies) do
			local row = findHygForObj(body)
			if row then return row end
		end
	end
end

--[[
would be nice to link this back to the float buffer that i'm using ...
how about instead of storing the constellation index, I store the hyg index? 
then for constellations look it up in the table?
honestly, nobody cares what constellation it is.  a better test would be to do a bounds check , and you can find the ra/dec grid boundaries of constellations if you really want it.
--]]
print('-- exoname\thyg')
local numFound = 0
for _,exoplanet in ipairs(exoplanets) do
	table.sort(exoplanet.names)
	exoplanet.name = exoplanet.names[1]
end
table.sort(exoplanets, function(a,b)
	return a.name < b.name
end)
for _,system in ipairs(exoplanets) do
	local hyg = findHygForObj(system)
	print('\t'..tolua{
		system.name,		-- luckily there's no duplicate system names, only duplicate body names
		hyg and hyg.id,
		--hyg and hyg.hip,
	}..',')
	if hyg then
		numFound = numFound + 1 
	end
end
print'}'
print('-- found '..numFound..' of '..#exoplanets..' systems')
