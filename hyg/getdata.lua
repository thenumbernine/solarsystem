#! /usr/bin/env luajit
require 'ext'
local json = require 'dkjson'
local ffi = require 'ffi'
local csv = require 'csv'

--[[ v2
local csvdata = csv.file'hygxyz.csv'
--]]
-- [[ v3
local csvdata = csv.file'hygdata_v3.csv'
--]]

local rowHeaders = csvdata.rows:remove(1)
for i=1,#rowHeaders do rowHeaders[i] = rowHeaders[i]:lower() end
csvdata:setColumnNames(rowHeaders)

-- sun is positioned relative to the earth, so offset everything for the sun to be the center
local sunRow = csvdata.rows[1]
local sunPos = {
	x = assert(tonumber(sunRow.x)),
	y = assert(tonumber(sunRow.y)),
	z = assert(tonumber(sunRow.z)),
}

local numRows = #csvdata.rows
local numElem = 9
local bufferType = 'float'
local buffer = ffi.new(bufferType..'[?]', numElem * numRows)	-- allocate space for data of each planet 

local maxAbs = 0
local namedStars = {}
local constellations = table()	-- list of unique names so far
local columnCounts = {}

local totalErrors = 0

local ranges = {}
ranges.mag = {}
ranges.colorIndex = {}
ranges.ra = {}
ranges.ra2 = {}
ranges.dec = {}
ranges.dist = {}
ranges.cons = {}

local numValidRows = 0
for i,row in ipairs(csvdata.rows) do
	for _,col in ipairs(csvdata.columns) do
		if #row[col] > 0 then
			columnCounts[col] = (columnCounts[col] or 0) + 1
		end
	end

	--[[ this assertion isn't required, I was just curious if it was enforced
	-- and as of the v3 release it's no longer true
	if numValidRows ~= tonumber(row.id) then
		error("zero-based index "..numValidRows.." not equal to row id "..row.id)
	end
	--]]

	if tonumber(row.dist) == 100000 then
		-- this is the magic number for 'idk'
		-- so skip these entries.  they're not even outside the milky way.  they just aren't filled in.
	else
		-- in parsecs
		local x = assert(tonumber(row.x)) - sunPos.x
		local y = assert(tonumber(row.y)) - sunPos.y
		local z = assert(tonumber(row.z)) - sunPos.z

		-- HYG is in equatorial coordinates
		-- rotate equatorial pos to ecliptic pos
		local epsilon = math.rad(23 + 1/60*(26 + 1/60*(21.4119)))	-- earth tilt
		local cosEps = math.cos(epsilon)
		local sinEps = math.sin(epsilon)
		y, z = cosEps * y + sinEps * z, -sinEps * y + cosEps * z
	
		-- reconstruct spherical to cartesian and see accurate the xyz is
		-- dec is ranged [-90,90], but ra is only [0,23]
-- TODO also why is ursa major pointing away from ursa minor?		
		local ra = math.rad(assert(tonumber(row.ra)))
		local dec = math.rad(assert(tonumber(row.dec)))
		local dist = assert(tonumber(row.dist))
		local rx = dist * math.cos(ra) * math.cos(dec) - sunPos.x
		local ry = dist * math.sin(ra) * math.cos(dec) - sunPos.y
		local rz = dist * math.sin(dec) - sunPos.z
		ry, rz = cosEps * ry + sinEps * rz, -sinEps * ry + cosEps * rz
		local rerr = math.sqrt((x-rx)^2 + (y-ry)^2 + (z-rz)^2)/math.max(1,dist)
		totalErrors = totalErrors + rerr


-- recalc ra, dec
ra = math.atan2(y, x)
dec = math.atan2(z, math.sqrt(x*x + y*y))
dist = math.sqrt(x*x + y*y + z*z)


		ranges.ra.min = math.min(ranges.ra.min or ra, ra)
		ranges.ra.max = math.max(ranges.ra.max or ra, ra)

local ra2 = ra % (2 * math.pi)
		ranges.ra2.min = math.min(ranges.ra2.min or ra2, ra2)
		ranges.ra2.max = math.max(ranges.ra2.max or ra2, ra2)
		
		ranges.dec.min = math.min(ranges.dec.min or dec, dec)
		ranges.dec.max = math.max(ranges.dec.max or dec, dec)
		
		ranges.dist.min = math.min(ranges.dist.min or dist, dist)
		ranges.dist.max = math.max(ranges.dist.max or dist, dist)

		local mag = assert(tonumber(row.absmag))	--absolute magnitude
		ranges.mag.min = math.min(ranges.mag.min or mag, mag)
		ranges.mag.max = math.max(ranges.mag.max or mag, mag)
		
		local colorIndex = tonumber(row.colorindex or row.ci) or .656	-- fill in blanks with something close to white
		ranges.colorIndex.min = math.min(ranges.colorIndex.min or colorIndex, colorIndex)
		ranges.colorIndex.max = math.max(ranges.colorIndex.max or colorIndex, colorIndex)

		-- in parsecs per year
		local vx = assert(tonumber(row.vx))
		local vy = assert(tonumber(row.vy))
		local vz = assert(tonumber(row.vz))


		local name = (#row.proper > 0 and row.proper)
				or (#row.bf > 0 and row.bf)
				or nil
		if name and #name > 0 then
			name = name:trim():gsub('%s+', ' ')
			table.insert(namedStars, {name=name, index=numValidRows})
		end
		
		local con = #row.con > 0 and row.con or nil
		local conInfo, constellationIndex
		if con and #con > 0 then
			constellationIndex, conInfo = constellations:find(nil, function(c) return c.name == con end)
		else
			constellationIndex = 1
			conInfo = constellations[constellationIndex]
		end
			
		if not conInfo then
			conInfo = {
				name = con,
				ra = {},
				ra2 = {},
				dec = {},
				dist = {},
				mag = {},
				indexes = {},	-- array of all indexes in the HYG data of all stars of the constellation 
			}
			constellations:insert(conInfo)
			constellationIndex = #constellations
		end

		conInfo.ra.min = math.min(conInfo.ra.min or ra, ra)
		conInfo.ra.max = math.max(conInfo.ra.max or ra, ra)
		
		conInfo.ra2.min = math.min(conInfo.ra2.min or ra2, ra2)
		conInfo.ra2.max = math.max(conInfo.ra2.max or ra2, ra2)
		
		conInfo.dec.min = math.min(conInfo.dec.min or dec, dec)
		conInfo.dec.max = math.max(conInfo.dec.max or dec, dec)
	
		conInfo.dist.min = math.min(conInfo.dist.min or dist, dist)
		conInfo.dist.max = math.max(conInfo.dist.max or dist, dist)
		
		conInfo.mag.min = math.min(conInfo.mag.min or mag, mag)
		conInfo.mag.max = math.max(conInfo.mag.max or mag, mag)

		-- without individual star info, constellsions.json is 22603 bytes
		-- with individual star info it is 678818 bytes
		table.insert(conInfo.indexes, numValidRows) 
		-- verify original row matches
		--table.insert(conInfo.indexes, {i, row.id, numValidRows})

--print(i, row.id, numValidRows, mag, row.proper)
		buffer[0 + numElem * numValidRows] = x
		buffer[1 + numElem * numValidRows] = y
		buffer[2 + numElem * numValidRows] = z
		buffer[3 + numElem * numValidRows] = vx
		buffer[4 + numElem * numValidRows] = vy
		buffer[5 + numElem * numValidRows] = vz
		buffer[6 + numElem * numValidRows] = mag
		buffer[7 + numElem * numValidRows] = colorIndex 
		buffer[8 + numElem * numValidRows] = constellationIndex-1
		
		maxAbs = math.max(maxAbs, x)
		maxAbs = math.max(maxAbs, y)
		maxAbs = math.max(maxAbs, z)

		numValidRows = numValidRows + 1
	end
end

-- now sort the constellation stars by their magnitude
for _,con in ipairs(constellations) do
	table.sort(con.indexes, function(a,b)
		-- sort by buffer magnitude
		return buffer[6 + numElem * a] < buffer[6 + numElem * b]
		-- sort by original row magnitude
		--return tonumber(csvdata.rows[a[1]].absmag) < tonumber(csvdata.rows[b[1]].absmag)
	end)
	--[=[ if you want to see those values:
	con.magnitudes = table.mapi(con.indexes, function(i)
		-- buffer magnitude
		return buffer[6 + numElem * i]
		-- original row magnitude
		--return tonumber(csvdata.rows[i[1]].absmag)
	end)
	--]=]
end

print('rows processed:', numRows)
print('valid entries:', numValidRows)

print('max reconstruction error',totalErrors/numValidRows)

print('abs max coordinate', maxAbs)
print('num columns provided:')
for _,name in ipairs(csvdata.columns) do
	print('', name, columnCounts[name])
end

print(tolua(ranges))
--[[
for k,v in pairs(ranges) do
	print(k,'min',ranges[k].min,'max',ranges[k].max)
end
--]]

for _,con in ipairs(constellations) do
	local dra = con.ra.max - con.ra.min
	local dra2 = con.ra2.max - con.ra2.min
	if dra2 < dra then
		con.ra = con.ra2
	end
	con.ra2 = nil
end

-- make name of sun consistent with NASA data
-- TODO call it Sol in the NASA data?
assert(namedStars[1].name == 'Sol')
namedStars[1].name = 'Sun'

-- write 
file['stardata.f32'] = ffi.string(buffer, ffi.sizeof(bufferType) * numElem * numValidRows)
file['namedStars.json'] = 'namedStars = ' .. json.encode(namedStars, {indent=true}) ..';'
file['constellations.json'] = 'constellations = '..json.encode(constellations, {indent=true}) .. ';'

--[[ plot dist density map
local f = io.open('dist-distribution.txt', 'w')
f:write'#dist_min\tdist_max\tcount\n'
local bins = 2000
local counts = range(bins):map(function(i) return 0 end)
for i,row in ipairs(csvdata.rows) do
	-- [=[ bin by dist
	local dist = assert(tonumber(row.dist))
	--]=]
	--[=[ bin by xyz dist
	local x = assert(tonumber(row.x)) - sunPos.x
	local y = assert(tonumber(row.y)) - sunPos.y
	local z = assert(tonumber(row.z)) - sunPos.z
	local dist = math.sqrt(x*x + y*y + z*z)
	--]=]

	local bin = 1 + math.floor(bins * (dist / ranges.dist.max))
	if bin >= 1 and bin <= bins then
		counts[bin] = counts[bin] + 1
	end
end
for bin,count in ipairs(counts) do
	local distMin = (bin-1) / bins * ranges.dist.max
	local distMax = bin / bins * ranges.dist.max
	f:write(distMin,'\t',distMax,'\t',count,'\n')
end
f:close()
--]]
