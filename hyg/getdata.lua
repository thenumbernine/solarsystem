#! /usr/bin/env luajit
require 'ext'
local json = require 'dkjson'
local ffi = require 'ffi'
local csv = require 'csv'

local csvdata = csv.file 'hygdata_v3.csv'
csvdata:setColumnNames(csvdata.rows:remove(1))

-- sun is positioned relative to the earth, so offset everything for the sun to be the center
local sunRow = csvdata.rows[1]
local sunPos = {
	x = assert(tonumber(sunRow.x)),
	y = assert(tonumber(sunRow.y)),
	z = assert(tonumber(sunRow.z)),
}

local numRows = #csvdata.rows
local numElem = 8
local buffer = ffi.new('float[?]', numElem * numRows)	-- allocate space for data of each planet 

local maxAbs = 0
local namedStars = {}
local columnCounts = {}

local totalErrors = 0

local zeroBasedIndex = 0
for i,row in ipairs(csvdata.rows) do
	for _,col in ipairs(csvdata.columns) do
		if #row[col] > 0 then
			columnCounts[col] = (columnCounts[col] or 0) + 1
		end
	end

	--[[ this assertion isn't required, I was just curious if it was enforced
	-- and as of the v3 release it's no longer true
	if zeroBasedIndex ~= tonumber(row.id) then
		error("zero-based index "..zeroBasedIndex.." not equal to row id "..row.id)
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
		local ra = assert(tonumber(row.ra))
		local dec = assert(tonumber(row.dec))
		local dist = assert(tonumber(row.dist))
		local rx = dist * math.cos(ra) * math.cos(dec) - sunPos.x
		local ry = dist * math.sin(ra) * math.cos(dec) - sunPos.y
		local rz = dist * math.sin(dec) - sunPos.z
		ry, rz = cosEps * ry + sinEps * rz, -sinEps * ry + cosEps * rz
		local rerr = math.sqrt((x-rx)^2 + (y-ry)^2 + (z-rz)^2)/math.max(1,dist)
		totalErrors = totalErrors + rerr

		local mag = assert(tonumber(row.absmag))	--apparent magnitude
		local colorIndex = tonumber(row.ci) or .656	-- fill in blanks with something close to white

		local vx = assert(tonumber(row.vx))
		local vy = assert(tonumber(row.vy))
		local vz = assert(tonumber(row.vz))
		
		buffer[0 + numElem * zeroBasedIndex] = x
		buffer[1 + numElem * zeroBasedIndex] = y
		buffer[2 + numElem * zeroBasedIndex] = z
		buffer[3 + numElem * zeroBasedIndex] = vx
		buffer[4 + numElem * zeroBasedIndex] = vy
		buffer[5 + numElem * zeroBasedIndex] = vz
		buffer[6 + numElem * zeroBasedIndex] = mag
		buffer[7 + numElem * zeroBasedIndex] = colorIndex 
		
		maxAbs = math.max(maxAbs, x)
		maxAbs = math.max(maxAbs, y)
		maxAbs = math.max(maxAbs, z)
		
		if #row.proper > 0 then
			table.insert(namedStars, {name=row.proper, index=zeroBasedIndex})
		end
		
		zeroBasedIndex = zeroBasedIndex + 1
	end
end

print('max reconstruction error',totalErrors/zeroBasedIndex)

print('abs max coordinate', maxAbs)
print('num columns provided:')
for _,name in ipairs(csvdata.columns) do
	print('', name, columnCounts[name])
end

-- make name of sun consistent with NASA data
assert(namedStars[1].name == 'Sol')
namedStars[1].name = 'Sun'

-- write 
file['stardata.f32'] = ffi.string(buffer, ffi.sizeof(buffer))
file['namedStars.json'] = 'namedStars = ' .. json.encode(namedStars, {indent=true}) ..';'

