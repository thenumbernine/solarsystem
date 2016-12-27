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
local numElem = 8
local bufferType = 'float'
local buffer = ffi.new(bufferType..'[?]', numElem * numRows)	-- allocate space for data of each planet 

local maxAbs = 0
local namedStars = {}
local columnCounts = {}

local totalErrors = 0

local range = {}
range.mag = {}
range.colorIndex = {}
range.dist = {}

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
		local ra = assert(tonumber(row.ra))
		local dec = assert(tonumber(row.dec))
		local dist = assert(tonumber(row.dist))
		local rx = dist * math.cos(ra) * math.cos(dec) - sunPos.x
		local ry = dist * math.sin(ra) * math.cos(dec) - sunPos.y
		local rz = dist * math.sin(dec) - sunPos.z
		ry, rz = cosEps * ry + sinEps * rz, -sinEps * ry + cosEps * rz
		local rerr = math.sqrt((x-rx)^2 + (y-ry)^2 + (z-rz)^2)/math.max(1,dist)
		totalErrors = totalErrors + rerr
		
		range.dist.min = math.min(range.dist.min or dist, dist)
		range.dist.max = math.max(range.dist.max or dist, dist)

		local mag = assert(tonumber(row.absmag))	--apparent magnitude
		range.mag.min = math.min(range.mag.min or mag, mag)
		range.mag.max = math.max(range.mag.max or mag, mag)
		
		local colorIndex = tonumber(row.colorindex or row.ci) or .656	-- fill in blanks with something close to white
		range.colorIndex.min = math.min(range.colorIndex.min or colorIndex, colorIndex)
		range.colorIndex.max = math.max(range.colorIndex.max or colorIndex, colorIndex)

		-- in parsecs per year
		local vx = assert(tonumber(row.vx))
		local vy = assert(tonumber(row.vy))
		local vz = assert(tonumber(row.vz))
		
		buffer[0 + numElem * numValidRows] = x
		buffer[1 + numElem * numValidRows] = y
		buffer[2 + numElem * numValidRows] = z
		buffer[3 + numElem * numValidRows] = vx
		buffer[4 + numElem * numValidRows] = vy
		buffer[5 + numElem * numValidRows] = vz
		buffer[6 + numElem * numValidRows] = mag
		buffer[7 + numElem * numValidRows] = colorIndex 
		
		maxAbs = math.max(maxAbs, x)
		maxAbs = math.max(maxAbs, y)
		maxAbs = math.max(maxAbs, z)

		local properName = row.proper
		if properName and #properName > 0 then
			table.insert(namedStars, {name=properName, index=numValidRows})
		end
		
		numValidRows = numValidRows + 1
	end
end

print('rows processed:', numRows)
print('valid entries:', numValidRows)

print('max reconstruction error',totalErrors/numValidRows)

print('abs max coordinate', maxAbs)
print('num columns provided:')
for _,name in ipairs(csvdata.columns) do
	print('', name, columnCounts[name])
end

for k,v in pairs(range) do
	print(k,'min',range[k].min,'max',range[k].max)
end

-- make name of sun consistent with NASA data
assert(namedStars[1].name == 'Sol')
namedStars[1].name = 'Sun'

-- write 
file['stardata.f32'] = ffi.string(buffer, ffi.sizeof(bufferType) * numElem * numValidRows)
file['namedStars.json'] = 'namedStars = ' .. json.encode(namedStars, {indent=true}) ..';'

-- [[ plot dist density map
local f = io.open('dist-distribution.txt', 'w')
f:write'#dist_min\tdist_max\tcount\n'
local bins = 2000
local counts = require 'ext.range'(bins):map(function(i) return 0 end)
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

	local bin = 1 + math.floor(bins * (dist / range.dist.max))
	if bin >= 1 and bin <= bins then
		counts[bin] = counts[bin] + 1
	end
end
for bin,count in ipairs(counts) do
	local distMin = (bin-1) / bins * range.dist.max
	local distMax = bin / bins * range.dist.max
	f:write(distMin,'\t',distMax,'\t',count,'\n')
end
f:close()
--]]
--[[ 
dist binned within 0.5 parsec:
0	0.4950495	1
0.4950495	0.990099	0
0.990099	1.4851485	3
1.4851485	1.980198	1
1.980198	2.4752475	1

|x,y,z| binned within 0.5 parsecs:
0	0.4950495	1
0.4950495	0.990099	0
0.990099	1.4851485	3
1.4851485	1.980198	1
1.980198	2.4752475	1

--]]
