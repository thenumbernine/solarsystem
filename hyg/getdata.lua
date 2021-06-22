#! /usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
local CSV = require 'csv'
local Stat = require 'stat'
local StatSet = require 'stat.set'

local csv = CSV.file'hygdata_v3.csv'

local outputLum = ... == 'lum'

local offsetFromSun = true
local rotateEquatorialToEcliptic = true

--[=[
M = -2.5 log(LStarOverLSun * LSunOverL0)/log(10)
LStarOverLSun = 10^(-M/2.5) / LSunOverL0
--]=]


local rowHeaders = csv.rows:remove(1)
for i=1,#rowHeaders do rowHeaders[i] = rowHeaders[i]:lower() end
csv:setColumnNames(rowHeaders)

-- sun is positioned relative to the earth, so offset everything for the sun to be the center
local sunRow = csv.rows[1]
local sunPos = {
	x = assert(tonumber(sunRow.x)),
	y = assert(tonumber(sunRow.y)),
	z = assert(tonumber(sunRow.z)),
}

local numRows = #csv.rows
local numElem = 9
local bufferType = 'float'
local buffer = ffi.new(bufferType..'[?]', numElem * numRows)	-- allocate space for data of each planet 

local maxAbsPos = 0
local namedStars = {}
local constellations = table()	-- list of unique names so far
local columnCounts = {}


local statset = StatSet(
	'absmag',	-- absolute magnitude
	'mag',		-- apparent magnitude
	'lum',		-- luminosity, in LSun (even though the sun isn't 1.0 ... )
	'colorIndex',
	'temp',
	'ra',
	'ra2',
	'dec',
	'dist',
	'posSphereError',
	'postEclipticDistError'
)

local numValidRows = 0
local numInvalidRows = 0
for i,row in ipairs(csv.rows) do
	for _,col in ipairs(csv.columns) do
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
		numInvalidRows = numInvalidRows + 1
	else
		-- in parsecs
		local x = assert(tonumber(row.x))
		local y = assert(tonumber(row.y))
		local z = assert(tonumber(row.z))
		
		-- in hour-angles
		local ra = assert(tonumber(row.ra))
		ra = ra * (2 * math.pi / 24)
		
		-- in degrees
		local dec = assert(tonumber(row.dec))
		dec = math.rad(dec)
		
		-- docs say parsecs
		local dist = assert(tonumber(row.dist))
		
		-- reconstruct spherical to cartesian and see accurate the xyz is
		local rx = dist * math.cos(ra) * math.cos(dec)
		local ry = dist * math.sin(ra) * math.cos(dec)
		local rz = dist * math.sin(dec)

		-- dist error in xyz <-> dist ra dec of the original data
		local posSphereError = math.sqrt(
			(x - rx) * (x - rx)
			+ (y - ry) * (y - ry)
			+ (z - rz) * (z - rz)
		)

		if offsetFromSun then
			x = x - sunPos.x
			y = y - sunPos.y
			z = z - sunPos.z
			rx = rx - sunPos.x
			ry = ry - sunPos.y
			rz = rz - sunPos.z
		end


		-- HYG is in equatorial coordinates
		-- rotate equatorial pos to ecliptic pos
		local postEclipticDistError = math.nan
		if rotateEquatorialToEcliptic then
			local epsilon = math.rad(23 + 1/60*(26 + 1/60*(21.4119)))	-- earth tilt
			local cosEps = math.cos(epsilon)
			local sinEps = math.sin(epsilon)
			y, z = cosEps * y + sinEps * z,
					-sinEps * y + cosEps * z
			ry, rz = cosEps * ry + sinEps * rz,
					-sinEps * ry + cosEps * rz
			
			postEclipticDistError = math.sqrt((x-rx)^2 + (y-ry)^2 + (z-rz)^2)
		end


--[[ recalc ra, dec
ra = math.atan2(y, x)
dec = math.atan2(z, math.sqrt(x*x + y*y))
dist = math.sqrt(x*x + y*y + z*z)
--]]

local ra2 = (ra + math.pi) % (2 * math.pi) - math.pi

		local mag = assert(tonumber(row.mag))		-- apparent magnitude
		local absmag = assert(tonumber(row.absmag))	-- absolute magnitude
		local lum = assert(tonumber(row.lum))		-- luminosity
		local colorIndex = tonumber(row.ci) or .656	-- fill in blanks with something close to white

		-- calc temp based on colorIndex BV using Newton method on the color temp BV function
		local temp = 4600 * (1 / (.92 * colorIndex + 1.7) + 1 / (.92 * colorIndex + .62))

		
		statset:accum(
			absmag,
			mag,
			lum,
			colorIndex,
			temp,
			ra,
			ra2,
			dec,
			dist,
			posSphereError,
			postEclipticDistError 
		)


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
				-- TODO use StatSet ?
				-- in fact, then just combine StatSet's at the end?
				ra = {},
				ra2 = {},
				dec = {},
				dist = {},
				mag = {},
				absmag = {},
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
		
		conInfo.absmag.min = math.min(conInfo.absmag.min or absmag, absmag)
		conInfo.absmag.max = math.max(conInfo.absmag.max or absmag, absmag)


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
	
		if outputLum then
			-- gaia is storing lum and temp
			buffer[6 + numElem * numValidRows] = lum
			buffer[7 + numElem * numValidRows] = colorIndex 
		else
			-- solarsystem is storing absmag and colorIndex (or should it be absmag?)
			buffer[6 + numElem * numValidRows] = absmag
			buffer[7 + numElem * numValidRows] = colorIndex 
		end
		buffer[8 + numElem * numValidRows] = constellationIndex-1
		
		maxAbsPos = math.max(maxAbsPos, x)
		maxAbsPos = math.max(maxAbsPos, y)
		maxAbsPos = math.max(maxAbsPos, z)

		numValidRows = numValidRows + 1
	end
end

-- now sort the constellation stars by their magnitude
for _,con in ipairs(constellations) do
	table.sort(con.indexes, function(a,b)
		-- sort by buffer abs mag / luminosity
		return buffer[6 + numElem * a] < buffer[6 + numElem * b]
		-- sort by original row magnitude
		--return tonumber(csv.rows[a[1]].absmag) < tonumber(csv.rows[b[1]].absmag)
	end)
	--[=[ if you want to see those values:
	con.magnitudes = table.mapi(con.indexes, function(i)
		-- buffer magnitude
		return buffer[6 + numElem * i]
		-- original row magnitude
		--return tonumber(csv.rows[i[1]].absmag)
	end)
	--]=]
end

print('rows processed:', numRows)
print('valid entries:', numValidRows)
print('invalid entries:', numInvalidRows)

print('max abs position coordinate:', maxAbsPos)

--[[ see what data is present:
print('num columns provided:')
for _,name in ipairs(csv.columns) do
	print('', name, columnCounts[name])
end
--]]

print(statset)

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

--[[ in json
local json = require 'dkjson'
file['namedStars.json'] = 'namedStars = ' .. json.encode(namedStars, {indent=true}) ..';'
file['constellations.json'] = 'constellations = '..json.encode(constellations, {indent=true}) .. ';'
--]]
-- [[ in lua
file['namedStars.lua'] = 'namedStars = ' .. tolua(namedStars)
file['constellations.lua'] = 'constellations = '..tolua(constellations)
--]]
print'done!'

--[[ plot dist density map
local f = io.open('dist-distribution.txt', 'w')
f:write'#dist_min\tdist_max\tcount\n'
local distbin = Bin(0, statset.dist.max, 2000)
for i,row in ipairs(csv.rows) do
	-- [=[ bin by dist
	local dist = assert(tonumber(row.dist))
	--]=]
	--[=[ bin by xyz dist
	local x = assert(tonumber(row.x)) - sunPos.x
	local y = assert(tonumber(row.y)) - sunPos.y
	local z = assert(tonumber(row.z)) - sunPos.z
	local dist = math.sqrt(x*x + y*y + z*z)
	--]=]
	distbin:accum(dist) 
end
for bin,count in ipairs(distbin) do
	local distMin, distMax = distbin:getBinBounds(bin)
	f:write(distMin,'\t',distMax,'\t',count,'\n')
end
f:close()
--]]
