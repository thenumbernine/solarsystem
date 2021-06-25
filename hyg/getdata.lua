#! /usr/bin/env luajit

-- should i move the data processing stuff to universe/offline?
-- should I put the whole universe/offline folder in a repo of its own for point processing?

require 'ext'
local ffi = require 'ffi'
local CSV = require 'csv'
local Stat = require 'stat'
local StatSet = require 'stat.set'

local csv = CSV.file'hygdata_v3.csv'

local outputLum = ... == 'lum'

local offsetFromSun = true
local rotateEquatorialToEcliptic = true

local compareAbsMagToLum = true
local compareAppMagToAbsMag = true
local compareTempToLuminosity = true

local rows = csv.rows

local rowHeaders = rows:remove(1)
for i=1,#rowHeaders do rowHeaders[i] = rowHeaders[i]:lower() end
csv:setColumnNames(rowHeaders)

--[[ debugging
for i,row in ipairs(rows) do
	if row.proper ~= '' then
		rows[row.proper] = row
	end
end
local function rowtotable(row)
	local t = {}
	for i,h in ipairs(rowHeaders) do
		t[h] = row[i]
	end
	return t
end
print(tolua(rowtotable(rows.Sol)))
print(tolua(rowtotable(rows.Polaris)))
print(tolua(rowtotable(rows.Merak)))
do return end
--]]

-- sun is positioned relative to the earth, so offset everything for the sun to be the center
local sunRow = rows[1]
local sunPos = {
	x = assert(tonumber(sunRow.x)),
	y = assert(tonumber(sunRow.y)),
	z = assert(tonumber(sunRow.z)),
}

local numRows = #rows
local numElem = 9
local bufferType = 'float'
local buffer = ffi.new(bufferType..'[?]', numElem * numRows)	-- allocate space for data of each planet 

local maxAbsPos = 0
local namedStars = {}
local constellations = table()	-- list of unique names so far
local columnCounts = {}


local statset = StatSet(
	'ra',
	'dec',
	'dist',
	'x', 'y', 'z',
	'ra2',
	
	'absmag',	-- absolute magnitude
	'mag',		-- apparent magnitude
	'lum',		-- luminosity, in LSun (even though the sun isn't 1.0 ... )
	'colorIndex',
	'temp',
	
	'posSphereError',
	'postEclipticDistError',
	'absMagToLumError',
	'appMagToAbsMagError'
)

local numValidRows = 0
local numInvalidRows = 0
for i,row in ipairs(rows) do
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
		-- https://en.wikipedia.org/wiki/Color_index 
		local temp = 4600 * (1 / (.92 * colorIndex + 1.7) + 1 / (.92 * colorIndex + .62))

		local absMagToLumError = 0
		if compareAbsMagToLum then
			local LStarOverLSun = lum
			
			-- https://www.omnicalculator.com/physics/luminosity 
			local LSun = 3.828e+26 	-- Watts
			local L0 = 3.0128e+28	-- Watts
			local LSunOverL0 = LSun / L0
			
			local M = (-2.5 / math.log(10)) * math.log(LStarOverLSun * LSunOverL0)
			-- max abs error is 3.5527136788005e-15
			absMagToLumError = M - absmag
		end

		local appMagToAbsMagError = 0
		if compareAppMagToAbsMag then
			local m = absmag - 5 + (5 / math.log(10)) * math.log(dist)
			-- I guess for the sun, which is at the center of our coordinate system, log(0) for the sun means inf apparent magnitude
			-- maybe that's why the earth is the original coordinate system origin?
			appMagToAbsMagError = m - mag
		end

		--[[ compare temp to luminosity -- this requires the star's radius to be known though ... hmm
		if compareTempToLuminosity then
			-- https://en.wikipedia.org/wiki/Luminosity#Luminosity_formulae
			local lumFromTemp = sigma * 4 * math.pi * radius * radius * temp * temp * temp * temp
			local dT = lum - lumFromTemp 
		end
		--]]
		
		statset:accum(
			ra,
			dec,
			dist,
			x,y,z,
			ra2,		
			
			absmag,
			mag,
			lum,
			colorIndex,
			temp,
			
			posSphereError,
			postEclipticDistError,
			absMagToLumError,
			appMagToAbsMagError 
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
				-- TODO don't use global statSet, only use constellation StatSets, and combine in the end?
				stats = StatSet(
					'ra',
					'ra2',
					'dec',
					'dist',
					'lum',
					'mag',
					'absmag'
				),
				-- array of all indexes in the HYG data of all stars of the constellation
				indexes = table(),
			}
			constellations:insert(conInfo)
			constellationIndex = #constellations
		end

		conInfo.stats:accum(
			ra,
			ra2,
			dec,
			dist,
			lum,
			mag,
			absmag
		)

		-- without individual star info, constellsions.json is 22603 bytes
		-- with individual star info it is 678818 bytes
		conInfo.indexes:insert(numValidRows) 
		-- verify original row matches
		--conInfo.indexes:insert{i, row.id, numValidRows}

--print(i, row.id, numValidRows, mag, row.proper)
		buffer[0 + numElem * numValidRows] = x
		buffer[1 + numElem * numValidRows] = y
		buffer[2 + numElem * numValidRows] = z
		buffer[3 + numElem * numValidRows] = vx
		buffer[4 + numElem * numValidRows] = vy
		buffer[5 + numElem * numValidRows] = vz
	
		-- gaia is storing lum and temp, and I think storing lum and calculating absmag can reduce us one log() call
		buffer[6 + numElem * numValidRows] = lum
		buffer[7 + numElem * numValidRows] = temp
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
		--return tonumber(rows[a[1]].absmag) < tonumber(rows[b[1]].absmag)
	end)
	--[=[ if you want to see those values:
	con.magnitudes = table.mapi(con.indexes, function(i)
		-- buffer magnitude
		return buffer[6 + numElem * i]
		-- original row magnitude
		--return tonumber(rows[i[1]].absmag)
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
	local dra = con.stats.ra.max - con.stats.ra.min
	local dra2 = con.stats.ra2.max - con.stats.ra2.min
	if dra2 < dra then
		con.stats.ra = con.stats.ra2
		-- don't forget the indexed reference to ra is still dangling
	end
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
file['constellations.lua'] = 'constellations = '..tolua(
	constellations:mapi(function(con) 
		local o = {
			name = con.name,
			indexes = table(con.indexes):setmetatable(nil),
		}
		for k,stat in pairs(con.stats) do
			if type(k) ~= 'number' then
				o[k] = {
					min = stat.min,
					max = stat.max,
				}
			end
		end
		return o
	end)
)
--]]
print'done!'

--[[ plot dist density map
local f = io.open('dist-distribution.txt', 'w')
f:write'#dist_min\tdist_max\tcount\n'
local distbin = Bin(0, statset.dist.max, 2000)
for i,row in ipairs(rows) do
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
