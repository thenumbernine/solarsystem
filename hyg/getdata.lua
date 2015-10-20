#! /usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
require 'ffi.c.stdio'	-- FILE and all its related functions

local verbose = false
if arg[1] == 'v' then
	verbose = true
end

local csv = require 'csv'.file 'hygdata_v3.csv'
csv:setColumnNames(csv.rows:remove(1))

local numRows = #csv.rows
local numElem = 8
local buffer = ffi.new('float[?]', numElem * numRows)	-- allocate space for data of each planet 

local maxAbs = 0
local namedStars = {}
local columnCounts = {}

for i,row in ipairs(csv.rows) do
	for _,col in ipairs(csv.columns) do
		if #row[col] > 0 then
			columnCounts[col] = (columnCounts[col] or 0) + 1
		end
	end
	local zeroBasedIndex = i-1
	
	--[[ this assertion isn't required, I was just curious if it was enforced
	-- and as of the v3 release it's no longer true
	if zeroBasedIndex ~= tonumber(row.id) then
		error("zero-based index "..zeroBasedIndex.." not equal to row id "..row.id)
	end
	--]]

	-- in parsecs
	local x = assert(tonumber(row.x))
	local y = assert(tonumber(row.y))
	local z = assert(tonumber(row.z))

	-- HYG is in equatorial coordinates
	-- rotate equatorial pos to ecliptic pos
	local epsilon = math.rad(23 + 1/60*(26 + 1/60*(21.4119)))	-- earth tilt
	local cosEps = math.cos(epsilon)
	local sinEps = math.sin(epsilon)
	y, z = cosEps * y + sinEps * z, -sinEps * y + cosEps * z
	
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
	if verbose then
		print(tolua(row))
		if io.read(1) == 'q' then break end
	end
end
print('abs max coordinate', maxAbs)
print('num columns provided:')
for _,name in ipairs(csv.columns) do
	print('', name, columnCounts[name])
end

-- make name of sun consistent with NASA data
assert(namedStars[1].name == 'Sol')
namedStars[1].name = 'Sun'

-- write 
local filename = 'stardata.f32'
local fh = ffi.C.fopen(filename, 'wb') or error("failed to open file "..filename.." for writing")
if ffi.C.fwrite(buffer, ffi.sizeof('float') * numElem * numRows, 1, fh) ~= 1 then error("failed to write data to "..filename) end
ffi.C.fclose(fh)

file['namedStars.json'] = 'namedStars = '
	.. require 'dkjson'.encode(namedStars, {indent=true}) ..';'
--[=[
..'[\n\t' ..
-- poor man's json encode 
	table(namedStars):map(function(namedStar)
		local kvs = table()
		for k,v in pairs(namedStar) do
			if type(v) == 'string' then
				v = ('%q'):format(v)
			else
				v = tostring(v)
			end
			kvs:insert(('%q'):format(k) .. ':' .. v)
		end
		return '{' .. kvs:concat(', ') .. '}'
	end):concat(',\n\t') .. '\n];'
--]=]

