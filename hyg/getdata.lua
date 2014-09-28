#! /usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
require 'ffi.c.stdio'	-- FILE and all its related functions

local verbose = false
if arg[1] == 'v' then
	verbose = true
end

local lines = io.readfile('hygxyz.csv'):trim():split('\n')
local header = lines:remove(1)
local columnNames = header:split(',')	-- columnNames[1] = 'StarID'
local columnForNames = columnNames:map(function(name,i) return i,name end)	-- columnForNames.StarID = 1

local numElem = 5
local buffer = ffi.new('float[?]', numElem * #lines)	-- allocate space for xyz of each planet 

local maxAbs = 0
local namedStars = {}
local columnCounts = {}
assert(lines[1]:split(',')[columnForNames.ProperName] == 'Sol')	-- assert that the first star is the sun, named "Sol" (and at xyz = 000)
for i=2,#lines do
	local line = lines[i]
	local parts = line:split(',')
	if #parts ~= #columnNames then error("line "..(i+2).." does not have enoguh columns\n"..line) end 
	local data = parts:map(function(part, i) 
		local columnName = columnNames[i]
		if #part > 0 then columnCounts[columnName] = (columnCounts[columnName] or 0) + 1 end
		return part, columnName
	end)
	local x = assert(tonumber(data.X))
	local y = assert(tonumber(data.Y))
	local z = assert(tonumber(data.Z))
	local mag = assert(tonumber(data.AbsMag))	--apparent magnitude
	local colorIndex = tonumber(data.ColorIndex) or .656	-- fill in blanks with something close to white
	local zeroBasedIndex = i-2
	buffer[0 + numElem * zeroBasedIndex] = x
	buffer[1 + numElem * zeroBasedIndex] = y
	buffer[2 + numElem * zeroBasedIndex] = z
	buffer[3 + numElem * zeroBasedIndex] = mag
	buffer[4 + numElem * zeroBasedIndex] = colorIndex 
	maxAbs = math.max(maxAbs, x)
	maxAbs = math.max(maxAbs, y)
	maxAbs = math.max(maxAbs, z)
	assert(zeroBasedIndex == tonumber(data.StarID)-1)	-- this assertion isn't required, I was just curious if it was enforced
	if #data.ProperName > 0 then
		table.insert(namedStars, {name=data.ProperName, index=zeroBasedIndex})
	end
	if verbose then
		print(toLua(data, ' ', ''))
		if io.read(1) == 'q' then break end
	end
end
print('abs max coordinate', maxAbs)
print('num columns provided:')
for _,name in ipairs(columnNames) do
	print('',name,columnCounts[name])
end

-- write 
local filename = 'stardata.f32'
local fh = ffi.C.fopen(filename, 'wb') or error("failed to open file "..filename.." for writing")
if ffi.C.fwrite(buffer, ffi.sizeof('float') * numElem * #lines, 1, fh) ~= 1 then error("failed to write data to "..filename) end
ffi.C.fclose(fh)

io.writefile('namedStars.json', 'namedStars = [\n\t' ..
-- poor man's json.encode ... I need luajit for the ffi, but my luarocks aren't compiled against luajit, so I can't use dkjson ...
	table(namedStars):map(function(namedStar)
		local kvs = table()
		for k,v in pairs(namedStar) do
			if type(v) == 'string' then
				v = ('%q'):format(v)
			else
				v = tostring(v)
			end
			kvs:insert('[' .. ('%q'):format(k) .. '] : ' .. v)
		end
		return '{' .. kvs:concat(', ') .. '}'
	end):concat(',\n\t') .. '\n];')
