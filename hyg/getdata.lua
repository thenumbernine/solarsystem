#! /usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
require 'ffi.c.stdio'	-- FILE and all its related functions

local lines = io.readfile('hygxyz.csv'):trim():split('\n')
local header = lines:remove(1)
local columnNames = header:split(',')	-- columnNames[1] = 'StarID'
local columnForNames = columnNames:map(function(name,i) return i,name end)	-- columnForNames.StarID = 1

local numElem = 4
local buffer = ffi.new('float[?]', numElem * #lines)	-- allocate space for xyz of each planet 

local maxAbs = 0
for i,line in ipairs(lines) do
	local parts = line:split(',')
	if #parts ~= #columnNames then error("row "..(i+1).." does not have enoguh columns\n"..line) end 
	local x = assert(tonumber(parts[columnForNames.X]))
	local y = assert(tonumber(parts[columnForNames.Y]))
	local z = assert(tonumber(parts[columnForNames.Z]))
	local mag = assert(tonumber(parts[columnForNames.AbsMag]))	--apparent magnitude
	buffer[0 + numElem * (i - 1)] = x
	buffer[1 + numElem * (i - 1)] = y
	buffer[2 + numElem * (i - 1)] = z
	buffer[3 + numElem * (i - 1)] = mag
	maxAbs = math.max(maxAbs, x)
	maxAbs = math.max(maxAbs, y)
	maxAbs = math.max(maxAbs, z)
end
print('abs max coordinate',maxAbs)

local filename = 'stardata.f32'
local fh = ffi.C.fopen(filename, 'wb') or error("failed to open file "..filename.." for writing")
if ffi.C.fwrite(buffer, ffi.sizeof('float') * numElem * #lines, 1, fh) ~= 1 then error("failed to write data to "..filename) end
ffi.C.fclose(fh)

