#! /usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
require 'ffi.c.stdio'	-- FILE and all its related functions

local lines = io.readfile('hygxyz.csv'):trim():split('\n')
local header = lines:remove(1)
local columnNames = header:split(',')	-- columnNames[1] = 'StarID'
local columnForNames = columnNames:map(function(name,i) return i,name end)	-- columnForNames.StarID = 1

local buffer = ffi.new('float[?]', 3*#lines)	-- allocate space for xyz of each planet 

local maxAbs = 0
for i,line in ipairs(lines) do
	local parts = line:split(',')
	if #parts ~= #columnNames then error("row "..(i+1).." does not have enoguh columns\n"..line) end 
	local x = assert(tonumber(parts[columnForNames.X]))
	local y = assert(tonumber(parts[columnForNames.Y]))
	local z = assert(tonumber(parts[columnForNames.Z]))
	buffer[0+3*(i-1)] = x
	buffer[1+3*(i-1)] = y
	buffer[2+3*(i-1)] = z
	maxAbs = math.max(maxAbs,x)
	maxAbs = math.max(maxAbs,y)
	maxAbs = math.max(maxAbs,z)
end
print('abs max coordinate',maxAbs)

local filename = 'starpositions.f32'
local fh = ffi.C.fopen(filename, 'wb') or error("failed to open file "..filename.." for writing")
if ffi.C.fwrite(buffer, ffi.sizeof('float') * 3 * #lines, 1, fh) ~= 1 then error("failed to write data to "..filename) end
ffi.C.fclose(fh)

