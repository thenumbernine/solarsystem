#! /usr/bin/env lua -lluarocks.require
require 'ext'
local json = require 'dkjson'

local verbose = false
if arg[1] == 'v' then
	verbose = true
end

local lines = io.readfile('exoplanet.eu_catalog.csv'):trim():split('\n')
local header = lines:remove(1):match('#(.*)')
local columnNames = header:split(','):map(string.trim)
local columnForNames = columnNames:map(function(name,i) return i,name end)	-- columnForNames.StarID = 1

local maxAbs = 0
local namedStars = {}
local columnCounts = {}
for i=1,#lines do
	local line = lines[i]
	local parts = line:split(',')		-- TODO a real parse, so quotes can contain commas
	if #parts ~= #columnNames then print("line "..(i+1).." has "..#parts.." columns when we expected "..#columnNames) end 
	local data = columnNames:map(function(part, i) 
		if i > #columnNames then return end
		local columnName = columnNames[i]
		if #part > 0 then columnCounts[columnName] = (columnCounts[columnName] or 0) + 1 end
		return part, columnName
	end)
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

-- make name of sun consistent with NASA data
assert(namedStars[1].name == 'Sol')
namedStars[1].name = 'Sun'

io.writefile('namedStars.json', 'namedStars = \n' .. json.encode(table, {indent=true}) .. ';')

