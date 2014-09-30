#! /usr/bin/env lua -lluarocks.require
require 'ext'
local json = require 'dkjson'
local csv = require 'csv'

local verbose = false
if arg[1] == 'v' then
	verbose = true
end

local data = csv.file('planets.csv')
local columnNames = table.remove(data.rows, 1)
setmetatable(columnNames, nil)
data:setColumnNames(columnNames)

local starNames = table()
for _,row in ipairs(data.rows) do
	starNames[row.pl_hostname] = true
end
starNames = starNames:keys():sort()
print('stars:')
print(starNames:concat('\n'))

