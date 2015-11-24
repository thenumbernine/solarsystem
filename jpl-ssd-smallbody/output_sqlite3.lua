local class = require 'ext.class'
local OutputMethod = require 'output'

local OutputToSQLite3 = class(OutputMethod)

--[[
I'm going to define columns up front, then sql merges will be easier
--]]
OutputToSQLite3.columns = require 'coldesc'

OutputToSQLite3.tableName = 'data'

function OutputToSQLite3:staticInit()
	local luasql = require 'luasql.sqlite3'
	self.env = assert(luasql.sqlite3())
	self.conn = assert(self.env:connect('database.sqlite3'))
	self.conn:execute('drop table '..self.tableName)	--might fail if the table is new

	local columnDescs = table()
	for _,column in ipairs(self.columns) do
		columnDescs:insert(column.name .. ' ' .. column.type)
	end
	local cmd = 'create table '..self.tableName..' ('..columnDescs:concat(', ')..')'
	-- don't create the table until the first body
	assert(self.conn:execute(cmd))
end


--[[
args:
	bodyType = 0 = comet, 1 = numbered asteroid, 2 = unnumbered asteroid
--]]
function OutputToSQLite3:init(args)
	self.bodyType = assert(args.bodyType)
end

function OutputToSQLite3:processBody(body)
	body.bodyType = self.bodyType
	local valueStrs = table()
	for _,column in ipairs(self.columns) do
		local value = body[column.name]
		if column.type == 'text' then
			if value then value = ('%q'):format(tostring(value)) end 	--preserve nils
		elseif column.type == 'number' then
			value = tonumber(value)
		end
		if value == nil then value = 'null' end
		valueStrs:insert(value)
	end
	local cmd = 'insert into '..self.tableName..' values ('..valueStrs:concat(',')..')'
	assert(self.conn:execute(cmd))
end

function OutputToSQLite3:staticDone()
	self.conn:close()
	self.env:close()
end

return OutputToSQLite3
