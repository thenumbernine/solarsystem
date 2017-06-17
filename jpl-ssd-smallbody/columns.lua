require 'ext'

--[[
holds list of {start, finish, name} for each column
--]]
local Columns = class()

function Columns:init(lines)
	local line1, line2 = lines[1], lines[2]
	self.columns = table()
	local current = 1
	while true do
		local start, finish = line2:find('%-+', current)
		if not start then break end
		assert(finish)
		self.columns:insert{
			start = start,
			finish = finish,
			name = line1:sub(start,finish):trim(),
		}
		current = finish + 1
	end
	self.columnsByName = self.columns:map(function(column)
		return column, column.name
	end)
end

local ColumnAccess = class()

function ColumnAccess:init(columns, line)
	self.columns = assert(columns)
	self.line = assert(line)
end
	
function ColumnAccess:__index(name)
	local m = getmetatable(self)
	if m[name] then return m[name] end
	local columns = assert(rawget(self, 'columns'))
	local line = assert(rawget(self, 'line'))
	local col = columns.columnsByName[name] or error("failed to find column "..tostring(name))
	local start, finish = col.start, col.finish
	if col == columns.columns:last() then 	-- allow last column to read to end of line
		return line:sub(start)
	else
		return line:sub(start, finish)
	end
end

function Columns:__call(line)
	return ColumnAccess(self, line)
end

return Columns
