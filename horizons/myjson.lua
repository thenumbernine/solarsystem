local table = require 'ext.table'
local tolua = require 'ext.tolua'
local dkjson = require 'dkjson'

local myjson = table(dkjson)

local function adjustMetatables(o)
	local m = getmetatable(o)
	if m ~= nil then
		if m == table then
			setmetatable(o, nil)
		else
			print("got an unknown metatable for " .. tolua(o):sub(1, 100))
			print("metatable: "..tolua(m):sub(1, 100))
			error'here'
		end
	end
	local keyorder = table.keys(o):mapi(tostring):sort()
	setmetatable(o, {
		__jsonorder = setmetatable(keyorder, nil),
	})
	for k,v in pairs(o) do
		if type(v) == 'table' then 
			adjustMetatables(v) 
		end
	end
end

-- I would like keys sorted, so I can diff files and quickly see differences
-- dkjson is great for everything except for doing this.  it only provides {keyorder=...} or the __jsonorder=... metafield.
function myjson.encode(...)
	local object = ...
	if type(object) == 'table' then
		adjustMetatables(object)
		return dkjson.encode(object, select(2, ...))
	else
		return dkjson.encode(...)
	end
end

local function removeMetatables(o)
	setmetatable(o, nil)
	for k,v in pairs(o) do
		if type(v) == 'table' then
			removeMetatables(v)
		end
	end
end

function myjson.decode(...)
	local result = table.pack(dkjson.decode(...))
	if result[1] ~= nil then
		if type(result[1]) == 'table' then
			removeMetatables(result[1])
		end
	end
	return table.unpack(result, 1, result.n)
end

return myjson
