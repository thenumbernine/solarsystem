module(..., package.seeall)
require 'ext'
require 'wsapi.request'
local json = require 'dkjson'

--[[
GET parameters:
	index = 0-based index of what star to retrieve information on
	TODO search = string to search for a star's name
--]]
function run(env)
	local headers = {['Content-type'] = 'text/javascript'}
	local req = wsapi.request.new(env)
	local get = req.GET

	local index
	if get then
		index = tonumber(get.index)
		if index then
			-- ...maybe I want this as a SQLite database as well?  I haven't checked access times of SQLite vs CSV ...
		end
	end

	local function text()
		if not index then
			coroutine.yield(json.encode{error='invalid index', index=index or 'nil'})
		end
		
		coroutine.yield'TODO FINISHME'
	end	
	return 200, headers, coroutine.wrap(text)
end
