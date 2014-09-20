--[[
remote request
--]]

module(..., package.seeall)

require 'ext'
require 'wsapi.request'

--[[
GET parameters:
	comets	= whether to search through comets
	numbered = whether to search through numbered asteroids
	unnumbered = whether to search through unnumbered asteroids
	text = what string to search for
	page = what page of results to look at
--]]
function run(env)
	local headers = {["Content-type"] = "text/javascript"}
	local req = wsapi.request.new(env)
	local get = req.GET

	local isComets, isNumbered, isUnnumbered, searchText, searchPage
	if not get then
		isComets = true
		isNumbered = true
		isUnnumbered = true
		searchText = ''
		searchPage = 1
	else
		isComets = get.comets
		isNumbered = get.numbered
		isUnnumbered = get.unnumbered
		searchText = get.text
		searchPage = math.max(1, tonumber(get.page) or 1)
	end

	if #searchText == 0 then
		searchText = '%'
	else
		searchText = '%'..searchText..'%'
	end

	local pageSize = 20
	local offset = (searchPage - 1) * pageSize

	local luasql = require 'luasql.sqlite3'
	local env, conn, cur
	xpcall(function()
		env = luasql.sqlite3()
		conn = assert(env:connect('database.sqlite3'))
		
		local bodyTypeCond = table()
		if isComets then bodyTypeCond:insert('bodyType == 0') end
		if isNumbered then bodyTypeCond:insert('bodyType == 1') end
		if isUnnumbered then bodyTypeCond:insert('bodyType == 2') end
	
		local cmd = 'select * from data where name like '..('%q'):format(searchText)..' and ('..bodyTypeCond:concat(' or ')..') limit '..offset..', '..pageSize
		cur = assert(con:execute(cmd))

		local results = {}

		-- TODO only collect what we need? convert to Lua types?
		local row = cur:fetch({}, 'a')
		while row do
			-- collect results
			local body = {}
			for k,v in pairs(row) do
				body[k] = row[k]
			end
			table.insert(results, body)
			
			row = cur:fetch(row, 'a')
		end

		local json = require 'dkjson'
		return 200, headers, function()
			coroutine.yield(json.encode(results)
		end
	end, function(err)
		if cur then cur:close() end
		if conn then conn:close() end
		if env then env:close() end
		error(err..'\n'..debug.traceback())
	end)
end

