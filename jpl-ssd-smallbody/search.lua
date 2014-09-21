module(..., package.seeall)
require 'ext'
require 'wsapi.request'

local bodyTypeForEnum = {
	[0] = 'comet',
	[1] = 'numbered asteroid',
	[2] = 'unnumbered asteroid',
}

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

	local isComet = false
	local isNumbered = false
	local isUnnumbered = false
	local searchText = ''
	local searchPage = 1
	if get then
		isComet = get.comet
		isNumbered = get.numbered
		isUnnumbered = get.unnumbered
		searchText = get.text or ''
		searchPage = math.max(1, tonumber(get.page) or 1)
	end

	if #searchText == 0 then
		searchText = '%'
	else
		searchText = '%'..searchText..'%'
	end

	local pageSize = 20
	local offset = (searchPage - 1) * pageSize	-- offset is 0-based

	local function text()
		local env, conn, cur
		if not isComet and not isNumbered and not isUnnumbered then
			coroutine.yield('please specify comet, numbered, or unnumbered')
			return
		end
		local results = select(2, xpcall(function()
			local luasql = require 'luasql.sqlite3'
			env = luasql.sqlite3()
			conn = assert(env:connect('database.sqlite3'))
			
			local bodyTypeCond = table()
			if isComet then bodyTypeCond:insert('bodyType == 0') end
			if isNumbered then bodyTypeCond:insert('bodyType == 1') end
			if isUnnumbered then bodyTypeCond:insert('bodyType == 2') end
	
			local fromStmt = 'data where name like '..('%q'):format(searchText)..' collate nocase and ('..bodyTypeCond:concat(' or ')..')'
			
			local cmd = 'select count() from '..fromStmt
			cur = assert(conn:execute(cmd))

			local row = cur:fetch({}, 'a')
			local count = row['count()']
			cur:close()

			local limitStmt = ' limit '..offset..', '..pageSize
			local cmd = 'select * from '..fromStmt..limitStmt
			cur = assert(conn:execute(cmd))

			local rows = {}

			-- TODO only collect what we need? convert to Lua types?
			local row = cur:fetch({}, 'a')
			while row do
				-- collect results
				local body = {}
				for k,v in pairs(row) do
					if k == 'pk' then	-- don't return PK
					elseif k == 'bodyType' then
						body[k] = bodyTypeForEnum[row[k]]
					else
						body[k] = row[k]
					end
				end
				table.insert(rows, body)
				
				row = cur:fetch(row, 'a')
			end
			cur:close()

			local json = require 'dkjson'
			local results = {rows=rows, count=count}
			return json.encode(results)
		end, function(err)
			if cur then cur:close() end
			if conn then conn:close() end
			if env then env:close() end
			return err..'\n'..debug.traceback()
		end))
		coroutine.yield(results)
	end
	return 200, headers, coroutine.wrap(text)
end
