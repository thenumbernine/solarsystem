#!/usr/bin/env wsapi.cgi
require 'ext'
local wsapi_request = require 'wsapi.request'

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
return {
run = function(env)
	local headers = {["Content-type"] = "text/javascript"}
	local req = wsapi_request.new(env)
	local get = req.GET

	local id
	local isComet = false
	local isNumbered = false
	local isUnnumbered = false
	local searchText = ''
	local searchPage = 1
	if get then
		-- either this
		id = get.id
		-- or this
		isComet = get.comet ~= '0'
		isNumbered = get.numbered ~= '0'
		isUnnumbered = get.unnumbered ~= '0'
		searchText = get.text or ''
		searchPage = math.max(1, tonumber(get.page) or 1)
	end

	searchText = searchText:lower()

	local pageSize = 20
	local offset = (searchPage - 1) * pageSize	-- offset is 0-based

	local function text()
		local env, conn, cur
		local json = require 'dkjson'
		if not isComet and not isNumbered and not isUnnumbered then
			coroutine.yield(json.encode{rows={}, count=0})
			return
		end
		local results = select(2, xpcall(function()

			local allNodeIDs = json.decode(file['octree.json']).nodes
			local nodes = {}

			local rows = {}
			for line in io.lines('node-dict.csv') do
				if line:sub(1,1) ~= '#' then
					local name, nodeIDStr, localIndexStr = line:match('^"([^"]*)",([^,]*),([^,]*)$')
					--assert(name and nodeIDStr and localIndexStr)
					local nodeID = tonumber(nodeIDStr)
					--assert(tostring(nodeID) == nodeIDStr)
					--assert(table.find(allNodeIDs, nodeID))
					local localIndex = tonumber(localIndexStr)
					--assert(tostring(localIndex) == localIndexStr)
					if name:lower():find(searchText,1,true) then
						local node = nodes[nodeID]
						if not node then
							node = json.decode(file['nodes/'..nodeID..'.json'])
							nodes[nodeID] = node
						end
						local row = assert(node[localIndex])
						local bodyType = row[17]
						if (isComet and bodyType == 0)
						or (isNumbered and bodyType == 1)
						or (isUnnumbered and bodyType == 2)
						then
							table.insert(rows, {
								pk = row[4],
								bodyType = bodyTypeForEnum[bodyType],
								idNumber = row[19],
								name = name,
								epoch = row[13],
								perihelionDistance = row[14],
								semiMajorAxis = row[5],
								eccentricity = row[9],
								inclination = row[8],
								argumentOfPeriapsis = row[7],
								longitudeOfAscendingNode = row[6],
								meanAnomalyAtEpoch = row[12],
								absoluteMagnitude =row[15],
								magnitudeOfSlopeParameter = row[16],
								timeOfPerihelionPassage = row[10],
								orbitSolutionReference = row[21],
							})
						end				
					end
				end
			end
			
			local count = #rows
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
}

