#!/usr/bin/env wsapi.cgi
require 'ext'
local wsapi_request = require 'wsapi.request'
local bit = require 'bit32' or require 'bit'
local body_t = require 'body_t_desc'
local convert = require 'convert'

local bodyTypeForEnum = {
	[0] = 'comet',
	[1] = 'numbered asteroid',
	[2] = 'unnumbered asteroid',
}

local orbitTypeForEnum = {
	[0] = 'elliptic',
	[1] = 'hyperbolic',
	[2] = 'parabolic',
}

-- fixing json's incomplete specs ...
local json_inf = setmetatable({}, { __tojson = function() return ('%q'):format'Infinity' end})
local json_ninf = setmetatable({}, { __tojson = function() return ('%q'):format'-Infinity' end})
local json_nan = setmetatable({}, {__tojson = function() return ('%q'):format'NaN' end})

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
		local startTime = os.clock()	-- TODO hires timer?
		
		local env, conn, cur
		local json = require 'dkjson'
		if not isComet and not isNumbered and not isUnnumbered then
			coroutine.yield(json.encode{rows={}, count=0})
			return
		end
		local results = select(2, xpcall(function()
			local rows = {}
			local f = io.open('alldata.raw', 'r')
			local filesize = f:seek'end'
			local count = filesize / body_t.size
			-- count should be 961405
			for i=0,count-1 do
				local nodepos = i * body_t.size
				f:seek('set', nodepos + body_t.fields.name.offset)
				local name = f:read(44):match'[^\0]*'
				if name:lower():find(searchText,1,true) then
					local row = {}
					
					local function read(fieldName)
						local field = body_t.fields[fieldName]
						f:seek('set', nodepos + field.offset)
						local d = f:read(field.size)
						assert(d)
						if field.type == 'char[44]'
						or field.type == 'char[13]'
						then
							return d:match'[^\0]*'
						elseif field.type == 'double[3]' then
							return {
								convert.double(d:sub(1,8)),
								convert.double(d:sub(9,16)),
								convert.double(d:sub(17,24)),
							}
						else
							-- json doesn't handle non-finite numbers...
							local x = convert[field.type](d)
							if x == math.huge then return json_inf end
							if x == -math.huge then return json_ninf end
							if x ~= x then return json_nan end
							return x
						end
					end
					
					
					local bodyType = read'bodyType'
					if (isComet and bodyType == 0)
					or (isNumbered and bodyType == 1)
					or (isUnnumbered and bodyType == 2)
					then
						table.insert(rows, {
							epoch = read'epoch',
							perihelionDistance = read'perihelionDistance',
							semiMajorAxis = read'semiMajorAxis',
							eccentricity = read'eccentricity',
							inclination = read'inclination',
							argumentOfPeriapsis = read'argumentOfPeriapsis',
							longitudeOfAscendingNode = read'longitudeOfAscendingNode',
							meanAnomalyAtEpoch = read'meanAnomalyAtEpoch',
							absoluteMagnitude = read'absoluteMagnitude',
							magnitudeSlopeParameter = read'magnitudeSlopeParameter',
							timeOfPerihelionPassage = read'timeOfPerihelionPassage',
							bodyType = bodyTypeForEnum[bodyType],
							horizonID = read'horizonID',
							name = name,
							orbitSolutionReference = read'orbitSolutionReference',
							index = read'index',
							pos = read'pos',
							vel = read'vel',
							A = read'A',
							B = read'B',
							eccentricAnomaly = read'eccentricAnomaly',
							timeOfPeriapsisCrossing = read'timeOfPeriapsisCrossing',
							meanAnomaly = read'meanAnomaly',
							orbitType = orbitTypeForEnum[read'orbitType'],
							orbitalPeriod = read'orbitalPeriod',
						})
					end				
				end
			end

			f:close()
			local endTime = os.clock()	-- TODO hires timer?

			local count = #rows
			local results = {rows=rows, count=count, time=endTime - startTime}
			return json.encode(results, {indent=true})
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
