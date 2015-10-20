#! /usr/bin/env lua
require 'ext'
local socket = require 'socket'
local json = require 'dkjson'
local julian = assert(loadfile('../lua/julian.lua'))()

-- try to match the astro-phys date on file

print('connecting...')
local conn = assert(socket.connect('horizons.jpl.nasa.gov', 6775))
conn:settimeout(0, 't')
print('connected!')

local f = io.open('horizons.txt', 'w')

local function readUntil(find)
	local data = ''
	
	-- gather all input until there's none left
	while true do
		if #socket.select({conn}, nil, 0) == 1 then
			local result, reason, partial = conn:receive(1)
			if reason == 'closed' then error 'closed' end
			result = result or partial
			if result then
				data = data .. result
				io.write(result)
				io.flush()
				f:write(result)
				f:flush()
				if data:find(find,nil,true) then break end
			end
		end
	end
	return data 
end

local function send(cmd)
	f:write('[[[my input:'..cmd..']]]\n')
	f:flush()
	conn:send(cmd..'\n')
	readUntil(cmd)
end

local function getEphemerisData(data)
	local lines = data:trim():split('\r\n')
	local i = assert(lines:find('$$SOE'))
	i = i + 1
	-- lines[i] has the current date
	i = i + 1
	local pos = lines[i]:trim():split('%s+'):map(function(w) return w:lower() end)	-- leave them as strings for now
	i = i + 1
	local vel = lines[i]:trim():split('%s+'):map(function(w) return w:lower() end)	-- no point in losing precision deserializing and reserializing data
	-- next line is the light time, range, and range-rate
	-- next line should be $$EOE
	return {pos=pos, vel=vel}
end

local function getMajorBodies(data)
	local lines = data:split('\r\n')
	local i = assert(lines:find(nil, function(line)
		local ws = line:trim():split('%s+')
		return ws[1] == 'ID#' and ws[2] == 'Name' and ws[3] == 'Designation'
	end))
	i = i + 1
	-- separator
	i = i + 1
	local bodies = table()
	while true do
		local line = lines[i]
		if line:match('^%s*$') then break end
		local id = tonumber(line:sub(3,9)) 
		local name = line:sub(12,45)
		local designation = line:sub(47,57)
		local alias = line:sub(60, 78)
		if id > 0 then
			bodies:insert{id=id, name=name:trim()}
		end
		i = i + 1
	end
	return bodies
end

local time = os.time()
local date = os.date('!*t', time)
local julianDate = julian.fromCalendar(date)
local startDate = os.date('%Y-%b-%d %H:%M', time)	--'2014-Sep-18 22:07'
local endDate = os.date('%Y-%b-%d %H:%M', time+61)	--'2014-Sep-18 22:08'

local entries = {}
local function run()

	socket.sleep(3)	-- stop two seconds for terminal negotiation ... or implement it
	readUntil('Horizons>')
	
	send('page')		-- get rid of interactive output
	readUntil('Horizons>')

	-- [[
	send('mb')	-- list major bodies
	local majorBodies = getMajorBodies(readUntil(': '))		-- list major bodies
	send('')			-- end final prompt
	readUntil('Horizons>')
	--]]
	--[[
	local majorBodies = {'10', '1'}	-- debug: sun and mercury
	--]]

-- NOTICE planet Hartley 2 has no orbit data past 2012-MAR-16 

	local first = true
	for _,body in ipairs(majorBodies) do
		-- first option
		send(body.id)		-- get planet data
		readUntil('<cr>:')
		send('e')		-- get its ephemeris data
		readUntil(' :')
		send('v')		-- read in vector format 
		readUntil(' :')
		send(first and '@0' or '')		-- centered at object 0 (solar barycentric coordinate center) 
		readUntil(' :')
		send('eclip')	-- ecliptical coordinates 
		readUntil(' :')
		send(startDate)	-- starting date (probably GMT)
		readUntil(' :')
		send(endDate)	-- ending date (can't be equal to starting date)
		readUntil(' :')
		send('')		-- keep output interval.  if it's bigger than the end date is from the start then we just get the starting coordinates
		readUntil(' :')
		--output table defaults - don't accept the first time through
		local ephemerisData
		if not first then
			send('')
		else
			send('n')		-- no, don't accept defaults
			readUntil(' :')
			send('')		-- keep J2000 coordinates
			readUntil(' :')
			send('')		-- keep no corrections
			readUntil(' :')
			send('3')		-- switch to km/day (from default of AU's)
			readUntil(' :')
			send('')		-- keep no CSV output format 
			readUntil(' :')
			send('')		-- keep no cartesian output labels
			readUntil(' :')
			send('')		-- keep output table type
		end
		ephemerisData = getEphemerisData(readUntil(', ? : '))		-- keep output table type to state vector (pos+vel) + extra stuff
		local ephemerisData = table(ephemerisData, body)
		setmetatable(ephemerisData, nil)
		table.insert(entries, ephemerisData)
		--io.writefile('horizons-results.json', json.encode(entries))
		-- done with output table stuff
		send('n')	-- exit sun ephemeris query
		readUntil('Horizons>')
		first = false
	end

	send('q')	-- quit
	readUntil('forever')
end
xpcall(run, function(err)
	if err == 'closed' then return end
	io.stderr:write(err..'\n'..debug.traceback()..'\n')
end)

f:close()
conn:close()

entries = {
	coords = entries,
	julianDate = julianDate,
}
io.writefile('dynamic-vars.json', 'horizonsDynamicData = '..json.encode(entries, {indent=true})..';')

