-- http://www.astro-phys.com/api

require 'ext'
local json = require 'json'
local entries = assert(json.decode(assert(io.readfile('earthquakes.json'))))

local planetnames = {
	'sun',
	'mercury',
	'venus',
	'earth',
	'moon',
	'mars',
	'jupiter',
	'saturn',
	'uranus',
	'neptune',
	'pluto',
}
local planetnamestr = table.concat(planetnames, ',')

local planetsForDates = {}

for _,entry in ipairs(entries) do
		
	if entry.year and entry.year > 800 then
		local year = entry.year
		local month = entry.month or 1
		local day = entry.day or 1
		local hour = entry.hour or 0
		local min = entry.min or 0
		local sec = entry.sec or 0
		
		-- some of the month/day's are 0 (for 'unknown')
		if month == 0 then month = 1 end
		if day == 0 then day = 1 end
		
		print('year',year)
		print('month',month)
		print('day',day)
		print('hour',hour)
		print('min',min)
		print('sec',sec)
		local datestr = ('%d-%d-%d %.02d:%.02d:%.02d'):format(year, month, day, hour, min, sec)
		print('datestr',datestr)
		
		
		-- [[
		local http = require 'socket.http'
		local url = require 'socket.url'
		d = assert(http.request('http://www.astro-phys.com/api/de406/states?date='..url.escape(datestr)..'&bodies='..planetnamestr))
		print(d)
		d = assert(json.decode(d))
		d.calendarDate = datestr
		table.insert(planetsForDates, d)
		--]]
	end
end

io.writefile('planets.json', json.encode(planetsForDates))