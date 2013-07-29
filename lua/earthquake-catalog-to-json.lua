-- http://earthquake.usgs.gov/research/data/scr_catalog.php

require 'ext'

local fn = 'catalog.txt'
local ls = assert(io.readfile(fn)):trim():split('[\r\n]')
local rows = table()
for _,l in ipairs(ls) do
	if not l:match('^%s') then	-- new entry
		rows:insert(l:trim():split('%s+'))
	else	-- append to last entry
		rows[#rows]:append(l:trim():split('%s+'))
	end
end

local entries = table()
for _,row in ipairs(rows) do
	print(row:concat(' | '))
	local year, month, day
	if row[1] then
		year, month, day = row[1]:match('^(%d%d%d%d?)(%d%d)(%d%d)$')
	else
		print("bad row: "..row:concat('\t'))
	end
	local hour, min, sec
	if row[5] then
		hour, min, sec = row[5]:match('^(%d%d?)(%d%d)(%d%d%.?%d*)$')
	end
	
	print(year,month,day,hour,min,sec)
	
	local entry = {
		year=tonumber(year),
		month=tonumber(month),
		day=tonumber(day),
		hour=tonumber(hour),
		min=tonumber(min),
		sec=tonumber(sec),
		lat=tonumber(row[3]),
		lon=tonumber(row[4]),
		depth=tonumber(row[7]),	-- km
		momentMagn=tonumber(row[8]),
		momentMagnVar=tonumber(row[9]),
	}
	
	if entry.year and entry.month and entry.day then
		entry.time = os.time(entry)
		if entry.time and entry.sec then entry.time = entry.time + math.fpart(entry.sec) end
	end
	
	entry.area = row[6]
	entries:insert(entry)
end

local json = require 'json'
io.writefile('earthquakes.json', json.encode(entries))