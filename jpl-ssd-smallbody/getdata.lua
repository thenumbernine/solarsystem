#!/usr/bin/env lua -lluarocks.require
require 'ext'
local socket = require 'socket'
local http = require 'socket.http'
local json = require 'dkjson'
local julian = assert(loadfile('../lua/julian.lua'))()	-- in ../lua/julian.lua

local function readCache(filename, url)
	if io.fileexists(filename) then return io.readfile(filename) end
	local results = {http.request(url)}
	print(unpack(results))
	local data = results[1]
	io.writefile(filename, data)
	return data
end

local numberedAsteroids = readCache('ELEMENTS.NUMBR', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR')
local unnumberedAsteriods = readCache('ELEMENTS.UNNUM', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM')
local comets = readCache('ELEMENTS.COMET', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET')

local auInM = 149597870700	-- 1AU in meters

-- process comets
local lines = io.readfile('ELEMENTS.COMET'):split('\n')
local cometData = {}
for i=3,#lines do	-- skip headers and ---- line
	local comet = {}
	local line = lines[i]
	if #line:trim() > 0 then
		xpcall(function()
			comet.idNumber = line:sub(1,4):trim()
			comet.name = line:sub(6,43):trim()
			comet.epoch = assert(tonumber(line:sub(45,51):trim()))
			comet.perihelionDistance = assert(tonumber(line:sub(53,63):trim())) * auInM
			comet.eccentricity = assert(tonumber(line:sub(65,74):trim()))
			comet.inclination = assert(tonumber(line:sub(76,84):trim()))	-- wrt J2000 ecliptic plane
			comet.augmentOfPerihelion = assert(tonumber(line:sub(86,94):trim()))
			comet.longitudeOfAscendingNode = assert(tonumber(line:sub(96,104):trim()))
		
			-- time of perihelion passage
			do
				local str = line:sub(106,119)
				local year = assert(tostring(str:sub(1,4)))
				local month = assert(tostring(str:sub(5,6)))
				local day = assert(tostring(str:sub(7)))
				
				-- TODO calendarToJulian() function.  this is a rough rough guess.
				comet.timeOfPerihelionPassage = julian.fromCalendar{year=year, month=month, day=day}
			end
			
			comet.orbitSolutionReference = line:sub(121):trim()
			table.insert(cometData, comet)
		end, function(err)
			io.stderr:write('failed on line '..i..'\n')
			io.stderr:write(err..'\n'..debug.traceback()..'\n')
			os.exit()
		end)
	end
end

io.writefile('comets.json', 'cometData = '..json.encode(cometData, {indent=true}))

