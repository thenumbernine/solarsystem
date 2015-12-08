#!/usr/bin/env lua -lluarocks.require
--[[
script to convert text files to another format
either JSON or SQLite3 at the moment 
--]]

require 'ext'

local julian = assert(loadfile('../horizons/julian.lua'))()
local json = require 'dkjson'

local function readCache(filename, url)
	if io.fileexists(filename) then return file[filename] end
	local http = require 'socket.http'
	local results = {http.request(url)}
	print(table.unpack(results))
	local data = results[1]
	file[filename] = data
	return data
end

local numberedAsteroids = readCache('ELEMENTS.NUMBR', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR')
local unnumberedAsteriods = readCache('ELEMENTS.UNNUM', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM')
readCache('ELEMENTS.COMET', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET')

local auInM = 149597870700	-- 1AU in meters

local Columns = require 'columns'

--[[
args:
	inputFile (text)
	processRow = function(line)
	outputMethod
--]]
local function processToFile(args)
	local inputFilename = assert(args.inputFilename)
	local processRow = assert(args.processRow)
	local outputMethod = assert(args.outputMethod)

	local lastTime = os.time()
	local lines = file[inputFilename]:split('\n')
	local cols = Columns(lines)
	for i=3,#lines do	-- skip headers and ---- line
		local line = lines[i]
		if #line:trim() > 0 then
			xpcall(function()
				outputMethod:processBody(assert(processRow(cols(line))))
			end, function(err)
				io.stderr:write('failed on file '..inputFilename..' line '..i..'\n')
				io.stderr:write(err..'\n'..debug.traceback()..'\n')
				os.exit(1)
			end)
		end
		local thisTime = os.time()
		if thisTime ~= lastTime then
			print(inputFilename..' '..math.floor(100*(i-2)/(#lines-2))..'% complete...') 
			lastTime = thisTime
		end
	end
	outputMethod:done()
end

--local outputMethod = require 'output_json' 
local outputMethod = require 'output_sqlite3'
--local outputMethod = require 'output_points'

outputMethod:staticInit()

-- process comets
processToFile{
	inputFilename = 'ELEMENTS.COMET',
	outputMethod = outputMethod{
		filename = 'comets.json',
		variableName = 'cometData',
		bodyType = 0,
	},
	processRow = function(row)
		local body = {}
		local numberAndName = row['Num  Name']
		body.idNumber = numberAndName:sub(1,4):trim()
		body.name = numberAndName:sub(6,43):trim()
		body.epoch = assert(tonumber(row.Epoch:trim()))
		body.perihelionDistance = assert(tonumber(row.q:trim())) * auInM
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclination = math.rad(assert(tonumber(row.i:trim())))	-- wrt J2000 ecliptic plane
		body.argumentOfPeriapsis = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
	
		-- time of perihelion passage
		do
			local str = row.Tp
			local year = assert(tonumber(str:sub(1,4)))
			local month = assert(tonumber(str:sub(5,6)))
			local day = assert(tonumber(str:sub(7)))
			-- TODO calendarToJulian() function.  this is a rough rough guess.
			body.timeOfPerihelionPassage = julian.fromCalendar{year=year, month=month, day=day}
		end
		
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}

-- process numbered bodies 
processToFile{
	inputFilename = 'ELEMENTS.NUMBR',
	outputMethod = outputMethod{
		filename = 'smallbodies-numbered.json',
		variableName = 'numberedSmallBodyData',
		bodyType = 1,
	},
	processRow = function(row)
		local body = {}
		body.idNumber = row.Num:trim()
		body.name = row.Name:trim()
		body.epoch = assert(tonumber(row.Epoch:trim()))
		body.semiMajorAxis = assert(tonumber(row.a:trim())) * auInM
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclination = math.rad(assert(tonumber(row.i:trim())))
		body.argumentOfPeriapsis = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
		body.meanAnomaly = math.rad(assert(tonumber(row.M:trim())))
		body.absoluteMagnitude = assert(tonumber(row.H:trim()))
		body.magnitudeSlopeParameter = assert(tonumber(row.G:trim()))
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}

-- process unnumbered bodies 
processToFile{
	inputFilename = 'ELEMENTS.UNNUM',
	outputMethod = outputMethod{
		filename = 'smallbodies-unnumbered.json',
		variableName = 'unnumberedSmallBodyData',
		bodyType = 2,
	},
	processRow = function(row)
		local body = {}
		body.name = row.Designation:trim()
		body.epoch = assert(tonumber(row.Epoch:trim()))
		body.semiMajorAxis = assert(tonumber(row.a:trim())) * auInM
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclination = math.rad(assert(tonumber(row.i:trim())))
		body.argumentOfPeriapsis = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
		body.meanAnomaly = math.rad(assert(tonumber(row.M:trim())))
		body.absoluteMagnitude = assert(tonumber(row.H:trim()))
		body.magnitudeSlopeParameter = assert(tonumber(row.G:trim()))
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}

outputMethod:staticDone()
