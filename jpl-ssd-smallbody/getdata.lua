#!/usr/bin/env lua -lluarocks.require
require 'ext'
local julian = assert(loadfile('../lua/julian.lua'))()	-- in ../lua/julian.lua
local json = require 'dkjson'

local function readCache(filename, url)
	if io.fileexists(filename) then return io.readfile(filename) end
	local http = require 'socket.http'
	local results = {http.request(url)}
	print(unpack(results))
	local data = results[1]
	io.writefile(filename, data)
	return data
end

local numberedAsteroids = readCache('ELEMENTS.NUMBR', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR')
local unnumberedAsteriods = readCache('ELEMENTS.UNNUM', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM')
readCache('ELEMENTS.COMET', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET')

local auInM = 149597870700	-- 1AU in meters

--[[
holds list of {start, finish, name} for each column
--]]
local Columns = class()

function Columns:init(lines)
	local line1, line2 = lines[1], lines[2]
	self.columns = table()
	local current = 1
	while true do
		local start, finish = line2:find('%-+', current)
		if not start then break end
		assert(finish)
		self.columns:insert{
			start = start,
			finish = finish,
			name = line1:sub(start,finish):trim(),
		}
		current = finish + 1
	end
	self.columnsByName = self.columns:map(function(column)
		return column, column.name
	end)
end

local ColumnAccess = class()

function ColumnAccess:init(columns, line)
	self.columns = assert(columns)
	self.line = assert(line)
end
	
function ColumnAccess:__index(name)
	local m = getmetatable(self)
	if m[name] then return m[name] end
	local columns = assert(rawget(self, 'columns'))
	local line = assert(rawget(self, 'line'))
	local col = columns.columnsByName[name] or error("failed to find column "..tostring(name))
	local start, finish = col.start, col.finish
	if col == columns.columns:last() then 	-- allow last column to read to end of line
		return line:sub(start)
	else
		return line:sub(start, finish)
	end
end

function Columns:__call(line)
	return ColumnAccess(self, line)
end

--[[
args:
	inputFile (text)
	outputFile (json)
--]]
local function processFile(args)
	local inputFilename = assert(args.inputFilename)
	local outputFilename = assert(args.outputFilename)
	local variableName = assert(args.variableName)
	local processRow = assert(args.processRow)

	local lastTime = os.time()
	local lines = io.readfile(inputFilename):split('\n')
	local cols = Columns(lines)
	local outputLines = table()
	outputLines:insert(variableName..' = [')
	for i=3,#lines do	-- skip headers and ---- line
		local line = lines[i]
		if #line:trim() > 0 then
			xpcall(function()
				local row = cols(line)
				local body = assert(processRow(row))
				outputLines:insert('\t'..json.encode(body)..',')
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
	outputLines:insert('];')
	io.writefile(outputFilename, outputLines:concat('\n'))

end

-- process comets
processFile{
	inputFilename = 'ELEMENTS.COMET',
	outputFilename = 'comets.json',
	variableName = 'cometData',
	processRow = function(row)
		local body = {}
		local numberAndName = row['Num  Name']
		body.idNumber = numberAndName:sub(1,4):trim()
		body.name = numberAndName:sub(6,43):trim()
		body.epoch = assert(tonumber(row.Epoch:trim()))
		body.perihelionDistance = assert(tonumber(row.q:trim())) * auInM
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclination = math.rad(assert(tonumber(row.i:trim())))	-- wrt J2000 ecliptic plane
		body.argumentOfPerihelion = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
	
		-- time of perihelion passage
		do
			local str = row.Tp
			local year = assert(tostring(str:sub(1,4)))
			local month = assert(tostring(str:sub(5,6)))
			local day = assert(tostring(str:sub(7)))
			
			-- TODO calendarToJulian() function.  this is a rough rough guess.
			body.timeOfPerihelionPassage = julian.fromCalendar{year=year, month=month, day=day}
		end
		
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}

-- process numbered bodies 
processFile{
	inputFilename = 'ELEMENTS.NUMBR',
	outputFilename = 'smallbodies-numbered.json',
	variableName = 'numberedSmallBodyData',
	processRow = function(row)
		local body = {}
		body.idNumber = row.Num:trim()
		body.name = row.Name:trim()
		body.epoch = assert(tonumber(row.Epoch:trim()))
		body.semiMajorAxis = assert(tonumber(row.a:trim()))
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclindation = math.rad(assert(tonumber(row.i:trim())))
		body.argumentOfPerihelion = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
		body.meanAnomaly = math.rad(assert(tonumber(row.M:trim())))
		body.absoluteMagnitude = assert(tonumber(row.H:trim()))
		body.magnitudeSlopeParameter = assert(tonumber(row.G:trim()))
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}

-- process unnumbered bodies 
processFile{
	inputFilename = 'ELEMENTS.UNNUM',
	outputFilename = 'smallbodies-unnumbered.json',
	variableName = 'unnumberedSmallBodyData',
	processRow = function(row)
		local body = {}
		body.name = row.Designation:trim()
		body.epoch = assert(tonumber(row.Epoch:trim()))
		body.semiMajorAxis = assert(tonumber(row.a:trim()))
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclindation = math.rad(assert(tonumber(row.i:trim())))
		body.argumentOfPerihelion = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
		body.meanAnomaly = math.rad(assert(tonumber(row.M:trim())))
		body.absoluteMagnitude = assert(tonumber(row.H:trim()))
		body.magnitudeSlopeParameter = assert(tonumber(row.G:trim()))
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}

