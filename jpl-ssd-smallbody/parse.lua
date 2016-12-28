#!/usr/bin/env luajit
-- script to convert text files to another format

require 'ext'

local ffi = require 'ffi'

local outputMethod = ... or 'points'
print('outputMethod = '..outputMethod)

--local outputClass = require 'output_sqlite3'
--local outputClass = require 'output_points'
local outputClass = require('output_'..outputMethod)

local julian = assert(loadfile('../horizons/julian.lua'))()
local mjdOffset = 2400000.5
local auInM = 149597870700	-- 1AU in meters

local Columns = require 'columns'

--[[
args:
	inputFile (text)
	processRow = function(line)
	outputObj
--]]
local function processToFile(args)
	local inputFilename = assert(args.inputFilename)
	local processRow = assert(args.processRow)
	local outputObj = assert(args.outputObj)

	local lastTime = os.time()
	local lines = assert(file[inputFilename], "failed to read file "..inputFilename):split('\n')
	local cols = Columns(lines)
	for i=3,#lines do	-- skip headers and ---- line
		local line = lines[i]
		if #line:trim() > 0 then
			xpcall(function()
				outputObj:processBody(assert(processRow(cols(line))))
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
	outputObj:done()
end

outputClass:staticInit()

-- these match the fields in coldesc.lua and body_write_t
local numberFields = table{ 
	'epoch',
	'perihelionDistance',		--comets
	'semiMajorAxis',		--asteroids
	'eccentricity',
	'inclination',
	'argumentOfPeriapsis',
	'longitudeOfAscendingNode',
	'meanAnomalyAtEpoch',		--asteroids
	'absoluteMagnitude',		--asteroids
	'magnitudeSlopeParameter',		--asteroids
	'timeOfPerihelionPassage',		--comets
}

if ffi then
	local template = require 'template'
	local code = template([[
typedef double real;

typedef struct {
<? for _,field in ipairs(numberFields) do
?>	real <?=field?>;
<? end
?>
	int bodyType;	// 0=comet, 1=numbered, 2=unnum
	char idNumber[4+1];	//+1 for null-term because I'm lazy and using string assignment in ffi
	char name[38+1];
	char orbitSolutionReference[10+1];

	//computed parameters:
	long index;
	real pos[3], A[3], B[3];
	real eccentricAnomaly;
	real timeOfPeriapsisCrossing;
	real meanAnomaly;
	int orbitType;	// 0=elliptic, 1=hyperbolic, 2=parabolic ... this can be derived from eccentricity
	real orbitalPeriod;
} body_t;
]], {
	numberFields = numberFields,
})
	ffi.cdef(code)
end

local bodyMT = not ffi and {
	__index = {
		getPos = function(self)
			return unpack(self.pos)
		end,
	},
} or {
	__index = {
		getPos = function(self)
			return self.pos[0], self.pos[1], self.pos[2]
		end,
	},
}

local bodyMetaType
if ffi then
	bodyMetaType = ffi.metatype('body_t', bodyMT) 
end

local function newBody()
	if not ffi then return setmetatable({}, bodyMT) end
	local body = bodyMetaType()
	for _,field in ipairs(numberFields) do
		body[field] = math.nan
	end
	body.absoluteMagnitude = math.huge
	return body
end

-- [[ process comets
processToFile{
	inputFilename = 'ELEMENTS.COMET',
	outputObj = outputClass{
		filename = 'comets.json',
		variableName = 'cometData',
		bodyType = 0,
	},
	processRow = function(row)
		local body = newBody()
		local numberAndName = row['Num  Name']
		body.idNumber = numberAndName:sub(1,4):trim()
		body.name = numberAndName:sub(6,43):trim()
		body.epoch = assert(tonumber(row.Epoch:trim())) + mjdOffset
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
--]]

-- [[ process numbered bodies 
processToFile{
	inputFilename = 'ELEMENTS.NUMBR',
	outputObj = outputClass{
		filename = 'smallbodies-numbered.json',
		variableName = 'numberedSmallBodyData',
		bodyType = 1,
	},
	processRow = function(row)
		local body = newBody()
		body.idNumber = row.Num:trim()
		body.name = row.Name:trim()
		body.epoch = assert(tonumber(row.Epoch:trim())) + mjdOffset
		body.semiMajorAxis = assert(tonumber(row.a:trim())) * auInM
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclination = math.rad(assert(tonumber(row.i:trim())))
		body.argumentOfPeriapsis = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
		body.meanAnomalyAtEpoch = math.rad(assert(tonumber(row.M:trim())))
		body.absoluteMagnitude = assert(tonumber(row.H:trim()))
		body.magnitudeSlopeParameter = assert(tonumber(row.G:trim()))
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}
--]]

-- [[ process unnumbered bodies 
processToFile{
	inputFilename = 'ELEMENTS.UNNUM',
	outputObj = outputClass{
		filename = 'smallbodies-unnumbered.json',
		variableName = 'unnumberedSmallBodyData',
		bodyType = 2,
	},
	processRow = function(row)
		local body = newBody()
		body.name = row.Designation:trim()
		body.epoch = assert(tonumber(row.Epoch:trim())) + mjdOffset
		body.semiMajorAxis = assert(tonumber(row.a:trim())) * auInM
		body.eccentricity = assert(tonumber(row.e:trim()))
		body.inclination = math.rad(assert(tonumber(row.i:trim())))
		body.argumentOfPeriapsis = math.rad(assert(tonumber(row.w:trim())))
		body.longitudeOfAscendingNode = math.rad(assert(tonumber(row.Node:trim())))
		body.meanAnomalyAtEpoch = math.rad(assert(tonumber(row.M:trim())))
		body.absoluteMagnitude = assert(tonumber(row.H:trim()))
		body.magnitudeSlopeParameter = assert(tonumber(row.G:trim()))
		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}
--]]

outputClass:staticDone()