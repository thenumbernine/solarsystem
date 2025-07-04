#!/usr/bin/env luajit
-- script to convert text files to another format

require 'ext'
local template = require 'template'
--local gcmem = require 'ext.gcmem'
local ffi = require 'ffi'

local OutputPoints = require 'output_points'

local Julian = require 'solarsystem.julian'
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

	local data = path(inputFilename):read()

	local lines = assert(data, "failed to read file "..inputFilename):split'\n'
	local cols = Columns(lines)

print(tolua(cols.columns))

	for i=3,#lines do	-- skip headers and ---- line
		local line = lines[i]
		if #line:trim() > 0 then
			assert(xpcall(function()
				outputObj:processBody(assert(processRow(cols(line))))
			end, function(err)
				return 'failed on file '..inputFilename..' line '..i..'\n'
					..err..'\n'
					..debug.traceback()
			end))
		end
		local thisTime = os.time()
		if thisTime ~= lastTime then
			print(inputFilename..' '..math.floor(100*(i-2)/(#lines-2))..'% complete...')
			lastTime = thisTime
		end
	end
end

OutputPoints:staticInit()

-- these match the fields in coldesc.lua and probably row-desc.json
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

local real
do
	real = 'double'
	local code = template([[
typedef <?=real?> real;

typedef struct {
<? for _,field in ipairs(numberFields) do
?>	real <?=field?>;
<? end
?>
	int bodyType;	// 0=comet, 1=numbered, 2=unnum
	int horizonID;	//for numbered asteroids only
	char name[43+1];
	char orbitSolutionReference[12+1];

	//computed parameters:
	long index;
	real pos[3], vel[3], A[3], B[3];
	real eccentricAnomaly;
	real timeOfPeriapsisCrossing;
	real meanAnomaly;
	int orbitType;	// 0=elliptic, 1=hyperbolic, 2=parabolic ... this can be derived from eccentricity
	real orbitalPeriod;
} body_t;
]], {
	real = real,
	numberFields = numberFields,
})
	ffi.cdef(code)

	local allFields = numberFields:mapi(function(name)
		return {name=name, type=real}
	end):append{
		{name='bodyType', type='int'},
		{name='horizonID', type='int'},
		{name='name', type='char[44]'},
		{name='orbitSolutionReference', type='char[13]'},
		{name='index', type='long'},
		{name='pos', type=real..'[3]'},
		{name='vel', type=real..'[3]'},
		{name='A', type=real..'[3]'},
		{name='B', type=real..'[3]'},
		{name='eccentricAnomaly', type=real},
		{name='timeOfPeriapsisCrossing', type=real},
		{name='meanAnomaly', type=real},
		{name='orbitType', type='int'},
		{name='orbitalPeriod', type=real},
	}

	local f = assert(path'body_t_desc.lua':open'w')
	f:write'return {\n'
	f:write('\tname = '..('%q'):format'body_t'..',\n')
	f:write('\tsize = '..ffi.sizeof'body_t'..',\n')
	f:write('\tjulianDay = '..OutputPoints.julianDay..',\n')
	f:write'\tfields = {\n'
	for _,field in ipairs(allFields) do
		f:write('\t\t'..field.name..' = {'
			..'type='..('%q'):format(field.type)..', '
			..'offset='..ffi.offsetof('body_t', field.name)..', '
			..'size='..ffi.sizeof(ffi.new(field.type))..', '
		..'},\n')
	end
	f:write'\t},\n'
	f:write'}\n'
end

local function newBody()
	local body = setmetatable({
		_ptr = ffi.new'body_t[1]', --gcmem.new('body_t',1),
		getPos = function(self)
			return self._ptr[0].pos[0], self._ptr[0].pos[1], self._ptr[0].pos[2]
		end,
	}, {
		-- can you use cdata as __index and __newindex values
		-- similar to how you can use tables?
		__index = function(self,k)
			return self._ptr[0][k]
		end,
		__newindex = function(self,k,v)
			self._ptr[0][k] = v
		end
	})
	for _,field in ipairs(numberFields) do
		body[field] = math.nan
	end
	body.horizonID = 0
	body.absoluteMagnitude = math.huge
	return body
end

-- [[ process comets
processToFile{
	inputFilename = 'ELEMENTS.COMET',
	outputObj = OutputPoints{
		filename = 'comets.json',
		variableName = 'cometData',
		bodyType = 0,
	},
	processRow = function(row)
		local body = newBody()
		local numberAndName = row['Num  Name']
		body.name = numberAndName:trim()
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
			body.timeOfPerihelionPassage = Julian.fromCalendar{year=year, month=month, day=day}
		end

		body.orbitSolutionReference = row.Ref:trim()
		return body
	end,
}
--]]

-- [[ process numbered bodies
processToFile{
	inputFilename = 'ELEMENTS.NUMBR',
	outputObj = OutputPoints{
		filename = 'smallbodies-numbered.json',
		variableName = 'numberedSmallBodyData',
		bodyType = 1,
	},
	processRow = function(row)
		local body = newBody()
		body.horizonID = assert(tonumber(row.Num:trim()))
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
	outputObj = OutputPoints{
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

OutputPoints:staticDone()
