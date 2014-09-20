#!/usr/bin/env lua -lluarocks.require
require 'ext'

-- TODO write out to sqlite3, see if that saves us any space

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

local OutputToJSON = class()

--[[
args:
	filename
	variableName
--]]
function OutputToJSON:init(args)
	self.filename = assert(args.filename)
	self.variableName = assert(args.variableName)

	self.outputLines = table()
	self.outputLines:insert(self.variableName..' = [')
end

function OutputToJSON:processBody(body)
	self.outputLines:insert('\t'..json.encode(body)..',')
end

function OutputToJSON:done()
	self.outputLines:insert('];')
	io.writefile(self.filename, self.outputLines:concat('\n'))
end

local OutputToSQLite3 = class()

--[[
args:
	variableName = the database filename, given a suffix of .sqlite3
--]]
function OutputToSQLite3:init(args)
	local databaseName = assert(args.variableName)

	local luasql = require 'luasql.sqlite3'
	self.env = assert(luasql.sqlite3())
	self.conn = assert(self.env:connect(databaseName..'.sqlite3'))
end

function OutputToSQLite3:processBody(body)
	if not self.createdTable then
		self.createdTable = true
		
		self.tableName = 'data'
		self.conn:execute('drop table '..self.tableName)	--might fail if the table is new

		self.columns = table()

		local luaToSqlite3Types = {
			number = 'real',
			string = 'text',
		}

		local columnDescs = table()
		for columnName,value in pairs(body) do
			local columnType = assert(luaToSqlite3Types[type(value)])
			self.columns:insert{
				name = columnName,
				type = columnType
			}
			columnDescs:insert(columnName .. ' ' .. columnType)
		end
		local cmd = 'create table '..self.tableName..'('..columnDescs:concat(', ')..')'
		print(cmd)	
		-- don't create the table until the first body
		assert(self.conn:execute(cmd))
	end

	local valueStrs = table()
	for _,column in ipairs(self.columns) do
		local value = body[column.name]
		if column.type == 'text' then
			value = ('%q'):format(value)
		end
		valueStrs:insert(value)
	end
	local cmd = 'insert into '..self.tableName..' values ('..valueStrs:concat(',')..')'
	print(cmd)	
	assert(self.conn:execute(cmd))
end

function OutputToSQLite3:done()
	self.conn:close()
	self.env:close()
end

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
	local lines = io.readfile(inputFilename):split('\n')
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

--local outputMethod = OutputToJSON
local outputMethod = OutputToSQLite3

-- process comets
processToFile{
	inputFilename = 'ELEMENTS.COMET',
	outputMethod = outputMethod{
		filename = 'comets.json',
		variableName = 'cometData',
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
processToFile{
	inputFilename = 'ELEMENTS.NUMBR',
	outputMethod = outputMethod{
		filename = 'smallbodies-numbered.json',
		variableName = 'numberedSmallBodyData',
	},
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
processToFile{
	inputFilename = 'ELEMENTS.UNNUM',
	outputMethod = outputMethod{
		filename = 'smallbodies-unnumbered.json',
		variableName = 'unnumberedSmallBodyData',
	},
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

