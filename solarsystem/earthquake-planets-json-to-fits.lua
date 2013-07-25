require 'ext'
local json = require 'json'

print('reading json ...')
local planetStates = assert(json.decode(assert(io.readfile('planets.json'))))

local ffi = require 'ffi'

local numRows = #planetStates
local planetNames = {'sun','mercury','venus','earth','moon','mars','jupiter','saturn','uranus','neptune','pluto'}

local planetDateStrs = table()

--[[
data for each entry:
	date
	per planet:
		pos (3)
		vel (3)
--]]

local numCols = #planetNames * 6 + 1
local data = ffi.new('double[?]', numCols * numRows)

math.nan = 0/0

print('generating buffer...')
for yplusone,planetState in ipairs(planetStates) do
	local y = yplusone-1
	data[y * numCols + 0] = planetState.date or math.nan	--julian date
	local planetResults = planetState.results
	for i,name in ipairs(planetNames) do
		local planet = planetResults and planetResults[name]
		local pos = planet and planet[1]
		local vel = planet and planet[2]
		for j=1,3 do
			data[y * numCols + (i-1)*6 + j] = pos and pos[j] or math.nan
			data[y * numCols + (i-1)*6 + j+3] = vel and vel[j] or math.nan
		end
	end
	planetDateStrs:insert(planetState.calendarDate or 'Unknown')
end

print('writing data...')
local fitsIO = require 'fits.io'
fitsIO.save('planets.fits', numCols, numRows, 1, data)

io.writefile('planets-fits-datestrs.json', json.encode(planetDateStrs))