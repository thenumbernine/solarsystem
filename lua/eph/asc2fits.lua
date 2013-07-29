--[[
utility function for converting the ascii data into fits data

based on ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/asc2bin.cc
--]]

local denum = 406

require 'ext'

local objNames = table{
	'Mercury', 'Venus', 'EM_Bary', 'Mars', 'Jupiter', 'Saturn',
	'Uranus', 'Neptune', 'Pluto',
	'GeoCMoon',
	'Sun',
	'Nutation',
	'Libration'
};
local objIndexForName = objNames:map(function(v,k) return k,v end)

local function ephToDec(x)
	return tonumber((x:gsub('D','e')))		-- eph uses aD+-b instead of ae+-b for its number format
end

-- read header file


-- I could make this a structure ... but I'm not using it concurrently so no worries

local lines
local lineIndex

function openLineParser(filename)
	lines = io.readfile(filename):split('[\r\n]')
	lineIndex = 0
end

local function nextLine()
	repeat
		lineIndex = lineIndex + 1
	until not lines[lineIndex] or lines[lineIndex]:trim() ~= ''
	return lines[lineIndex]
end

local function nextLineWords()
	return nextLine():trim():split('%s+')
end

local function nextLineNumbers()
	return nextLineWords():map(ephToDec)
end

local function nextGroup()
	local line
	repeat
		line = nextLine()
	until line:trim() ~= ''
	local group = line:match('GROUP%s+(%d+)')
	return assert(tonumber(group))
end


openLineParser(denum..'/header.'..denum)

local hdr = table()

local ksize
ksize, hdr.numCoeffs = nextLine():match('KSIZE=%s*(%d+)%s+NCOEFF=%s*(%d+)')
hdr.numCoeffs = tonumber(hdr.numCoeffs)
assert(hdr.numCoeffs > 0)

assert(nextGroup() == 1010)

-- then comes three title lines
hdr.title = ''
do
	local sep = ''
	for i=1,3 do
		hdr.title = hdr.title .. sep .. nextLine():trim()
		sep = '\n'
	end
end

assert(nextGroup() == 1030)

hdr.epoch1, hdr.epoch2, hdr.interval = unpack(nextLineNumbers())

assert(nextGroup() == 1040)

local numConsts = assert(tonumber(nextLine():trim()))

local constNames = table()
for i=0,numConsts-1,10 do
	constNames:append(nextLineWords())
end
assert(#constNames == numConsts)

assert(nextGroup() == 1041)
assert(tonumber(nextLine():trim()) == numConsts)

local constValues = table()
for i=0,numConsts-1,3 do
	constValues:append(nextLineNumbers())
end
assert(#constValues == numConsts)

hdr.vars = constNames:map(function(name, i) return constValues[i], name end)
hdr.au = assert(hdr.vars.AU)
hdr.emrat = assert(hdr.vars.EMRAT)
hdr.DEnumber = assert(hdr.vars.DENUM)

assert(nextGroup() == 1050)

do
	local offsets = nextLineNumbers()
	local numCoeffs = nextLineNumbers()
	local numSubIntervals = nextLineNumbers()

	local numObj = #offsets
	assert(numObj == #numCoeffs)
	assert(numObj == #numSubIntervals)

	hdr.objs = table()
	for i=1,numObj do
		hdr.objs[i] = table{
			name = objNames[i],	-- lowercase?
			offset = offsets[i],
			numCoeffs = numCoeffs[i],
			numSubIntervals = numSubIntervals[i],
			numComponents = i == objIndexForName.Nutation and 2 or 3,
		}
	end
end

assert(nextGroup() == 1070)

local json = require 'json'
io.writefile(denum..'/header.json', json.encode(hdr))	-- pretty print?


--[[
read ascii / write raw
--]]

local ffi = require 'ffi'
local lfs = require 'lfs'
--local fitsIO = require 'fits.io'

--lfs.mkdir(denum..'/fits/')
lfs.mkdir(denum..'/f64/')

require 'ffi.c.stdio'

-- one lump file
local file = ffi.C.fopen(denum..'/f64/de'..denum..'.f64.raw', 'wb')
assert(file ~= nil)

for year = -3000,2900,100 do
	local yearName = (year < 0 and 'm' or 'p') .. ('%04d'):format(math.abs(year))
	local srcFilename = denum..'/asc' .. yearName .. '.'..denum
	print(srcFilename)
	openLineParser(srcFilename)

	local records = table()	-- store in lua before converting to C
	while true do
		local line = nextLine()
		if not line then break end
		local recordNum, numCoeffs = unpack(line:trim():split('%s+'):map(ephToDec))
		numCoeffs = tonumber(numCoeffs)
		assert(numCoeffs == hdr.numCoeffs)

		local record = table()
		for i=0,hdr.numCoeffs-1,3 do
			record:append(nextLineNumbers())
		end
		assert(#record >= hdr.numCoeffs and #record < hdr.numCoeffs+3)	
		records:insert(record)
	end

	local data = ffi.new('double[?]', #records * hdr.numCoeffs)
	for y=0,#records-1 do
		local record = records[y+1]
		for x=0,hdr.numCoeffs-1 do
			data[x+hdr.numCoeffs*y] = record[x+1]
		end
	end

	--[[
	convert to FITS files
		width = numCoeffs
		height = numRecords
		channels = 1
		
	the buffer passed into the fitsIO rw is column-major with the origin at the upper-left
	(to keep my fits io rw buffers compatible with glTexImage2D)
	(so in the fits file the data will be transposed and flipped)
	--]]
	--[[
	fitsIO.save(denum..'/fits/' .. yearName .. '.fits', hdr.numCoeffs, #records, 1, data)
	--]]
	
	-- write the raw data out too ... for fast seeking
	-- separate files:
	--local file = assert(ffi.C.fopen(denum..'/f64/'..yearName..'.f64.raw', 'wb'))
	-- one lump file
	local num = ffi.C.fwrite(data, ffi.sizeof('double'), hdr.numCoeffs * #records, file)
	assert(num == hdr.numCoeffs * #records)
	--ffi.C.fclose(file)
end
ffi.C.fclose(file)
