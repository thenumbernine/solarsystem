#!/usr/bin/env luajit
--[[
extract GAIA star data
--]]

-- TODO maybe universal FITS -> raw / json / lua / csv tool ?

local ffi = require 'ffi'
local fits = require 'ffi.req' 'fitsio'
local stdio = require 'ffi.req' 'c.stdio'
local table = require 'ext.table'

-- add 'last' to the last element of the parameter pack
local function addlast(last, ...)
	if select('#', ...) == 0 then return last end
	return (...), addlast(last, select(2, ...))
end

-- TODO put in ffi/fitsio.lua also?
local function fitsStatus(f, ...)
	local status = ffi.new('int[1]', 0)
	local res = f(addlast(status, ...))
	return res, status[0]
end

local function fitsError(status)
	fits.fits_report_error(stdio.stderr, status)
	local buffer = ffi.new('char[256]')	-- TODO size?
	ffi.fill(buffer, 0, ffi.sizeof(buffer))
	fits_get_errstatus(status, buffer)
	-- https://heasarc.gsfc.nasa.gov/fitsio/c/c_user/node34.html
	-- this says the buffer is 80 characters
	-- also says to call repeatedly until it returns 0 ...
	-- but it's a void functin
	-- maybe they meant until ... hmm there's no error-code-getter function either ...
	-- *shrug*
	error("FITS error "..status .. ': ' ..ffi.string(buffer))
end

local function fitsSafe(f, ...)
	local res, status = fitsStatus(f, ...)
	if status ~= 0 then
		fitsError(status)
	end
	return res
end

-- [[ reverse-mapping for the T* FITS type codes
local fitsNameForType = {}
for _,name in ipairs{'TBIT', 'TBYTE', 'TSBYTE', 'TLOGICAL', 'TSTRING', 'TUSHORT', 'TSHORT', 'TUINT', 'TINT', 'TULONG', 'TLONG', 'TINT32BIT', 'TFLOAT', 'TULONGLONG', 'TLONGLONG', 'TDOUBLE', 'TCOMPLEX', 'TDBLCOMPLEX'} do
	local val = fits[name]
	assert(val)	-- if fits wasn't wrapped then ffi would throw an error if the key wasn't present
	fitsNameForType[val] = name
end

local fitsTypeForCType = {
	bool = 'TBIT',
	int8_t = 'TSBYTE',
	uint8_t = 'TBYTE',
	--TLOGICAL,
	--TSTRING,
	uint16_t = 'TUSHORT',
	int16_t = 'TSHORT',
	-- FITS standards here and my gcc compiler start to disagree so 
	--TUINT,	-- is this a short or long?  depends? does FITS think INT is 16-bit?
	--TINT,
	uint32_t = 'TULONG',	-- fits seems to say LONG is 32-bit ... olllld C standard
	int32_t = 'TLONG',
	--TINT32BIT,
	uint64_t = 'TULONGLONG',
	int64_t = 'TLONGLONG',
	float = 'TFLOAT',
	double = 'TDOUBLE',
	['complex double'] = 'TCOMPLEX',
	--TDBLCOMPLEX,
}
local ctypeForFitsType = table.map(fitsTypeForCType, function(v,k) return k,v end):setmetatable(nil)

local filename = '1512497275987O-result.fits'

local filep = ffi.new('fitsfile*[1]', nil)
fitsSafe(fits.fits_open_table, filep, filename, fits.READONLY)
--file = file[1]	-- will this refcount and free early? yes.
local file = filep[0]		-- so give it a dif name

local numRows = ffi.new('long[1]', 0)
fitsSafe(fits.fits_get_num_rows, file, numRows)
numRows = numRows[0]

local numCols = ffi.new('int[1]', 0)
fitsSafe(fits.fits_get_num_cols, file, numCols)
numCols = numCols[0]

print('numCols:', numCols)
print('numRows:', numRows)

-- get col names:
-- reuse status between fits_get_colname calls for enumeration
do
	local status = ffi.new('int[1]', 0)
	while true do
		local colName = ffi.new'char[256]'
		ffi.fill(colName, 0)
		local colNum = ffi.new('int[1]', 0)
		colNum[0] = 0
		fits.fits_get_colname(file, fits.CASESEN, ffi.cast('char*', '*'), colName, colNum, status)
		if status[0] == fits.COL_NOT_FOUND then break end
		if status[0] ~= 0 
		and status[0] ~= fits.COL_NOT_UNIQUE 
		then
			fitsError(status[0])
		end
		colName = ffi.string(colName)
		
		local fitsType = ffi.new('int[1]', 0)
		-- FITS own macros use "long" for int32 and "long long" for int64, but the modern standard is "int" for int32, "long" and "long long" for int64
		--  so how come I have a funny feeling there might be memory errors in FITS with primitives improperly cast / overlapping writes in memory? 
		local colRepeat = ffi.new('long[1]', 0)
		local colWidth = ffi.new('long[1]', 0)
		fitsSafe(fits.fits_get_coltype, file, colNum[0], fitsType, colRepeat, colWidth)

		local fitsTypeName = assert(fitsNameForType[fitsType[0]], "failed to find FITS type name for FITS type "..tostring(fitsType[0]))
		local ctypeName = assert(ctypeForFitsType[fitsTypeName], "failed to find ctype name for FITS type name "..tostring(fitsTypeName))
--print(colNum[0], colName, fitsTypeName, colRepeat[0], colWidth[0])
		assert(colWidth[0] == ffi.sizeof(ctypeName))	-- should always be true? except for arrays?  can you have cols of arrays?
		assert(colRepeat[0] == 1)
		print(colNum[0], ctypeName, colName) 
	end
end
