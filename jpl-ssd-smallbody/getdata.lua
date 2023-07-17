#!/usr/bin/env lua
require 'ext'
local https = require 'ssl.https'
local ltn12 = require 'ltn12'

local function download(filename, url)
	if path(filename):exists() then
		print('already have file '..filename..', so skipping the download and using the cached version')
		return
	end
	print('downloading url '..url..' ...')
	local data = table()
	assert(https.request{
		url = url,
		sink = ltn12.sink.table(data),
		protocol = 'tlsv1',
	})
	data = data:concat()
	print('writing file '..filename..' with this much data: '..#data)
	path(filename):write(data)
end

--[[ "no protocols available ...
download('ELEMENTS.NUMBR', 'https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR')
download('ELEMENTS.UNNUM', 'https://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM')
download('ELEMENTS.COMET', 'https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET')
--]]
-- [[
print(os.execute'wget https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR')
print(os.execute'wget https://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM')
print(os.execute'wget https://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET')
--]]
