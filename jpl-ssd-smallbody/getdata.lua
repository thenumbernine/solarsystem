#!/usr/bin/env lua
local http = require 'socket.http'
require 'ext'

local function download(filename, url)
	if io.fileexists(filename) then
		print('already have file '..filename..', so skipping the download and using the cached version')
		return
	end
	print('downloading url '..url..' ...')
	local results = {assert(http.request(url))}
	local data = assert(results[1], "didn't get any data back from our request")
	print('writing file '..filename..' with this much data: '..#data)
	file[filename] = data
end

download('ELEMENTS.NUMBR', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR')
download('ELEMENTS.UNNUM', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM')
download('ELEMENTS.COMET', 'http://ssd.jpl.nasa.gov/dat/ELEMENTS.COMET')
