#!/usr/bin/env lua
--[[
I made this solar system visualization app way back in 2013 (or earlier?)
But looks like finally in 2022 or so NASA switched from telnet to http.
Good job NASA, welcome to the information age of HTTP.  Taxpayer money hard at work there.
So here is me getting around to updating this app.

https://ssd-api.jpl.nasa.gov/doc/horizons.html

...
... what?
... the JSON is just the text telnet output wrapped in a JSON array?
... 
...
... smh
--]]
local table = require 'ext.table'
local URL = require 'socket.url'
local parsed = URL.parse[[https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='499'&OBJ_DATA='YES'&MAKE_EPHEM='YES'&EPHEM_TYPE='OBSERVER'&CENTER='500@399'&START_TIME='2006-01-01'&STOP_TIME='2006-01-20'&STEP_SIZE='1%20d'&QUANTITIES='1,9,20,23,24,29']]
print(require 'ext.tolua'(parsed))
local query = table{
	command = '499',
	obj_data = 'YES',
	make_ephem = 'YES',
	ephem_type = 'OBSERVER',
	center = '500@399',
	start_time = '2006-01-01',
	stop_time = '2006-01-20',
	step_size = '1 d',
	quantities = '1,9,20,23,24,29',
}
local parsed = {
	scheme = 'https',
	host = 'ssd.jpl.nasa.gov',
	path = '/api/horizons.api',
	query = 'format=text&'..query:map(function(v,k,t)
		-- I'm escaping single-quotes from ' to \' since in the example all values are wrapped in '' quotes.  Not sure if this is necessary or how Horizons internally handles single-quotes in values.
		return k:upper().."='"..URL.escape(v:gsub("'", "\\'")).."'", #t+1
	end):concat'&',
}
print('rebuilt:')
print(URL.build(parsed))
