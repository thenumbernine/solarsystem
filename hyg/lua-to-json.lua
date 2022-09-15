-- if you have to run luajit separate of lua ...
require 'ext'
local json = require 'dkjson'
-- assigns globals:
local namedStars = fromlua(file'namedStars.lua':read())
local constellations = fromlua(file'constellations.lua':read())
-- writes files:
file'namedStars.json':write('namedStars = ' .. json.encode(namedStars, {indent=true}) ..';')
file'constellations.json':write('constellations = '..json.encode(constellations, {indent=true}) .. ';')
