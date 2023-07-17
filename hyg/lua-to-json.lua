-- if you have to run luajit separate of lua ...
require 'ext'
local json = require 'dkjson'
-- assigns globals:
local namedStars = fromlua(path'namedStars.lua':read())
local constellations = fromlua(path'constellations.lua':read())
-- writes files:
path'namedStars.json':write('namedStars = ' .. json.encode(namedStars, {indent=true}) ..';')
path'constellations.json':write('constellations = '..json.encode(constellations, {indent=true}) .. ';')
