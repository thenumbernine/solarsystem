-- if you have to run luajit separate of lua ...
require 'ext'
local json = require 'dkjson'
-- assigns globals:
local namedStars = fromlua(file['namedStars.lua'])
local constellations = fromlua(file['constellations.lua'])
-- writes files:
file['namedStars.json'] = 'namedStars = ' .. json.encode(namedStars, {indent=true}) ..';'
file['constellations.json'] = 'constellations = '..json.encode(constellations, {indent=true}) .. ';'
