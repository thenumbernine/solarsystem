#!/usr/bin/env lua
require 'ext'
-- this just removes .vars from each planet
-- it writes results to stdout
local json = require 'myjson'
local staticIntro = 'horizonsStaticData = '
local filename = ...
local static = json.decode(assert(file(filename):read()):match(staticIntro..'(.*)'))
for _,p in ipairs(static) do
	p.vars = nil
end
print(staticIntro..json.encode(static, {indent=true})..';')
