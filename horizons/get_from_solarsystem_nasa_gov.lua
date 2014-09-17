--[[
based on http://solarsystem.nasa.gov/planets/profile.cfm?Object=Jup_Lysithea
this isn't horizons.  i am keeping the master copy in this folder.  i should move it elsewhere and list these as sources.
run this after the fact to fill in the missing blanks on mass and radius of objects we've already got position and velocity info for
--]]

require 'ext'
local json = require 'dkjson'
local socket = require 'socket'
local http = require 'socket.http'
require 'htmlparser.htmlparser'
require 'htmlparser.common'
require 'htmlparser.xpath'

local planetsFilename = 'full-horizons-extra-data.json'
local planetsData = io.readfile(planetsFilename)
-- strip off variable wrapper
local planetsLines = planetsData:trim():split('\n')
planetsLines[1] = assert(planetsLines[1]:match('horizonsExtraData = (.*)'))
print('last line',planetsLines[#planetsLines])
planetsLines[#planetsLines] = assert(planetsLines[#planetsLines]:match('(.*);'))
planetsData = planetsLines:concat('\n')
local planets = assert(json.decode(planetsData))

local function upperCase(name)
	-- separators?
	return name:sub(1,1):upper()..name:sub(2):lower()
end

local parentPrefixes = {
	Sun = '',
	Jupiter = 'Jup_',
	Saturn = 'Sat_',
	Uranus = 'Ura_',
	Neptune = 'Nep_',
	Pluto = 'Plu_',
}

local function getInfo(planet)
	local prefix = parentPrefixes[planet.parent]
	if not prefix then 
		print("Couldn't find prefix for parent "..tostring(planet.parent))
		return
	end
	local urlname = upperCase(prefix)..upperCase(planet.name)
	local url = [[http://solarsystem.nasa.gov/planets/profile.cfm?Object=]]..urlname
	local page = assert(http.request(url))
	local tree = htmlparser.new(page):parse()
	local spans = htmlparser.xpath(tree, '//span')
	print(#spans,'spans')
	spans = spans:map(function(span) return flattenText(span) end)
	for _,span in ipairs(spans) do print(span) end
	do
		local i = 1
		local function readValue(expectedUnits)
			i = i + 1
			local valueWords = spans[i]:gsub(',',''):split('%s+')
			local value = valueWords[1]
			local units = valueWords[#valueWords]
			if units ~= expectedUnits then 
				error("expected units "..expectedUnits.." but found "..units.." for field "..value)
			end
			local valueNumber = tonumber(value)
			if not valueNumber then
				error("failed to parse to number: "..value)
			end
			return valueNumber
		end
		--local volume
		while i <= #spans do
			local key = spans[i]
			if i == #spans then break end	-- already at the end
			if key:lower():find('radius') then
				planet.radius = readValue('km') * 1000
			elseif key:lower():find('mass') then
				planet.mass = readValue('kg')
			--elseif key:lower():find('volume') then
			--	volume = readValue('km3')
			end
			i = i + 1
		end
		--if not planet.mass and planet.radius and volume then
			--approx via sphere
			-- .. but then we need the density ... in the "More Facts" section ...
		--end
	end
	--io.writefile('tmp.html', page)
end

for _,planet in ipairs(planets) do
	planet.name = planet.name:trim()
	if planet.name:match('Barycenter') then	--skip barycenters
	else
		if not planet.mass or not planet.radius then
			local oldMass = planet.mass
			local oldRadius = planet.radius
			getInfo(planet)
			print(planet.parent, planet.name, 'old mass',oldMass,'old radius',oldRadius, 'mass', planet.mass, 'radius', planet.radius)
			if not planet.mass or not planet.radius then
				print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
				print('!!!! failed to find mass for planet '..planet.name..' !!!')
				print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
			end
		end
	end
end

local planetsData = 'horizonsExtraData = '..json.encode(planets, {indent=true})..';'
io.writefile(planetsFilename, planetsData)
