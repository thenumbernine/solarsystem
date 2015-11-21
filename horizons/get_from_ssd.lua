require 'ext'
local json = require 'dkjson'
local socket = require 'socket'
local http = require 'socket.http'
local htmlparser = require 'htmlparser.htmlparser'
require 'htmlparser.common'
require 'htmlparser.xpath'

-- [[
local planetsFilename = 'static-vars.json'
local planetsData = file[planetsFilename]
-- strip off variable wrapper
local planetsLines = planetsData:trim():split('\n')
planetsLines[1] = assert(planetsLines[1]:match('horizonsStaticData = (.*)'))
print('last line',planetsLines[#planetsLines])
planetsLines[#planetsLines] = assert(planetsLines[#planetsLines]:match('(.*);'))
planetsData = planetsLines:concat('\n')
local planets = assert(json.decode(planetsData))
--]]

local G = 6.67259e-20
local parsecInKM = 3.08567758e+13
local auInKM = 149597871

-- http://idahoptv.org/ntti/nttilessons/lessons2000/lau1.html
local distFromEarthToMoonInKM = 384400
local distFromSunInKM = {
	mercury = 57910000,
	venus = 108200000,
	earth = 149600000,
	mars = 227940000,
	jupiter = 778330000,
	saturn = 1424600000,
	uranus = 2873550000,
	neptune = 4501000000,
	pluto = 5945900000,
}

local function parseFloat(s)
	s = s:trim()
	if s == '?' then return 0 end
	s = s:match('[%d%.e%+%-]*')
	return tonumber(s)
end

local vars = {}

local cacheFileName = 'ssd_sat_phys_par.html'
local page = file[cacheFileName]
if not page then
	local url = 'http://ssd.jpl.nasa.gov/?sat_phys_par'
	page = assert(http.request(url))
	file[cacheFileName] = page
end
page = page:gsub('\r\n', '\n')
-- fix nasa's crappy html
page = page:gsub('(D>\n)<TR', '%1</TR><TR')
local tree = htmlparser.parse(page)
local ps = htmlparser.xpath(tree, '//p')
for _,p in ipairs(ps) do
	if p.child then
		local tabl = findchild(p, 'table')
		local font = findchild(p, 'font')
		if tabl and tabl.child and font then
			local parentPlanetName = flattenText(font):match('%w+'):lower()
			print(parentPlanetName)
			parentPlanetName = ({
				martian = 'mars',
				jovian = 'jupiter',
				saturnian = 'saturn',
				uranian = 'uranus',
				neptunian = 'neptune',
			})[parentPlanetName] or parentPlanetName
			local trs = findchilds(tabl, 'tr')
			for i=2,#trs do	-- trs[1] is the header column ...
				local tr = trs[i]
				local tds = findchilds(tr, 'td')
				for i=#tds,1,-1 do
					local colspan = tonumber(findattr(tds[i],'colspan') or 1)
					for i=1,colspan-1 do
						table.insert(tds,i+1,{})	-- fill empty columns so common column cells are at matching offsets 
					end
				end
				assert(#tds == 10, "found "..#tds.." cols")
		
				local name = flattenText(tds[1])
				local GM = parseFloat(flattenText(tds[2]))
				local mass = GM / G
				local radius = 1e+3 * parseFloat(flattenText(tds[4])) -- radius in m
				local density = parseFloat(flattenText(tds[6])) * 1e+12 -- density in kg/m^3 = kg/g (cm/m)^3 g/cm^3 = g/cm^3 * 1e+3
				local meanOppositionMagnitude = parseFloat(flattenText(tds[7]))	-- V0 or R (if has a R suffix)
				local albedo = parseFloat(flattenText(tds[9]))
				
				if mass == 0 and density ~= 0 then
					mass = density * (4/3 * math.pi * radius^3)
					print('calculating mass of ',name,' to be',mass)
				end
				
				-- "mean opposition magnitude" is the magnitude when the planet and earth are in opposition and at mean distances from the sun 
				-- so this is the apparent magnitude with distance of this avg distance from sun plus earth's avg distance from sun 
				
				local distInKM
				if parentPlanetName == 'earth' then 	-- because moon orbits earth where observations are from, use this distance:
					distInKM = distFromEarthToMoonInKM
				else
					distInKM = distFromSunInKM[parentPlanetName]
				end
				local absMag = meanOppositionMagnitude - 5 * (math.log(distInKM / parsecInKM, 10) - 1)

				print(name, mass, radius, density, absMag, albedo)
				name = name:trim():lower()
				name = ({
					herse = 'herse (2003j17)',
					prospero = 'prospero (1999u3)',
					setebos = 'setebos (1999u1)',
					stephano = 'stephano (1999u2)',
					['s/2003 j15']	= '2003j15',
					['s/2003 j18']	= '2003j18',
					['s/2003 j19']	= '2003j19',
					['s/2011 j1']	= '2011j1',
					['s/2006 s3']	= '2006s3',
					['s/2000 j11']	= 'dia (2000j11)',
					['s/2003 j12']	= '2003j12',
					['s/2003 j3']	= '2003j3',
					['s/2010 j1']	= '2010j1',
					['s/2003 j10']	= '2003j10',
					['s/2010 j2']	= '2010j2',
					['s/2011 j2']	= '2011j2',
					['s/2004 n1']	= 's 2004 n1',
					['s/2003 j5']	= '2003j5',
					['s/2003 j16']	= '2003j16',
					['s/2006 s1']	= '2006s01',
					['s/2007 s2']	= '2007s2',
					['s/2007 s3']	= '2007s3',
					['s/2004 s13']	= '2004s13',
					['s/2004 s12']	= '2004s12',
					['s/2003 j23']	= '2003j23',
					['s/2003 j9']	= '2003j9',
					['s/2003 j4']	= '2003j4',
					['s/2004 s17']	= '2004s17',
					['s/2004 s7']	= '2004s7',
					['s/2003 j2']	= '2003j2',
				})[name] or name
				
				vars[name] = {
					mass = mass,
					radius = radius,
					density = density,
					magnitude = absMag,
					albedo = albedo,
				}
			end
		end
	end
end

-- [[
local cacheFileName = 'ssd_planet_phys_par.html'
local page = file[cacheFileName]
if not page then
	local url = 'http://ssd.jpl.nasa.gov/?planet_phys_par'
	page = assert(http.request(url))
	file[cacheFileName] = page
end
page = page:gsub('\r\n', '\n')
-- fix nasa's crappy html
page = page:gsub('<tr>', '</tr>')

local tree = htmlparser.parse(page)
local tds = htmlparser.xpath(tree, '//td')
print('#tds',#tds)
for _,td in ipairs(tds) do
	if td.child then
		local tabls = findchilds(td, 'table')
		for _,tabl in ipairs(tabls) do
			local trs = findchilds(tabl, 'tr')
			if trs[1] then
				local td = findchild(trs[1], 'td')
				if td and flattenText(td) == 'Planet' then
					for i=3,#trs do	-- trs[1] is the header column ...
						local tr = trs[i]
						local tds = findchilds(tr, 'td')
						local function stupidNasaHTML(s)
							return flattenText(s):match('(%S+)%.')
						end
						local name = stupidNasaHTML(tds[1]) 
						local equatorialRadius = 1e+3 * parseFloat(stupidNasaHTML(tds[2]))
						local meanRadius = 1e+3 * parseFloat(stupidNasaHTML(tds[3]))
						local mass = 1e+24 * parseFloat(stupidNasaHTML(tds[4]))
						local bulkDensity = 1e+3 * parseFloat(stupidNasaHTML(tds[5]))
						local rotationPeriod = parseFloat(stupidNasaHTML(tds[6]))
						local orbitPeriod = parseFloat(stupidNasaHTML(tds[7]))
						local magnitude = parseFloat(stupidNasaHTML(tds[8]))
						local albedo = parseFloat(stupidNasaHTML(tds[9]))
						local equatorialGravity = parseFloat(stupidNasaHTML(tds[10]))
						local escapeVelocity = parseFloat(stupidNasaHTML(tds[11]))
						
						name = name:lower()
					
						-- convert from V(1,0) magnitude to abs magnitude
						-- magnitude is V(1,0) which I've found in one source to be the magnitude when planet is opposite the sun of earth, with distance measured in AU
						local distToObserverInKM = distFromSunInKM[name] + distFromSunInKM.earth
						magnitude = magnitude - 5 * (math.log(distToObserverInKM / parsecInKM, 10) - 1)
						
						print(name, equatorialRadius, meanRadius, mass, bulkDensity, rotationPeriod, orbitPeriod, magnitude, albedo, equatorialGravity, escapeVelocity) 
						
						name = ({
							pluto = '134340 pluto',
						})[name] or name
						
						vars[name] = {
							equatorialRadius = equatorialRadius,
							radius = meanRadius,
							mass = mass,
							density = bulkDensity,
							magnitude = magnitude,
						}
					end
				end
			end
		end
	end
end
--]]

-- [[
for _,planet in ipairs(planets) do
	planet.name = planet.name:trim()
	if planet.name:match('Barycenter') then
	else
		local v = vars[planet.name:lower()]
		if v then
			local fields = {'mass', 'radius', 'density', 'magnitude', 'albedo'}
			for _,field in ipairs(fields) do
				print('replacing',planet.name,field,'from',planet[field],'to',v[field])
				planet[field] = v[field]
			end
			v.read = true
		end
	end
end
for name,var in pairs(vars) do
	if not var.read then print('didnt read from ',name) end
end
--]]

-- [[
local planetsData = 'horizonsStaticData = '..json.encode(planets, {indent=true})..';'
file[planetsFilename] = planetsData
--]]
