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

local vars = {}
for _,p in ipairs(ps) do
	if p.child then
		local tabl = findchild(p, 'table')
		local font = findchild(p, 'font')
		if tabl and tabl.child and font then
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
		
				local function parseFloat(s)
					s = s:trim()
					if s == '?' then return 0 end
					s = s:match('[%d%.e%+%-]*')
					return tonumber(s)
				end
				local G = 6.67259e-20
				local name = flattenText(tds[1])
				local GM = parseFloat(flattenText(tds[2]))
				local mass = GM / G
				local radius = parseFloat(flattenText(tds[4]))
				local density = parseFloat(flattenText(tds[6]))
				local magnitude = parseFloat(flattenText(tds[7]))	-- V0 or R (if has a R suffix)
				if mass == 0 and density ~= 0 then
					-- density in kg/km^3 = kg/g (cm/km)^3 g/cm^3 = g/cm^3 * 1e+12
					mass = 1e+12 * density * (4/3 * math.pi * radius^3)
					print('calculating mass of ',name,' to be',mass)
				end
				-- from here on out, radius in meters
				radius = radius * 1000
				local albedo = parseFloat(flattenText(tds[9]))
				print(name, mass, radius, density, magnitude, albedo)
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
					magnitude = magnitude,
					albedo = albedo,
				}
			end
		end
	end
end

-- [[
for _,planet in ipairs(planets) do
	planet.name = planet.name:trim()
	if planet.name:match('Barycenter') then
	else
		local v = vars[planet.name:lower()]
		if v then
			print('replacing',planet.name,'mass from',planet.mass,'to',v.mass)
			planet.mass = v.mass
			print('replacing',planet.name,'radius from',planet.radius,'to',v.radius)
			planet.radius = v.radius
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
