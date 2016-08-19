#!/usr/bin/env luajit
require 'ext'
local json = require 'dkjson'
local htmlparser = require 'htmlparser'
local htmlparser_common = require 'htmlparser.common'
local findattr = htmlparser_common.findattr
local findtags = htmlparser_common.findtags
local findchilds = htmlparser_common.findchilds
require 'htmlparser.xpath'
htmlparser.htmlnonclosing = {}
local tree = htmlparser.parse(file['systems.xml'])

local function getText(child)
	if not (#child.child == 1) then error("tried to parse text of a child with multiple children: " ..tolua(child)) end
	assert(type(child.child[1]) == 'string')
	return child.child[1]
end

local function getNumber(child)
	local upperlimit = tonumber(findattr(child, 'upperlimit') or nil)
	local lowerlimit = tonumber(findattr(child, 'lowerlimit') or nil)
	local value = child.child and tonumber(getText(child))
	if value then return value end
	if upperlimit and lowerlimit then return .5 * (upperlimit + lowerlimit) end
	return upperlimit or lowerlimit or error("failed to find number value")
end

local auInM = 149597870700
local parsecsToMeters = 3.08567758e+16
local jupiterMassInKg = 1.89813e+27
local jupiterRadiusInM = 6.9911e+7
local sunMassInKg = 1.9891e+30
local sunRadiusInM = 6.960e+8

--[[
results = {
	{	-- system
		name = string,
		{ -- body: planet/star/binary barycenter
			name = string
			parent = <- some way to identify the system and body,
			...
		}, ...
	}, ...
}
--]]

local processStar
local processPlanet
local processBinary
local processSystem

local function setField(t, k, v)
	if t[k] then error("tried to overwrite body field "..k.." prior v "..t[k].." with v "..v) end
	t[k] = v
end

-- accumulate multiple names
local nameSeparator = ' / '
local function setName(obj, name)
	if not obj.names then obj.names = {} end
	table.insert(obj.names, name)
	if obj.name then 
		obj.name = obj.name .. nameSeparator .. name
	else
		obj.name = name
	end
end

--[[
local function getSeparation(node)
	-- separation looks to be equivalent of the radius of the orbit
	-- this should be a fallback for when semimajoraxis is not provided
	if findattr(node, 'unit') == 'arcsec' then	-- then process separtion in arcseconds?
		--error('how do we handle arcsecond separation?')
	elseif findattr(node, 'unit') == 'AU' then
		return getNumber(node) * auInM
	else
		-- sometimes units aren't given ... in that case what do we do?
		return getNumber(node) * auInM
		--error('unknown separation type '..tolua(node))
	end
end
--]]

local function setSeparation(obj, node)
	-- separation looks to be equivalent of the radius of the orbit
	-- this should be a fallback for when semimajoraxis is not provided
	local unit = findattr(node, 'unit')
	if not unit then
		io.stderr:write("ignoring separation with no specified units: ",tolua(node),'\n')
		return
	elseif unit == 'arcsec' then	-- then process separtion in arcseconds?
	elseif unit == 'AU' then
		setField(obj, 'separation', getNumber(node) * auInM)
	else
		error('unknown separation type '..tolua(node.attrs))
	end
end

function processPlanet(node, resultSystem)
	local resultBody = {type='planet'}
	table.insert(resultSystem.bodies, resultBody)
	for _,child in ipairs(node.child) do
		local tag = child.tag:lower()
		if tag == 'name' then
			setName(resultBody, getText(child))
		elseif tag == 'semimajoraxis' then
			setField(resultBody, 'semiMajorAxis', getNumber(child) * auInM)
		elseif tag == 'separation' then
			setSeparation(resultBody, child)
		elseif tag == 'eccentricity' then
			setField(resultBody, 'eccentricity', getNumber(child)) 	-- I'm guessing eccentricity of orbit, and not of spheroid
		elseif tag == 'periastron' then
			setField(resultBody, 'longitudeOfPeriapsis', math.rad(getNumber(child))) 	-- degrees -> radians
		elseif tag == 'periastrontime' then
			setField(resultBody, 'timeOfPerihelionPassage', math.rad(getNumber(child)))
		elseif tag == 'longitude' then
			setField(resultBody, 'meanLongitude', math.rad(getNumber(child)))		-- degrees -> radians
		elseif tag == 'ascendingnode' then
			setField(resultBody, 'longitudeOfAscendingNode', math.rad(getNumber(child)))	-- degrees -> radians
		elseif tag == 'inclination' then
			setField(resultBody, 'inclination', math.rad(getNumber(child)))		-- degrees -> radians
		elseif tag == 'period' then
			setField(resultBody, 'orbitalPeriod', getNumber(child))	-- days
		elseif tag == 'transittime' then	-- who cares?
		elseif tag == 'mass' then
			setField(resultBody, 'mass', getNumber(child) * jupiterMassInKg)
		elseif tag == 'radius' then
			setField(resultBody, 'radius', getNumber(child) * jupiterRadiusInM)
		elseif tag == 'temperature' then
			setField(resultBody, 'temperature', getNumber(child))
		elseif tag == 'age' then
			setField(resultBody, 'age', getNumber(child)) 		-- in Gyr
		elseif tag == 'spectraltype' then
			setField(resultBody, 'spectralType', getText(child))
		elseif tag == 'magb' then
			setField(resultBody, 'BMagnitude', getNumber(child))
		elseif tag == 'magv' then
			setField(resultBody, 'visualMagnitude', getNumber(child))
		elseif tag == 'magr' then
			setField(resultBody, 'RMagnitude', getNumber(child))
		elseif tag == 'magi' then
			setField(resultBody, 'IMagnitude', getNumber(child))
		elseif tag == 'magj' then
			setField(resultBody, 'JMagnitude', getNumber(child))
		elseif tag == 'magh' then
			setField(resultBody, 'HMagnitude', getNumber(child))
		elseif tag == 'magk' then
			setField(resultBody, 'KMagnitude', getNumber(child))
		elseif tag == 'meananomaly' then
			setField(resultBody, 'meanAnomaly', getNumber(child))
		elseif tag == 'discoverymethod' then	-- who cares?
		elseif tag == 'istransiting' then
		elseif tag == 'description' then
		elseif tag == 'discoveryyear' then
		elseif tag == 'lastupdate' then
		elseif tag == 'spinorbitalignment' then
			setField(resultBody, 'spinOrbitAlignment', math.rad(getNumber(child)))	-- degrees -> radians
		elseif tag == 'list' then
			resultBody.lists = resultBody.lists or {}
			table.insert(resultBody.lists, getText(child))
		elseif tag == 'positionangle' then
			setField(resultBody, 'positionAngle', math.rad(getNumber(child)))	-- no documentation, guessing degrees -> radians
		elseif tag == 'image' then
		elseif tag == 'imagedescription' then
		elseif tag == 'maximumrvtime' then
		elseif tag == 'impactparameter' then
		elseif tag == 'metallicity' then
		elseif tag == 'new' then
		else
			error('unknown parameter '..tag..' for planet')
		end
	end
	return resultBody
end

function processStar(node, resultSystem)
	local resultBody = {type='star'}
	table.insert(resultSystem.bodies, resultBody)
	for _,child in ipairs(node.child) do
		local tag = child.tag:lower()
		if tag == 'planet' then
			local body = processPlanet(child, resultSystem)
			body.parent = resultBody
		elseif tag == 'name' then
			setName(resultBody, getText(child))
		elseif tag == 'mass' then
			setField(resultBody, 'mass', getNumber(child) * sunMassInKg)
		elseif tag == 'radius' then
			setField(resultBody, 'radius', getNumber(child) * sunRadiusInM)
		elseif tag == 'temperature' then
			setField(resultBody, 'temperature', getNumber(child))
		elseif tag == 'age' then
			setField(resultBody, 'age', getNumber(child)) 		-- in Gyr
		elseif tag == 'metallicity' then
			setField(resultBody, 'metallicity', getText(child))
		elseif tag == 'spectraltype' then
			setField(resultBody, 'spectralType', getText(child))
		elseif tag == 'magb' then
			setField(resultBody, 'BMagnitude', getNumber(child))
		elseif tag == 'magv' then
			setField(resultBody, 'visualMagnitude', getNumber(child))
		elseif tag == 'magr' then
			setField(resultBody, 'RMagnitude', getNumber(child))
		elseif tag == 'magi' then
			setField(resultBody, 'IMagnitude', getNumber(child))
		elseif tag == 'magj' then
			setField(resultBody, 'JMagnitude', getNumber(child))
		elseif tag == 'magh' then
			setField(resultBody, 'HMagnitude', getNumber(child))
		elseif tag == 'magk' then
			setField(resultBody, 'KMagnitude', getNumber(child))
		-- undocumented
		elseif tag == 'positionangle' then
			setField(resultBody, 'positionAngle', math.rad(getNumber(child)))	-- no documentation, guessing degrees -> radians
		elseif tag == 'magu' then
		else
			error('unknown parameter '..tag..' for star')
		end
	end
	return resultBody
end

function processBinary(node, resultSystem)
	local resultBody = {type='barycenter'}
	table.insert(resultSystem.bodies, resultBody)
	for _,child in ipairs(node.child) do
		local tag = child.tag:lower()
		if tag == 'planet' then
			local body = processPlanet(child, resultSystem)
			body.parent = resultBody
		elseif tag == 'star' then
			local body = processStar(child, resultSystem)
			body.parent = resultBody
		elseif tag == 'binary' then
			local body = processBinary(child, resultSystem)
			body.parent = resultBody
		elseif tag == 'name' then
			setName(resultBody, getText(child))
		elseif tag == 'semimajoraxis' then
			setField(resultBody, 'semiMajorAxis', getNumber(child) * auInM)
		elseif tag == 'separation' then
			setSeparation(resultBody, child)
		elseif tag == 'eccentricity' then
			setField(resultBody, 'eccentricity', getNumber(child)) 	-- I'm guessing eccentricity of orbit, and not of spheroid
		elseif tag == 'periastron' then
			setField(resultBody, 'longitudeOfPeriapsis', math.rad(getNumber(child))) 	-- degrees -> radians
		elseif tag == 'periastrontime' then
			setField(resultBody, 'timeOfPerihelionPassage', math.rad(getNumber(child)))
		elseif tag == 'longitude' then
			setField(resultBody, 'meanLongitude', math.rad(getNumber(child)))		-- degrees -> radians
		elseif tag == 'ascendingnode' then
			setField(resultBody, 'longitudeOfAscendingNode', math.rad(getNumber(child)))	-- degrees -> radians
		elseif tag == 'inclination' then
			setField(resultBody, 'inclination', math.rad(getNumber(child)))		-- degrees -> radians
		elseif tag == 'period' then
			setField(resultBody, 'orbitalPeriod', getNumber(child))	-- days
		elseif tag == 'transittime' then	-- who cares?
		elseif tag == 'magb' then
			setField(resultBody, 'BMagnitude', getNumber(child))
		elseif tag == 'magv' then
			setField(resultBody, 'visualMagnitude', getNumber(child))
		elseif tag == 'magr' then
			setField(resultBody, 'RMagnitude', getNumber(child))
		elseif tag == 'magi' then
			setField(resultBody, 'IMagnitude', getNumber(child))
		elseif tag == 'magj' then
			setField(resultBody, 'JMagnitude', getNumber(child))
		elseif tag == 'magh' then
			setField(resultBody, 'HMagnitude', getNumber(child))
		elseif tag == 'magk' then
			setField(resultBody, 'KMagnitude', getNumber(child))
		--undocumented:
		elseif tag == 'positionangle' then
			setField(resultBody, 'positionAngle', math.rad(getNumber(child)))	-- no documentation, guessing degrees -> radians
		elseif tag == 'temperature' then
			setField(resultBody, 'temperature', getNumber(child))
		else
			error('unknown parameter '..tag..' for binary')
		end
	end
	return resultBody
end

function processSystem(node)
	local resultSystem = {
		bodies = {}
	}
	for _,child in ipairs(node.child) do
		local tag = child.tag:lower()
		if tag == 'star' then
			local body = processStar(child, resultSystem)
		elseif tag == 'planet' then
			local body = processPlanet(child, resultSystem)
		elseif tag == 'binary' then
			local body = processBinary(child, resultSystem)
		elseif tag == 'declination' then
			local parts = getText(child):split(' ')
			assert(#parts == 3)
			local deg = assert(tonumber(parts[1]))
			local min = assert(tonumber(parts[2]))
			local sec = assert(tonumber(parts[3]))
			setField(resultSystem, 'declination', math.rad(deg + 1/60 * (min + 1/60 * sec)))	--degree minute second -> radians
		elseif tag == 'rightascension' then
			local parts = getText(child):split(' ')
			assert(#parts == 3)
			local hours = assert(tonumber(parts[1]))
			local min = assert(tonumber(parts[2]))
			local sec = assert(tonumber(parts[3]))
			setField(resultSystem, 'rightAscension', math.rad(360/24 * (hours + 1/60 * (min + 1/60 * sec))))	-- hour minute second-> radians
		elseif tag == 'distance' then
			setField(resultSystem, 'distance', getNumber(child) * parsecsToMeters)
		elseif tag == 'name' then
			setName(resultSystem, getText(child))
		elseif tag == 'epoch' then
			setField(resultSystem, 'epoch', getText(child))
		--undocumented
		elseif tag == 'videolink' then
		elseif tag == 'positionangle' then
			setField(resultSystem, 'positionAngle', math.rad(getNumber(child)))	-- no documentation, guessing degrees -> radians
		elseif tag == 'spectraltype' then
			setField(resultSystem, 'spectralType', getText(child))
		elseif tag == 'magj' then
		elseif tag == 'magh' then
		elseif tag == 'magk' then
		else
			error('unknown parameter '..tag..' for system')
		end
	end
	return resultSystem
end

local resultSystems = {}
for _,systems in ipairs(findtags(tree, 'systems')) do
	for _,node in ipairs(findchilds(systems, 'system')) do
		local resultSystem = processSystem(node)
		table.insert(resultSystems, resultSystem)
	end
end

do
	local found = 0
	for i=#resultSystems,1,-1 do
		if resultSystems[i].name == 'Sun' then
			found = found + 1 
			table.remove(resultSystems, i)
		end
	end
	if found ~= 1 then error("failed to find and remove our own solar system") end
end

local greekLettersForAbbrevs = {
	alf = 'alpha',
	bet = 'beta',
	gam = 'gamma',
	delta = 'delta',
	eps = 'epsilon',
	zeta = 'zeta',
	eta = 'eta',
	theta = 'theta',
	iota = 'iota',
	kap = 'kappa',
	lambda = 'lambda',
	mu = 'mu',
	nu = 'nu',
	xi = 'xi',
	ksi = 'ksi',	-- Cyrillic variation of Xi
	omi = 'omicron',
	pi = 'pi',
	rho = 'rho',
	sigma = 'sigma',
	tau = 'tau',	-- notice 'tau' is also an abbreviation for 'Taurus'
	ups = 'upsilon',
	phi = 'phi',
	chi = 'chi',
	psi = 'psi',
	ome = 'omega',
}
local greekAbbrevsForLetters = setmetatable(table.map(greekLettersForAbbrevs, function(v,k)
	return k,v
end), nil)

-- TODO map names from abbreviations to full names?
-- got translation from here: http://stars.astro.illinois.edu/sow/sowlist.html
-- replace 'abbrev' with 'star', group by 'constellation'
local nameInfo = table{
	{abbrev='and', star='Andromedae', constellation='Andromeda'},
	{abbrev='ant', star='Antliae', constellation='Antlia'},
	{abbrev='aps', star='Apodis', constellation='Apus'},
	{abbrev='aql', star='Aquilae', constellation='Aquila'},
	{abbrev='aqr', star='Aquarii', constellation='Aquarius'},
	{abbrev='ara', star='Arae', constellation='Ara'},
	{abbrev='ari', star='Arietis', constellation='Aries'},
	{abbrev='aur', star='Aurigae', constellation='Auriga'},
	{abbrev='boo', star='Bootis', constellation='Bootes'},
	{abbrev='boo', star='Boötis', constellation='Bootes'},	-- for the accents on the spelling ...
	{abbrev='cae', star='Caeli', constellation='Caelum'},
	{abbrev='cam', star='Camelopardalis', constellation='Camelopardalis'},
	{abbrev='cap', star='Capricorni', constellation='Capricornus'},
	{abbrev='car', star='Carinae', constellation='Carina'},
	{abbrev='cas', star='Cassiopeiae', constellation='Cassiopeia'},
	{abbrev='cen', star='Centauri', constellation='Centaurus'},
	{abbrev='cep', star='Cephei', constellation='Cepheus'},
	{abbrev='cet', star='Ceti', constellation='Cetus'},
	{abbrev='cha', star='Chamaeleontis', constellation='Chameleon'},
	{abbrev='cir', star='Circini', constellation='Circinus'},
	{abbrev='cma', star='Canis Majoris', constellation='Canis Major'},
	{abbrev='cmi', star='Canis Minoris', constellation='Canis Minor'},
	{abbrev='cnc', star='Cancri', constellation='Cancer'},
	{abbrev='col', star='Columbae', constellation='Columba'},
	{abbrev='com', star='Comae Berenices', constellation='Coma Berenices'},
	{abbrev='cra', star='Coronae Australis', constellation='Corona Australis'},
	{abbrev='crb', star='Coronae Borealis', constellation='Corona Borealis'},
	{abbrev='crb', star='Corona Borealis', constellation='Corona Borealis'},	-- another alternative spelling for the star name ...
	{abbrev='crv', star='Corvi', constellation='Corvus'},
	{abbrev='crt', star='Crateris', constellation='Crater'},
	{abbrev='cru', star='Crucis', constellation='Crux'},
	{abbrev='cyg', star='Cygni', constellation='Cygnus'},
	{abbrev='cvn', star='Canum Venaticorum', constellation='Canes Venatici'},
	{abbrev='del', star='Delphini', constellation='Delphinus'},
	{abbrev='dor', star='Doradus', constellation='Dorado'},
	{abbrev='dra', star='Draconis', constellation='Draco'},
	{abbrev='equ', star='Equulei', constellation='Equuleus'},
	{abbrev='eri', star='Eridani', constellation='Eridanus'},
	{abbrev='for', star='Fornacis', constellation='Fornax'},
	{abbrev='gem', star='Geminorum', constellation='Gemini'},
	{abbrev='gru', star='Gruis', constellation='Grus'},
	{abbrev='her', star='Herculis', constellation='Hercules'},
	{abbrev='hor', star='Horologii', constellation='Horologium'},
	{abbrev='hya', star='Hydrae', constellation='Hydra'},
	{abbrev='hyi', star='Hydri', constellation='Hydrus'},
	{abbrev='ind', star='Indi', constellation='Indus'},
	{abbrev='lac', star='Lacertae', constellation='Lacerta'},
	{abbrev='leo', star='Leonis', constellation='Leo'},
	{abbrev='lep', star='Leporis', constellation='Lepus'},
	{abbrev='lib', star='Librae', constellation='Libra'},
	{abbrev='lmi', star='Leonis Minoris', constellation='Leo Minor'},
	{abbrev='lup', star='Lupi', constellation='Lupus'},
	{abbrev='lyn', star='Lyncis', constellation='Lynx'},
	{abbrev='lyr', star='Lyrae', constellation='Lyra'},
	{abbrev='men', star='Mensae', constellation='Mensa'},
	{abbrev='mic', star='Microscopii', constellation='Microscopium'},
	{abbrev='mon', star='Monocerotis', constellation='Monoceros'},
	{abbrev='mus', star='Muscae', constellation='Musca'},
	{abbrev='nor', star='Normae', constellation='Norma'},
	{abbrev='oct', star='Octantis', constellation='Octans'},
	{abbrev='oph', star='Ophiuchi', constellation='Ophiuchus'},
	{abbrev='ori', star='Orionis', constellation='Orion'},
	{abbrev='pav', star='Pavonis', constellation='Pavo'},
	{abbrev='peg', star='Pegasi', constellation='Pegasus'},
	{abbrev='per', star='Persei', constellation='Perseus'},
	{abbrev='phe', star='Phoenicis', constellation='Phoenix'},
	{abbrev='pic', star='Pictoris',	constellation='Pictor'},
	{abbrev='psa', star='Piscis Austrini', constellation='Piscis Austrinus'},
	{abbrev='psc', star='Piscium', constellation='Pisces'},
	{abbrev='pup', star='Puppis', constellation='Puppis'},
	{abbrev='pyx', star='Pyxidis', constellation='Pyxis'},
	{abbrev='ret', star='Reticuli', constellation='Reticulum'},
	{abbrev='scl', star='Sculptoris', constellation='Sculptor'},
	{abbrev='sco', star='Scorpii', constellation='Scorpius'},
	{abbrev='sct', star='Scuti', constellation='Scutum'},
	{abbrev='ser', star='Serpentis', constellation='Serpens'},
	{abbrev='sex', star='Sextantis', constellation='Sextans'},
	{abbrev='sge', star='Sagittae', constellation='Sagitta'},
	{abbrev='sgr', star='Sagittarii', constellation='Sagittarius'},
	{abbrev='tau', star='Tauri', constellation='Taurus'},
	{abbrev='tel', star='Telescopii', constellation='Telescopium'},
	{abbrev='tra', star='Trianguli Australis', constellation='Triangulum Australe'},
	{abbrev='tri', star='Trianguli', constellation='Triangulum'},
	{abbrev='tuc', star='Tucanae', constellation='Tucana'},
	{abbrev='uma', star='Ursae Majoris', constellation='Ursa Major'},
	{abbrev='umi', star='Ursae Minoris', constellation='Ursa Minor'},
	{abbrev='vel', star='Velorum', constellation='Vela'},
	{abbrev='vir', star='Virginis', constellation='Virgo'},
	{abbrev='vol', star='Volantis', constellation='Volans'},
	{abbrev='vul', star='Vulpeculae', constellation='Vulpecula'},
	-- extras:
	-- ['PH-(%d)'] = "Planet Hunters %1"
	-- ['TrES-(%d)'] = "Trans-Atlantic Exoplanet Survey %1"
	-- "96 G. Piscium" is also "HD 4628"
}
local nameInfoForAbbrev = setmetatable(nameInfo:map(function(info)
	return info, info.abbrev
end), nil)
local nameInfoForStar = setmetatable(nameInfo:map(function(info)
	return info, info.star:lower()
end), nil)
local nameInfoForConstellation = setmetatable(nameInfo:map(function(info)
	return info, info.constellation:lower()
end), nil)

--[[
- find prefix associated with constellation / survey
- ... how to store so that the names can be sorted / recalled easily?
	names = {
		[surveyName] = survey's name for star / planet
	}
- return value+key of fixed name (abbreviations removed) + associated survey
--]]
local function parseName(name)
	local index = 1
	local lasttoken = ''
	local function matches(pattern)
		local result = name:sub(index):match(pattern)
		if not result then return end
		lasttoken = result
		return result
	end
	local function matchesKey(t)
		for k,v in pairs(t) do
			-- if the key matches
			if name:sub(index, index+#k-1):lower() == k:lower() then
				--print('matched '..('%q'):format(name:sub(index, index+#k-1)))
				--print('#name ('..#name..') vs index ('..index..') + #k ('..#k..')')
				-- and the end of the name
				if (#name == index-1 + #k)
				-- .. or it's either before a space or a parenthesis (this only applies to matchesKey() calls on 2nd tokens
				or ({
					[' ']=1,
					['(']=1,
				})[name:sub(index+#k,index+#k)]
				then
					--print('... to either end or space') 
					lasttoken = k
					--print('got token '..k)
					return k
				end
			end
		end
	end
	local function nexttoken()
		index = index + #lasttoken
	end

	-- prefix for name of survey
	if false then
	elseif matches('^1RXS') then	-- only "1RXS1609" which is in constellation Upper Scorpius
		return name, 'Scorpius'
	elseif matches('^1SWASP') then	-- only "1SWASP J022837.22-070338.4" which is a synonym for "WASP-77"
		return name, '1SWASP'
	elseif matches('^2MASS ') then	-- "2 Micron All-Sky Survey"
		return name, '2MASS'
	elseif matches('^BPM ') then	-- ???
		return name, 'BPM'
	elseif matches('^GJ ') then		-- "Gliese Catalogue of Nearby Stars" numbers 1000-1294 and 2001-2159, GJ stands for "W. Gliese and H. Jahreiß"
		return name, 'GJ'
	elseif matches('^Gj ') then		-- same as above?  why the lower case? 
		return name, 'GJ'
	elseif matches('^Gliese ') then	-- "Gliese Catalogue of Nearby Stars"
		return name, 'Gliese'
	elseif matches('^Gl ') then		-- same as above?  why the lower case? 
		return name, 'Gliese'
	elseif matches('^GSC ') then	-- "Guide Star Catalog"
		return name, 'GSC'
	elseif matches('^HD') then		-- "Henry Draper"
		return name, 'HD'
	elseif matches('^HIP') then	-- "Hipparcos Catalogue"
		return name, 'HIP'
	elseif matches('^HR') then		-- "Harvard Research"
		return name, 'HR'
	elseif matches('^KIC ')  then	-- "Kepler Input Catalog" 
		return name, 'KIC'
	elseif matches('^KIC-') then	-- same as above.  some have dashes, some have spaces ...
		return name, 'KIC'
	elseif matches('^NGC ')  then	-- "New General Catalogue" / "New General Catalogue of Nebulae and Clusters of Stars"
		return name, 'NGC'
	elseif matches('^ROXs ') then	-- only "ROXs 12" and "ROXs 42 B" ... maybe "Rho Ophiuchi cloud complex"
		return name, 'ROXs'
	elseif matches('^SR ') then		-- only "SR 12 AB" which is ... where? 
		return name, 'SR'
	elseif matches('^TYC ') then	-- "Tycho Catalogue"
		return name, 'TYC'
	elseif matches('^WD ') then		-- only "WD 0806-661" is maybe for "White Dwarf"
		return name, 'WD'
	elseif matches('^2M') then		-- only shows up as a shorthand for "2MASS"
		return name, '2M'
	elseif matches('^BD') then		-- "Bonner Durchmusterung"
		return name, 'BD'
	elseif matches('^CFBDS') then	-- "CF" could be for the "Canada-France-Hawaii Telescope", "BDS" for "Brown Dwarf Star", (optionally) "IR" for "using infra-red camera NIRC2"
		return name, 'CFBDS'
	elseif matches('^CoRoT-') then	-- "Convection Rotation and Planetary Transits"
		return name, 'CoRoT'
	elseif matches('^HAT-') then	-- "Hungarian Automated Telescope Network", "HAT-P-#" for north, "HATS-#" for south
		return name, 'HAT'
	elseif matches('^HATS-') then	-- "Hungarian Automated Telescope Network", "HAT-P-#" for north, "HATS-#" for south
		return name, 'HAT'
	elseif matches('^HIC') then		-- only "HIC 78407" aka "BD+15 2940".  Not sure what for.
		return name, 'HIC'
	elseif matches('^KELT-') then	-- "Kilodegree Extremely Little Telescope"
		return name, 'KELT'
	elseif matches('^Kepler-') then	-- "Kepler" spacecraft
		return name, 'Kepler'
	elseif matches('^KOI-') then	-- "Kepler Objects of Interest"
		return name, 'KOI'
	elseif matches('^MOA-') then	-- "Microlensing Observations in Astrophysics"
		return name, 'MOA'
	elseif matches('^OGLE') then	-- "Optical Gravitational Lensing Experiment"
		return name, 'OGLE'
	elseif matches('^PH-') then		-- "Planet Hunters"
		return name, 'PH'
	elseif matches('^Qatar-') then	-- only "Qatar-1" or "Qatar-2", discovered by Dr. Khalid Al Subai from Qatar
		return name, 'Qatar'
	elseif matches('^SAND') then	-- only "SAND364"
		return name, 'SAND'
	elseif matches('^SAO ') then	-- ???
		return  name, 'SAO'
	elseif matches('^SDSS') then	-- "Sloan Digital Sky Survey"
		return name, 'SDSS'
	elseif matches('^SWEEPS-') then	-- "Sagittarius Window Eclipsing Extrasolar Planet Search"
		return name, 'SWEEPS'
	elseif matches('^TrES-') then	-- Trans-Atlantic Exoplanet Survey
		return name, 'TrES'
	elseif matches('^UGA-') then	-- only "UGA-1785", probably for "University of Georgia"
		return name, 'UGA'
	elseif matches('^WASP-') then	-- "Wide Angle Search for Planets"
		return name, 'WASP'
	elseif matches('^WISE') then	-- ???
		return name, 'WISE'
	elseif matches('^WTS-') then	-- only "WTS-1" and "WTS-2"
		return name, 'WTS'
	elseif matches('^XO-') then		-- "XO Telescope" on Haleakala, Maui
		return name, 'XO'
	elseif matches('^YBP') then		-- only "YBP1194" and "YBP1514", both discovered by "European Southern Observatory"
		return name, 'YBP'
	-- prefix is number, rest is name of survey
	-- all these capture all but the last space
	-- so be sure to read it afterwards
	elseif matches('^([%a%d]) ') 
	or matches('^([%a%d][%a%d]) ') 
	or matches('^(V%d+) ')					-- some are of the form "V#### Peg", etc
	or matchesKey(greekLettersForAbbrevs)	-- finds a matching key with no regards to suffix
	or matchesKey(greekAbbrevsForLetters)
	then
		local prefix = lasttoken
		if greekLettersForAbbrevs[prefix:lower()] then
			prefix = greekLettersForAbbrevs[prefix:lower()]
			prefix = prefix:sub(1,1):upper() .. prefix:sub(2)
		elseif greekAbbrevsForLetters[prefix:lower()] then
			prefix = prefix:sub(1,1):upper() .. prefix:sub(2)
		end
		nexttoken()	-- skip the prefix
		
		if not matches('^ ') then 
			error("didn't expect just one part for name "..('%q'):format(name).." with remaining "..('%q'):format(name:sub(index))) 
		end
		nexttoken()	-- skip the space
		
		-- now there could be an optional '1' or 'G.' ...
		-- TODO this isn't picking up ...
		if matches('^1') then 
			prefix = prefix .. ' 1'
			nexttoken()
		end
		if matches('^G%.') then
			prefix = prefix .. ' G.'
			nexttoken()
		end
		
		if matchesKey(nameInfoForAbbrev) then
			local part = assert(nameInfoForAbbrev[lasttoken:lower()].star)
			nexttoken()
			local rest = index < #name and name:sub(index) or nil
			if rest and rest:sub(1,1) == ' ' then rest = rest:sub(2) end
			return table{prefix, part, rest}:concat' ', 
					assert(nameInfoForAbbrev[lasttoken:lower()].constellation)
		elseif matchesKey(nameInfoForStar) then
			local part = assert(nameInfoForStar[lasttoken:lower()].star)
			nexttoken()
			local rest = index < #name and name:sub(index) or nil
			if rest and rest:sub(1,1) == ' ' then rest = rest:sub(2) end
			return table{prefix, part, rest}:concat' ', 
					assert(nameInfoForStar[lasttoken:lower()].constellation)
		else
			print("couldn't classify name with prefix "..name)
			return name, ''
		end
	elseif matchesKey(nameInfoForAbbrev) then
		local part = assert(nameInfoForAbbrev[lasttoken:lower()].star)
		nexttoken()
		local rest = index < #name and name:sub(index) or nil
		if rest and rest:sub(1,1) == ' ' then rest = rest:sub(2) end
		return table{part, rest}:concat' ', 
				assert(nameInfoForAbbrev[lasttoken:lower()].constellation)
	-- some more classificaitons:
	-- CHXR 73 ... I'm betting the "CH" stands for "Chameaeleontis", because this is in the "Chameleon" constellation
	-- Lupus-TR 3 ... is in fact in "Lupus" even though it doesn't have the term "Lupi" in the name, and even though it has the hyphen after the constellation "Lupus"
	-- Ross 458 ... now "Ross 458 AB" is also "DT Virginis", but this has no "AB" suffix ... 
	-- UScoCTIO 108 ... USco stands for "Upper Scorpius"
	-- V391 Peg ... is in fact "V391 Pegasi", aka "HS 2201+2610"
	elseif matches('^CHXR ') then
		return name, 'Chameleon'
	elseif matches('^Lupus-') then
		return name, 'Lupus'
	elseif matches('^UScoCTIO ') then
		return name, 'Scorpius'
	else
		print("couldn't classify name "..name)
		return name, ''
	end
end

-- fixes a batch of names
-- also returns the new combination name
local function fixnames(names)
	local newnames = {}	-- new name table is [survey] => [newname]
	local newvalues = {}
	for i=1,#names do
		local oldname = names[i]
		local name, survey = parseName(oldname)
		--if oldname ~= name then print(oldname..' => '..name) end
		-- re-add name here
		newnames[survey] = name
		table.insert(newvalues, name)
	end
	-- remove duplicates from newvalues
	for i=#newvalues-1,1,-1 do
		for j=#newvalues,i+1,-1 do
			if newvalues[i] == newvalues[j] then
				table.remove(newvalues, j)
			end
		end
	end
	return newnames, table.concat(newvalues, nameSeparator)
end

-- fixes a single name
local function fixname(oldname)
	local name, survey = parseName(oldname)
	--if oldname ~= name then print(oldname..' => '..name) end
	-- re-add name here
	return name
end

for _,system in ipairs(resultSystems) do
	if system.name == '' then system.name = nil end
	if system.names then
		system.names, system.name = fixnames(system.names)
	else
		system.name = fixname(system.name)
	end
	
	for i=#system.bodies,1,-1 do
		local body = system.bodies[i]
		if body.name == '' then body.name = nil end
		if body.names then
			body.names, body.name = fixnames(body.names)
		elseif body.name then
			body.name = fixname(body.name)
		-- else the object is a barycenter and its name will be created in the next block of code 
		end
		if body.name == '' then error("body has empty name for star "..('%q'):format(system.name).." #".._) end 
	end
end

-- TODO sort by (1) no-group (i.e. common names) (2) constellation (3) survey
local function getSortForSystem(s)
	local names = s.names
	if not names then 
		if not s.name then return '0' end
		names = {s.name} 
	end
	for k,info in pairs(nameInfoForConstellation) do
		if names[info.constellation] then return '1'..k:lower() end
	end
	if names[''] then return '2'..names['']:lower() end
	return '3'..next(s):lower()
end
table.sort(resultSystems, function(sa, sb)
	return getSortForSystem(sa) < getSortForSystem(sb)
end)
for _,system in ipairs(resultSystems) do
	table.sort(system.bodies, function(ba, bb)
		return getSortForSystem(ba) < getSortForSystem(bb)
	end)
end

-- now that all names have been accumulated,
for _,system in ipairs(resultSystems) do
	for i=#system.bodies,1,-1 do
		local body = system.bodies[i]
		-- some barycenters have no names ...
		if not body.name then
			-- give them the name of the solar system?  we don't worry about duplicate names colliding with the systems, since lookup is unique to each system
			if body.type == 'barycenter' then
				local parts = {}
				for j=1,#system.bodies do
					local otherBody = system.bodies[j]
					if otherBody.parent == body then
						if not otherBody.name then
							error("failed to find name of child of barycenter for system "..system.name)
						end
						table.insert(parts, assert(otherBody.name))
					end
				end
				assert(#parts > 0)
				body.name = table.concat(parts, ' + ')
			end
		end
	end
end


-- assert that all names are unique
local allNames = {}
for _,system in ipairs(resultSystems) do
	for _,body in ipairs(system.bodies) do
		if allNames[body.name] then 
			if body.type == 'star' then
			-- multiple named stars, give them #2, #3 etc suffix
				if allNames[body.name].type == 'star' then
					local nextName
					local counter = 2
					while true do
						nextName = body.name..' #'..counter
						if not allNames[nextName] then break end
						counter = counter + 1
					end
					body.name = assert(nextName)
				else
					error('found duplicate name '..body.name) 
				end
			else
				error('found duplicate name '..('%q'):format(body.name))
			end
		end
		allNames[body.name] = body 
	end
end

-- recursive change .parent field from the object itself to the object's name
for _,system in ipairs(resultSystems) do
	for _,body in ipairs(system.bodies) do
		if body.parent then 
			if not body.parent.name then
				error('body parent has no name: '.. tolua(body.parent))
			end
			body.parent = assert(body.parent.name) 
		end
	end
end

--[[
I feel like I should be doing the data pre-processing here in this file, to make it more than just a conversion from xml's 2 forms of keys (attrs and children) to json's 1 form of keys (object fields)
--]]

local results = {systems=resultSystems}
file['openExoplanetCatalog.json'] = json.encode(results, {indent=true})
