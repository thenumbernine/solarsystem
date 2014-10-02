#!/usr/bin/env lua -lluarocks.require
require 'ext'
local json = require 'dkjson'
local htmlparser = require 'htmlparser'
require 'htmlparser.xpath'
htmlparser.htmlnonclosing = {}
local tree = htmlparser.parse(io.readfile('systems.xml'))

local function getText(child)
	if not (#child.child == 1) then error("tried to parse text of a child with multiple children: " ..toLua(child)) end
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

local auToM = 149597870700
local parsecsToMeters = 3.08567758e+16
local jupiterMassInKg = 1.89813e+27
local jupiterRadiusInM = 6.9911e+7
local sunMassInKg = 1.9891e+30
local sunRadiusInM = 6.960e+8

local function getSeparation(child)
	-- separation looks to be equivalent of the radius of the orbit
	-- this should be a fallback for when semimajoraxis is not provided
	if findattr(child, 'unit') == 'arcsec' then	-- then process separtion in arcseconds?
		--error('how do we handle arcsecond separation?')
	elseif findattr(child, 'unit') == 'AU' then
		return getNumber(child) * auToM
	else
		-- sometimes units aren't given ... in that case what do we do?
		return getNumber(child) * auToM
		--error('unknown separation type '..toLua(child))
	end
end



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
local function setName(t, v)
	if t.name then 
		t.name = t.name .. ' / ' .. v
	else
		t.name = v
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
			setField(resultBody, 'semiMajorAxis', getNumber(child) * auToM)
		elseif tag == 'separation' then
			setField(resultBody, 'separation', getSeparation(child))
		elseif tag == 'eccentricity' then
			setField(resultBody, 'eccentricity', getNumber(child)) 	-- I'm guessing eccentricity of orbit, and not of spheroid
		elseif tag == 'periastron' then
			setField(resultBody, 'longitudeOfPeriapsis', math.rad(getNumber(child))) 	-- degrees -> radians
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
			setField(resultBody, 'semiMajorAxis', getNumber(child) * auToM)
		elseif tag == 'separation' then
			-- separation looks to be equivalent of the radius of the orbit
			-- this should be a fallback for when semimajoraxis is not provided
			if findattr(child, 'unit') == 'arcsec' then	-- then process separtion in arcseconds?
			elseif findattr(child, 'unit') == 'AU' then
				setField(resultBody, 'separation', getNumber(child) * auToM)
			else
				error('unknown separation type '..toLua(child.attrs))
			end
		elseif tag == 'eccentricity' then
			setField(resultBody, 'eccentricity', getNumber(child)) 	-- I'm guessing eccentricity of orbit, and not of spheroid
		elseif tag == 'periastron' then
			setField(resultBody, 'longitudeOfPeriapsis', math.rad(getNumber(child))) 	-- degrees -> radians
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
				error('found duplicate name '..body.name) 
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
				error('body parent has no name: '.. toLua(body.parent))
			end
			body.parent = assert(body.parent.name) 
		end
	end
end

local results = {systems=resultSystems}
io.writefile('openExoplanetCatalog.json', 'exoplanetCatalogueResults = \n' ..json.encode(results, {indent=true}) .. '\n;')
