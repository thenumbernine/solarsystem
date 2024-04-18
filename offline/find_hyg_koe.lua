#!/usr/bin/env luajit
--[[
ok when they did this for the local supercluster, they used a 'Wiener map', which is something about following vector field flow or something.
maybe I can do that ...
maybe I can build a poisson function based on all the stars and look at each of their gravitational acceleration ...
but that still won't tell me orbits ...

so how to determine orbits of planets?
especially one orbit vs another ...
for example, take the moons of jupiter or saturn etc
how do we ensure they are orbitting their parent planet and not orbiting one another?
https://en.wikipedia.org/wiki/Standard_gravitational_parameter#Small_body_orbiting_a_central_body
Kepler's 1-2-3 law
	... for circular and elliptic:
	μ = (2 π / T)^2 a^3
	... for circular orbits:
	μ = r v^2
	μ = r^3 ω^2
	v ≈ r / ω
	ω = 2 π / T
	r ≈ a
	... for parabolic:
	μ = 1/2 r v^2

1) read HYG data
	and filter out approximate unaffected ... ?
	... how to tell when sun vs other star is unaffected without considering all ours and their neighbors?
	how about just filter by distance at first?

2) calc mass vs distance between all bodies ... see which is the strongest for each?

... I might just merge this with visualize-stars ...
--]]

-- [[ test run on the horizons solar system bodies
local path = require 'ext.path'
local table = require 'ext.table'
local range = require 'ext.range'
local tolua = require 'ext.tolua'
local string = require 'ext.string'
local vec3d = require 'vec-ffi.vec3d'
local json = require 'dkjson'

local function jsonstrtovec3(t)
	return vec3d(table.mapi(t, function(x) return tonumber(string.trim(x)) end):unpack())
end

local staticVarsFileName = '../horizons/static-vars.json'
local staticVarsJSON = assert(path(staticVarsFileName):read())	-- really is js, has a var = at the front
local bodies = table((assert(json.decode((staticVarsJSON:match'^.-=(.*)$')))))

local dynamicVarsFileName = '../horizons/dynamic-vars.json'
local dynamicVarsJSON = assert(path(dynamicVarsFileName):read())	-- really is js, has a var = at the front
local dynamicVars = assert(json.decode((dynamicVarsJSON:match'^.-=(.*)$')))
assert(#bodies == #dynamicVars.coords)


-- uhh should I be doing this in horizons/getdata.lua ?
-- do this before dynObj / parentIndex calculation
local pluto = assert(select(2, bodies:find(nil, function(body) return body.name == '134340 Pluto' end)))
pluto.name = 'Pluto'	-- fix for barycenters

-- add mass=0 to all bodies without mass
for _,body in ipairs(bodies) do
	body.mass = body.mass or 0
	body.radius = body.radius or 0
end

for i,dynObj in ipairs(dynamicVars.coords) do
	local body = bodies[i]
	assert(body.id == dynObj.id)
	for k,v in pairs(dynObj) do
		-- hmm how come the dynamic vars have dif names than the static vars?
		-- yeah dynamic vars use 'Barycenter' even for non-barycenter planets (without any moons)
		-- and for earth, dynamic vars uses 'Earth-Moon Barycenter'
		-- so I'll just use the static var name
		if body[k] then
			if k ~= 'name'
			and body[k] ~= v
			then
				error("don't match "..tolua{
					id = body.id,
					field = k,
					lhs = body[k],
					rhs = v,
				})
			end
		else
			body[k] = v
		end
	end
	-- convert km to m
	body.pos = jsonstrtovec3(body.pos) * 1000
	-- convert km/day to m/day
	body.vel = jsonstrtovec3(body.vel) * 1000
	if body.parent then
		-- [[ only add 'Barycenter' if the parent has moons
		local parentName = ({
			Earth = true,
			Mars = true,
			Jupiter = true,
			Saturn = true,
			Uranus = true,
			Neptune = true,
			Pluto = true,
		})[body.parent] and body.parent..' Barycenter' or body.parent
		--]]
		--[[ always add Barycenter, since there's an entry there even if a planet has no moons
		-- but shoudl they be ?  why have a seprate barycenter of a single-body system?
		local parentName = body.parent == 'Sun' and body.parent or body.parent..' Barycenter'
		--]]
		body.parentIndex = assert((bodies:find(nil, function(body2)
			return body2.name == parentName
		end)), "couldn't find parent named "..tolua(parentName).." for body named "..tolua(body.name))
	end
end

-- TODO accumulate mass into parents for barycenters ...

-- rename the barycenters of planets without moons?
-- or just remove them?
-- for some reason the barycenters do exist, but their name has "barycenter" removed....
bodies:remove(
	assert((
		bodies:find(nil, function(body)
			return body.id == 1
		end)
	), "couldn't find Mercury Barycenter to remove it")
)
bodies:remove(
	assert((
		bodies:find(nil, function(body)
			return body.id == 2
		end)
	), "couldn't find Venus Barycenter to remove it")
)

do -- no unique names?
	local namesUsed = {}
	for _,body in ipairs(bodies) do
		if namesUsed[body.name] then
			error("duplicate name "..body.name)
		end
		namesUsed[body.name] = true
	end
end


local bodiesWithMass = table(bodies):sort(function(a,b)
	return a.mass > b.mass	-- largest mass first
end)

-- force of B on A
local function gravForceMag(a, b)
	return b.mass / (a.pos - b.pos):lenSq()
end

local G = 6.67408e-11	-- m^3 / (kg s^2)

--[[
for all pairs of bodies, look at COM
look at the if the COM is within either body
if it's within one body then that should be the parent binary
but what if the small body has a more massive and closer?
	cuz we can get false positives: mercury can seem to be the parent of jupiters moons if we ignore jupiter.
we could sort them from smallest to biggest mass
and for grav attraction we need to consider inv sq dist?
--]]
--print(o.name, so.parent, so.mass, jsonstrtovec3(o.pos))
local Image = require 'image'
local img = Image(#bodiesWithMass, #bodiesWithMass, 1, 'double')
for ia=#bodiesWithMass,1,-1 do
	local a = bodiesWithMass[ia]
	local ibs = range(#bodiesWithMass)
	ibs:remove(ia)
	ibs:sort(function(i,j)
		return gravForceMag(a, bodiesWithMass[i]) > gravForceMag(a, bodiesWithMass[j])
	end)
	--[[
	so from here we have all other planets sorted by their gravitational force
	any two could be binary orbits
	how do we know when the force of two pairs is comparable such that it needs to be modeled as a 3-body problem?
	ex. the Sun exerts a force on 2003J12 of 3821175.3127415 m/day
	while Jupiter exerts a force on 2003J12 of 2806107.5098147 m/s
	the sun is stronger
	however Jupiter is the parent planet of the moon
	does that mean this moon is realy orbiting the Sun and not Jupiter?
	technically I think the answer is 'no' because the barycenter will lie well within Jupiter (so the pair is not a binary)

	so now I think I have to look at gravitational potential isobars ...
	and this is only convenient when I'm in a plane and it is isobars and not isosurfaces ...
	but for a generic N-body problem (like the stellar neighborhood) it's going to be isosurfaces ...
	hmm and at this point, we're talking about the isobars of the function of gravitation from *all* bodies, so that's a big sum at every single point in space...
	--]]
	--[=[
	if bodiesWithMass[ibs[1]].name ~= a.parent then
		print(a.name..' has parent '..tostring(a.parent)..' but most significant gravitation force from '..bodiesWithMass[ibs[1]].name)
	end
	--]=]
	-- [[
	local Gpot = 0
	--for _,ib in ipairs(ibs:sub(1,10)) do
	for _,ib in ipairs(ibs) do
		local b = bodiesWithMass[ib]

		local dsq = (a.pos - b.pos):lenSq()

		-- calc grav pull
		local F = b.mass / dsq

		-- [=[
		local d = math.sqrt(dsq)
		-- calc COM / grav param between a & b
		local com = (a.pos * a.mass + b.pos * b.mass) / (a.mass + b.mass)
		local da = (a.pos - com):length()
		local relda = da / d
		-- da < planet.radius mean this is a planet
		--local db = (b.pos - com):length()
		print(F, a.name, b.name,
			da < a.radius and "..."..a.name.." is the parent of "..b.name.." (should be "..b.parent..")"
			or (
				(d - da) < b.radius and "..."..b.name.." is the parent of "..a.name.." (should be "..a.parent..")"
				or ''
			)
		)
		--]=]
		img.buffer[(ia-1) + #bodiesWithMass * (ib-1)] = F

		-- negative potential
		Gpot = Gpot + (G * b.mass) / d
	end
	--]]
	a.Gpot = Gpot
end
img:map(function(x) return math.log(x+1) end)
:normalize()
:save'gravbetweenbodies.png'
--]]
print()

print('gravitational potential')
-- TODO I could do this for all bodies, even those without masses, since Gpot doesn't consider its own body
for _,body in ipairs(table(bodiesWithMass):sort(function(a,b)
	--return a.Gpot > b.Gpot
	return a.mass / a.Gpot > b.mass / b.Gpot
end)) do
	print(body.Gpot, body.mass, body.mass / body.Gpot,
		body.name)
end

