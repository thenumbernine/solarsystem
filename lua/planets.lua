local vec3 = require 'vec.vec3'
local class = require 'ext.class'
require 'ext.serialize'
require 'ext.table'	--map

--[[
planet:
	pos	- vec3 in kilometers
	vel - vec3 in kilometers per julian day
	mass - number in 
--]]
local Planet = class()

function Planet:init(args)
	self.pos = vec3(args and args.pos)
	self.vel = vec3(args and args.vel)
end

function Planet.__add(a,b)
	assert(getmetatable(a) == getmetatable(b))
	local p = getmetatable(a)()
	p.pos = a.pos + b.pos
	p.vel = a.vel + b.vel
	return p
end

function Planet.__mul(a,b)
	assert(a:isa(Planet))
	assert(type(b) == 'number')
	local p = getmetatable(a)()
	p.pos = a.pos * b
	p.vel = a.vel * b
	return p
end

function Planet:__tostring() return toLua(self, '', ' ') end
function Planet.__concat(a,b) return tostring(a) .. tostring(b) end

function Planet:geodeticPosition(lat, lon, height)
	local phi = math.rad(lat)
	local lambda = math.rad(lon)
	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)
	
	local equatorialRadius = self.equatorialRadius or self.radius
	if not equatorialRadius then
		error("don't know how to calculate this planet's surface"..self)
	end
	local inverseFlattening = self.inverseFlattening
	if inverseFlattening then
		local eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening)
		local sinPhiSquared = sinPhi * sinPhi
		local N = equatorialRadius / math.sqrt(1 - eccentricitySquared * sinPhiSquared)
		local NPlusH = N + height
		local x = NPlusH * cosPhi * math.cos(lambda) 
		local y = NPlusH * cosPhi * math.sin(lambda)
		local z = (N * (1 - eccentricitySquared) + height) * math.sin(phi)
		return x, y, z
	else
		local NPlusH = equatorialRadius + height
		local x = NPlusH * cosPhi * math.cos(lambda)
		local y = NPlusH * cosPhi * math.sin(lambda)
		local z = NPlusH * math.sin(phi)
		return x, y, z
	end
end

function Planet:geodeticNormal(lat, lon)
	local phi = math.rad(lat)
	local lambda = math.rad(lon)
	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)

	local equatorialRadius = self.equatorialRadius or self.radius
	if not equatorialRadius then
		error("don't know how to calculate this planet's surface"..self)
	end
	local inverseFlattening = self.inverseFlattening
	if inverseFlattening then
		local eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening)
		local sinPhiSquared = sinPhi * sinPhi
		local N = equatorialRadius / math.sqrt(1 - eccentricitySquared * sinPhiSquared)
		local oneMinusEccSq = 1 - eccentricitySquared
		local NPlusH = N + height
		local NPrime = sinPhi * eccentricitySquared / (equatorialRadius * equatorialRadius) * N * N * N
		local sinPhi = math.sin(phi)
		local cosPhi = math.cos(phi)
		local sinLambda = math.sin(lambda)
		local cosLambda = math.cos(lambda)
		local nx = oneMinusEccSq * NPlusH * cosPhi * cosLambda * (NPrime * sinPhi + NPlusH * cosPhi)
		local nx = oneMinusEccSq * NPlusH * cosPhi * sinLambda * (NPrime * sinPhi + NPlusH * cosPhi)
		local nz = -NPlusH * cosPhi * (NPrime * cosPhi - NPlusH * sinPhi)
		local l = math.sqrt(nx*nx + ny*ny + nz*nz)
		return nx/l, ny/l, nz/l
		
		--[[
		N = r (1 - e^2 sin(phi)^2)^-1/2
		dN/dphi = sin(phi) e^2 r (1 - e^2 sin(phi)^2)^-3/2 = sin(phi) N^3 e^2 / r^2
		x = [(N + height) cos(phi) cos(lambda), (N + height) cos(phi) sin(lambda), (N (1 - e^2) + height) sin(phi)]
		
		d/dlambda[x] = [-(N + height) cos(phi) sin(lambda), (N + height) cos(phi) cos(lambda), 0]
		d/dphi[x] = [N' cos(phi) cos(lambda) - (N + height) sin(phi) cos(lambda), N' cos(phi) sin(lambda) - (N + height) sin(phi) sin(lambda), N' (1 - e^2) sin(phi) + (N (1 - e^2) + height) cos(phi)]

		d/dlambda x d/dphi
		
		= [-(N + height) cos(phi) sin(lambda), (N + height) cos(phi) cos(lambda), 0]
		x [N' cos(phi) cos(lambda) - (N + height) sin(phi) cos(lambda), N' cos(phi) sin(lambda) - (N + height) sin(phi) sin(lambda), N' (1 - e^2) sin(phi) + (N (1 - e^2) + height) cos(phi)]

		= [
			(N+height) cos(phi) cos(lambda) (N' (1 - e^2) sin(phi) + (N (1 - e^2) + height) cos(phi)),
			(N' (1 - e^2) sin(phi) + (N (1 - e^2) + height) cos(phi)) (N + height) cos(phi) sin(lambda),
			-(N + height) cos(phi) sin(lambda) (N' cos(phi) sin(lambda) - (N + height) sin(phi) sin(lambda)) - (N' cos(phi) cos(lambda) - (N + height) sin(phi) cos(lambda)) (N + height) cos(phi) cos(lambda)
		]

		= [
			(1-e^2) (N+height) cos(phi) cos(lambda) (N' sin(phi) + (N+height) cos(phi)),
			(1-e^2) (N' sin(phi) + (N+height) cos(phi)) (N + height) cos(phi) sin(lambda),
			-(N+height) cos(phi) (N' cos(phi) - (N + height) sin(phi))
		]

		= [
			(1-e^2) (N+height) cos(phi) cos(lambda) (N' sin(phi) + (N+height) cos(phi)),
			(1-e^2) (N+height) cos(phi) sin(lambda) (N' sin(phi) + (N+height) cos(phi)),
			-(N+height) cos(phi) (N' cos(phi) - (N + height) sin(phi))
		]
		
		... now in terms of x,y,z ...
		
		= [
			(1-e^2) x (N' sin(phi) + (N+height) cos(phi)),
			(1-e^2) y (N' sin(phi) + (N+height) cos(phi)),
			
		]
		
		--]]
	else
		local x = cosPhi * math.cos(lambda)
		local y = cosPhi * math.sin(lambda)
		local z = math.sin(phi)
		return x, y, z
		
		--[[
		x = [(N+height) cos(phi) cos(lambda), (N+height) cos(phi) sin(lambda), (N+height) sin(phi)]
		d/dphi[x] = [-(N+height) sin(phi) cos(lambda), -(N+height) sin(phi) sin(lambda), (N+height) cos(phi)]
		d/dlambda[x] = [-(N+height) cos(phi) sin(lambda), (N+height) cos(phi) cos(lambda), 0]
		
		d/dlambda x d/dphi

		= 	[-(N+height) cos(phi) sin(lambda), (N+height) cos(phi) cos(lambda), 0]
		x	[-(N+height) sin(phi) cos(lambda), -(N+height) sin(phi) sin(lambda), (N+height) cos(phi)]

		= (N+height)^2 * [
			cos(phi)^2 cos(lambda),
			cos(phi)^2 sin(lambda),
			sin(phi) cos(phi) cos(lambda)^2 + sin(phi) cos(phi) sin(lambda)^2 = sin(phi) cos(phi)
		]

		= (N+height)^2 cos(phi) [cos(phi) cos(lambda), cos(phi) sin(lambda), sin(phi)]
		= (N+height)^2 cos(phi) [x,y,z]
		= normalize([x,y,z])
		--]]
	end
end



-- integrate all-at-once
local Planets = class()
Planets.Planet = Planet
--[[
classes of each individual planet
	- name
	- mass (kg)
	- inverseFlattening (?)
	- radius (m)	sun: photosphere radius; jupiter, saturn, uranus, neptune: volumetric mean radius; moon, pluto: only radius; all else: mean radius
	we have a few options for the radius
	- equatorial radius (m)
	- radius (m)			<- for the sun: the photosphere radius
upon construction they have the following instance vars:
	- pos (m)
	- vel (m/julian day)
--]]
Planets.planetClasses = {
	class(Planet, {name='sun', mass=1.9891e+30, radius=6.960e+8}),
	class(Planet, {name='mercury', mass=3.302e+23, radius=2440e+3, equatorialRadius=2440e+3}),
	class(Planet, {name='venus', mass=4.8685e+24, radius=6051.8e+3, equatorialRadius=6051.893e+3}),
	class(Planet, {name='earth', mass=5.9736e+24, radius=6371.01e+3, equatorialRadius=6378.136e+3, inverseFlattening=298.257223563}),
	class(Planet, {name='moon', mass=7.349e+22, radius=1737.53e+3}),
	class(Planet, {name='mars', mass=6.4185e+23, radius=3389.9e+3, equatorialRadius=3397e+3, inverseFlattening=154.409}),
	class(Planet, {name='jupiter', mass=1.89813e+27, radius=69911e+3, equatorialRadius=714192e+3, inverseFlattening=1/0.06487}),
	class(Planet, {name='saturn', mass=5.68319e+26, radius=58232e+3, equatorialRadius=60268e+3, inverseFlattening=1/0.09796}),
	class(Planet, {name='uranus', mass=8.68103e+25, radius=25362e+3, equatorialRadius=25559e+3, inverseFlattening=1/0.02293}),
	class(Planet, {name='neptune', mass=1.0241e+26, radius=24624e+3, equatorialRadius=24766e+3, inverseFlattening=1/0.0171}),
	class(Planet, {name='pluto', mass=1.314e+22, radius=1151e+3}),
}

Planets.indexes = table.map(Planets.planetClasses, function(cl, i) return i, cl.name end)

function Planets:init()
	for i=1,#self.planetClasses do
		self[i] = self.planetClasses[i]()
		self[i].index = i	-- does anyone use this anymore?
	end
end

-- convert from json query object from http://www.astro-phys.com
function Planets.fromAstroPhys(astroPhysPlanets)
	local planets = Planets()
	for i=1,#planets do
		local planet = planets[i]
		local astroPhysPlanet = astroPhysPlanets[name]
		planet.pos = vec3(unpack(astroPhysPlanet[1])) * 1000	-- in m
		planet.vel = vec3(unpack(astroPhysPlanet[2])) * 1000	-- in m/s
	end
	return planets
end

-- convert an ephemeris dataset
-- date is julian date
--	dir is only used if solarsys.eph hasn't been initialized
function Planets.fromEphemeris(date, denum, dir)
	local eph = require 'solarsystem.eph'
	if not eph.hasInitialized() then
		eph.init(denum, dir)
	end
	local planets = Planets()
	for i=1,#planets do
		local planet = planets[i]
		local pos, vel = eph[planet.name](date)	-- returns pos & vel in terms of km and km/julian day
		planet.pos = pos * 1000
		planet.vel = vel * 1000
	end
	return planets
end

--[[
astro-phys vs ephemeris 406 pos & vel magnitude differences at julian date 2456049.41235 (in meters)
3.3981182490492e-06	2.8836300808522e-09
0.00057943229815819	3.6608122786031e-05
0.0012150050330084	2.8597560700262e-05
0.0030637220461267	5.0634628114541e-05
0.0036651042771833	4.6430602743782e-05
0.0024137630710962	4.0148053453085e-06
0.005456164046066	6.8293502747548e-06
0.038529200782488	9.0394919443745e-07
0.019099597605569	6.536673376946e-06
0.044723237426297	4.7332243194799e-06
0.037458086432356	3.0577894415391e-06
--]]

function Planets.__add(a,b)	-- vector addition
	local newPlanets = Planets()
	for i=1,#a do
		local pa = a[i]
		local pb = b[i]
		newPlanets[i] = pa + pb
	end
	return newPlanets
end
function Planets.__mul(a,b)	-- scalar multiplication
	local newPlanets = Planets()
	for i=1,#a do
		newPlanets[i] = a[i] * b
	end
	return newPlanets
end
function Planets.__sub(a,b)	-- rkf45 uses vector subtraction
	return a + (b*-1)	-- quick and inefficient implementation
end

function Planets:length()
	local l = 0
	for i=1,#self do
		local p = self[i]
		l = l + p.pos:length() + p.vel:length()
	end
	return l
end

function Planets:__tostring() return toLua(self, '', ' ') end
function Planets.__concat(a,b) return tostring(a) .. tostring(b) end

return Planets
