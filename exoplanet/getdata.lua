#! /usr/bin/env lua -lluarocks.require
require 'ext'
local json = require 'dkjson'
local csv = require 'csv'

local verbose = false
if arg[1] == 'v' then
	verbose = true
end

local data = csv.file('planets.csv')
local columnNames = table.remove(data.rows, 1)
setmetatable(columnNames, nil)
data:setColumnNames(columnNames)

if arg[1] == 'v' then	-- verbose print all
	local columnNames = {'ra', 'st_rah', 'dec', 'st_glon', 'st_glat', 'st_elon', 'st_elat'}
	for i,row in ipairs(data.rows) do
		print(i)
		for _,name in ipairs(columnNames) do
			local v = row[name]
			if v == '' then v = nil end
			if v then print('', name, v) end
		end
		if io.read(1) == 'q' then break end
	end
	os.exit()
end

local auInM = 149597870700			-- 1AU, in meters
local parsecInM = 3.08567758e+16	-- 1 parsec, in meters
local earthMassInKg = 5.9736e+24	-- mass of Earth, in kg
local earthRadiusInM = 6.37101e+6	-- radius of Earth, in m
local jupiterMassInKg = 1.89813e+27	-- mass of Jupiter, in kg
local jupiterRadiusInM = 6.9911e+7	-- radius of Jupiter, in m
local sunMassInKg = 1.9891e+30		-- mass of Sun, in kg
local sunRadiusInM = 6.960e+8		-- radius of Sun, in m

local results = {
	stars = {}
}
for _,row in ipairs(data.rows) do
	local starName = assert(row.pl_hostname)
	if not results.stars[starName] then
	
		-- "right ascension of the planetary system" / "declination of the planetary system"
		-- are these the ra/dec of the star wrt the earth, or of the planet wrt its host star? 
		local rightAscension = math.rad(assert(row.ra))			-- right ascension, in degress -> radians
		local declination = math.rad(assert(row.dec))			-- declination, in degrees -> radians	
		-- there is also row.st_rah, the right ascension in hours, and galactic and ecliptic coordinates
		local rightAscensionHours = math.rad(assert(row.st_rah)) * math.pi * 2 / 24	-- right ascension, in hours -> radians
		local galacticLongitude = math.rad(row.st_glon)
		local galacticLatitude = math.rad(row.st_glat)
		local eclipticLongitude = math.rad(row.st_elon)
		local eclipticLatitude = math.rad(row.st_elat)

		local vec3 = require 'vec.vec3'
		local quat = require 'vec.quat'
		
		local eclipticCartesian = 
			Quat():fromAngleAxis(0,0,1,math.deg(eclipticLongitude)):rotate(	-- right ascension / longitude / eastward angle
				Quat():fromAngleAxis(0,-1,0,math.deg(eclipticLatitude)):rotate(	-- declination / latitude / northward angle
					vec3(1,0,0)	-- starting point, sun @ vernal equinox
				)
			)
		print('ecliptic', eclipticCartesian)
--[[		
		local eclipticInEquatorial = 
			Quat():fromAngleAxis(-1,0,0,math.rad(23.4)):rotate(
				eclipticCartesian
			)
		print('ecliptic in equatorial', eclipticInEquatorial)
--]]
		local galacticCenterInEcliptic = 
			Quat():fromAngleAxis(0,0,1, math.pi/12*(17 + 1/60*(45 + 1/60*(37.19910)))):rotate(
				Quat():fromAngleAxis(0,-1,0, math.pi/180*(-28 + 1/60*(56 + 1/60*(10.2207)))):rotate(
					vec3(1,0,0)
				)
			)
		local galacticCartesian = 
			Quat():fromAngleAxis(0,0,1,math.deg(galacticLongitude)):rotate(	-- right ascension / longitude / eastward angle
				Quat():fromAngleAxis(0,-1,0,math.deg(galacticLatitude)):rotate(	-- declination / latitude / northward angle
					vec3(1,0,0)
				)
			)

		
		-- lets see how well I got these ...
		-- compare galactic vs ecliptic coordinate reconstruction
		
		
		local distance = assert(row.st_dist) * parsecInM		-- distance, in parsecs -> meters
		-- use some of the above to get the x,y,z coordinates, in parsecs

		local temperature = row.st_teff							-- the temperature, in K.  one way of measuring color index
		local mass = row.st_mass and row.st_mass * sunMassInKg	-- mass, in solar masses -> kg
		local radius = row.st_rad and row.st_rad * sunRadiusInM	-- radius, in solar radii -> m
		local density = row.st_dens and row.st_dens * 1000		-- density, in g/cm^3 -> kg/m^3
		-- surface gravity?
		-- luminosity?
		-- metalicity?
		-- rotational velocity?
		results.stars[starName] = {
			name = starName,
			temperature = temperature,
			mass = mass,
			density = density,
			planets = {},
		}
	end
	local star = results.stars[starName]

	-- inclination is "angular distance of the orbital plane from the line of sight"
	local inclination = math.rad(assert(row.pl_orbincl))		-- inclination, in degrees -> radians

	local radius		-- radius, in m
	if row.pl_rade then
		radius = row.pl_rade * earthMassInKg	-- Earth radii 
	elseif row.pl_radj then
		radius = row.pl_radj * jupiterRadiusInM	-- Jupiter radii 
	elseif row.pl_rads then
		radius = row.pl_rads * sunRadiusInM -- Sun radii
	end

	local density
	if row.pl_dens then
		density = row.pl_dens * 1000		-- density, in g/cm^3 -> kg/m^3
	end

	local mass
	if row.pl_masse then
		mass = row.pl_masse * earthMassInKg
	elseif row.pl_massj then
		mass = row.pl_massj * jupiterMassInKg
	elseif row.pl_msinie then	-- minimal mass
		mass = row.pl_msinie * earthMassInKg
	elseif row.pl_msinij then	-- minimal mass
		mass = row.pl_msinij * jupiterMassInKg
	elseif radius and density then
		mass = 4/3 * math.pi * radius^3
	end

	local planet = {
		name = (row.pl_name),								-- should be pl_hostname + pl_letter
		orbitalPeriod = assert(row.pl_orbper),				-- orbital period, in days
		semiMajorAxis = assert(row.pl_orbsmax) * auInM,		-- semi-major axis, in AU -> m
		eccentricity = assert(row.pl_orbeccen),				-- eccentricity (unitless)
		mass = mass,
		radius = radius,
		density = density,
		inclination = inclination,
		timeOfPeriastronCrossing = assert(row.pl_orbtper),	-- time of periastron crossing, in days
		argumentOfPericenter = math.rad(assert(row.pl_orblper)),	-- longitude of periastron <=> argument of periapsis / pericenter, in degrees -> radians
		temperature = assert(row.pl_eqt),					-- equilibrium temperature, in K
		numberOfMoons = row.pl_mnum,
	}
	table.insert(star.planets, planet)
end


