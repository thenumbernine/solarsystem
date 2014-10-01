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

-- replace empties with nil
for _,row in ipairs(data.rows) do
	assert(#row == #columnNames)
	for i=1,#columnNames do
		if row[i] == '' then row[i] = nil end
	end
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
local maxAngleError = 0	-- angle error of galactic (to equatorial) to ecliptical coordinates
local foundDistance = 0
local didntFindDistance = 0
for _,row in ipairs(data.rows) do
	local starName = assert(row.pl_hostname)
	if not results.stars[starName] then
	
		-- "right ascension of the planetary system" / "declination of the planetary system"
		-- are these the ra/dec of the star wrt the earth, or of the planet wrt its host star? 
		local rightAscension = math.rad(assert(tonumber(row.ra)))			-- right ascension, in degress -> radians
		local declination = math.rad(assert(tonumber(row.dec)))			-- declination, in degrees -> radians	
		-- there is also row.st_rah, the right ascension in hours, and galactic and ecliptic coordinates
		local rightAscensionHours = math.rad(assert(tonumber(row.st_rah))) * math.pi * 2 / 24	-- right ascension, in hours -> radians
		local galacticLongitude = math.rad(tonumber(row.st_glon))
		local galacticLatitude = math.rad(tonumber(row.st_glat))
		local eclipticLongitude = math.rad(tonumber(row.st_elon))
		local eclipticLatitude = math.rad(tonumber(row.st_elat))

		local vec3 = require 'vec.vec3'
		local quat = require 'vec.quat'
		
		local eclipticCartesian = vec3(
			math.cos(eclipticLatitude) * math.cos(eclipticLongitude),
			math.cos(eclipticLatitude) * math.sin(eclipticLongitude),
			math.sin(eclipticLatitude))
			--[[
			quat():fromAngleAxis(0,0,-1,math.deg(eclipticLongitude)):rotate(	-- right ascension / longitude / eastward angle
				quat():fromAngleAxis(0,1,0,math.deg(eclipticLatitude)):rotate(	-- declination / latitude / northward angle
					vec3(1,0,0)	-- starting point, sun @ vernal equinox
				)
			)
			--]]
		--print('ecliptic', eclipticCartesian, 'length', eclipticCartesian:length())

		local galacticCartesian = vec3(
			math.cos(galacticLatitude) * math.cos(galacticLongitude),
			math.cos(galacticLatitude) * math.sin(galacticLongitude),
			math.sin(galacticLatitude))		
			--[[
			quat():fromAngleAxis(0,0,-1,math.deg(galacticLongitude)):rotate(	-- right ascension / longitude / eastward angle
				quat():fromAngleAxis(0,1,0,math.deg(galacticLatitude)):rotate(	-- declination / latitude / northward angle
					vec3(1,0,0)
				)
			)
			--]]
		--print('galacticCartesian', galacticCartesian, 'length', galacticCartesian:length())

		local epsilon = math.rad(23.4)
		local cosEps = math.cos(epsilon)
		local sinEps = math.sin(epsilon)
		local equatorialCartesian = vec3(
			eclipticCartesian[1],
			eclipticCartesian[2] * cosEps - eclipticCartesian[3] * sinEps,
			eclipticCartesian[2] * sinEps + eclipticCartesian[3] * cosEps)
		--print('equatorialCartesian', equatorialCartesian, 'length', equatorialCartesian:length())

		-- now use the "Reconsidering the galactic coordinate system" paper to convert galactic coordinates to ecliptic coordinates ...
		-- [[ eqn 9, J2000 equatorial to galactic
		local m = {	-- row major, so m[i][j] = m_ij
			{-0.054875539390, -0.873437104725, -0.483834991775},
			{0.494109453633, -0.444829594298, 0.746982248696},
			{-0.867666135681, -0.198076389622, 0.455983794523},	
		}
		local eclipticToGalacticCartesian = vec3(
			equatorialCartesian[1] * m[1][1] + equatorialCartesian[2] * m[1][2] + equatorialCartesian[3] * m[1][3],
			equatorialCartesian[1] * m[2][1] + equatorialCartesian[2] * m[2][2] + equatorialCartesian[3] * m[2][3],
			equatorialCartesian[1] * m[3][1] + equatorialCartesian[2] * m[3][2] + equatorialCartesian[3] * m[3][3])

		--print('eclipticToGalacticCartesian', eclipticToGalacticCartesian, 'length', eclipticToGalacticCartesian:length())
		
		local err = math.acos(vec3.dot(galacticCartesian, eclipticToGalacticCartesian))
		maxAngleError = maxAngleError and math.max(maxAngleError, err) or err
		--print('error (degrees)', math.deg(err))
		--print('error (delta)', galacticCartesian - eclipticToGalacticCartesian)

		-- lets see how well I got these ...
		-- compare galactic vs ecliptic coordinate reconstruction
		
		
		local distance
		if row.st_dist then
			distance = assert(tonumber(row.st_dist)) * parsecInM		-- distance, in parsecs -> meters
		end
		local x, y, z
		if distance then
			x = eclipticCartesian[1] * distance
			y = eclipticCartesian[2] * distance
			z = eclipticCartesian[3] * distance
			foundDistance = foundDistance + 1
		else	
			didntFindDistance = didntFindDistance + 1
			print('star '..starName..' missing distance!')
		end
		
		-- use some of the above to get the x,y,z coordinates, in parsecs

		local temperature = row.st_teff							-- the temperature, in K.  one way of measuring color index
		local mass = row.st_mass and row.st_mass * sunMassInKg	-- mass, in solar masses -> kg
		local radius = row.st_rad and row.st_rad * sunRadiusInM	-- radius, in solar radii -> m
		local density
		if row.st_dens then
			density = (tonumber(row.st_dens) or error('failed to deduce density from '..tostring(row.st_dens))) * 1000		-- density, in g/cm^3 -> kg/m^3
		end
		-- surface gravity?
		-- luminosity?
		-- metalicity?
		-- rotational velocity?
		results.stars[starName] = {
			name = starName,
			temperature = temperature,
			mass = mass,
			density = density,
			x = x, y = y, z = z,	-- meters
			planets = {},
		}
	end
	local star = results.stars[starName]

	local orbitalPeriod
	if row.pl_orbper then
		orbitalPeriod = assert(tonumber(row.pl_orbper))				-- orbital period, in days
	end
	
	local semiMajorAxis
	if row.pl_orbsmax then
		semiMajorAxis = assert(tonumber(row.pl_orbsmax)) * auInM		-- semi-major axis, in AU -> m
	end
		
	local eccentricity
	if row.pl_orbeccen then
		eccentricity = assert(tonumber(row.pl_orbeccen))				-- eccentricity (unitless)
	end

	-- inclination is "angular distance of the orbital plane from the line of sight"
	-- does this mean we should add (or subtract) the star coordinate (wrt the sun)'s declination?
	local inclination
	if row.pl_orbincl then
		inclination = math.rad(assert(tonumber(row.pl_orbincl)))		-- inclination, in degrees -> radians
	end

	local timeOfPeriastronCrossing 	-- time of periastron crossing, in days
	if row.pl_orbtper then
		timeOfPeriastronCrossing = assert(tonumber(row.pl_orbtper))
	end

	local argumentOfPericenter 
	if row.pl_orblper then
		argumentOfPericenter = math.rad(assert(row.pl_orblper))		-- longitude of periastron <=> argument of periapsis / pericenter, in degrees -> radians
	end

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
		orbitalPeriod = orbitalPeriod,
		semiMajorAxis = semiMajorAxis,
		eccentricity = eccentricity,
		mass = mass,
		radius = radius,
		density = density,
		inclination = inclination,
		timeOfPeriastronCrossing = timeOfPeriastronCrossing,
		argumentOfPericenter = argumentOfPericenter,
		temperature = tonumber(row.pl_eqt),					-- equilibrium temperature, in K
		numberOfMoons = row.pl_mnum,
	}
	table.insert(star.planets, planet)

	if arg[1] == 'v' then
		if io.read(1) == 'q' then break end
	end
end

print('max angle error (deg)', math.deg(maxAngleError))
print("didn't find distance for "..didntFindDistance.." of "..(foundDistance + didntFindDistance).." stars")

io.writefile('exoplanetResults.json', 'exoplanetInfo = \n' .. json.encode(results, {indent=true}) .. '\n;')
