local class = require 'ext.class'
local OutputMethod = require 'output'

local sunMassInKg = 1.9891e+30	-- kg
local gravitationalConstant = 6.6738480e-11		-- m^3 / (kg * s^2)
local julian = assert(loadfile('../horizons/julian.lua'))()

local julianDate = julian.fromCalendar(os.date'!*t')

-- and for just outputting a cloud of points ...
local OutputToPoints = class(OutputMethod)

function OutputToPoints:init(args)
	self.currentBodyType = args.bodyType
end

function OutputToPoints:staticInit()
	self.pts = table()
end

function math.cosh(x)
	return (math.exp(x) + math.exp(-x)) / 2
end
function math.sinh(x)
	return (math.exp(x) - math.exp(-x)) / 2
end

function OutputToPoints:processBody(body)
	local isComet = self.currentBodyType == 0
	local isAsteroid = self.currentBodyType > 0

	local semiMajorAxis
	local eccentricity
	local eccentricAnomaly
	local orbitType
	local meanAnomaly
	local A, B

	local eccentricity = body.eccentricity

	--[[
	local parabolicEccentricityEpsilon = 1e-7
	if math.abs(eccentricity - 1) < parabolicEccentricityEpsilon then
		orbitType = 'parabolic'
	else
	--]]
	if eccentricity > 1 then
		orbitType = 'hyperbolic'
	else
		orbitType = 'elliptic'
	end

	local pericenterDistance
	if isComet then
		pericenterDistance = assert(body.perihelionDistance)
		if orbitType ~= 'parabolic' then
			semiMajorAxis = pericenterDistance / (1 - eccentricity)
		else		--otherwise, if it is parabolic, we don't get the semi-major axis ...
			semiMajorAxis = 0/0
		end
	elseif isAsteroid then
		semiMajorAxis = assert(body.semiMajorAxis)
		pericenterDistance = semiMajorAxis * (1 - eccentricity)
	end

	local gravitationalParameter = gravitationalConstant * sunMassInKg	--assuming the comet mass is negligible, since the comet mass is not provided
	local semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis
	
	--orbital period is only defined for circular and elliptical orbits (not parabolic or hyperbolic)
	local orbitalPeriod
	if orbitType == 'elliptic' then
		orbitalPeriod = 2 * math.pi * math.sqrt(semiMajorAxisCubed / gravitationalParameter) / (60*60*24)	--julian day
	end

	local longitudeOfAscendingNode = assert(body.longitudeOfAscendingNode)
	local cosAscending = math.cos(longitudeOfAscendingNode)
	local sinAscending = math.sin(longitudeOfAscendingNode)

	local argumentOfPeriapsis = assert(body.argumentOfPeriapsis)
	local cosPericenter = math.cos(argumentOfPeriapsis)
	local sinPericenter = math.sin(argumentOfPeriapsis)

	local inclination = assert(body.inclination)
	local cosInclination = math.cos(inclination)
	local sinInclination = math.sin(inclination)

	local oneMinusEccentricitySquared = 1 - eccentricity * eccentricity
	--magnitude of A is a 
	A = {semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
		 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
		 semiMajorAxis * sinPericenter * sinInclination}
	--magnitude of B is a * sqrt(|1 - e^2|)
	B = {semiMajorAxis * math.sqrt(math.abs(oneMinusEccentricitySquared)) * -(cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
		 semiMajorAxis * math.sqrt(math.abs(oneMinusEccentricitySquared)) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
		 semiMajorAxis * math.sqrt(math.abs(oneMinusEccentricitySquared)) * cosPericenter * sinInclination}
	--inner product: A dot B = 0

	local timeOfPeriapsisCrossing
	if isComet then
		timeOfPeriapsisCrossing = body.timeOfPerihelionPassage	--julian day
	end

	if orbitType == 'parabolic' then
		eccentricAnomaly = math.tan(argumentOfPeriapsis / 2)
		meanAnomaly = eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3 
	elseif orbitType == 'hyperbolic' then
		assert(timeOfPeriapsisCrossing)	--only comets are hyperbolic, and all comets have timeOfPeriapsisCrossing defined
		meanAnomaly = math.sqrt(-gravitationalParameter / semiMajorAxisCubed) * timeOfPeriapsisCrossing * 60*60*24	--in seconds
	elseif orbitType == 'elliptic' then
		--in theory I can say 
		--eccentricAnomaly = math.acos((eccentricity + math.cos(argumentOfPeriapsis)) / (1 + eccentricity * math.cos(argumentOfPeriapsis)))
		-- ... but math.acos has a limited range ...

		if isComet then
			timeOfPeriapsisCrossing = body.timeOfPerihelionPassage	--julian day
			local timeSinceLastPeriapsisCrossing = julianDate - timeOfPeriapsisCrossing
			meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * math.pi / orbitalPeriod
		elseif isAsteroid then
			meanAnomaly = body.meanAnomaly
		else
		--	error 'here'
		end
	else
		error 'here'
	end
	
	--solve Newton Rhapson
	--for elliptical orbits:
	--	f(E) = M - E + e sin E = 0
	--	f'(E) = -1 + e cos E
	--for parabolic oribts:
	--	f(E) = M - E - E^3 / 3
	--	f'(E) = -1 - E^2
	--for hyperbolic orbits:
	--	f(E) = M - e sinh(E) - E
	--	f'(E) = -1 - e cosh(E)
	eccentricAnomaly = meanAnomaly
	for i=1,10 do
		local func, deriv
		if orbitType == 'parabolic' then	--parabolic
			func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3
			deriv = -1 - eccentricAnomaly * eccentricAnomaly
		elseif orbitType == 'elliptic' then 	--elliptical
			func = meanAnomaly - eccentricAnomaly + eccentricity * math.sin(eccentricAnomaly)
			deriv = -1 + eccentricity * math.cos(eccentricAnomaly)	--has zeroes ...
		elseif orbitType == 'hyperbolic' then	--hyperbolic
			func = meanAnomaly + eccentricAnomaly - eccentricity  * math.sinh(eccentricAnomaly)
			deriv = 1 - eccentricity * math.cosh(eccentricAnomaly)
		else
			error 'here'
		end

		local delta = func / deriv
		if math.abs(delta) < 1e-15 then break end
		eccentricAnomaly = eccentricAnomaly - delta
	end

	local sinEccentricAnomaly = math.sin(eccentricAnomaly)
	local cosEccentricAnomaly = math.cos(eccentricAnomaly)
	
	--parabolas and hyperbolas don't define orbitalPeriod
	-- so no need to recalculate it
	if orbitalPeriod and meanAnomaly then
		timeOfPeriapsisCrossing = meanAnomaly * orbitalPeriod / (2 * math.pi) --if it is a comet then we're just reversing the calculation above ...
	end

	local dt_dE
	if orbitType == 'parabolic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 + eccentricAnomaly * eccentricAnomaly)
	elseif orbitType == 'elliptic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 - eccentricity * math.cos(eccentricAnomaly))
	elseif orbitType == 'hyperbolic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (eccentricity * math.cosh(eccentricAnomaly) - 1)
	end
	local dE_dt = 1/dt_dE
	--finally using http://en.wikipedia.org/wiki/Kepler_orbit like I should've in the first place ...
	local coeffA, coeffB
	local coeffDerivA, coeffDerivB
	if orbitType == 'parabolic' then
		--...?
	elseif orbitType == 'elliptic' then
		coeffA = math.cos(eccentricAnomaly) - eccentricity
		coeffB = math.sin(eccentricAnomaly)
		coeffDerivA = -math.sin(eccentricAnomaly) * dE_dt
		coeffDerivB = math.cos(eccentricAnomaly) * dE_dt
	elseif orbitType == 'hyperbolic' then
		coeffA = eccentricity - math.cosh(eccentricAnomaly)
		coeffB = math.sinh(eccentricAnomaly)
		coeffDerivA = -math.sinh(eccentricAnomaly) * dE_dt
		coeffDerivB = math.cosh(eccentricAnomaly) * dE_dt
	elseif orbitType == 'parabolic' then
		coeffA = 0/0
		coeffB = 0/0
	end
	
	local posX = A[1] * coeffA + B[1] * coeffB
	local posY = A[2] * coeffA + B[2] * coeffB
	local posZ = A[3] * coeffA + B[3] * coeffB
	local r = math.sqrt(posX^2 + posY^2 + posZ^2)
	OutputToPoints.maxR = math.max(OutputToPoints.maxR or r, r)

	self.pts:insert(posX)
	self.pts:insert(posY)
	self.pts:insert(posZ)
end
function OutputToPoints:staticDone()
	local ffi = require 'ffi'
	require 'ffi.c.stdio'
	local fh = ffi.C.fopen('output.f32', 'wb')
	ffi.C.fwrite(ffi.new('float[?]', #self.pts, self.pts), #self.pts * 4, 1, fh)
	ffi.C.fclose(fh)
	
	print('max radius:',self.maxR)
end

return OutputToPoints
