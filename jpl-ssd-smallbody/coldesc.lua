return {
	{name='pk', type='integer primary key'},
	{name='bodyType', type='integer'},	-- 0 = comet, 1 = numbered asteroid, 2 = unnumbered asteroid
	{name='idNumber', type='text'},	-- comets, asteroids
	{name='name', type='text'},
	{name='epoch', type='number'},
	{name='perihelionDistance', type='number'},	-- comets
	{name='semiMajorAxis', type='number'},	-- asteroids
	{name='eccentricity', type='number'},
	{name='inclination', type='number'},
	{name='argumentOfPeriapsis', type='number'},
	{name='longitudeOfAscendingNode', type='number'},
	{name='meanAnomalyAtEpoch', type='number'},	-- asteroids
	{name='absoluteMagnitude', type='number'},	-- asteroids
	{name='magnitudeSlopeParameter', type='number'},	-- asteroids
	{name='timeOfPerihelionPassage', type='number'},	-- comets
	{name='orbitSolutionReference', type='text'},
}

--[[
if I were to put this in one giant texture
and update comet positions on the GPU, all at once ...
what's the minimum that we need?
A, B, eccentricity, eccentricAnomaly, dE_dt 
dE_dt is based on semiMajorAxis, eccentricity, eccentricAnomaly
A, B are based on semiMajorAxis, longitudeOfAscendingNode, argumentOfPeriapsis, inclination


Ax Ay Az 
Bx By Bz 
a e E 

Ax Ay Az Bx By Bz  <- come from 
	semiMajorAxis
	longitudeOfAscendingNode
	argumentOfPeriapsis
	inclination

coeffA, coeffB <- come from
	eccentricity
	eccentricAnomaly <- comes from
		for parabolic orbits ...
			argumentOfPeriapsis
		for hyperbolic orbits ...
			meanAnomaly <- comes from
				gravitationalParameter	<- G * sun mass in Kg
				semiMajorAxis
				timeOfPeriapsisCrossing (i.e. timeOfPerihelionPassage, for comets only)
		for elliptic orbits ...
			current time
			for comets ...
				meanAnomaly <- comes from
					timeOfPeriapsisCrossing (i.e. timeOfPerihelionPassage, for comets only)
					orbitalPeriod
			for asteroids ...
				meanAnomaly <- comes from
					meanAnomalyAtEpoch
					orbitalPeriod
					epoch
	dE/dt <- comes from
		semiMajorAxis
		gravitationalParameter	<- G * sun mass in Kg
		eccentricity
		eccentricAnomaly

parameters required:
	semiMajorAxis
	longitudeOfAscendingNode
	argumentOfPeriapsis
	inclination
	eccentricity
								-- nothing extra for parabolic orbits
	timeOfPerihelionPassage		-- for hyperbolic comets and elliptic comets
								-- I've got no calculations for hyperbolic asteroids
	orbitalPeriod				-- for elliptic comets and elliptic asteroids
	meanAnomalyAtEpoch			-- for elliptic asteroids
	epoch						-- for elliptic asteroids

9 values per body


--]]
