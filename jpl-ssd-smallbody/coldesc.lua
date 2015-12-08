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
