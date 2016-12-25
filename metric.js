var Metric = makeClass();

//low gravity newtonian approximation
var NewtonApproximateMetric = makeClass({
	super : Metric,

	/*
	geodesic calculation:
	x''^u = -Conn^u_ab x'^a x'^b

	Newtonian approximation:
	Phi = G M / |x|
	Conn^j_tt = dPhi/dx^j = G M x^j / |x|^3
	so x''^j = G M x^j / |x|^3
	*/
	calcGravity : function(accel, pos) {
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.mass === undefined) continue;
			var x = pos[0] - planet.pos[0];
			var y = pos[1] - planet.pos[1];
			var z = pos[2] - planet.pos[2];
			r = Math.sqrt(x * x + y * y + z * z);
			var r2 = r * r;
			var accelMagn = gravitationalConstant * planet.mass / r2;
			accel[0] -= x/r * accelMagn;
			accel[1] -= y/r * accelMagn;
			accel[2] -= z/r * accelMagn;
		}
	},

	/*
	geodesic deviation:
	x''^i = R^i_jkl u^j dx^k v^l

	x''^i = R^i_ttt u^t dx^t v^t	= 0 by Riemann anitsymmetry
		  + R^i_jtt u^j dx^t v^t	= 0 by Riemann antisymmetry
		  + R^i_tkt u^t dx^k v^t
		  + R^i_ttl u^t dx^t v^l
		  + R^i_jkt u^j dx^k v^t
		  + R^i_tkl u^t dx^k v^l
		  + R^i_jtl u^j dx^t v^l
		  + R^i_jkl u^j dx^k v^l

	Newtonian approximation:
	R^i_tjt = -R^i_ttj = Phi_,ij

	substitute:
	x''^i = u^t Phi_,ij (dx^j v^t - dx^t v^j)
	Phi_,ij = G M (3 xi xj / |x|^5 - delta_ij / |x|^3)
	x''^i = u^t (dx^j v^t - dx^t v^j) G M (3 xi xj / |x|^5 - delta_ij / |x|^3)
	
	u^t = 1 or else everything would be 0
	either...
		v^j = n^j
		dx^t = 1 or else we'd lose normal component
		either v^t = 0 or dx^j = 0 or both
	or...
		dx^j = n^j
		v^t = 1
		v^j or dx^t = 0 or both

	x''^i = G M (3 xi (x dot n) / |x|^5 - n / |x|^3)

	But what if phi changes wrt time? then phi_,tt is nonzero, right? How does our Riemann metric change?
	*/
	calcTidal : function(accel, pos, srcPlanet) {
		var nx = pos[0] - srcPlanet.pos[0];
		var ny = pos[1] - srcPlanet.pos[1];
		var nz = pos[2] - srcPlanet.pos[2];

		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.index === srcPlanet.index) continue;
			if (planet.mass === undefined) continue;

			var x = pos[0] - planet.pos[0];
			var y = pos[1] - planet.pos[1];
			var z = pos[2] - planet.pos[2];
			var r2 = x*x + y*y + z*z; 
			var r = Math.sqrt(r2);
			var r3 = r * r2;
			var r5 = r2 * r3;
			var xDotN = x * nx + y * ny + z * nz;

			accel[0] += gravitationalConstant * planet.mass * (3 * xDotN * x / r5 - nx / r3);
			accel[1] += gravitationalConstant * planet.mass * (3 * xDotN * y / r5 - ny / r3);
			accel[2] += gravitationalConstant * planet.mass * (3 * xDotN * z / r5 - nz / r3);
		}
	}
});

//rotation-less spherical body
var SchwarzschildMetric = makeClass({
	super : Metric,

	/*
	geodesic calculation:
	x''^u = -Conn^u_ab x'^a x'^b
	
	Schwarzschild:
	Conn^r_tt = R(r-R)/(2r^3)
	*/
	calcGravity : function(accel, pos) {
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.mass === undefined) continue;
			
			var x = pos[0] - planet.pos[0];
			var y = pos[1] - planet.pos[1];
			var z = pos[2] - planet.pos[2];
			var r2 = x*x + y*y + z*z;
			var r = Math.sqrt(r2);
			var r3 = r * r2;
			var R = 2. * planet.mass * kilogramsPerMeter;
			var accelMagn = R * (r - R) / (2 * r3);
			accelMagn *= metersPerUnits.ls * metersPerUnits.ls;
			
			accel[0] -= x/r * accelMagn;
			accel[1] -= y/r * accelMagn;
			accel[2] -= z/r * accelMagn;
		}
	},
	
	/*
	geodesic deviation:
	x''^i = R^i_jkl u^j dx^k v^l

	Schwarzschild:
	u = v = t, dx = n 
	x''^i = R^i_tkt n^k

	r = r, h = theta, p = phi
	x''^r = R^r_trt n^r + R^r_tht n^h + R^r_tpt n^p
	x''^h = R^h_trt n^r + R^h_tht n^h + R^h_tpt n^p
	x''^p = R^p_trt n^r + R^p_tht n^h + R^p_tpt n^p
	
	x''^r = -R(R-r)/r^4 n^r + 0 n^h + 0 n^p
	x''^h = 0 n^r + R(R-r)/(2r^4) n^h + 0 n^p
	x''^p = 0 n^r + 0 n^h + R(R-r)/(2r^4) n^p
	
	x''^r = -R(R-r)/r^4 n^r
	x''^h = R(R-r)/(2r^4) n^h
	x''^p = R(R-r)/(2r^4) n^p
	
	x''^r = R(R-r)/(2r^4) n^r * (1 - 3)
	x''^h = R(R-r)/(2r^4) n^h
	x''^p = R(R-r)/(2r^4) n^p

	x'' = R(R-r)/(2r^4) (n - 3 n^r e_r)

	difference between this and Newton is on the order of 1e-6 for Earth & Moon
	*/
	calcTidal : function(accel, pos, srcPlanet) {
		var nx = pos[0] - srcPlanet.pos[0];
		var ny = pos[1] - srcPlanet.pos[1];
		var nz = pos[2] - srcPlanet.pos[2];
		
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.index === srcPlanet.index) continue;
			if (planet.mass === undefined) continue;
			
			var x = pos[0] - planet.pos[0];
			var y = pos[1] - planet.pos[1];
			var z = pos[2] - planet.pos[2];
			
			var r2 = x*x + y*y + z*z;
			var r = Math.sqrt(r2);

			var R = 2. * planet.mass * kilogramsPerMeter;
			var scale = R * (R - r) / (2 * r2 * r2);
			scale *= metersPerUnits.ls * metersPerUnits.ls;
		
			var nDotR = nx * x + ny * y + nz * z;
			accel[0] += scale * (nx - 3 * nDotR * x / r2);
			accel[1] += scale * (ny - 3 * nDotR * y / r2);
			accel[2] += scale * (nz - 3 * nDotR * z / r2);
		}
	}
});

/*
uncharged, rotating spherical body
using cartesian implementation from http://arxiv.org/pdf/0706.0622.pdf
ds^2 = -dt^2 + dx^2 + dy^2 + dz^2 + (2mr^3)/(r^4 + a^2 z^2) (
	dt + r(x dx + y dy) / (a^2 + r^2) + a(y dx - x dy)/(a^2 + r^2) + z/r dz
)^2
for x^2 + y^2 + z^2 = r^2 + a^2(1 - z^2/r^2)
r^4 + r^2(a^2 - x^2 - y^2 - z^2) - a^2 z^2 = 0
r^2 = 1/2( -(a^2 - x^2 - y^2 - z^2) +- sqrt((a^2 - x^2 - y^2 - z^2)^2 - 4(-a^2 z^2)))


theta = pi/2.										# polar angle
angularVelocity = 2. * pi / (60. * 60. * 24.) / c	# angular velocity, in m^-1
inertia = 2. / 5. * M * R**2						# moment of inertia about a sphere, in m^3
angularMomentum = inertia * angularVelocity			# angular momentum in m^2
a = angularMomentum / M			
Delta(r) = r**2 - 2.*m(r) * r + a**2
Sigma(r) = r**2 + a**2 * cos(theta)**2
kerr_gravity(r) = -2.*m(r) * Delta(r) * (r**2 - a**2 * cos(theta)**2) / (2 * Sigma(r)**3) * c**2

*/
var KerrMetric = makeClass({
	super : Metric,
	calcGravity : function(accel, pos) {
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.mass === undefined) continue;

			var x = pos[0] - planet.pos[0];
			var y = pos[1] - planet.pos[1];
			var z = pos[2] - planet.pos[2];
			var r2DSq = x*x + y*y;
			var r2D = Math.sqrt(r2DSq);
			var rSq = r2DSq + z*z;
			var r = Math.sqrt(rSq);
			var cosTheta = r2D / r;

			var M = planet.mass * gravitationalConstant / metersPerUnits.ls / metersPerUnits.ls;	//mass, in m

			var rotationPeriod = planet.rotationPeriod !== undefined ? planet.rotationPeriod : 0;
			var angularVelocity = 2. * Math.PI / ( rotationPeriod * (60. * 60. * 24.)) / metersPerUnits.ls;	// angular velocity, in 1/m
			var inertia = 2./5. * M * planet.radius * planet.radius;					// moment of inertia about a sphere, in m^3
			var angularMomentum = inertia * angularVelocity;						// angular momentum in m^2
			var a = M == 0 ? 0 : angularMomentum / M;											// m

			var Delta = r * r - 2 * M * r + a * a;
			var Sigma = r * r + a * a * cosTheta * cosTheta;
			var accelMagn = 2 * M * Delta * (r * r - a * a * cosTheta * cosTheta) / (2 * Sigma * Sigma * Sigma);
			accelMagn *= metersPerUnits.ls * metersPerUnits.ls;

			accel[0] -= x/r * accelMagn;
			accel[1] -= y/r * accelMagn;
			accel[2] -= z/r * accelMagn;
		}
	},
	
	//haven't got this yet
	calcTidal : function(accel, pos, srcPlanet) {}
});

var metricInfos = [
	{name:'Newtonian', classObj:NewtonApproximateMetric},
	{name:'Schwarzschild', classObj:SchwarzschildMetric},
	{name:'Kerr', classObj:KerrMetric}
];

//var metric = new NewtonApproximateMetric();
var metric = new SchwarzschildMetric();

var calcMetricForce;
(function(){
	var accel = [];
	var norm = [];
	calcMetricForce = function(x, planet) {
		accel[0] = accel[1] = accel[2] = 0;
		switch (displayMethod) {
		case 'Tangent Tidal':
		case 'Normal Tidal':
		case 'Total Tidal':
			metric.calcTidal(accel, x, planet);
			break;
		case 'Tangent Gravitational':
		case 'Normal Gravitational':
		case 'Total Gravitational':
			metric.calcGravity(accel, x);
			break;
		case 'Tangent Total':
		case 'Normal Total':
		case 'Total':
			metric.calcTidal(accel, x, planet);
			metric.calcGravity(accel, x);
			break;
		}

		norm[0] = x[0] - planet.pos[0];
		norm[1] = x[1] - planet.pos[1];
		norm[2] = x[2] - planet.pos[2];
		var normLen = Math.sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
		norm[0] /= normLen;
		norm[1] /= normLen;
		norm[2] /= normLen;

		//var toTheMoon = (solarSystem.planets[solarSystem.indexes.Moon].pos - x):normalize()
		//var normCrossMoon = norm:cross(toTheMoon)	//points upwards, tangent, right angle to norm and moon
		//var tangentTowardsMoon = normCrossMoon:cross(norm)
		//var tidalAccel = accel:dot(tangentTowardsMoon)
		//var tidalAccel = accel:dot(toTheMoon)	// moonward component

		var t;
		switch (displayMethod) {
		case 'Tangent Tidal':
		case 'Tangent Gravitational':
		case 'Tangent Total':
			var dot = accel[0] * norm[0] + accel[1] * norm[1] + accel[2] * norm[2];
			var dx = accel[0] - norm[0] * dot;
			var dy = accel[1] - norm[1] * dot;
			var dz = accel[2] - norm[2] * dot;
			t = Math.sqrt(dx * dx + dy * dy + dz * dz);
			break;
		case 'Normal Tidal':
		case 'Normal Gravitational':
		case 'Normal Total':
			t = accel[0] * norm[0] + accel[1] * norm[1] + accel[2] * norm[2];
			break;
		case 'Total Tidal':
		case 'Total Gravitational':
		case 'Total':
			t = Math.sqrt(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]);
			break;
		}
		return t;
	};
})();
