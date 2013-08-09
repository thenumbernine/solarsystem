// with a little help from
// http://www.cv.nrao.edu/~rfisher/Ephemerides/ephem_descr.html
// ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/jpl_eph.cc
// http://www.astro-phys.com/js/astro/api.js
var EphemerisData = makeClass({

	//parsed from header.406
	hdr : {
		"numCoeffs" : 728, 
		"epoch1" : 625360.5,
		"epoch2" : 2816912.5,
		"title" : "JPL Planetary Ephemeris DE406\/DE406\nStart Epoch :  JED=   625360.5-3000 FEB 23 00 : 00 : 00\nFinal Epoch :  JED=  2816912.5 3000 MAY 06 00 : 00 : 00",
		"au" : 149597870.691,
		"DEnumber" : 406,
		"interval" : 64,
		"emrat" : 81.30056,
		"vars" : {
			"XD4" : 0.014481653059735,
			"YD4" : 0.0002424631177603,
			"ZD4" : -0.00028152073424801,
			"AM" : 1738,
			"LBET" : 0.00063161213421947,
			"Z7" : -0.25010348393738,
			"GM2" : 7.2434524861627e-10,
			"Y9" : -0.8724783554958,
			"YM" : -0.0019946300016204,
			"XM" : -0.00080817732791148,
			"ZM" : -0.0010872626608381,
			"YD5" : -0.0065181165657927,
			"XD5" : 0.0010920115430089,
			"LGAM" : 0.00022785831399776,
			"GM5" : 2.8253459095242e-07,
			"ZS" : 0.00026605680517703,
			"XS" : 0.0045025081562339,
			"ZD6" : 0.0019286465656538,
			"YD6" : 0.0043358098589553,
			"XD6" : -0.0032175552393033,
			"S41M" : 2.94743374914e-06,
			"AU" : 149597870.691,
			"CLIGHT" : 299792.458,
			"J4E" : -1.616e-06,
			"JDEPOC" : 2440400.5,
			"YD7" : -0.0037624759328466,
			"ZD7" : -0.001651014703068,
			"S43M" : -7.88967312839e-07,
			"GM7" : 1.292024916782e-08,
			"XD8" : 0.0026427710433669,
			"THTC" : 0.4088443,
			"YD8" : -0.0014983144553592,
			"ZD8" : -0.00067904190301796,
			"MGMIS" : 1,
			"YD9" : -0.003143570302152,
			"TDATEF" : 0,
			"S44M" : 5.6404155572e-08,
			"C32M" : 4.8798073627551e-06,
			"THT" : 0.38239065587686,
			"ASUN" : 696000,
			"GDOT" : 0,
			"ZD9" : -0.0010779488297393,
			"GM9" : 2.188699765426e-12,
			"RAD1" : 2439.76,
			"K2E2" : 0.3,
			"ZB" : -0.40154022080181,
			"XB" : 0.1205174172954,
			"YB" : -0.92583847893945,
			"S31M" : 4.2593286300303e-06,
			"Z2" : -0.19527828898023,
			"XDS" : -3.5174820964519e-07,
			"YDS" : 5.1776253995848e-06,
			"ZDS" : 2.2291018543917e-06,
			"ROTEX" : 0,
			"C41M" : -7.17780149806e-06,
			"X4" : -0.11018607428286,
			"PHI" : 0.0051299597051581,
			"Z4" : -0.6058891326142,
			"Y6" : 4.5964778016269,
			"K2M" : 0.029922116659705,
			"Z8" : -9.4001567235488,
			"J2M" : 0.00020431200665465,
			"Z6" : 1.5586975735303,
			"Y8" : -23.942181216179,
			"X6" : 7.8943924419791,
			"C43M" : -8.54788154819e-08,
			"OMGCX" : 0,
			"GMAST2" : 1.2774811891041e-14,
			"PSIDOT" : 0,
			"BETA" : 1,
			"YDB" : 0.0017483093133816,
			"GAMMA" : 1,
			"XDB" : 0.016811268303993,
			"ZDB" : 0.00075820286899054,
			"GMB" : 8.9970113467125e-10,
			"S42M" : -2.8843721272e-06,
			"OMEGAX" : 4.5247044990228e-05,
			"J4M" : -1.45383007072e-07,
			"DROTEY" : -0.001193,
			"ROTEY" : 0,
			"J2SUN" : 2e-07,
			"YDM" : -0.00016744546061515,
			"OMEGAZ" : 0.22994485870137,
			"YD1" : 0.024894520467649,
			"C31M" : 3.0803809783429e-05,
			"OMGCY" : -1.5816707e-06,
			"PSIC" : -1.714509,
			"PHIC" : -0.0042595183,
			"XD9" : 0.00032221044772359,
			"Z9" : 8.9115630409899,
			"GMAST3" : 3.3340587729603e-15,
			"GMAST1" : 6.4668254338426e-14,
			"TAUM" : 0.16671655584928,
			"DROTEX" : 0.000244,
			"TAUE2" : 0.0069417855840523,
			"J3E" : -2.533e-06,
			"TAUE0" : 0,
			"Y7" : -1.1619445055181,
			"ZD2" : 0.0063311055536011,
			"K2E0" : 0.34,
			"LENUM" : 406,
			"TAUE1" : 0.01290895939156,
			"C33M" : 1.7701764624348e-06,
			"J2E" : 0.001082626,
			"YS" : 0.00076707470093238,
			"RAD4" : 3397.515,
			"X8" : -16.055042583768,
			"C44M" : -1.5490389313e-07,
			"XD7" : 0.00022118841741777,
			"Y4" : -1.3275994561326,
			"XD1" : 0.0033674939139841,
			"MAD1" : 1.8,
			"ZD1" : 0.01294630068865,
			"GM4" : 9.5495351057793e-11,
			"X2" : 0.61275194134184,
			"GM1" : 4.9125474514508e-11,
			"MA0002" : 2.9591220828559e-14,
			"AE" : 6378.137,
			"GMS" : 0.00029591220828559,
			"Y2" : -0.34836536849497,
			"Y1" : -0.090781967729586,
			"KVC" : 0,
			"MA0001" : 1.3907873789423e-13,
			"CENTER" : 0,
			"OMGCZ" : 0.229888,
			"TDATEB" : 1.1997060422594e+16,
			"K2E1" : 0.3,
			"RAD2" : 6052.3,
			"XD2" : 0.010952068361699,
			"GM8" : 1.5243589007843e-08,
			"Z5" : -0.22482870022837,
			"X7" : -18.265398306822,
			"C22M" : 2.2517824391662e-05,
			"X1" : 0.36176271460351,
			"YD2" : 0.015617684365262,
			"Z1" : -0.085714983181763,
			"RE" : 6378.137,
			"MAD3" : 5,
			"X9" : -30.483319603999,
			"XDM" : 0.00060108481665913,
			"MAD2" : 2.4,
			"ZDM" : -8.5562144973986e-05,
			"IFAC" : 0.0003,
			"J3M" : 8.7854695077931e-06,
			"MA0004" : 3.8468587077127e-14,
			"PSI" : 1.2941422241103,
			"S33M" : -2.7097009665669e-07,
			"OMEGAY" : -2.2309276319874e-06,
			"EMRAT" : 81.30056,
			"GM6" : 8.4597151856807e-08,
			"C42M" : -1.43951838385e-06,
			"DENUM" : 406,
			"S32M" : 1.6955162680004e-06,
			"Y5" : -0.83048058146015,
			"X5" : -5.3797068988359,
			"ZD5" : -0.002820783165365
		},
		"objs" : [
			{"numSubIntervals" : 4,"numComponents" : 3,"offset" : 3,"name" : "Mercury","numCoeffs" : 14},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 171,"name" : "Venus","numCoeffs" : 12},
			{"numSubIntervals" : 2,"numComponents" : 3,"offset" : 207,"name" : "EM_Bary","numCoeffs" : 9},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 261,"name" : "Mars","numCoeffs" : 10},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 291,"name" : "Jupiter","numCoeffs" : 6},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 309,"name" : "Saturn","numCoeffs" : 6},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 327,"name" : "Uranus","numCoeffs" : 6},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 345,"name" : "Neptune","numCoeffs" : 6},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 363,"name" : "Pluto","numCoeffs" : 6},
			{"numSubIntervals" : 8,"numComponents" : 3,"offset" : 381,"name" : "GeoCMoon","numCoeffs" : 13},
			{"numSubIntervals" : 1,"numComponents" : 3,"offset" : 693,"name" : "Sun","numCoeffs" : 12},
			{"numSubIntervals" : 0,"numComponents" : 2,"offset" : 729,"name" : "Nutation","numCoeffs" : 0},
			{"numSubIntervals" : 0,"numComponents" : 3,"offset" : 729,"name" : "Libration","numCoeffs" : 0}
		]
	},

	interp : function(coeff, numCoeffs, intervalLength, time) {
		var pc = new Float64Array(18);//ffi.new('double[18]')
		var vc = new Float64Array(18);//ffi.new('double[18]')
	
		if (time < 0 || time > 1) throw 'out of bounds';
		// tc is the normalized chebyshev time (-1 <= tc <= 1)
		var tc = 2 * time - 1;
		
		pc[0] = 1;
		pc[1] = tc;

		var twot = tc + tc
		for (var i = 2; i < numCoeffs; ++i) {
			pc[i] = twot * pc[i-1] - pc[i-2];
		}
		
		var pos = 0;
		for (var i = numCoeffs-1; i >= 0; --i) {
			pos += pc[i] * coeff[i];
		}
		
		vc[0] = 0;
		vc[1] = 1;
		vc[2] = twot + twot;
		for (var i=3; i < numCoeffs; ++i) {
			vc[i] = twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
		}
		var vel = 0;
		for (var i = numCoeffs-1; i >= 1; --i) {
			vel += vc[i] * coeff[i];
		}
		vel *= 2 / intervalLength;
		return [pos, vel];
	},

	posVel : function(planetIndex, timeOrigin, timeOffset) {
		timeOffset = +timeOffset;
		var coeffBuffer = this.getCoeffBuffer(timeOrigin, timeOffset);
		var coeff, numCoeffs, subintervalLength, subintervalFrac = getCoeffSubinterval(planetIndex, coeffBuffer, timeOrigin, timeOffset);
		var pos = [];
		var vel = [];
		var tmp;
		for (var i = 0; i < 3; ++i) {
			tmp = interp(coeff, numCoeffs, subintervalLength, subintervalFrac);
			pos[i] = tmp[0];
			vel[i] = tmp[1];
			coeff = coeff + numCoeffs
		}
		return [pos, vel];
	},

	mercury : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Mercury, timeOrigin, timeOffset); },
	venus : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Venus, timeOrigin, timeOffset); },
	mars : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Mars, timeOrigin, timeOffset); },
	jupiter : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Jupiter, timeOrigin, timeOffset); },
	saturn : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Saturn, timeOrigin, timeOffset); },
	uranus : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Uranus, timeOrigin, timeOffset); },
	neptune : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Neptune, timeOrigin, timeOffset); },
	pluto : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Pluto, timeOrigin, timeOffset); },
	sun : function(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Sun, timeOrigin, timeOffset); },

	earth : function(timeOrigin, timeOffset) {
		var tmp = this.posVel(this.objIndexForName.EM_Bary, timeOrigin, timeOffset);
		var earthMoonPos = tmp[0];
		var earthMoonVel = tmp[1] 
		var tmp = this.posVel(this.objIndexForName.GeoCMoon, timeOrigin, timeOffset);
		var geoMoonPos = tmp[0];
		var geoMoonVel = tmp[1];
		var scale = 1 / (1 + this.hdr.emrat);
		
		var earthPos = [
			earthMoonPos[0] - geoMoonPos[0] * scale,
			earthMoonPos[1] - geoMoonPos[1] * scale,
			earthMoonPos[2] - geoMoonPos[2] * scale
		];
		var earthVel = [
			earthMoonVel[0] - geoMoonVel[0] * scale,
			earthMoonVel[1] - geoMoonVel[1] * scale,
			earthMoonVel[2] - geoMoonVel[2] * scale
		];

		return [earthPos, earthVel];
	},

	moon : function(timeOrigin, timeOffset) {
		var tmp = this.posVel(this.objIndexForName.EM_Bary, timeOrigin, timeOffset);
		var earthMoonPos = tmp[0];
		var earthMoonVel = tmp[1];
		var tmp = this.posVel(this.objIndexForName.GeoCMoon, timeOrigin, timeOffset);
		var geoMoonPos = tmp[0];
		var geoMoonVel = tmp[1];
		var scale = 1 / (1 + this.hdr.emrat);
		
		var earthPos = [
			earthMoonPos[0] - geoMoonPos[0] * scale,
			earthMoonPos[1] - geoMoonPos[1] * scale,
			earthMoonPos[2] - geoMoonPos[2] * scale
		];
		var earthVel = [
			earthMoonVel[0] - geoMoonVel[0] * scale,
			earthMoonVel[1] - geoMoonVel[1] * scale,
			earthMoonVel[2] - geoMoonVel[2] * scale
		];
		var moonPos = [
			geoMoonPos[0] + earthPos[0],
			geoMoonPos[1] + earthPos[1],
			geoMoonPos[2] + earthPos[2]
		];
		var moonVel = [
			geoMoonVel[0] + earthVel[0],
			geoMoonVel[1] + earthVel[1],
			geoMoonVel[2] + earthVel[2]
		];
		
		return [moonPos, moonVel];
	},

	objNames : [
		'Mercury',
		'Venus',
		'EM_Bary',
		'Mars',
		'Jupiter',
		'Saturn',
		'Uranus',
		'Neptune',
		'Pluto',
		'GeoCMoon',
		'Sun',
		'Nutation',
		'Libration'
	]
});
EphemerisData.prototype.objIndexForName = {};
$.each(EphemerisData.prototype.objNames, function(objIndex,objName) {
	EphemerisData.prototype.objIndexForName[objName] = objIndex;
});
var ephemerisData = new EphemerisData();

var Planet = makeClass({

	init : function(args) {
		if (args !== undefined) {
			if (args.pos !== undefined) {
				this.pos = [args.pos[0], args.pos[1], args.pos[2]];
			}
			if (args.vel !== undefined) {
				this.vel = [args.vel[0], args.vel[1], args.vel[2]];
			}
		}
		if (this.pos === undefined) {
			this.pos = [0,0,0];
		}
		if (this.vel === undefined) {
			this.vel = [0,0,0];
		}
	},

	add : function(other) {
		var a = this;
		var b = other;
		assert(getmetatable(a) == getmetatable(b))
		var p = new a.init();	//init is the ctor is the class
		
		p.pos[0] = a.pos[0] + b.pos[0];
		p.pos[1] = a.pos[1] + b.pos[1];
		p.pos[2] = a.pos[2] + b.pos[2];
		p.vel[0] = a.vel[0] + b.vel[0];
		p.vel[1] = a.vel[1] + b.vel[1];
		p.vel[2] = a.vel[2] + b.vel[2];
		
		return p;
	},

	//scalar multiply
	mul : function(other) {
		var a = this;
		var b = other;
		var p = new a.init();
		
		p.pos[0] = a.pos[0] * b;
		p.pos[1] = a.pos[1] * b;
		p.pos[2] = a.pos[2] * b;
		p.vel[0] = a.vel[0] * b;
		p.vel[1] = a.vel[1] * b;
		p.vel[2] = a.vel[2] * b;
		
		return p;
	},

	toString : function() {
		return 'Planet('+JSON.stringify(this)+')';
	},
	
	geodeticPosition : function(destX, lat, lon, height) {
		var phi = Math.rad(lat);
		var lambda = Math.rad(lon);
		var cosPhi = Math.cos(phi);
		var sinPhi = Math.sin(phi);
		
		var equatorialRadius = this.equatorialRadius !== undefined 
			? this.equatorialRadius
			: this.radius;
		if (equatorialRadius === undefined) {
			throw "don't know how to calculate this planet's surface "+this;
		}
		var inverseFlattening = this.inverseFlattening
		if (inverseFlattening !== undefined) { 
			var eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening);
			var sinPhiSquared = sinPhi * sinPhi;
			var N = equatorialRadius / Math.sqrt(1 - eccentricitySquared * sinPhiSquared);
			var NPlusH = N + height;
			destX[0] = NPlusH * cosPhi * Math.cos(lambda);
			destX[1] = NPlusH * cosPhi * Math.sin(lambda);
			destX[2] = (N * (1 - eccentricitySquared) + height) * Math.sin(phi);
		} else {
			var NPlusH = equatorialRadius + height;
			destX[0] = NPlusH * cosPhi * Math.cos(lambda);
			destX[1] = NPlusH * cosPhi * Math.sin(lambda);
			destX[2] = NPlusH * Math.sin(phi);
		}
	},

	geodeticNormal : function(lat, lon) {
		var phi = Math.rad(lat)
		var lambda = Math.rad(lon)
		var cosPhi = Math.cos(phi)
		var sinPhi = Math.sin(phi)

		var equatorialRadius = this.equatorialRadius !== undefined
			? this.equatorialRadius
			: this.radius;
		if (equatorialRadius === undefined) {
			throw "don't know how to calculate this planet's surface "+this;
		}
		var inverseFlattening = this.inverseFlattening;
		if (inverseFlattening !== undefined) {
			var eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening);
			var sinPhiSquared = sinPhi * sinPhi;
			var N = equatorialRadius / Math.sqrt(1 - eccentricitySquared * sinPhiSquared);
			var oneMinusEccSq = 1 - eccentricitySquared;
			var NPlusH = N + height;
			var NPrime = sinPhi * eccentricitySquared / (equatorialRadius * equatorialRadius) * N * N * N;
			var sinPhi = Math.sin(phi);
			var cosPhi = Math.cos(phi);
			var sinLambda = Math.sin(lambda);
			var cosLambda = Math.cos(lambda);
			var nx = oneMinusEccSq * NPlusH * cosPhi * cosLambda * (NPrime * sinPhi + NPlusH * cosPhi);
			var nx = oneMinusEccSq * NPlusH * cosPhi * sinLambda * (NPrime * sinPhi + NPlusH * cosPhi);
			var nz = -NPlusH * cosPhi * (NPrime * cosPhi - NPlusH * sinPhi);
			var l = Math.sqrt(nx*nx + ny*ny + nz*nz);
			return [nx/l, ny/l, nz/l];
		} else {
			var x = cosPhi * Math.cos(lambda);
			var y = cosPhi * Math.sin(lambda);
			var z = Math.sin(phi);
			return [x, y, z]
		}
	}
});

var Planets;
Planets = makeClass({
	Planet : Planet,
	/*
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
	*/
	planetClasses : [
		makeClass({super:Planet, name:'sun', mass:1.9891e+30, radius:6.960e+8}),
		makeClass({super:Planet, name:'mercury', mass:3.302e+23, radius:2440e+3, equatorialRadius:2440e+3}),
		makeClass({super:Planet, name:'venus', mass:4.8685e+24, radius:6051.8e+3, equatorialRadius:6051.893e+3}),
		makeClass({super:Planet, name:'earth', mass:5.9736e+24, radius:6371.01e+3, equatorialRadius:6378.136e+3, inverseFlattening:298.257223563}),
		makeClass({super:Planet, name:'moon', mass:7.349e+22, radius:1737.53e+3}),
		makeClass({super:Planet, name:'mars', mass:6.4185e+23, radius:3389.9e+3, equatorialRadius:3397e+3, inverseFlattening:154.409}),
		makeClass({super:Planet, name:'jupiter', mass:1.89813e+27, radius:69911e+3, equatorialRadius:71492e+3, inverseFlattening:1/0.06487}),
		makeClass({super:Planet, name:'saturn', mass:5.68319e+26, radius:58232e+3, equatorialRadius:60268e+3, inverseFlattening:1/0.09796}),
		makeClass({super:Planet, name:'uranus', mass:8.68103e+25, radius:25362e+3, equatorialRadius:25559e+3, inverseFlattening:1/0.02293}),
		makeClass({super:Planet, name:'neptune', mass:1.0241e+26, radius:24624e+3, equatorialRadius:24766e+3, inverseFlattening:1/0.0171}),
		makeClass({super:Planet, name:'pluto', mass:1.314e+22, radius:1151e+3}),
	],

	init : function() {
		for (var i = 0; i < this.planetClasses.length; ++i) {
			this[i] = new (this.planetClasses[i])();
			this[i].index = i;
		}
		this.length = this.planetClasses.length;
	},

	// convert from json query object from http://www.astro-phys.com
	fromAstroPhysState : function(astroPhysPlanets) {
		var planets = new Planets();
		for (var i = 0; i < planets.length; ++i) {
			var planet = planets[i];
			var name = planet.name;
			var astroPhysPlanet = astroPhysPlanets[name];
			
			planet.pos[0] = astroPhysPlanet[0][0] * 1000;
			planet.pos[1] = astroPhysPlanet[0][1] * 1000;
			planet.pos[2] = astroPhysPlanet[0][2] * 1000;
			planet.vel[0] = astroPhysPlanet[1][0] * 1000;
			planet.vel[1] = astroPhysPlanet[1][1] * 1000;
			planet.vel[2] = astroPhysPlanet[1][2] * 1000;
		}
		return planets;
	},

//ephemeris data is a few hundred megs
//so i'm not putting it on my serverspace
//so don't use fromEphemeris
	// convert an ephemeris dataset
	// date is julian date
	//	dir is only used if solarsys.ephemerisData hasn't been initialized
	fromEphemeris : function(date, denum, dir) {
		throw "not implemented";
		if (!ephemerisData.hasInitialized()) {
			ephemerisData.initialize(denum, dir);
		}
		var planets = new Planets();
		for (var i=0; i < planets.length; ++i) {
			var planet = planets[i];
			var tmp = ephemerisData[planet.name](date);	// returns pos & vel in terms of km and km/julian day
			var pos = tmp[0];
			var vel = tmp[1];
			planet.pos[0] = pos[0] * 1000;
			planet.pos[1] = pos[1] * 1000;
			planet.pos[2] = pos[2] * 1000;
			planet.vel[0] = vel[0] * 1000;
			planet.vel[1] = vel[1] * 1000;
			planet.vel[2] = vel[2] * 1000;
		}
		return planets;
	},

	add : function(other) {
		var a = this;
		var b = other;
		var newPlanets = new Planets();
		for (var i=0; i < a.length; ++i) {
			var pa = a[i];
			var pb = b[i];
			newPlanets[i] = pa.add(pb);
		}
		return newPlanets;
	},

	mul : function(other) { // scalar multiplication
		var a = this;
		var b = other;
		var newPlanets = new Planets();
		for (var i=0; i < a.length; ++i) {
			newPlanets[i] = a[i].mul(b);
		}
		return newPlanets;
	},
	
	sub : function(other) {	// rkf45 uses vector subtraction
		var a = this;
		var b = other;
		return a.add(b.mul(-1));	// quick and inefficient implementation
	},

	length : function() {
		var l = 0;
		for (var i=0; i < this.length; ++i) {
			var p = this[i];
			l += Math.sqrt(p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1] + p.pos[2] * p.pos[2]);
			l += Math.sqrt(p.vel[0] * p.vel[0] + p.vel[1] * p.vel[1] + p.vel[2] * p.vel[2]);
		}
		return l;
	},

	toString : function() {
		return 'Planets('+JSON.stringify(this)+')';
	}
});
mergeInto(Planets.prototype, Array.prototype);
Planets.prototype.indexes = {};
for (var i = 0; i < Planets.prototype.planetClasses.length; ++i) {
	Planets.prototype.indexes[Planets.prototype.planetClasses[i].prototype.name] = i;
}


var planets = undefined;
var julianDate = 0;
var targetJulianDate = undefined;
var dateTime = Date.now();
var dateStr = new Date('%Y-%m-%d %H:%M:%S', dateTime);

// track ball motion variables
var orbitPlanetIndex;
var orbitGeodeticLocation;
var orbitDistance;
var orbitTargetDistance;
var orbitZoomFactor = .0003;	// upon mousewheel

var mouse;
var mouseDir;

var gravitationConstant = 6.6738480e-11;	// m^3 / (kg * s^20

//this does not zero accel beforehand!
var calcTidalForce;
(function(){
	x = [];
	srcPlanetToPos = [];

	calcTidalForce = function(accel, pos, srcPlanet) {
		srcPlanetToPos[0] = pos[0] - srcPlanet.pos[0];
		srcPlanetToPos[1] = pos[1] - srcPlanet.pos[1];
		srcPlanetToPos[2] = pos[2] - srcPlanet.pos[2];
		
		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			var planet = planets[planetIndex];
			if (!planetInfluences[planetIndex]) continue;
			if (planet.index === srcPlanet.index) continue;
			
			x[0] = pos[0] - planet.pos[0];
			x[1] = pos[1] - planet.pos[1];
			x[2] = pos[2] - planet.pos[2];
			var xLength = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
			var xToTheSecond = xLength * xLength;
			var xToTheThird = xLength * xToTheSecond;
			var xToTheFourth = xLength * xToTheThird;
			var xToTheFifth = xLength * xToTheFourth;
			
			// a^i = -R^i_jkl dt^j dx^k dt^l = -R^i_tjt = R^i_ttj = -phi_,ij
			// looks like dt^j = [1,0,0,0] so that we only get the t part of R (which is zero)
			// but what if phi changes wrt time? then phi_,tt is nonzero, and how does our Riemann metric change?
			for (i = 0; i < 3; ++i) {
				for (j = 0; j < 3; ++j) { 
					if (i == j) {
						phi_ij = gravitationConstant * planet.mass * (3 * x[i] * x[i] / xToTheFifth - 1 / xToTheThird);
					} else {
						phi_ij = gravitationConstant * planet.mass * (3 * x[i] * x[j]) / xToTheFifth;
					}
					accel[i] -= phi_ij * srcPlanetToPos[j];
				}
			}
		}
	};
})();

//this does not zero accel beforehand!
function calcGravitationForce(accel, pos) {
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		if (!planetInfluences[planetIndex]) continue;
		var planet = planets[planetIndex];
		var x = pos[0] - planet.pos[0];
		var y = pos[1] - planet.pos[1];
		var z = pos[2] - planet.pos[2];
		xLength = Math.sqrt(x * x + y * y + z * z);
		var xToTheSecond = xLength * xLength;
		var xToTheThird = xLength * xToTheSecond;
		accel[0] -= x * (gravitationConstant * planet.mass / xToTheThird);
		accel[1] -= y * (gravitationConstant * planet.mass / xToTheThird);
		accel[2] -= z * (gravitationConstant * planet.mass / xToTheThird);
	}
}

var quatMul;
(function(){
	quatMul = function(res, q, r) {
		var a = (q[3] + q[0]) * (r[3] + r[0]);
		var b = (q[2] - q[1]) * (r[1] - r[2]);
		var c = (q[0] - q[3]) * (r[1] + r[2]);
		var d = (q[1] + q[2]) * (r[0] - r[3]);
		var e = (q[0] + q[2]) * (r[0] + r[1]);
		var f = (q[0] - q[2]) * (r[0] - r[1]);
		var g = (q[3] + q[1]) * (r[3] - r[2]);
		var h = (q[3] - q[1]) * (r[3] + r[2]);

		res[0] =  a - .5 * ( e + f + g + h);
		res[1] = -c + .5 * ( e - f + g - h);
		res[2] = -d + .5 * ( e - f - g + h);
		res[3] =  b + .5 * (-e - f + g + h);
	};
})();

function vec3TransformQuat(dest, src, q) {
	var vx = src[0];
	var vy = src[1];
	var vz = src[2];
	var x = q[0];
	var y = q[1];
	var z = q[2];
	var w = q[3];
	var xx = x*2;
	var yy = y*2;
	var zz = z*2
	var wxx = w*xx;
	var wyy = w*yy;
	var wzz = w*zz
	var xxx = x*xx;
	var xyy = x*yy;
	var xzz = x*zz
	var yyy = y*yy;
	var yzz = y*zz;
	var zzz = z*zz;
	
	dest[0] = (1.0-(yyy+zzz)) * vx + (xyy-wzz) * vy + (xzz+wyy) * vz;
	dest[1] = (xyy+wzz) * vx + (1.0-(xxx+zzz)) * vy + (yzz-wxx) * vz;
	dest[2] = (xzz-wyy) * vx + (yzz+wxx) * vy + (1.0-(xxx+yyy)) * vz;
}

function planetCartesianToSolarSystemBarycentric(destX, srcX, planet) {
	// right now planet angles aren't stored in the planet state (for adaptive integration's sake -- how to weight and combine time + space measurements)
	vec3TransformQuat(destX, srcX, planet.angle);
	
	// now rotate by axial tilt (along
	var tiltAngle = planet.tiltAngle;
	if (tiltAngle !== undefined) {
		vec3TransformQuat(destX, destX, tiltAngle);
	}
	
	destX[0] += planet.pos[0];
	destX[1] += planet.pos[1];
	destX[2] += planet.pos[2];
}

function planetGeodeticToSolarSystemBarycentric(destX, planet, lat, lon, height) {
	planet.geodeticPosition(destX, lat, lon, height);		// position relative to the planet center
	planetCartesianToSolarSystemBarycentric(destX, destX, planet);
}

function mouseRay() {
	var aspectRatio = canvas.width / canvas.height;
	// ray intersect
	var v = [2 * mouse.xf - 1, 1 - 2 * mouse.yf, -1]; 
	// I want to inverse transform the matrices ... especially the infinite projection matrix ...
	v[0] = aspectRatio * v[0] / GL.projMat[0];
	v[1] = aspectRatio * v[1] / GL.projMat[5];
	vec3TransformQuat(v, v, GL.view.angle);
	var l = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	v[0] /= l;
	v[1] /= l;
	v[2] /= l;
	return v;
}

function chooseNewPlanet(mouseDir,doChoose) {
	var bestDot = 0;
	var bestPlanet = undefined; 
	for (var i = 0; i < planets.length; ++i) {
		var planet = planets[i];
		var deltaX = planet.pos[0] - GL.view.pos[0];
		var deltaY = planet.pos[1] - GL.view.pos[1];
		var deltaZ = planet.pos[2] - GL.view.pos[2];
		var deltaLength = Math.sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
		deltaX /= deltaLength;
		deltaY /= deltaLength;
		deltaZ /= deltaLength;
		var dot = deltaX * mouseDir[0] + deltaY * mouseDir[1] + deltaZ * mouseDir[2];
		if (dot > bestDot) {
			bestDot = dot;
			bestPlanet = planet;
		}
	}
	if (bestPlanet !== undefined) {
		hoverPlanetText.text(bestPlanet.name);
		if (bestPlanet.index !== orbitPlanetIndex && doChoose) {
			orbitPlanetIndex = bestPlanet.index;
			var deltaX = GL.view.pos[0] - bestPlanet.pos[0];
			var deltaY = GL.view.pos[1] - bestPlanet.pos[1];
			var deltaZ = GL.view.pos[2] - bestPlanet.pos[2];
			orbitDistance = Math.sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
			orbitTargetDistance = 2 * bestPlanet.radius;
			refreshOrbitTargetDistanceText();
			orbitPlanetText.text(bestPlanet.name);
			refreshMeasureText();
		}
	}
}

function refreshMeasureText() {
	var planet = planets[orbitPlanetIndex];
	measureMinText.text(planet.measureMin === undefined ? '' : (planet.measureMin.toExponential() + ' m/s^2'));
	measureMaxText.text(planet.measureMax === undefined ? '' : (planet.measureMax.toExponential() + ' m/s^2'));
}

// routines for calculating the earth's angle ...
// ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/astrtime.cc
function getLeapSec(mjd, utc) {
	// oh great there is a leap-second file...
	return 0;
}
function tai2utc(taimjd, tai) {
	var mjd = taimjd;
	var ls = getLeapSec(taimjd, tai);
	var utc = tai - ls / 86400;
	if (utc < 0) {
		++utc;
		--mjd;
	}
	if (ls > getLeapSec(mjd, utc) + .00001) {
		utc += 1/86400;
		if (utc > 1) {
			--utc;
			++mjd;
		}
	}
	return [mjd, utc];
}
function tt2utc(ttmjd, tt) {
	var tai = tt - 32.184 / 86400;
	var taimjd = ttmjd;
	if (tai < 0) {
		++tai;
		--taimjd;
	}
	return tai2utc(taimjd, tai);
}
function getUt1Offset(mjd, utc) {
	// yet another absurd time offset file
	return 0;
}
function utc2ut1(mjd, utc) {
	var ut1mjd = mjd;
	var ut1 = utc + getUt1Offset(mjd, utc) / 86400;
	if (ut1 > 1) {
		--ut1;
		++ut1mjd;
	}
	if (ut1 < 0) {
		++ut1;
		--ut1mjd;
	}
	return [ut1mjd, ut1];
}
function iauGmst82() { return 0; } // sofa
function utc2gmst(mjd, utc) {
	var tmp = utc2ut1(mjd, utc);
	var ut1mjd = tmp[0];
	var ut1 = tmp[1];
	return iauGmst82(2400000.5, ut1mjd + ut1) / (2 * Math.pi);
}
function utc2tai(mjd, utc) {
	var taimjd = mjd;
	var tai = utc;
	tai += getLeapSec(mjd, utc) / 86400;
	if (tai > 1) { 
		--tai;
		++taimjd;
	}
	return [taimjd, tai];
}
function utc2tt(mjd, utc) {
	return utc2tai(mjd, utc);
}
function iauEqeq94() { return 0; }	// sofa
function utc2gast(mjd, utc) {
	var gmst = utc2gmst(mjd, utc);
	var tmp = utc2tt(mjd, utc);
	var ttmjd = tmp[0];
	var tt = tmp[1];
	return gmst + iauEqeq94(2400000.5, ttmjd + tt) / (2 * Math.pi);
}
function utc2last(mjd,utc,lon) {
	var lon = 0;	// longitude of our observatory ... from space!
	var last = utc2gast(mjd, utc) + lon / (2 * Math.pi);
	if (last < 0) ++last;	//modulo?
	if (last > 1) --last;
	return last;
}
function getEarthAngle(jd) {
	var jdofs = julianDate - 2400000.5;
	var ttmjd = Math.floor(jdofs);
	var tt = jdofs - ttmjd;
	var tmp = tt2utc(ttmjd, tt);
	var mjd = tmp[0];
	var utc = tmp[1];
	var arg = utc2last(mjd,utc) * Math.pi * 2;
	return deg;
}

var fragmentPrecision = 'precision mediump float;\n';
var vertexPrecision = '';
var gl;
var canvas;
var panel;
var panelContent;
var measureMinText;
var measureMaxText;
var orbitTargetDistanceText;

var hsvTex;
var colorShader;
var texShader;
var hsvShader;
var displayMethods = [
	'Plain',
	'Tangent Tidal',
	'Normal Tidal',
	'Total Tidal',
	'Tangent Gravitational',
	'Normal Gravitational',
	'Total Gravitational',
	'Tangent Total',
	'Normal Total',
	'Total'
];
var displayMethod = 'Plain';
var planetInfluences = [];
var showLinesToOtherPlanets = false;
var showLatAndLonLines = false;
var hoverPlanetText;
var orbitPlanetText;

var heatAlpha = .5;
var colorBarHSVRange = 2/3;	// how much of the rainbow to use

var updatePlanetClassSceneObj;
(function(){
	var x = [];
	var accel = [];
	var norm = [];
	var proj = [];
	updatePlanetClassSceneObj = function(planet) {
		var planetClassPrototype = planet.init.prototype;

		//update tide attributes
		var showTide = displayMethod != 'Plain';
		if (!showTide) {
			if (planetClassPrototype.tex === undefined) {
				if (planetClassPrototype.sceneObj.shader !== colorShader) {
					planetClassPrototype.sceneObj.shader = colorShader;
					planetClassPrototype.sceneObj.texs = [];
				}
			} else {
				if (planetClassPrototype.sceneObj.shader !== texShader) {
					planetClassPrototype.sceneObj.shader = texShader;
					planetClassPrototype.sceneObj.texs = [planetClassPrototype.tex];
				}
			}
		} else {
			if (planetClassPrototype.sceneObj.shader !== hsvShader) {
				planetClassPrototype.sceneObj.shader = hsvShader;
				planetClassPrototype.sceneObj.texs = [planetClassPrototype.tex, hsvTex];
			}
			//and update calculated variable if it is out of date ...
			if (planetClassPrototype.lastMeasureCalcDate !== julianDate) {
				planetClassPrototype.lastMeasureCalcDate = julianDate;
				var measureMin = undefined;
				var measureMax = undefined;
				var vertexIndex = 0;
				for (var tideIndex = 0; tideIndex < planetClassPrototype.sceneObj.attrs.tide.data.length; ++tideIndex) {
					x[0] = planetClassPrototype.sceneObj.attrs.vertex.data[vertexIndex++];
					x[1] = planetClassPrototype.sceneObj.attrs.vertex.data[vertexIndex++];
					x[2] = planetClassPrototype.sceneObj.attrs.vertex.data[vertexIndex++];
					planetCartesianToSolarSystemBarycentric(x, x, planet);
				
					accel[0] = accel[1] = accel[2] = 0;
					switch (displayMethod) {
					case 'Tangent Tidal':
					case 'Normal Tidal':
					case 'Total Tidal':
						calcTidalForce(accel, x, planet);
						break;
					case 'Tangent Gravitational':
					case 'Normal Gravitational':
					case 'Total Gravitational':
						calcGravitationForce(accel, x);
						break;	
					case 'Tangent Total':
					case 'Normal Total':
					case 'Total':
						calcTidalForce(accel, x, planet);
						calcGravitationForce(accel, x);
						break;
					}
					
					norm[0] = x[0] - planet.pos[0];
					norm[1] = x[1] - planet.pos[1];
					norm[2] = x[2] - planet.pos[2];
					var normLen = Math.sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
					norm[0] /= normLen;
					norm[1] /= normLen;
					norm[2] /= normLen;
					
					//var toTheMoon = (planets[planets.indexes.moon].pos - x):normalize()
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
					
					if (measureMin === undefined || t < measureMin) measureMin = t;
					if (measureMax === undefined || t > measureMax) measureMax = t;
					planetClassPrototype.sceneObj.attrs.tide.data[tideIndex] = t;
				}
				planetClassPrototype.measureMin = measureMin;
				planetClassPrototype.measureMax = measureMax;
				for (var i = 0; i < planetClassPrototype.sceneObj.attrs.tide.data.length; ++i) {
					planetClassPrototype.sceneObj.attrs.tide.data[i] = 
						(255/256 - (planetClassPrototype.sceneObj.attrs.tide.data[i] - measureMin) 
							/ (measureMax - measureMin) * 254/256) * colorBarHSVRange;
				}
				planetClassPrototype.sceneObj.attrs.tide.updateData();
				//if it updated...
				if (planet.index == orbitPlanetIndex) {
					refreshMeasureText();
				}
			}
		}
	};
})();

var drawScene;
(function(){
	var delta = vec3.create();
	var viewAngleInv = quat.create();
	var invRotMat = mat4.create();
	var viewPosInv = vec3.create();
	
	drawScene = function() {
		mat4.identity(GL.mvMat);
		
		quat.conjugate(viewAngleInv, GL.view.angle);
		mat4.fromQuat(invRotMat, viewAngleInv);
		mat4.multiply(GL.mvMat, GL.mvMat, invRotMat);
		
		vec3.negate(viewPosInv, GL.view.pos);
		mat4.translate(GL.mvMat, GL.mvMat, viewPosInv);

		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			var planet = planets[planetIndex];
			var planetClassPrototype = planet.init.prototype;

			if (showLinesToOtherPlanets) {
				if (orbitPlanetIndex !== undefined) {
					var orbitPlanet = planets[orbitPlanetIndex];
					if (orbitPlanet !== planet) {
						//while here, update lines
						
						delta[0] = planet.pos[0] - orbitPlanet.pos[0];
						delta[1] = planet.pos[1] - orbitPlanet.pos[1];
						delta[2] = planet.pos[2] - orbitPlanet.pos[2];
						//var dist = Math.sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
						//delta[0] /= dist;
						//delta[1] /= dist;
						//delta[2] /= dist;
						
						//vec3.sub(delta, planet.pos, orbitPlanet.pos);
						//var dist = vec3.length(delta);
						//vec3.scale(delta, delta, 1/dist);
						
						/*	
						planet.lineObj.attrs.vertex.data[0] = orbitPlanet.pos[0] + delta[0] * orbitPlanet.radius;
						planet.lineObj.attrs.vertex.data[1] = orbitPlanet.pos[1] + delta[1] * orbitPlanet.radius;
						planet.lineObj.attrs.vertex.data[2] = orbitPlanet.pos[2] + delta[2] * orbitPlanet.radius;
						planet.lineObj.attrs.vertex.data[4] = orbitPlanet.pos[0] + delta[0] * 100 * orbitPlanet.radius;
						planet.lineObj.attrs.vertex.data[5] = orbitPlanet.pos[1] + delta[1] * 100 * orbitPlanet.radius;
						planet.lineObj.attrs.vertex.data[6] = orbitPlanet.pos[2] + delta[2] * 100 * orbitPlanet.radius;
						*/	
						planet.lineObj.attrs.vertex.data[0] = orbitPlanet.pos[0];
						planet.lineObj.attrs.vertex.data[1] = orbitPlanet.pos[1];
						planet.lineObj.attrs.vertex.data[2] = orbitPlanet.pos[2];
						planet.lineObj.attrs.vertex.data[3] = orbitPlanet.pos[0] + delta[0];
						planet.lineObj.attrs.vertex.data[4] = orbitPlanet.pos[1] + delta[1];
						planet.lineObj.attrs.vertex.data[5] = orbitPlanet.pos[2] + delta[2];
						
						planet.lineObj.attrs.vertex.updateData();
						planet.lineObj.draw();
					}
				}
			}
			
			//update vis ratio
			var dx = planet.pos[0] - GL.view.pos[0];
			var dy = planet.pos[1] - GL.view.pos[1];
			var dz = planet.pos[2] - GL.view.pos[2];
			planet.visRatio = planet.radius / Math.sqrt(dx * dx + dy * dy + dz * dz); 

			//update scene object
			planet.sceneObj.pos[0] = planet.pos[0];
			planet.sceneObj.pos[1] = planet.pos[1];
			planet.sceneObj.pos[2] = planet.pos[2];
			if (planet.tiltAngle) {
				//quat.multiply(planet.sceneObj.angle, planet.tiltAngle, planet.angle);
				quatMul(planet.sceneObj.angle, planet.tiltAngle, planet.angle);
			} else {
				//quat.copy(planet.sceneObj.angle, planet.angle);
				planet.sceneObj.angle[0] = planet.angle[0];
				planet.sceneObj.angle[1] = planet.angle[1];
				planet.sceneObj.angle[2] = planet.angle[2];
				planet.sceneObj.angle[3] = planet.angle[3];
			}

			if (planet.visRatio >= .005) {
				updatePlanetClassSceneObj(planet);
				planet.sceneObj.draw();
			} else {
				var push = mat4.clone(GL.mvMat);
				mat4.translate(GL.mvMat, GL.mvMat, planet.pos);
				
				pointObj.draw({
					uniforms : {
						color : planetClassPrototype.color
					}
				});
			
				mat4.copy(GL.mvMat, push);
			}
			if (showLatAndLonLines) {
				planet.latLonObj.pos[0] = planet.sceneObj.pos[0];
				planet.latLonObj.pos[1] = planet.sceneObj.pos[1];
				planet.latLonObj.pos[2] = planet.sceneObj.pos[2];
				planet.latLonObj.angle[0] = planet.sceneObj.angle[0];
				planet.latLonObj.angle[1] = planet.sceneObj.angle[1];
				planet.latLonObj.angle[2] = planet.sceneObj.angle[2];
				planet.latLonObj.angle[3] = planet.sceneObj.angle[3];
				planet.latLonObj.draw();
			}
		}
	};
})();

//quad of tris 
var quad = [[0,0],[0,1],[1,1],[1,1],[1,0],[0,0]];
var latMin = -90;
var latMax = 90;
var latStep = 5;
var lonMin = -180;
var lonMax = 180;
var lonStep = 5;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	panel.css('height', window.innerHeight);
	GL.resize();

	/**/
	//setup infinite projection matrix
	//http://www.terathon.com/gdc07_lengyel.pdf
	//http://www.gamasutra.com/view/feature/131351/the_mechanics_of_robust_stencil_.php?page=2
	var aspectRatio = canvas.width / canvas.height;
	var nearHeight = Math.tan(GL.view.fovY * Math.PI / 360);
	var nearWidth = aspectRatio * nearHeight;
	var epsilon = 1e-1;
	mat4.identity(GL.projMat);
	GL.projMat[0] = GL.view.zNear / nearWidth;
	GL.projMat[5] = GL.view.zNear / nearHeight;
	GL.projMat[10] = epsilon-1;
	GL.projMat[11] = (epsilon - 2) * GL.view.zNear;
	GL.projMat[14] = -1;
	GL.projMat[15] = 0;
/**/
}

function invalidateForces() {
	//invalidate all
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		var planet = planets[planetIndex];
		var planetClassPrototype = planet.init.prototype;
		planetClassPrototype.lastMeasureCalcDate = undefined;
	}
}

function refreshOrbitTargetDistanceText() {
	orbitTargetDistanceText.text(orbitTargetDistance.toExponential()+' m');
}

$(document).ready(init1);

function init1() {
	panel = $('#panel');
	panelContent = $('#content');
	$('#menu').click(function() {
		if (panelContent.css('display') == 'none') {
			panelContent.show();
		} else {
			panelContent.hide();
		}
	});
	measureMinText = $('#measureMin');
	measureMaxText = $('#measureMax');
	orbitTargetDistanceText = $('#orbitTargetDistance');
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute',
			background : 'red'
		}
	}).prependTo(document.body).get(0);

	$(canvas).disableSelection();

	try {
		gl = GL.init(canvas);
	} catch (e) {
		panel.remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}

	var vtxhigh = gl.getShaderPrecisionFormat(gl.VERTEX_SHADER, gl.HIGH_FLOAT)
	console.log('vertex high? '+vtxhigh.rangeMin+' '+vtxhigh.rangeMax+' '+vtxhigh.precision);
	if (vtxhigh.rangeMin !== 0 && vtxhigh.rangeMax !== 0 && vtxhigh.precision !== 0) {
		console.log('vertex using high precision float');
		vertexPrecision = 'precision highp float;\n';
	}
	var fraghigh = gl.getShaderPrecisionFormat(gl.FRAGMENT_SHADER, gl.HIGH_FLOAT)
	console.log('fragment high? '+fraghigh.rangeMin+' '+fraghigh.rangeMax+' '+fraghigh.precision);
	if (fraghigh.rangeMin !== 0 && fraghigh.rangeMax !== 0 && fraghigh.precision !== 0) {
		console.log('fragment using high precision float');
		fragmentPrecision = 'precision highp float;\n';
	}

	GL.view.zNear = 1;
	GL.view.zFar = 1e+7;//Infinity;

	$('<span>', {text:'Display Method:'}).appendTo(panelContent);
	$('<br>').appendTo(panelContent);
	$.each(displayMethods, function(displayMethodIndex,thisDisplayMethod) {
		var radio = $('<input>', {
			type : 'radio',
			value : displayMethodIndex,
			click : function() {
				displayMethod = thisDisplayMethod;
				invalidateForces();
			}
		})
			.attr('name', 'display')
			.appendTo(panelContent);
		if (thisDisplayMethod == displayMethod) radio.attr('checked', 'checked');
		$('<span>', {text:thisDisplayMethod}).appendTo(panelContent);
		$('<br>').appendTo(panelContent);
	});
	$('<br>').appendTo(panelContent);
	$('<span>', {text:'Influencing Planets:'}).appendTo(panelContent);
	$('<br>').appendTo(panelContent);
	$.each(Planets.prototype.planetClasses, function(planetIndex,planetClass) {
		planetInfluences[planetIndex] = true;
		var planetName = planetClass.prototype.name;
		var checkbox;
		checkbox = $('<input>', {
			type : 'checkbox',
			value : '1',
			change : function() {
				planetInfluences[planetIndex] = checkbox.is(':checked');
				invalidateForces();
			}
		})
			.attr('checked', 'checked')
			.appendTo(panelContent);
		$('<span>', {text:planetName}).appendTo(panelContent);
		$('<br>').appendTo(panelContent);
	});
	hoverPlanetText = $('#hoverPlanetText');
	orbitPlanetText = $('#orbitPlanetText');

	(function(){
		var checkbox = $('#showLinesToOtherPlanets');
		checkbox.change(function() {
			showLinesToOtherPlanets = checkbox.is(':checked');
		});
	})();

	(function(){
		var checkbox = $('#showLatAndLonLines ');
		checkbox.change(function() {
			showLatAndLonLines = checkbox.is(':checked');
		});
	})();

	$(window).resize(resize);
	resize();

	if (true) { //getting single sets of states ...
		//returning a single state
		var planetNames = [];
		for (var planetIndex = 0; planetIndex < Planets.prototype.planetClasses.length; ++planetIndex) {
			planetClass = Planets.prototype.planetClasses[planetIndex];
			planetNames.push(planetClass.prototype.name);
		}

		if (true) {
			var url = 'http://www.astro-phys.com/api/de406/states?bodies=' + planetNames.join(',');
			//var url = 'astro-phys-state.json';
			console.log('reading from '+url);
			$.ajax({
				url : url,
				dataType : 'jsonp'
			}).done(function(d) {
				planets = Planets.prototype.fromAstroPhysState(d.results);
				julianDate = d.date;
				init2();
			});
		}

		if (false) { // astro.getStates not working as well as my own manual call ...
			console.log('planetNames '+JSON.stringify(planetNames));
			console.log('julianDate '+julianDate);
			astro.getStates(planetNames, julianDate, function(results, date) {
				console.log('results '+JSON.stringify(results));
				console.log('date '+date);
				planets = Planets.prototype.fromAstroPhysState(results);
				julianDate = date;
				init2();
			});
		}
	}

	if (false) { //getting records ...
		var url = 'http://www.astro-phys.com/api/de406/records';
		//var url = 'astro-phys-records.json';
		console.log('reading from '+url);
		$.ajax({
			url : url,
			dataType : 'jsonp'
		}).done(function(data) {
			//copied from astro-api.  i need to change the apicall but keep the function the same
			var results = data.results;
			var start = data.start;
			var end = data.end;
			for(var body in results) {
				var nchunks = results[body].length
				var c_step = (end - start) / nchunks;
				for(var i = 0; i < nchunks; ++i) {
					var c_start = start + i * c_step;
					var c_end = c_start + c_step;
					var coeffs = results[body][i];
					results[body][i] = new astro.CoeffSet(coeffs, c_start, c_end);
				}
			}
			
			var record = new astro.Record(results, start, end);
			julianDate = data.date;
			
			//record.getStates is giving back the whole coefficients ... ?	
			var positions = record.getPositions(julianDate);
			planets = new Planets();
			for (var i = 0; i < planets.length; ++i) {
				var planet = planets[i];
				var name = planet.name;
				planet.pos[0] = positions[name][0] * 1000;
				planet.pos[1] = positions[name][1] * 1000;
				planet.pos[2] = positions[name][2] * 1000;
				planet.vel[0] = 0;
				planet.vel[1] = 0;
				planet.vel[2] = 0;
			}
			init2();
		});
	}
}

function init2() {
	
	var imgs = [];
	for (var i = 0; i < planets.length; ++i) {
		imgs.push('textures/'+planets[i].name+'.png');
	}
	console.log('loading '+imgs.join(', '));
	$(imgs).preload(function() {
		$('#loadingDiv').hide();
		$('#menu').show();
		init3();
	}, function(percent){
		$('#loading').attr('value', parseInt(100*percent));
	});
}

function init3() {
	hsvTex = new GL.HSVTexture(256);

	colorShader = new GL.ShaderProgram({
		vertexCode : vertexPrecision,
		vertexCodeID : 'color-vsh',
		fragmentCode : fragmentPrecision,
		fragmentCodeID : 'color-fsh',
		uniforms : {
			color : [1,1,1,1]
		}
	});

	texShader = new GL.ShaderProgram({
		vertexCode : vertexPrecision,
		vertexCodeID : 'tex-vsh',
		fragmentCode : fragmentPrecision,
		fragmentCodeID : 'tex-fsh',
		uniforms : {
			tex : 0
		}
	});

	hsvShader = new GL.ShaderProgram({
		vertexCode : vertexPrecision,
		vertexCodeID : 'heat-vsh',
		fragmentCode : fragmentPrecision,
		fragmentCodeID : 'heat-fsh',
		uniforms : {
			tex : 0,
			hsvTex : 1
		}
	});
	
	$('#overlay-slider').slider({
		range : 'max',
		width : '200px',
		min : 0,
		max : 100,
		value : 100 * heatAlpha,
		slide : function(event, ui) {
			heatAlpha = ui.value / 100;
			gl.useProgram(hsvShader.obj);
			gl.uniform1f(hsvShader.uniforms.heatAlpha.loc, heatAlpha);
			gl.useProgram(null);
		}
	});
	gl.useProgram(hsvShader.obj);
	gl.uniform1f(hsvShader.uniforms.heatAlpha.loc, heatAlpha);
	gl.useProgram(null);


	pointObj = new GL.SceneObject({
		mode : gl.POINTS,
		shader : colorShader,
		attrs : {
			vertex : new GL.ArrayBuffer({data:[0,0,0]})
		},
		pos : [0,0,0],
		parent : null,
		static : false
	});

	var colors = {
		sun:[1,1,0],
		mercury:[.7,0,.2],
		venus:[0,1,0],
		earth:[0,0,1],
		moon:[.6,.6,.6],
		mars:[1,0,0],
		jupiter:[1,.5,0],
		saturn:[1,0,.5],
		uranus:[0,1,1],
		neptune:[1,0,1],
		pluto:[0,.5,1]
	};
	var planetsDone = 0;

	for (var planetIndex_ = 0; planetIndex_ < planets.length; ++planetIndex_) { (function(){
		var planetIndex = planetIndex_;
		var planet = planets[planetIndex];
		var planetClassPrototype = planet.init.prototype;
		var color = colors[planet.name];
		planetClassPrototype.color = [color[0], color[1], color[2], 1];
		
		planetClassPrototype.angle = quat.create();			// rotation ... only used for earth at the moment

		var triIndexArray = [];
		var latLonIndexArray = [];
		var vertexArray = [];
		var texCoordArray = [];
		var tideArray = [];

		var latdiv = Math.floor((latMax-latMin)/latStep);
		var londiv = Math.floor((lonMax-lonMin)/lonStep);
		var vtx = [0,0,0];
		for (var loni=0; loni <= londiv; ++loni) {
			var lon = lonMin + loni * lonStep;
			for (var lati=0; lati <= latdiv; ++lati) {
				var lat = latMin + lati * latStep;
				
				planetClassPrototype.geodeticPosition(vtx, lat, lon, 0);
				for (var j = 0; j < 3; ++j) {
					vertexArray.push(vtx[j]);
				}
				
				texCoordArray.push(lon / 360 + .5);
				texCoordArray.push(-lat / 180 + .5);

				tideArray.push(0);

				if (loni < londiv && lati < latdiv) {
					for (var j = 0; j < quad.length; ++j) {
						var ofs = quad[j];
						var index = (lati + ofs[0]) + (latdiv + 1) * (loni + ofs[1]);
						triIndexArray.push(index);
					}
					//if we're using 5 div step then every 6 will be 30 degrees
					if ((loni * lonStep) % 30 == 0) {
						latLonIndexArray.push(lati + (latdiv+1) * loni);
						latLonIndexArray.push(lati+1 + (latdiv+1) * loni);
					}
					if ((lati * latStep) % 30 == 0) {
						latLonIndexArray.push(lati + (latdiv+1) * loni);
						latLonIndexArray.push(lati + (latdiv+1) * (loni+1));
					}
				}
			}
		}
		
		var checkDone = function() {
			var vertexBuffer = new GL.ArrayBuffer({data:vertexArray, keep:true});
			planetClassPrototype.sceneObj = new GL.SceneObject({
				mode : gl.TRIANGLES,
				indexes : new GL.ElementArrayBuffer({data:triIndexArray}),
				attrs : {
					vertex : vertexBuffer,
					texCoord : new GL.ArrayBuffer({dim:2, data:texCoordArray}),
					tide : new GL.ArrayBuffer({dim:1, data:tideArray, keep:true, usage:gl.DYNAMIC_DRAW})
				},
				uniforms : {
					color : planetClassPrototype.color
				},
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null,
				static : false
			});
			planetClassPrototype.latLonObj = new GL.SceneObject({
				mode : gl.LINES,
				indexes : new GL.ElementArrayBuffer({data:latLonIndexArray}),
				shader : colorShader,
				attrs : {
					vertex : vertexBuffer
				},
				uniforms : {
					color : [1,1,1,.2]
				},
				useDepth : false,
				blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null,
				static : false
			});
			planetClassPrototype.lineObj = new GL.SceneObject({
				mode : gl.LINES,
				shader : colorShader,
				attrs : {
					vertex : new GL.ArrayBuffer({
						data : new Float32Array(6),
						keep : true, 
						usage : gl.DYNAMIC_DRAW
					}),
				},
				uniforms : {
					color : planetClassPrototype.color
				},
				parent : null,
				static : true 
			});

			++planetsDone;
			if (planetsDone == planets.length) {
				init4();
			}
		};

		// load texture
		var img = new Image();
		img.onload = function() {
			planetClassPrototype.tex = new GL.Texture2D({
				//flipY : true,
				data : img,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				generateMipmap : true
			});
			checkDone();
		};
		img.onerror = function() {
			console.log('failed to find texture for planet '+planet.name);
			checkDone();
		};
		img.src = 'textures/'+planet.name+'.png';


	})(); }
}

function init4() {
	//gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL);
	gl.clearColor(0,0,0,0);

	var trackPlanet = planets[planets.indexes.earth];
	orbitPlanetIndex = trackPlanet.index;
	orbitTargetDistance = 2. * trackPlanet.radius;
	refreshOrbitTargetDistanceText();
	orbitDistance = orbitTargetDistance;
	orbitPlanetText.text(trackPlanet.name);

	var dragging = false;
	var tmpQ = quat.create();	
	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
			dragging = false;
		},
		move : function(dx,dy) {
			dragging = true;
			var rotAngle = Math.PI / 180 * .01 * Math.sqrt(dx*dx + dy*dy);
			quat.setAxisAngle(tmpQ, [-dy, -dx, 0], rotAngle);
			quat.multiply(GL.view.angle, GL.view.angle, tmpQ);
			quat.normalize(GL.view.angle, GL.view.angle);
		},
		passiveMove : function() {
			mouseDir = mouseRay();
			chooseNewPlanet(mouseDir, false);
		},
		zoom : function(zoomChange) {
			dragging = true;
			var scale = Math.exp(-orbitZoomFactor * zoomChange);
			orbitTargetDistance *= scale;
			refreshOrbitTargetDistanceText();
		},
		click : function() {
			if (dragging) return;
			mouseDir = mouseRay();
			chooseNewPlanet(mouseDir, true);
		}
	});

	update();
}

function update() {
	if (mouse.leftClick) {
		if (mouseOverEvent) {
			targetJulianDate = mouseOverEvent.julianDate;
		}
	} else if (mouse.leftDragging) {
		var magn = Math.sqrt(mouse.deltaX * mouse.deltaY + mouse.deltaY * mouse.deltaY) * 1000;
		if (magn > 0) {
			var normDelta = vec2.fromValues(mouse.deltaX / magn, mouse.deltaY / magn);
			var r = quat.create();
			quat.fromAxisAngle(r, [-normDelta[2], normDelta[1], 0], -magn);
			quat.mul(GL.view.angle, GL.view.angle, r);
			quat.normalize(GL.view.angle, GL.view.angle);
		}
	}

	// track ball orbit

	var orbitCenter;
	if (orbitGeodeticLocation !== undefined) {
		planetGeodeticToSolarSystemBarycentric(
			orbitCenter,
			planets[orbitPlanetIndex], 
			orbitGeodeticLocation.lat, 
			orbitGeodeticLocation.lon, 
			orbitGeodeticLocation.height);
	} else {
		orbitCenter = planets[orbitPlanetIndex].pos;
	}
	var viewAngleZAxisX = 2 * (GL.view.angle[0] * GL.view.angle[2] + GL.view.angle[3] * GL.view.angle[1]); 
	var viewAngleZAxisY = 2 * (GL.view.angle[1] * GL.view.angle[2] - GL.view.angle[3] * GL.view.angle[0]); 
	var viewAngleZAxisZ = 1 - 2 * (GL.view.angle[0] * GL.view.angle[0] + GL.view.angle[1] * GL.view.angle[1]); 
	GL.view.pos[0] = orbitCenter[0] + viewAngleZAxisX * orbitDistance;
	GL.view.pos[1] = orbitCenter[1] + viewAngleZAxisY * orbitDistance;
	GL.view.pos[2] = orbitCenter[2] + viewAngleZAxisZ * orbitDistance;
	{
		var logDist = Math.log(orbitDistance);
		var logTarget = Math.log(orbitTargetDistance);
		var coeff = .05;
		var newLogDist = (1 - coeff) * logDist + coeff * logTarget;
		orbitDistance = Math.exp(newLogDist);
	}

/*
	var lastJulianDate = julianDate;
	if (targetJulianDate !== undefined) {
		var logDate = Math.log(julianDate);
		var logTarget = Math.log(targetJulianDate);
		var coeff = .5;
		var deltaLog = logTarget - logDate;
		if (deltaLog < 0) {
			var newLogDate = logDate + coeff * deltaLog;
			julianDate = Math.exp(newLogDate);
		} else {
			julianDate = targetJulianDate;
			targetJulianDate = undefined; 
		}
		planets = Planets.prototype.fromEphemeris(julianDate);
	} else
*/
/*
		if (integrateTimeStep !== undefined) {
			// if we're integrating ...
			planets = integrate(julianDate, planets, integrateTimeStep, integrateFunction, integrationMethod, integrationArgs);
			julianDate += integrateTimeStep;
		}
		datetable = julian.toCalendar(julianDate);
		
		if (lastJulianDate !== julianDate) {
			var deltaJulianDate = julianDate - lastJulianDate;
			var deltaAngle = quat.create();
			quat.fromAxisAngle(deltaAngle, [0,0,1], deltaJulianDate * 2 * Math.PI);
			var deltaPos = vec3.create();
			vec3.sub(deltaPos, GL.view.pos, orbitCenter);
			//deltaAngle[4] = -deltaAngle[4]
			vec3.transformQuat(deltaPos, deltaPos, deltaAngle);
			//deltaAngle[4] = -deltaAngle[4]
			vec3.add(GL.view.pos, GL.view.pos, orbitCenter);
			quat.multiply(viewAngle = deltaAngle * viewAngle
		end
	
*/	
	
	// set the angle for the time
	//2455389.287573 is aligned with the eclipse of 2455389.315718
	//so scale the angle by (2455389.287573 - 2455034.608623) / (2455389.315718 - 2455034.608623)
	{
		//Nowhere in any of this do I seem to be taking tiltAngle into account ...
		var angleOffset = .37;
		var angleScale = 1;//(2455389.287573 - 2455034.608623) / (2455389.315718 - 2455034.608623);
		var qdst = planets[planets.indexes.earth].angle;
		quat.identity(qdst);
		quat.rotateZ(qdst, qdst, ((angleScale * julianDate + angleOffset) % 1) * 2 * Math.PI);
	}


	GL.setupMatrices();
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	drawScene();
	GL.clearAlpha();

	requestAnimFrame(update);
}

