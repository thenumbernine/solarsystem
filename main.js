
//provided z, calculate x and y such that x,y,z form a basis
var calcBasis = function(x,y,z) {
	var cxx = 0;
	var cxy = z[2];
	var cxz = -z[1];
	var cyx = -z[2];
	var cyy = 0;
	var cyz = z[0];
	var czx = z[1];
	var czy = -z[0];
	var czz = 0;
	var lx = Math.sqrt(cxy * cxy + cxz * cxz);
	var ly = Math.sqrt(cyx * cyx + cyz * cyz);
	var lz = Math.sqrt(czx * czx + czy * czy);
	if (lx < ly) {
		if (lx < lz) {	//x is smallest
			x[0] = cyx;
			x[1] = cyy;
			x[2] = cyz;
			y[0] = czx;
			y[1] = czy;
			y[2] = czz;
		} else {		//z is smallest
			x[0] = cxx;
			x[1] = cxy;
			x[2] = cxz;
			y[0] = cyx;
			y[1] = cyy;
			y[2] = cyz;
		}
	} else {
		if (ly < lz) {	//y is smallest
			x[0] = czx;
			x[1] = czy;
			x[2] = czz;
			y[0] = cxx;
			y[1] = cxy;
			y[2] = cxz;
		} else {		//z is smallest
			x[0] = cxx;
			x[1] = cxy;
			x[2] = cxz;
			y[0] = cyx;
			y[1] = cyy;
			y[2] = cyz;
		}
	}
	var xLen = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	x[0] /= xLen; x[1] /= xLen; x[2] /= xLen;
	var yLen = Math.sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
	y[0] /= yLen; y[1] /= yLen; y[2] /= yLen;
};

//out = a^T * v
vec3.transformMat3Tr = function(out, v, a) {
	var x = a[0] * v[0] + a[1] * v[1] + a[2] * v[2];
	var y = a[3] * v[0] + a[4] * v[1] + a[5] * v[2];
	out[2] = a[6] * v[0] + a[7] * v[1] + a[8] * v[2];
	out[0] = x;
	out[1] = y;
}

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

	clone : function() {
		var p = new this.init();
		p.pos[0] = this.pos[0];
		p.pos[1] = this.pos[1];
		p.pos[2] = this.pos[2];
		p.vel[0] = this.vel[0];
		p.vel[1] = this.vel[1];
		p.vel[2] = this.vel[2];
		return p;
	},

	add : function(other) {
		var a = this;
		var b = other;
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
			"don't know how to calculate this planet's surface "+this.name;
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
			//mind you this is the planet shape eccentricity, not to be confused with the orbit path eccentricity
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
		- radius (m)
			sun: photosphere radius
			jupiter, saturn, uranus, neptune: volumetric mean radius
			moon, pluto: radius
			all else: mean radius
		we have a few options for the radius
		- equatorial radius (m)
		- radius (m)			<- for the sun: the photosphere radius
	upon construction they have the following instance vars:
		- pos (m)
		- vel (m/julian day)
	*/
	planetClasses : [
		makeClass({super:Planet, id:10, name:'Sun', mass:1.9891e+30, radius:6.960e+8}),
		makeClass({super:Planet, id:199, name:'Mercury', parent:'Sun', mass:3.302e+23, radius:2.440e+6, equatorialRadius:2440e+3}),
		makeClass({super:Planet, id:299, name:'Venus', parent:'Sun', mass:4.8685e+24, radius:6.0518e+6, equatorialRadius:6051.893e+3}),
		makeClass({super:Planet, id:399, name:'Earth', parent:'Sun', mass:5.9736e+24, radius:6.37101e+6, equatorialRadius:6378.136e+3, inverseFlattening:298.257223563}),
		makeClass({super:Planet, id:301, name:'Moon', parent:'Earth', mass:7.349e+22, radius:1.73753e+6}),
		makeClass({super:Planet, id:499, name:'Mars', parent:'Sun', mass:6.4185e+23, radius:3.3899e+6, equatorialRadius:3397e+3, inverseFlattening:154.409}),
		makeClass({super:Planet, id:599, name:'Jupiter', parent:'Sun', mass:1.89813e+27, radius:6.9911e+7, equatorialRadius:71492e+3, inverseFlattening:1/0.06487}),
		makeClass({super:Planet, id:699, name:'Saturn', parent:'Sun', mass:5.68319e+26, radius:5.8232e+7, equatorialRadius:60268e+3, inverseFlattening:1/0.09796}),
		makeClass({super:Planet, id:799, name:'Uranus', parent:'Sun', mass:8.68103e+25, radius:2.5362e+7, equatorialRadius:25559e+3, inverseFlattening:1/0.02293}),
		makeClass({super:Planet, id:899, name:'Neptune', parent:'Sun', mass:1.0241e+26, radius:2.4624e+7, equatorialRadius:24766e+3, inverseFlattening:1/0.0171}),
		makeClass({super:Planet, id:999, name:'Pluto', parent:'Sun', mass:1.314e+22, radius:1.151e+6}),
	],

	init : function() {
		for (var i = 0; i < this.planetClasses.length; ++i) {
			this[i] = new (this.planetClasses[i])();
		}
		this.length = this.planetClasses.length;
	},

	// convert from json query object from http://www.astro-phys.com
	fromAstroPhysState : function(astroPhysPlanets) {
		var planets = new Planets();
		for (var i = 0; i < planets.length; ++i) {
			var planet = planets[i];
			var astroPhysPlanet = astroPhysPlanets[planet.name.toLowerCase()];
		
			//convert km to m
			planet.pos[0] = astroPhysPlanet[0][0] * 1000;
			planet.pos[1] = astroPhysPlanet[0][1] * 1000;
			planet.pos[2] = astroPhysPlanet[0][2] * 1000;
			planet.vel[0] = astroPhysPlanet[1][0] * 1000;
			planet.vel[1] = astroPhysPlanet[1][1] * 1000;
			planet.vel[2] = astroPhysPlanet[1][2] * 1000;
		}
		return planets;
	},

	//static function for adding extra planets
	addExtraHorizonData : function() {
		//start out with the current planets
		var currentIDs = {};
		for (var i = 0; i < Planets.prototype.planetClasses.length; ++i) {
			currentIDs[Planets.prototype.planetClasses[i].prototype.id] = true;
		}

		//I need to fix up my export script ...
		if (horizonsData.coords.length != horizonsExtraData.length) throw 'static to dynamic data lengths differ: dynamic has '+horizonsData.coords.length+' while static has '+horizonsExtraData.length;
		for (var i = 0; i < horizonsData.coords.length; ++i) {
			var data = horizonsData.coords[i];
			var extra = horizonsExtraData[i];
			assert(data.id == extra.id, "got some bad data");
			if (data.name.match(/^L[1-5]/)) {
				//391-395 are Lagrangian points - skip
			} else if (data.name.match(/Barycenter$/)) {
				//1-9 are Barycenter points - skip
			} else {
				if (currentIDs[data.id]) {
					//if the id already exists then skip it
				} else {
					//... finally, add the planet
					/*
					Things we are going to want:
					mass
					radius
					equatorialRadius (all but sun, pluto, moon)
					inverseFlattening (all but sun, mercury, venus, moon, pluto)
					*/				
					Planets.prototype.planetClasses.push(
						makeClass({
							super:Planet,
							id:data.id,
							name:data.name,
							mass:extra.mass,
							radius:extra.radius,
							equatorialRadius:extra.equatorialRadius,
							inverseFlattening:extra.inverseFlattening,
							parent:extra.parent
						})
					);
					currentIDs[data.id] = true;
				}
			}
		}
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
			var tmp = ephemerisData[planet.name.toLowerCase()](date);	// returns pos & vel in terms of km and km/julian day
			var pos = tmp[0];
			var vel = tmp[1];
			
			//convert km to m
			planet.pos[0] = pos[0] * 1000;
			planet.pos[1] = pos[1] * 1000;
			planet.pos[2] = pos[2] * 1000;
			planet.vel[0] = vel[0] * 1000;
			planet.vel[1] = vel[1] * 1000;
			planet.vel[2] = vel[2] * 1000;
		}
		return planets;
	},

	clone : function() {
		var newPlanets = new Planets();
		for (var i = 0; i < this.length; ++i) {
			newPlanets[i] = this[i].clone();
		}
		return newPlanets;
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
	},

	//assume class will be removed separately, just remove the instance
	remove : function(index) {
		if (!(index >= 0 && index < this.length)) {
			throw "tried to remove OOB planet "+index;
		}
		for (var i = index; i < this.length-1; ++i) {
			this[i] = this[i+1];
		}
		this[this.length-1] = undefined;
		--this.length;
	}
});
mergeInto(Planets.prototype, Array.prototype);

Planets.prototype.addExtraHorizonData();	//add the 180 other satellites out there...

Planets.prototype.indexes = {};
Planets.prototype.planetClassForHorizonID = {};
for (var i = 0; i < Planets.prototype.planetClasses.length; ++i) {
	//I am suspicious that, because 'init' is also considered the class itself, it doesn't inherit by default -- and always gets overridden
	var oldPrototype = Planets.prototype.planetClasses[i].prototype;
	Planets.prototype.planetClasses[i] = function() {
		Planet.apply(this, arguments);
	};
	Planets.prototype.planetClasses[i].prototype = oldPrototype;
	Planets.prototype.planetClasses[i].prototype.index = i;
	var planetClass = Planets.prototype.planetClasses[i];
		
	Planets.prototype.indexes[planetClass.prototype.name] = i;
	Planets.prototype.planetClassForHorizonID[planetClass.prototype.id] = planetClass;
}
//convert parent from name to class (or undefined if no such name exists)
for (var i = 0; i < Planets.prototype.planetClasses.length; ++i) {
	var planetClass = Planets.prototype.planetClasses[i];
	if (planetClass.prototype.parent !== undefined) {
		planetClass.prototype.parent = Planets.prototype.indexes[planetClass.prototype.parent];
	}
}

var planets = undefined;
var julianDate = 0;
var lastJulianDate = 0;
var initPlanets = undefined;
var initJulianDate = 0;
var targetJulianDate = undefined;
var dateTime = Date.now();

// track ball motion variables
var mouseOverPlanetIndex;
var orbitPlanetIndex;
var orbitGeodeticLocation;
var orbitDistance;
var orbitOffset = [0,0,0];
var orbitTargetDistance;
var orbitZoomFactor = .0003;	// upon mousewheel

var mouse;
var mouseDir;
var skyCubeObj;
var pointObj;
var starsObj;

var skyTexFilenames = [
	'textures/sky-visible-cube-xp.png',
	'textures/sky-visible-cube-xn.png',
	'textures/sky-visible-cube-yp.png',
	'textures/sky-visible-cube-yn.png',
	'textures/sky-visible-cube-zp.png',
	'textures/sky-visible-cube-zn.png'
];

var glMaxCubeMapTextureSize;

var speedOfLight = 299792458;	// m/s
var gravitationalConstant = 6.6738480e-11;	// m^3 / (kg * s^2)

//1 = c m/s <-> c = s/m
var secondsPerMeter = speedOfLight;

//1 = G m^3 / (kg s^2) <=> G m^3 / (c m)^2 = kg <=> G/c^2 = kg/m
var kilogramsPerMeter = gravitationalConstant / (speedOfLight * speedOfLight);

//this does not zero accel beforehand!
var calcTidalForce;
(function(){
	x = [];
	n = [];

	calcTidalForce = function(accel, pos, srcPlanet) {
		n[0] = pos[0] - srcPlanet.pos[0];
		n[1] = pos[1] - srcPlanet.pos[1];
		n[2] = pos[2] - srcPlanet.pos[2];
		
		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = planets[planetIndex];
			if (planet.index === srcPlanet.index) continue;
			if (planet.mass === undefined) continue;
			
			x[0] = pos[0] - planet.pos[0];
			x[1] = pos[1] - planet.pos[1];
			x[2] = pos[2] - planet.pos[2];
			var xLength = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
			var xToTheSecond = xLength * xLength;
			var xToTheThird = xLength * xToTheSecond;
			var xToTheFourth = xLength * xToTheThird;
			var xToTheFifth = xLength * xToTheFourth;

			var xDotN = x[0] * n[0] + x[1] * n[1] + x[2] * n[2];
			
			// a^i = R^i_jkl n^j dx^k n^l
			// R^i_tjt = R^i_ttj = phi_,ij
			// but what if phi changes wrt time? then phi_,tt is nonzero, and how does our Riemann metric change?
			for (i = 0; i < 3; ++i) {
				accel[i] = gravitationalConstant * planet.mass * (3 * xDotN * x[i] / xToTheFifth - n[i] / xToTheThird);
			}
		}
	};
})();

/*
this does not zero accel beforehand!

geodesic calculation:

x''^u = -G^u_ab x'^a x'^b
*/
function calcGravitationForce(accel, pos) {
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		if (!planetInfluences[planetIndex]) continue;
		var planet = planets[planetIndex];
		if (planet.mass === undefined) continue;
		var x = pos[0] - planet.pos[0];
		var y = pos[1] - planet.pos[1];
		var z = pos[2] - planet.pos[2];
		xLength = Math.sqrt(x * x + y * y + z * z);
		var xToTheSecond = xLength * xLength;
		var xToTheThird = xLength * xToTheSecond;
		accel[0] -= x * (gravitationalConstant * planet.mass / xToTheThird);
		accel[1] -= y * (gravitationalConstant * planet.mass / xToTheThird);
		accel[2] -= z * (gravitationalConstant * planet.mass / xToTheThird);
	
/*
 		//this code is in natural units
		//so either convert xyz into natural units
		//or insert c's into here where they belong ...
		//also don't forget that G is in m^3 / (kg s^2) whereas our velocity is in m/day	
		
		x''^u = -conn^u_ab x'^a x'^b

		t'' = -rs / (r(r-rs)) t' r'
		r'' = -[c^2 rs (r - rs) / (2r^3) t'^2 - rs/(2r(r-rs)) r'^2 - (r-rs) (theta'^2 + sin^2theta phi'^2)]
	...and so on ...
		
		var r = Math.sqrt(x*x + y*y + z*z);
		var oneMinus2MOverR = 1 - 2*planet.mass/r;
		var posDotVel = x * oldVx + y * oldVy + z * oldVz;
		var velDotVel = oldVx * oldVx + oldVy * oldVy + oldVz * oldVz;
		var r2 = r * r;
		var invR2M = 1 / (r * oneMinus2MOverR);
		var rMinus2MOverR2 = oneMinus2MOverR / r;
		var MOverR2 = planet.mass / r2;
		accel[0] = -deltaLambda * MOverR2 * (rMinus2MOverR2 * x * oldVt * oldVt + invR2M * (x * velDotVel - 2 * oldVx * posDotVel));
		accel[1] = -deltaLambda * MOverR2 * (rMinus2MOverR2 * y * oldVt * oldVt + invR2M * (y * velDotVel - 2 * oldVy * posDotVel));
		accel[2] = -deltaLambda * MOverR2 * (rMinus2MOverR2 * z * oldVt * oldVt + invR2M * (z * velDotVel - 2 * oldVz * posDotVel));
		//accel[3] = deltaLambda * 2 * MOverR2 * invR2M * posDotVel * oldVt;
*/
	}
}

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

var defaultIntegrateTimeStep = 1/(24*60);	//1 min/frame
var integrateTimeStep = defaultIntegrateTimeStep;
var integrationPaused = true;
/*
local integrationMethod = 'rkf45'
local integrationArgs = {accuracy=1, iterations=100, norm=Planets.length}
*/
var integrationMethod = 'rk4';
var integrationArgs = undefined; 

function integrateFunction(time, planets) {
	var secondsPerJulianDay = 86400;
	// m^3 / (kg * julianday^2) = 6.6738480e-11 m^3 / (kg * s^2) * (second/julianday)^2
	var G = gravitationalConstant * secondsPerJulianDay * secondsPerJulianDay	// gravitationalConstant is proprotional to s^-2 and our time is in julian days, so convert
	// f'' = G m1 m2 / r^2 = rdir G m1 m2 / r^3
	// x1'' = rhat G m2 / r^3
	var deltaPlanets = new Planets();
	var accel = [];
	var delta = [];
	for (var i = 0; i < planets.length; ++i) {
		var pi = planets[i];
		accel[0] = accel[1] = accel[2] = 0;
		for (var j = 0; j < planets.length; ++j) {
			var pj = planets[j];
			if (i !== j && pj.mass !== undefined) {	//looks like we're not integrating against objects with no mass ...
				delta[0] = pj.pos[0] - pi.pos[0];
				delta[1] = pj.pos[1] - pi.pos[1];
				delta[2] = pj.pos[2] - pi.pos[2];
				var deltaLen = Math.sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
				var scalar = G * pj.mass / (deltaLen * deltaLen * deltaLen);
				accel[0] += delta[0] * scalar;
				accel[1] += delta[1] * scalar;
				accel[2] += delta[2] * scalar;
			}
		}
		deltaPlanets[i].pos[0] = pi.vel[0];
		deltaPlanets[i].pos[1] = pi.vel[1];
		deltaPlanets[i].pos[2] = pi.vel[2];
		deltaPlanets[i].vel[0] = accel[0];
		deltaPlanets[i].vel[1] = accel[1];
		deltaPlanets[i].vel[2] = accel[2];
	}
	return deltaPlanets;
}

function planetCartesianToSolarSystemBarycentric(destX, srcX, planet) {
	// right now planet angles aren't stored in the planet state (for adaptive integration's sake -- how to weight and combine time + space measurements)
	vec3TransformQuat(destX, srcX, planet.angle);
	
	// now rotate by axial tilt
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

var glutil;

//TODO use glutil.mouseDir?
function mouseRay() {
	var viewX = glutil.view.pos[0];
	var viewY = glutil.view.pos[1];
	var viewZ = glutil.view.pos[2];
	var viewFwdX = -2 * (glutil.view.angle[0] * glutil.view.angle[2] + glutil.view.angle[3] * glutil.view.angle[1]); 
	var viewFwdY = -2 * (glutil.view.angle[1] * glutil.view.angle[2] - glutil.view.angle[3] * glutil.view.angle[0]); 
	var viewFwdZ = -(1 - 2 * (glutil.view.angle[0] * glutil.view.angle[0] + glutil.view.angle[1] * glutil.view.angle[1])); 
	var viewRightX = 1 - 2 * (glutil.view.angle[1] * glutil.view.angle[1] + glutil.view.angle[2] * glutil.view.angle[2]); 
	var viewRightY = 2 * (glutil.view.angle[0] * glutil.view.angle[1] + glutil.view.angle[2] * glutil.view.angle[3]); 
	var viewRightZ = 2 * (glutil.view.angle[0] * glutil.view.angle[2] - glutil.view.angle[3] * glutil.view.angle[1]); 
	var viewUpX = 2 * (glutil.view.angle[0] * glutil.view.angle[1] - glutil.view.angle[3] * glutil.view.angle[2]);
	var viewUpY = 1 - 2 * (glutil.view.angle[0] * glutil.view.angle[0] + glutil.view.angle[2] * glutil.view.angle[2]);
	var viewUpZ = 2 * (glutil.view.angle[1] * glutil.view.angle[2] + glutil.view.angle[3] * glutil.view.angle[0]);
	var aspectRatio = canvas.width / canvas.height;
	var mxf = mouse.xf * 2 - 1;
	var myf = 1 - mouse.yf * 2;
	var tanFovY = Math.tan(glutil.view.fovY * Math.PI / 360);
	var mouseDirX = viewFwdX + tanFovY * (viewRightX * mxf * aspectRatio + viewUpX * myf);
	var mouseDirY = viewFwdY + tanFovY * (viewRightY * mxf * aspectRatio + viewUpY * myf);
	var mouseDirZ = viewFwdZ + tanFovY * (viewRightZ * mxf * aspectRatio + viewUpZ * myf);
	var mouseDirLength = Math.sqrt(mouseDirX * mouseDirX + mouseDirY * mouseDirY + mouseDirZ * mouseDirZ);
	return [mouseDirX/mouseDirLength, mouseDirY/mouseDirLength, mouseDirZ/mouseDirLength];
}

function chooseNewPlanet(mouseDir,doChoose) {
	var bestDot = 0;
	var bestPlanet = undefined;
	var bestPlanetIndex = undefined;
	var currentOrbitPlanet = planets[orbitPlanetIndex];
	for (var i = 0; i < planets.length; ++i) {
		var planet = planets[i];
		if (planet.hide) continue;
		var deltaX = planet.pos[0] - glutil.view.pos[0] - currentOrbitPlanet.pos[0];
		var deltaY = planet.pos[1] - glutil.view.pos[1] - currentOrbitPlanet.pos[1];
		var deltaZ = planet.pos[2] - glutil.view.pos[2] - currentOrbitPlanet.pos[2];
		var deltaLength = Math.sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
		deltaX /= deltaLength;
		deltaY /= deltaLength;
		deltaZ /= deltaLength;
		var dot = deltaX * mouseDir[0] + deltaY * mouseDir[1] + deltaZ * mouseDir[2];
		if (dot > bestDot) {
			bestDot = dot;
			bestPlanet = planet;
			bestPlanetIndex  = i;
		}
	}
	mouseOverPlanetIndex = undefined;
	if (bestPlanet !== undefined) {
		$('#hoverPlanetText').text(bestPlanet.name);
		mouseOverPlanetIndex = bestPlanetIndex;
		if (bestPlanet.index !== orbitPlanetIndex && doChoose) {
			vec3.sub(orbitOffset, planets[orbitPlanetIndex].pos, planets[bestPlanetIndex].pos);
			orbitPlanetIndex = bestPlanet.index;
			$('#orbitPlanetText').text(bestPlanet.name);
			refreshMeasureText();
		}
	}
}

function refreshMeasureText() {
	var planet = planets[orbitPlanetIndex];
	$('#measureMin').text(planet.measureMin === undefined ? '' : (planet.measureMin.toExponential() + ' m/s^2'));
	$('#measureMax').text(planet.measureMax === undefined ? '' : (planet.measureMax.toExponential() + ' m/s^2'));
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

var gl;
var canvas;

var hsvTex;
var colorShader;
var planetColorShader;
var planetTexShader;
var hsvShader;
var displayMethods = [
	'None',
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
var displayMethod = 'None';
var planetInfluences = [];

var showLinesToOtherPlanets = false;
var showVelocityVectors = false;
var velocityVectorScale = 30;
var showRotationAxis = false;
var showOrbitAxis = false;
var showLatAndLonLines = false;
var showGravityWell = false;
var showDistantPoints = true;
var showOrbits = true;
var showStars = false;
var gravityWellScaleNormalized = true;
var gravityWellScaleFixed = false;
var gravityWellScaleFixedValue = 2000;
var gravityWellRadialMinLog100 = -1;
var gravityWellRadialMaxLog100 = 2;

var heatAlpha = .5;
var colorBarHSVRange = 2/3;	// how much of the rainbow to use

var depthConstant = 1e-6;//2 / Math.log(1e+7 + 1);


var updatePlanetClassSceneObj;
(function(){
	var x = [];
	var accel = [];
	var norm = [];
	var proj = [];
	updatePlanetClassSceneObj = function(planet) {
		var planetClassPrototype = planet.init.prototype;

		//update tide attributes
		var showTide = displayMethod != 'None';
		if (!showTide) {
			if (planetClassPrototype.tex === undefined) {
				if (planetClassPrototype.sceneObj.shader !== planetColorShader) {
					planetClassPrototype.sceneObj.shader = planetColorShader;
					planetClassPrototype.sceneObj.texs = [];
				}
			} else {
				if (planetClassPrototype.sceneObj.shader !== planetTexShader) {
					planetClassPrototype.sceneObj.shader = planetTexShader;
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
					
					//var toTheMoon = (planets[planets.indexes.Moon].pos - x):normalize()
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
var gravityWellZScale = 1;
var gravityWellTargetZScale = 1;
var planetPointVisRatio = .001;
(function(){
	var delta = vec3.create();//[];
	var viewAngleInv = quat.create();
	var invRotMat = mat4.create();
	var viewPosInv = vec3.create();
	
	drawScene = function() {
		var orbitPlanet = planets[orbitPlanetIndex];
		
		mat4.identity(glutil.scene.mvMat);
		
		quat.conjugate(viewAngleInv, glutil.view.angle);
		mat4.fromQuat(invRotMat, viewAngleInv);
		mat4.multiply(glutil.scene.mvMat, glutil.scene.mvMat, invRotMat);

		//TODO pull from matrix
		var viewfwd = vec3.create();
		vec3.quatZAxis(viewfwd, glutil.view.angle);
		vec3.scale(viewfwd, viewfwd, -1);

		if (skyCubeObj) {
			gl.disable(gl.DEPTH_TEST);
			skyCubeObj.draw();
			gl.enable(gl.DEPTH_TEST);
		}
	
		vec3.scale(viewPosInv, glutil.view.pos, -1);
		mat4.translate(glutil.scene.mvMat, glutil.scene.mvMat, viewPosInv);

		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			var planet = planets[planetIndex];
			var planetClassPrototype = planet.init.prototype;
			if (planet.hide) continue;
			if (planet.pos === undefined) continue;	//comets don't have pos yet, but I'm working on that
			if (orbitPlanet.pos === undefined) continue;
			if (planet.lineObj === undefined) continue;

			if (showLinesToOtherPlanets) {
				if (orbitPlanet !== planet) {
					//while here, update lines
				
					vec3.sub(delta, planet.pos, orbitPlanet.pos);
					
					var dist = vec3.length(delta);
					vec3.scale(delta, delta, 1/dist);
					
					planet.lineObj.attrs.vertex.data[0] = delta[0] * orbitPlanet.radius;
					planet.lineObj.attrs.vertex.data[1] = delta[1] * orbitPlanet.radius;
					planet.lineObj.attrs.vertex.data[2] = delta[2] * orbitPlanet.radius;
					planet.lineObj.attrs.vertex.data[3] = delta[0] * orbitDistance;
					planet.lineObj.attrs.vertex.data[4] = delta[1] * orbitDistance;
					planet.lineObj.attrs.vertex.data[5] = delta[2] * orbitDistance;
					
					planet.lineObj.attrs.vertex.updateData();
					planet.lineObj.draw();
				}
			}

			if (showVelocityVectors) {
				vec3.sub(delta, planet.pos, orbitPlanet.pos);
				planet.lineObj.attrs.vertex.data[0] = delta[0];
				planet.lineObj.attrs.vertex.data[1] = delta[1];
				planet.lineObj.attrs.vertex.data[2] = delta[2];
				planet.lineObj.attrs.vertex.data[3] = delta[0] + planet.vel[0] * velocityVectorScale;
				planet.lineObj.attrs.vertex.data[4] = delta[1] + planet.vel[1] * velocityVectorScale;
				planet.lineObj.attrs.vertex.data[5] = delta[2] + planet.vel[2] * velocityVectorScale;
				planet.lineObj.attrs.vertex.updateData();
				planet.lineObj.draw();
			}
		
			if (showRotationAxis) {
				vec3.sub(delta, planet.pos, orbitPlanet.pos);
				var axis = [0,0,1];
				if (planet.tiltAngle) {
					vec3.quatZAxis(axis, planet.tiltAngle);
				}
				planet.lineObj.attrs.vertex.data[0] = delta[0] + axis[0] * 2 * planet.radius;
				planet.lineObj.attrs.vertex.data[1] = delta[1] + axis[1] * 2 * planet.radius;
				planet.lineObj.attrs.vertex.data[2] = delta[2] + axis[2] * 2 * planet.radius;
				planet.lineObj.attrs.vertex.data[3] = delta[0] + axis[0] * -2 * planet.radius;
				planet.lineObj.attrs.vertex.data[4] = delta[1] + axis[1] * -2 * planet.radius;
				planet.lineObj.attrs.vertex.data[5] = delta[2] + axis[2] * -2 * planet.radius;
				planet.lineObj.attrs.vertex.updateData();
				planet.lineObj.draw();
			}
			
			if (showOrbitAxis) {
				vec3.sub(delta, planet.pos, orbitPlanet.pos);
				planet.lineObj.attrs.vertex.data[0] = delta[0] + planet.orbitAxis[0] * 2 * planet.radius;
				planet.lineObj.attrs.vertex.data[1] = delta[1] + planet.orbitAxis[1] * 2 * planet.radius;
				planet.lineObj.attrs.vertex.data[2] = delta[2] + planet.orbitAxis[2] * 2 * planet.radius;
				planet.lineObj.attrs.vertex.data[3] = delta[0] + planet.orbitAxis[0] * -2 * planet.radius;
				planet.lineObj.attrs.vertex.data[4] = delta[1] + planet.orbitAxis[1] * -2 * planet.radius;
				planet.lineObj.attrs.vertex.data[5] = delta[2] + planet.orbitAxis[2] * -2 * planet.radius;
				planet.lineObj.attrs.vertex.updateData();
				planet.lineObj.draw();
			}
	
		}

		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			var planet = planets[planetIndex];
			var planetClassPrototype = planet.init.prototype;
			if (planet.hide) continue;
			
			//update vis ratio
			var dx = planet.pos[0] - glutil.view.pos[0] - orbitPlanet.pos[0];
			var dy = planet.pos[1] - glutil.view.pos[1] - orbitPlanet.pos[1];
			var dz = planet.pos[2] - glutil.view.pos[2] - orbitPlanet.pos[2];
			planet.visRatio = planet.radius / Math.sqrt(dx * dx + dy * dy + dz * dz); 

			//some planets have no radii ... so no object
			if (planet.sceneObj !== undefined) {
			
				//update scene object
				vec3.sub(planet.sceneObj.uniforms.pos, planet.pos, orbitPlanet.pos);
				if (planet.tiltAngle) {
					quat.multiply(planet.sceneObj.uniforms.angle, planet.tiltAngle, planet.angle);
				} else {
					quat.copy(planet.sceneObj.uniforms.angle, planet.angle);
				}
				
				//only allow ring obj if there's a radii and a sphere 
				//... this way i can just copy over the pos and angle	
				if (planet.ringObj !== undefined) {
					vec3.sub(planet.ringObj.pos, planet.pos, orbitPlanet.pos);
					quat.copy(planet.ringObj.angle, planet.sceneObj.uniforms.angle);	
				}
			}
		}

		//draw sphere planets
		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			var planet = planets[planetIndex];
			var planetClassPrototype = planet.init.prototype;
			if (planet.hide) continue;

			if (planet.sceneObj && (planet.visRatio >= planetPointVisRatio)) {
				updatePlanetClassSceneObj(planet);
			
				var Sun = planets[planets.indexes.Sun];
				vec3.sub(planet.sceneObj.uniforms.sunPos, Sun.pos, planet.pos);
				
				planet.sceneObj.draw();
			
				if (showLatAndLonLines) {
					vec3.copy(planet.latLonObj.pos, planet.sceneObj.uniforms.pos);
					quat.copy(planet.latLonObj.angle, planet.sceneObj.uniforms.angle);
					planet.latLonObj.draw();
				}
			
				if (planet.ringObj !== undefined) {
					gl.disable(gl.CULL_FACE);
					planet.ringObj.draw();
					gl.enable(gl.CULL_FACE);
				}
			}
		}
		
		//draw point planets
		// goes slow when comets are included
		//TODO make a buffer with a vertex per-planet rather than changing a single vertex
		for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
			var planet = planets[planetIndex];
			var planetClassPrototype = planet.init.prototype;
			if (planet.hide) continue;
			
			if ((!planet.sceneObj || planet.visRatio < planetPointVisRatio) && showDistantPoints) {
				vec3.sub(pointObj.attrs.vertex.data, planet.pos, orbitPlanet.pos);
				pointObj.attrs.vertex.updateData();
				pointObj.draw({
					uniforms : {
						color : planetClassPrototype.color,
						pointSize : 4
					}
				});
			}
		}

		if (mouseOverPlanetIndex !== undefined) {
			var planet = planets[mouseOverPlanetIndex];
			if (planet !== undefined) {
				gl.disable(gl.DEPTH_TEST);
				pointObj.attrs.vertex.data[0] = planet.pos[0] - orbitPlanet.pos[0];
				pointObj.attrs.vertex.data[1] = planet.pos[1] - orbitPlanet.pos[1];
				pointObj.attrs.vertex.data[2] = planet.pos[2] - orbitPlanet.pos[2];
				pointObj.attrs.vertex.updateData();
				pointObj.draw({
					uniforms : {
						color : [0,1,0,1],
						pointSize : 8
					},
					blend : [gl.SRC_ALPHA, gl.ONE]
				});
				gl.enable(gl.DEPTH_TEST);
			}
		}

		if (showOrbits) {
			window.orbitPathsDrawn = 0;
			for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
				var planet = planets[planetIndex];
				var planetClassPrototype = planet.init.prototype;
				if (planet.hide) continue;
				
				if (planet.orbitPathObj) {

					var semiMajorAxis = planet.keplerianOrbitalElements.semiMajorAxis;
					var eccentricity = planet.keplerianOrbitalElements.eccentricity;
					var distPeriapsis = semiMajorAxis * (1 + eccentricity);	//largest distance from the parent planet
				
					//vector from view to parent planet
					var parentPlanet = planets[planet.parent];
					vec3.sub(delta, parentPlanet.pos, orbitPlanet.pos);
					vec3.sub(delta, delta, glutil.view.pos);
					var deltaLength = vec3.length(delta); 

					var visRatio = distPeriapsis / deltaLength;

					var cosAngleToParent = vec3.dot(delta, viewfwd) / deltaLength;
					var angleToParent = Math.acos(Math.clamp(cosAngleToParent, -1, 1));

					//TODO cone/sphere intersection.  sphere center is delta, radius is distPeriapsis, cone center is origin, angle is angle
					var distToParentPlanet = vec3.dot(delta, viewfwd);
					//if (distToParentPlanet > 0) 
					{
						var tanAngleOfParentInScreen = distPeriapsis / distToParentPlanet;
						var angleOfParentInScreen = Math.abs(Math.atan(tanAngleOfParentInScreen));
						//if (distToParentPlanet > 0 && (angleToParent < glutil.view.fovY * Math.PI / 180 + angleOfParentInScreen) ||
						//	deltaLength < distPeriapsis)
						{
						
							//if sphere around parent planet of size max orbit size 
							// intersects with the view frustum then don't draw it
							//test whether the view pos is in the sphere, or if the view dir dot delta
							//if (deltaLength < distPeriapsis ||	//if we're within the periapsis
							//	cosAngle > 0)
							{
								//recenter around orbitting planet
								vec3.sub(planet.orbitPathObj.pos, planet.pos, orbitPlanet.pos);
								//global positions
								//less accurate, but can store all as one buffer
								//you see the orbit paths break down as far out as Saturn 
								//vec3.scale(planet.orbitPathObj.pos, orbitPlanet.pos, -1);
								planet.orbitPathObj.draw();
								++orbitPathsDrawn;
							}
						}
					}
				}
			}
		}
		
		if (showGravityWell) {
			//do transformation math in double
			//TODO just give gl-matrix a type param in its init
			var glMvMat = [];
			var viewAngleInvd = [];
			quat.conjugate(viewAngleInvd, glutil.view.angle);
			quat.normalize(viewAngleInvd, viewAngleInvd);	//normalize in double precision
			mat4.fromQuat(glMvMat, viewAngleInvd);
			var viewPosInvd = [];
			vec3.negate(viewPosInvd, glutil.view.pos);
			mat4.translate(glMvMat, glMvMat, viewPosInvd);
			
			var orbitBasis = [];
			var orbitBasisInv = [];
			var mvMat = [];
			for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
				var planet = planets[planetIndex];
				var planetClassPrototype = planet.init.prototype;
				if (planet.hide) continue;
				if (planet.radius === undefined) continue;	
				//max radial dist is R * Math.pow(100, gravityWellRadialMaxLog100)
				if (planet.visRatio * Math.pow(100, gravityWellRadialMaxLog100) < .001) continue;
				
				mat4.identity(orbitBasis);
				mat4.identity(orbitBasisInv);
				for (var i = 0; i < 16; ++i) {
					mvMat[i] = glMvMat[i];
				}
				for (var i = 0; i < 3; ++i) {
					orbitBasis[i] = planetClassPrototype.orbitBasis[0][i];
					orbitBasis[4+i] = planetClassPrototype.orbitBasis[1][i];
					orbitBasis[8+i] = planetClassPrototype.orbitBasis[2][i];
				}
				//gravWellObjBasis = orbitBasis * gravityWellZScale * zOffsetByRadius * orbitBasis^-1 * planetPosBasis
			
				//calc this for non-sceneObj planets.  maybe I should store it as a member variable?
				var relPos = [
					planet.pos[0] - orbitPlanet.pos[0],
					planet.pos[1] - orbitPlanet.pos[1],
					planet.pos[2] - orbitPlanet.pos[2]
				];
				
				//translate to planet center
				mat4.translate(mvMat, mvMat, relPos);
				
				//align gravity well with orbit plane
				mat4.multiply(mvMat, mvMat, orbitBasis);
				
				//align base gravity well with base of planet
				mat4.translate(mvMat, mvMat, [0,0,-planet.radius]);
				
				//scale for visual effect
				//too flat for earth, too sharp for sun			
				if (gravityWellScaleFixed) {
					gravityWellTargetZScale = gravityWellScaleFixedValue;
				} else if (gravityWellScaleNormalized) {
					//causes larger wells to be much smaller and sharper ...
					//var gravityWellTargetZScale = 2000000 * Math.log(1 + z);
					//normalized visually per-planet.  scale is not 1-1 across planets
					gravityWellTargetZScale = 1 / orbitPlanet.gravityWellScalar;
				}
				gravityWellZScale += .01 * (gravityWellTargetZScale - gravityWellZScale); 
				if (gravityWellZScale !== gravityWellZScale) gravityWellZScale = 1;
				mat4.scale(mvMat, mvMat, [1,1,gravityWellZScale * planet.gravityWellScalar]);

				//apply orbit basis inverse
				mat4.transpose(orbitBasisInv, orbitBasis);
				mat4.multiply(mvMat, mvMat, orbitBasisInv);
			
				for (var i = 0; i < 16; ++i) {
					planet.gravWellObj.uniforms.mvMat[i] = mvMat[i];
				}

				planet.gravWellObj.draw();
			}
		}
	
		if (showStars) {
			starsObj.draw();
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
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('height', window.innerHeight);
	});
	
	glutil.resize();

	//glutil.view.fovY = canvas.height / canvas.width * 90;

	/*
	var aspectRatio = canvas.width / canvas.height;
	var nearHeight = Math.tan(glutil.view.fovY * Math.PI / 360);
	var nearWidth = aspectRatio * nearHeight;
	*/
	
	/** /
	//setup infinite projection matrix
	//http://www.terathon.com/gdc07_lengyel.pdf
	//http://www.gamasutra.com/view/feature/131351/the_mechanics_of_robust_stencil_.php?page=2
	var epsilon = 0;
	mat4.identity(glutil.scene.projMat);
	glutil.scene.projMat[0] = glutil.view.zNear / nearWidth;
	glutil.scene.projMat[5] = glutil.view.zNear / nearHeight;
	glutil.scene.projMat[10] = epsilon-1;
	glutil.scene.projMat[11] = (epsilon - 2) * glutil.view.zNear;
	glutil.scene.projMat[14] = -1;
	glutil.scene.projMat[15] = 0;
	/**/

	/* or we could just linearly ramp the depth over the range we want, since the rest is going to full anyways * /
	//(-z - zNear) / (zFar - zNear) * 2 - 1
	//z * (-2 / (zFar - zNear)) + ( - 2 * zNear / (zFar - zNear) - 1)
	var resZNear = 1e+4;
	var resZFar = 1e+8;
	mat4.identity(glutil.scene.projMat);
	glutil.scene.projMat[0] = glutil.view.zNear / nearWidth;
	glutil.scene.projMat[5] = glutil.view.zNear / nearHeight;
	glutil.scene.projMat[10] = -2 * (resZFar - resZNear);
	glutil.scene.projMat[11] = -(1 + 2 * resZNear / (resZFar - resZNear));
	glutil.scene.projMat[14] = -1;
	glutil.scene.projMat[15] = 0;
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

var unitsPerM = [
	{name:'Mpc',	value:1000000*648000/Math.PI*149597870700},
	{name:'lyr', 	value:9460730472580800},
	{name:'AU', 	value:149597870700},
	{name:'km',		value:1000},
	{name:'m',		value:1}
];
function refreshOrbitTargetDistanceText() {
	var units = 'm';
	var dist = orbitTargetDistance; 
	for (var i = 0; i < unitsPerM.length; ++i) {
		var ratio = orbitTargetDistance / unitsPerM[i].value;
		if (ratio > .1 || i == unitsPerM[i].length-1) {
			units = unitsPerM[i].name;
			dist = ratio;
			break;
		}
	}
	$('#orbitTargetDistance').text(dist.toFixed(4)+' '+units);
}

function refreshCurrentTimeText() {
	$('#currentTimeText').text(astro.julianToCalendar(julianDate));
}

var primaryPlanetHorizonIDs = [10, 199, 299, 301, 399, 499, 599, 699, 799, 899, 999];

var ModifiedDepthShaderProgram;

$(document).ready(init1);

var slideDuration = 500;
var slideWidth = 300;
var currentOpenSidePanelID = undefined;
var allSidePanelIDs = [];

function showSidePanel(sidePanelID) {
	if (currentOpenSidePanelID !== undefined) {
		hideSidePanel(currentOpenSidePanelID, true);
	}
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('z-index', 1);
	});
	var sidePanel = $('#'+sidePanelID);
	sidePanel.css('z-index', 2);
	sidePanel.show();
	$('#menu')
	.animate({left:slideWidth}, {duration:slideDuration})
	.attr('src', 'reverse.png');
	sidePanel.animate({left:0}, {duration:slideDuration});
	currentOpenSidePanelID = sidePanelID;
}

function hideSidePanel(sidePanelID, dontMoveOpenButton) {
	if (sidePanelID === currentOpenSidePanelID) currentOpenSidePanelID = undefined;
	var sidePanel = $('#'+sidePanelID);
	if (!dontMoveOpenButton) {
		$('#menu')
		.animate({left:0}, {duration:slideDuration})
		.attr('src', 'play.png');
	}
	sidePanel.animate({left:-slideWidth}, {
		duration:slideDuration,
		complete:function() {
			sidePanel.hide();
		}
	});
}

function calendarToJulian(d) {
	var year = d.getUTCFullYear();
	var month = d.getUTCMonth();
	var day = d.getUTCDate();
	var hour = d.getUTCHours();
	var min = d.getUTCMinutes();
	var sec = d.getUTCSeconds();
//http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
//testing against results in Wikipedia page: http://en.wikipedia.org/wiki/Julian_day
// (notice I couldn't recreate these results with the algorithm on the same page)
//for UTC date 2000 jan 1 12:00:00 this gives 2451545 - right on
//for UTC date 2013 jan 1 00:30:00 this gives 2456294.0208333335 - right on
	var oneBasedMonth = 1 + month;
	var astronomicalYear = year < 0 ? year + 1 : year;
	if (oneBasedMonth <= 2) {	//jan or feb
		--astronomicalYear;
		oneBasedMonth += 12;
	}
	var a = Math.floor(astronomicalYear / 100);
	var b = Math.floor(a / 4);
	var c = 2 - a + b;
	var e = Math.floor(365.25 * (astronomicalYear + 4716));
	var f = Math.floor(30.6001 * (oneBasedMonth + 1));
	var jdn = c + day + e + f - 1524.5;
	jdn += (hour + (min + sec / 60) / 60) / 24;
	return jdn;
}

function init1() {
	allSidePanelIDs = [
		'mainSidePanel',
		'overlaySidePanel',
		'displayOptionsSidePanel',
		'celestialBodiesSidePanel'
	];
	mainSidePanel = $('#mainSidePanel');
	
	//keep track of what menu is open
	//if none are open, open the main menu
	//if any are ... if it's not main 
	$('#menu').click(function() {
		if (currentOpenSidePanelID === undefined) {
			showSidePanel('mainSidePanel');
		} else if (currentOpenSidePanelID == 'mainSidePanel') {
			hideSidePanel(currentOpenSidePanelID);
		} else {
			showSidePanel('mainSidePanel');
		}
	});

	$.each([
		{buttonID:'mainButtonOverlay', divID:'overlaySidePanel'},
		{buttonID:'mainButtonDisplayOptions', divID:'displayOptionsSidePanel'},
		{buttonID:'mainButtonCelestialBodies', divID:'celestialBodiesSidePanel'},
	], function(i, info) {
		$('#'+info.buttonID).click(function() {
			showSidePanel(info.divID);
		});
	});
	
	$('#reset').click(function() {
		integrationPaused = true;
		integrateTimeStep = defaultIntegrateTimeStep; 
		planets = initPlanets;
		julianDate = initJulianDate;
		refreshCurrentTimeText();
	});
	$('#play').click(function() {
		integrationPaused = false;
		integrateTimeStep = Math.abs(integrateTimeStep);
	});
	$('#reverse').click(function() {
		integrationPaused = false;
		integrateTimeStep = -Math.abs(integrateTimeStep);
	});
	$('#pause').click(function() {
		integrationPaused = true;
	});
	//fast/slow vs ffwd/rewind?
	$('#ffwd').click(function() {
		integrateTimeStep *= 2;
	});
	$('#rewind').click(function() {
		integrateTimeStep /= 2;
	});
	
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
		glutil = new GLUtil({canvas:canvas});
		gl = glutil.context;
	} catch (e) {
		$('#menu').remove();
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}

	ModifiedDepthShaderProgram = makeClass({
		super : glutil.ShaderProgram,
		init : function(args) {
			
			// maximizing depth range: http://outerra.blogspot.com/2012/11/maximizing-depth-buffer-range-and.html
			args.vertexCode = mlstr(function(){/*
	uniform float zNear, zFar, depthConstant;
	float depthfunction(vec4 v) {
		//return (log(v.w + 1.) * depthConstant - 1.) * v.w;
		//return (2.0 * log(v.w / zNear) / log(zFar / zNear) - 1.) * v.w; 
		return v.z;
	}
	*/}) + (args.vertexCode || '');
			if (args.uniforms === undefined) args.uniforms = {};
			args.uniforms.zNear = glutil.view.zNear;
			args.uniforms.zFar = glutil.view.zFar;
			args.uniforms.depthConstant = depthConstant;
			args.vertexPrecision = 'best';
			args.fragmentPrecision = 'best';
			ModifiedDepthShaderProgram.super.call(this, args);
		}
	});

	glMaxCubeMapTextureSize = gl.getParameter(gl.MAX_CUBE_MAP_TEXTURE_SIZE);

	/*glutil.view.angle[0] = -0.4693271591372717;
	glutil.view.angle[1] = 0.7157221264895661;
	glutil.view.angle[2] = -0.4298784661116332;
	glutil.view.angle[3] = 0.28753844912098436;*/
	glutil.view.zNear = 1e+3;
	glutil.view.zFar = 1e+25;
	
	
	// overlay side panel 


	var overlaySidePanel = $('#overlaySidePanel');
	$('<span>', {text:'Overlay:'}).appendTo(overlaySidePanel);
	$('<br>').appendTo(overlaySidePanel);
	$.each(displayMethods, function(displayMethodIndex,thisDisplayMethod) {
		var radio = $('<input>', {
			type : 'radio',
			name : 'displayMethods',
			value : displayMethodIndex,
			click : function() {
				displayMethod = thisDisplayMethod;
				invalidateForces();
			}
		})
			.attr('name', 'display')
			.appendTo(overlaySidePanel);
		if (thisDisplayMethod == displayMethod) radio.attr('checked', 'checked');
		$('<span>', {text:thisDisplayMethod}).appendTo(overlaySidePanel);
		$('<br>').appendTo(overlaySidePanel);
	});
	$('<br>').appendTo(overlaySidePanel);
	$('<span>', {text:'Influencing Planets:'}).appendTo(overlaySidePanel);
	$('<br>').appendTo(overlaySidePanel);

	//add radio buttons hierarchically ...
	var overlayControlsForPlanets = {};

	var HierarchicalCheckboxControl = makeClass({
		/*
		args:
			title	<- used to identify this checkbox
			change		<- callback upon changing the checkbox value
			isChecked
			... and anything else that can be referenced through this.args
		*/
		init : function(args) {
			this.args = args;
			this.childControls = [];
			this.div = $('<div>', {
				css : {paddingLeft:'5px'}
			});
	
			var thiz = this;
			this.checkbox = $('<input>', {
				type : 'checkbox',
				change : function() {
					//refresh all parent controls' moon checkboxes -- to whiteout or grey them
					for (var c = thiz.parentControls; c; c = c.parentControls) {
						c.recomputeMoonCheckbox();
					}
					
					args.change.call(thiz);
				}
			})
				.prop('checked', args.isChecked)
				.appendTo(this.div);
			
			$('<span>', {
				text : args.title
			}).appendTo(this.div);
		
			this.toggleChildDiv = $('<span>', {
				css : {
					paddingLeft : '10px',
					cursor : 'pointer'
				}
			}).appendTo(this.div);

			this.moonCheckbox = $('<input>', {
				type : 'checkbox',
				change : function() {
					//select / deselect all children
					thiz.setAllChildren($(this).prop('checked'));
				}
			})
				.prop('checked', 1)
				.appendTo(this.toggleChildDiv);
	
			$('<span>', {
				css : {
					cursor : 'pointer'
				},
				text : '...',
				click : function() {
					if (thiz.childDiv.css('display') == 'none') {
						thiz.childDiv.show();
					} else {
						thiz.childDiv.hide();
					}
				}	
			}).appendTo(this.toggleChildDiv);
	
			$('<br>').appendTo(this.div);
		
			this.childDiv = $('<div>').appendTo(this.div);

		},
		addChild : function(childControl) {
			childControl.div.appendTo(this.childDiv);
			childControl.parentControls = this;
			this.childControls.push(childControl);
		},
		setAllChildren : function(checked) {
			for (var i = 0; i < this.childControls.length; ++i) {
				var ch = this.childControls[i].checkbox;
				if (checked) {
					if (!ch.prop('checked')) { ch.prop('checked', 1); ch.trigger('change'); }
				} else {
					if (ch.prop('checked')) { ch.prop('checked', 0); ch.trigger('change'); }
				}
				this.childControls[i].setAllChildren(checked);
			}
		},
		recomputeMoonCheckbox : function() {
			var numChecked = 0;
			var total = 0;
			for (var i = 0; i < this.childControls.length; ++i) {
				++total;
				//check the child
				if (this.childControls[i].checkbox.prop('checked')) {
					++numChecked;
				}
				//check the child's children if they exist
				if (this.childControls[i].childControls.length > 0) {
					if (this.childControls[i].moonCheckbox.prop('checked')) {
						if (this.childControls[i].moonCheckbox.prop('indeterminate')) {
							numChecked += .5;
						} else {
							++numChecked;
						}
					}
					++total;
				}
			}
			if (numChecked == 0) {
				this.moonCheckbox.prop('checked', 0)
					.prop('indeterminate', 0);
			} else if (numChecked == total) {
				this.moonCheckbox.prop('checked', 1)
					.prop('indeterminate', 0);
			} else {
				this.moonCheckbox.prop('checked', 1)
					.prop('indeterminate', 1);
			}
		}
	});

	$.each(Planets.prototype.planetClasses, function(planetIndex,planetClass) {
		//if any other planet doesn't have recorded mass then skip it
		if (planetClass.prototype.mass === undefined) return;

		var parentPlanetIndex = planetClass.prototype.parent;
		if (parentPlanetIndex !== undefined && parentPlanetIndex >= planetIndex) throw "parent index should be < planet index or undefined";

		var controls = new HierarchicalCheckboxControl({
			title : planetClass.prototype.name,
			isChecked : true,
			change : function() {			
				planetInfluences[this.args.planetIndex] = this.checkbox.is(':checked');
				invalidateForces();
			},
			planetIndex : planetIndex
		});

		//add to parent or to the side panel
		if (parentPlanetIndex === undefined) {
			controls.div.appendTo(overlaySidePanel);
		} else {
			overlayControlsForPlanets[parentPlanetIndex].addChild(controls);
		}
		
		planetInfluences[planetIndex] = true;
		
		overlayControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
	});

	for (var planetIndex in overlayControlsForPlanets) {
		var planetIndex = +planetIndex;
		var controls = overlayControlsForPlanets[planetIndex];
	
		controls.recomputeMoonCheckbox();

		//if a planet didn't get any children, hide its 'toggle child div'
		if (controls.childDiv.children().length == 0) {
			controls.toggleChildDiv.hide();
			continue;
		}

		//if a planet got children and it's not the sun or earth then hide by default
		if (planetIndex !== Planets.prototype.indexes.Sun && 
			planetIndex !== Planets.prototype.indexes.Earth)
		{
			controls.childDiv.hide();
		}
	}
	
	$('<br>').appendTo(overlaySidePanel);
	
	
	// display options side panel
	

	$.each([
		'showLinesToOtherPlanets',
		'showVelocityVectors',
		'showRotationAxis',
		'showOrbitAxis',
		'showLatAndLonLines',
		'showGravityWell',
		'showOrbits',
		'showDistantPoints',
		'gravityWellScaleNormalized',
		'gravityWellScaleFixed'
	], function(_, toggle) {
		(function(){
			var checkbox = $('#'+toggle);
			if (window[toggle]) checkbox.attr('checked', 'checked');
			checkbox.change(function() {
				window[toggle] = checkbox.is(':checked');
				//TODO generalize for radius?
				if (toggle == 'gravityWellScaleNormalized' && gravityWellScaleNormalized) gravityWellScaleFixed = false;
				if (toggle == 'gravityWellScaleFixed' && gravityWellScaleFixed) gravityWellScaleNormalized = false;
			});
		})();
	});

	$.each([
		'gravityWellScaleFixedValue',
	], function(_, toggle) {
		(function(){
			var textfield = $('#'+toggle);
			textfield.val(window[toggle]);
			textfield.change(function() {
				window[toggle] = textfield.val();
			});
		})();
	});


	// celestial bodies side panel


	var celestialBodiesSidePanel = $('#celestialBodiesSidePanel');

	//add radio buttons hierarchically ...
	var celestialBodiesControlsForPlanets = {};

	var cometParent;

	$.each(Planets.prototype.planetClasses, function(planetIndex,planetClass) {

		var parentPlanetIndex = planetClass.prototype.parent;
		if (parentPlanetIndex !== undefined && parentPlanetIndex >= planetIndex) throw "parent index should be < planet index or undefined";

		var controls = new HierarchicalCheckboxControl({
			title : planetClass.prototype.name,
			isChecked : !planetClass.prototype.hide,
			change : function() {
				Planets.prototype.planetClasses[this.args.planetIndex].prototype.hide = !this.checkbox.is(':checked');		
			},
			planetIndex : planetIndex
		});
	
		if (parentPlanetIndex === undefined) {
			controls.div.appendTo($('#celestialBodiesVisibleBodies'));
		} else {
			celestialBodiesControlsForPlanets[parentPlanetIndex].addChild(controls);
		}

		celestialBodiesControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
	});

	if (cometParent) cometParent.recomputeMoonCheckbox();
	for (var planetIndex in celestialBodiesControlsForPlanets) {
		var planetIndex = +planetIndex;
		var controls = celestialBodiesControlsForPlanets[planetIndex];
		
		controls.recomputeMoonCheckbox();
		
		//if a planet didn't get any children, hide its 'toggle child div'
		if (controls.childDiv.children().length == 0) {
			controls.toggleChildDiv.hide();
			continue;
		}

		//if a planet got children and it's not the sun or earth then hide by default
		if (planetIndex !== Planets.prototype.indexes.Sun && 
			planetIndex !== Planets.prototype.indexes.Earth)
		{
			controls.childDiv.hide();
		}
	}
	
	$('<br>').appendTo($('#celestialBodiesVisibleBodies'));

	//these are added to the end of the result
	//they should get greyed upon new query (search, prev, next click)
	//they should be regenerated upon new content
	//they should be re-enabled upon error
	var nextButton = undefined;
	var prevButton = undefined;

	var searchResults = [];

	var processSearch = function(pageIndex) {
		var button = $('#celestialBodiesSearch');
		var searchText = $('#celestialBodiesSearchText');
		var searchStr = searchText.val();
		button.prop('disabled', 1);
		searchText.val('searching...');
		searchText.prop('disabled', 1);
			
		if (prevButton) prevButton.prop('disabled', 1);
		if (nextButton) nextButton.prop('disabled', 1);
					
		$.ajax({
			url : '/solarsystem/jpl-ssd-smallbody/search.lua',
			dataType : 'json',
			data : {
				comet : $('#celestialBodiesSearchComets').prop('checked')?1:0,
				numbered : $('#celestialBodiesSearchNumbered').prop('checked')?1:0,
				unnumbered : $('#celestialBodiesSearchUnnumbered').prop('checked')?1:0,
				text : searchStr,
				page : pageIndex+1	//1-based
			},
			timeout : 5000
		}).error(function() {
			searchText.val(searchStr);
			searchText.prop('disabled', 0);
			button.prop('disabled', 0);
			//TODO animate background color of search text
			if (prevButton) prevButton.prop('disabled', 0);
			if (nextButton) nextButton.prop('disabled', 0);
			
			var warning = $('<div>', {text:'Connection Failed!', css:{color:'red'}});
			$('#celestialBodiesSearchComets').before(warning);
			setTimeout(function() {
				//after five seconds, fade away
				warning.animate({
					height : 0,
					opacity : 0
				}, {
					duration : 500,
					complete : function() {
						warning.remove();
					}
				});
			}, 3000);
		}).done(function(results) {
			searchText.val(searchStr);
			searchText.prop('disabled', 0);
			button.prop('disabled', 0);
			
			var pageSize = 20;	//fixed atm
			var pageMax = Math.floor((results.count-1) / pageSize);

			searchResults = [];
			
			var resultsDiv = $('#celestialBodiesSearchResults');
			resultsDiv.empty();
			$.each(results.rows, function(i,row) {
				var rowDiv = $('<div>');
				rowDiv.appendTo(resultsDiv);

				var name = row.name;
				if (row.idNumber) {
					name = row.idNumber+'/'+name;
				}
				
				var checkbox = $('<input>', {
					type : 'checkbox',
					change : function() {
					
						if (!$(this).is(':checked')) {	//uncheck checkbox => remove planet

							//only remove if it's already there 
							if (planets.indexes[name] !== undefined) {
								var index = planets.indexes[name];
								planets.remove(index);
								initPlanets.remove(index);
								Planets.prototype.planetClasses.splice(index, 1);
								//now remap indexes
								for (var i = index; i < Planets.prototype.planetClasses.length; ++i) {
									Planets.prototype.planetClasses[index].prototype.index = i;
								}
								if (orbitPlanetIndex == index) orbitPlanetIndex = planets.indexes.Sun;
								if (orbitPlanetIndex > index) --orbitPlanetIndex;
								//TODO destruct WebGL geometry?  or is it gc'd automatically?
								//now rebuild indexes
								Planets.prototype.indexes = {};
								for (var i = 0; i < Planets.prototype.planetClasses.length; ++i) {
									Planets.prototype.indexes[Planets.prototype.planetClasses[i].prototype.name] = i;
								}
							}

						} else {	//check checkbox => add planet
						
							//only add if it's not there
							if (planets.indexes[name] === undefined) {	

								//add the row to the bodies
							
								var index = Planets.prototype.planetClasses.length;
								var newPlanetClass = makeClass({
									super : Planet,
									name : name,
									isComet : row.bodyType == 'comet',
									isAsteroid : row.bodyType == 'numbered asteroid' || row.bodyType == 'unnumbered asteroid',
									orbitData : row,
									parent : planets.indexes.Sun,
									index : index
								});
								
								Planets.prototype.planetClasses.push(newPlanetClass);
								
								//instanciate it
								var planet = new newPlanetClass();
								planets[index] = planet;
								planets.length = index+1;
								
								//add copy to initPlanets for when we get reset() working
								initPlanets[index] = planet.clone();
								initPlanets.length = planets.length;

								planets.indexes[newPlanetClass.prototype.name] = index;

								initPlanetColorSchRadiusAngle(planet);
								initPlanetSceneLatLonLineObjs(planet);
								initPlanetOrbitPathObj(planet);

							}
						}
					}
				})
					.prop('checked', planets.indexes[name] !== undefined)
					.appendTo(rowDiv);

				$('<span>', {text:name}).appendTo(rowDiv);
				//TODO put an 'add' button next to each
				//on clicking it, add the body to the planet list 
				//and repopulate the 'extra' div
			
				searchResults.push({
					checkbox : checkbox,
					div : rowDiv,
					data : row,
					name : name
				});
				
			});

			//and add prev/next/pages if there is 
			if (pageMax > 0) {
				var changePage = function(dir) {
					//TODO remove or grey out results as well?
					if (nextButton) nextButton.prop('disabled', 1);
					if (prevButton) prevButton.prop('disabled', 1);
					processSearch(pageIndex+dir);
				};
				if (pageIndex > 0) {
					prevButton = $('<button>', {
						text : 'Prev',
						click : function() {
							changePage(-1);
						}
					}).appendTo($('#celestialBodiesSearchResults'));
				}
				if (pageIndex < pageMax) {
					nextButton = $('<button>', {
						text : 'Next',
						click : function() {
							changePage(1);
						}
					}).appendTo($('#celestialBodiesSearchResults'));
				}
				$('<span>', {text:(pageIndex+1)+' of '+pageMax}).appendTo($('#celestialBodiesSearchResults'));
			}
		});
	};
	
	$('#celestialBodiesSearchText').keydown(function(e){
		if (e.keyCode == 13) {
			$('#celestialBodiesSearch').trigger('click');
		}
	});

	//change a check box, immediately update search results
	$.each([$('#celestialBodiesSearchComets'), $('#celestialBodiesSearchNumbered'), $('#celestialBodiesSearchUnnumbered')], function(_,checkbox) {
		checkbox.change(function() {
			$('#celestialBodiesSearch').trigger('click');
		});
	});

	$('#celestialBodiesSearchToggleAll').click(function() {
		var checked = $(this).is(':checked');
		$.each(searchResults, function(i,result) {
			if (result.checkbox.prop('checked') != checked) {
				result.checkbox.trigger('click');
			}
		});
	});

	$('#celestialBodiesSearch').click(function() {
		processSearch(0);
	});
	$('#celestialBodiesSearch').trigger('click');	//fire one off 

	// rest of the init


	$(window).resize(resize);
	resize();

	if (true) {
		//extra horizon data
		//has to be telnetted out of jpl and that takes about 5 mins for everything
		//some dude must have typed in every one of the 200 info cards by hand
		// because there's nothing consistent about formatting or variable names
		//so I'm thinking cron job to update and then integrate to extrapolate for later time.
		
		planets = new Planets();
		for (var i = 0; i < horizonsData.coords.length; ++i) {
			var data = horizonsData.coords[i];
			var planetClass = Planets.prototype.planetClassForHorizonID[data.id];
			if (planetClass) {	//excluding the BCC and Ln points
				var planet = planets[planetClass.prototype.index];
				for (var j = 0; j < 3; ++j) {
					//convert km to m
					planet.pos[j] = data.pos[j] * 1000;
					planet.vel[j] = data.vel[j] * 1000;
				}
			}	
		}

		
		julianDate = horizonsData.julianDate;	//get most recent recorded date from planets
		//TODO get current julian date from client
		// integrate forward by the small timestep between the two
		
		
		initJulianDate = julianDate;
		refreshCurrentTimeText();
		
		initPlanets = planets.clone();
		
		init2();
	}

	if (false) { //getting single sets of states ...
		//returning a single state
		var planetNames = [];
		$.each(primaryPlanetHorizonIDs, function(_,horizonID) {
			var planetClassPrototype = Planets.prototype.planetClassForHorizonID[horizonID].prototype;
			planetNames.push(planetClassPrototype.name.toLowerCase());
		});

		//works, but not for the extra horizon data
		if (false) {
			var url = 'http://www.astro-phys.com/api/de406/states?bodies=' + planetNames.join(',');
			//var url = 'astro-phys-state.json';
			console.log('reading from '+url);
			$.ajax({
				url : url,
				dataType : 'jsonp'
			}).done(function(d) {
				planets = Planets.prototype.fromAstroPhysState(d.results);
				julianDate = d.date;
				
				initPlanets = planets.clone();
				initJulianDate = julianDate;
				refreshCurrentTimeText();
				
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
				
				initPlanets = planets.clone();
				initJulianDate = julianDate;
				refreshCurrentTimeText();
				
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
			initPlanets = planets.clone();
			initJulianDate = julianDate;
			refreshCurrentTimeText();
			
			//record.getStates is giving back the whole coefficients ... ?	
			var positions = record.getPositions(julianDate);
			planets = new Planets();
			for (var i = 0; i < planets.length; ++i) {
				var planet = planets[i];
				//convert km to m
				planet.pos[0] = positions[planet.name.toLowerCase()][0] * 1000;
				planet.pos[1] = positions[planet.name.toLowerCase()][1] * 1000;
				planet.pos[2] = positions[planet.name.toLowerCase()][2] * 1000;
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

	$.each(primaryPlanetHorizonIDs, function(_,horizonID) {
		var planetClassPrototype = Planets.prototype.planetClassForHorizonID[horizonID].prototype;
		imgs.push('textures/'+planetClassPrototype.name.toLowerCase()+'.png');
	});
	imgs.push('textures/saturn-rings.png');
	for (var i = 0; i < skyTexFilenames.length; ++i) {
		imgs.push(skyTexFilenames[i]);
	}
	console.log('loading '+imgs.join(', '));
	$(imgs).preload(function(){
		$('#loadingDiv').hide();
	}, function(percent){
		$('#loading').attr('value', parseInt(100*percent));
	});
	
	$('#menu').show();
	$('#timeDiv').show();
	$('#infoDiv').show();

	init3();
}

function initStars() {
	var xhr = new XMLHttpRequest();
	xhr.open('GET', 'hyg/starpositions.f32', true);
	xhr.responseType = 'arraybuffer';
	/* if we want a progress bar ...
	xhr.onprogress = function(e) {
		if (e.total) {
			progress.attr('value', parseInt(e.loaded / e.total * 100));
		}
	};
	*/
	xhr.onload = function(e) {
		console.log('loaded star data!');
		var arrayBuffer = this.response;
		var data = new DataView(arrayBuffer);
		
		var floatBuffer = new Float32Array(data.byteLength / Float32Array.BYTES_PER_ELEMENT);
		var len = floatBuffer.length;
		for (var j = 0; j < len; ++j) {
			var x = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);
		
			//convert from parsec coordinates to meters ... max float is 10^38, so let's hope (/warn) if an incoming value is close to 10^22
			if (Math.abs(x) > 1e+20) {
				console.log('star '+Math.floor(j/3)+' has coordinate that position exceeds fp resolution'); 
			}
			x *= 3.08567758e+16;

			floatBuffer[j] = x;
		}

		//now that we have the float buffer ...
		starsObj = new glutil.SceneObject({
			mode : gl.POINTS,
			shader : colorShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({data : floatBuffer})
			},
			uniforms : {
				pointSize : 1,
				color : [1,1,1,1],
			},
			parent : null
		});

	};
	xhr.send();
}

function initPlanetColorSchRadiusAngle(planet) {
	var planetClassPrototype = planet.init.prototype;
	var colors = {
		Sun:[1,1,0],
		Mercury:[.7,0,.2],
		Venus:[0,1,0],
		Earth:[0,0,1],
		Moon:[.6,.6,.6],
		Mars:[1,0,0],
		Jupiter:[1,.5,0],
		Saturn:[1,0,.5],
		Uranus:[0,1,1],
		Neptune:[1,0,1],
		Pluto:[0,.5,1]
	};
	var color = colors[planet.name];
	if (!color) {
		//console.log("failed to find color for "+planet.name);
		planetClassPrototype.color = [Math.random(), Math.random(), Math.random(), 1];
		vec3.normalize(planetClassPrototype.color, planetClassPrototype.color);
	} else {
		planetClassPrototype.color = [color[0], color[1], color[2], 1];
	}
	planetClassPrototype.schwarzschildRadius = 2 * planetClassPrototype.mass * kilogramsPerMeter; 
	planetClassPrototype.angle = [0,0,0,1];			// rotation ... only used for earth at the moment
}

function initPlanetSceneLatLonLineObjs(planet) {
	var planetClassPrototype = planet.init.prototype;
	if (planet.radius === undefined) {
		//only/always use a point/basis/etc?
	} else {
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
				texCoordArray.push(lat / 180 + .5);

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
		
		var vertexBuffer = new glutil.ArrayBuffer({
			data : vertexArray
		});
		planetClassPrototype.sceneObj = new glutil.SceneObject({
			mode : gl.TRIANGLES,
			indexes : new glutil.ElementArrayBuffer({
				data : triIndexArray
			}),
			attrs : {
				vertex : vertexBuffer,
				texCoord : new glutil.ArrayBuffer({
					dim : 2,
					data : texCoordArray
				}),
				tide : new glutil.ArrayBuffer({
					dim : 1,
					data : tideArray,
					usage : gl.DYNAMIC_DRAW
				})
			},
			uniforms : {
				color : planetClassPrototype.color,
				pos : [0,0,0],
				angle : [0,0,0,1],
				sunPos : [0,0,0],
				ambient : planetClassPrototype.name == 'Sun' ? 1 : .3
			},
			parent : null,
			static : true 
		});
		planetClassPrototype.latLonObj = new glutil.SceneObject({
			mode : gl.LINES,
			indexes : new glutil.ElementArrayBuffer({
				data : latLonIndexArray
			}),
			shader : colorShader,
			attrs : {
				vertex : vertexBuffer
			},
			uniforms : {
				color : [1,1,1,.2]
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			pos : [0,0,0],
			angle : [0,0,0,1],
			parent : null
		});
	}
	
	if (!(planetClassPrototype.isComet || planetClassPrototype.isAsteroid) &&
		(planetClassPrototype.parent === planets.indexes.Sun ||
		planet === planets[planets.indexes.Sun] ||
		planet === planets[planets.indexes.Moon]))
	{
		planetClassPrototype.lineObj = new glutil.SceneObject({
			mode : gl.LINES,
			shader : colorShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({
					count : 2,
					usage : gl.DYNAMIC_DRAW
				})
			},
			uniforms : {
				color : planetClassPrototype.color
			},
			parent : null,
			static : true 
		});
	}
}

function init3() {
	hsvTex = new glutil.HSVTexture(256);
	
	colorShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex, 1.);
	gl_Position = projMat * vtx4;
	gl_PointSize = pointSize;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
void main() {
	gl_FragColor = color;
}
*/}),
		uniforms : {
			color : [1,1,1,1],
			pointSize : 4
		}
	});
	
	planetColorShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform vec3 sunPos;
uniform float pointSize;
varying vec3 lightDir;
varying vec3 normal;
vec3 quatRotate(vec4 q, vec3 v){ 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	vec3 vtx3 = quatRotate(angle, vertex) + pos;
	normal = quatRotate(angle, normalize(vertex));
	lightDir = normalize(sunPos - vtx3);
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_PointSize = pointSize;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
varying vec3 lightDir;
varying vec3 normal;
uniform float ambient;
void main() {
	gl_FragColor = color * max(ambient, dot(lightDir, normal));
}
*/}),
		uniforms : {
			color : [1,1,1,1],
			pointSize : 4
		}
	});

	planetTexShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
attribute vec2 texCoord;
varying vec2 texCoordv;
uniform vec3 pos;
uniform vec4 angle;
uniform vec3 sunPos;
uniform mat4 mvMat;
uniform mat4 projMat;
varying vec3 lightDir;
varying vec3 normal;
vec3 quatRotate(vec4 q, vec3 v){ 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	texCoordv = texCoord;
	vec3 vtx3 = quatRotate(angle, vertex) + pos;
	normal = quatRotate(angle, normalize(vertex));
	lightDir = normalize(sunPos - vtx3);
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
varying vec2 texCoordv;
varying vec3 lightDir;
varying vec3 normal;
uniform sampler2D tex;
uniform float ambient;
void main() {
	gl_FragColor = texture2D(tex, texCoordv) * max(ambient, dot(lightDir, normal));
}
*/}),
		uniforms : {
			tex : 0
		}
	});

	hsvShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
attribute vec2 texCoord;
attribute float tide;
uniform vec3 pos;
uniform vec4 angle;
uniform mat4 mvMat;
uniform mat4 projMat;
varying float tidev;
varying vec2 texCoordv;
vec3 quatRotate(vec4 q, vec3 v){ 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	tidev = tide;
	texCoordv = texCoord;
	vec3 vtx3 = quatRotate(angle, vertex) + pos;
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
varying float tidev;
varying vec2 texCoordv;
uniform sampler2D tex;
uniform sampler2D hsvTex;
uniform float heatAlpha;
void main() {
	vec4 hsvColor = texture2D(hsvTex, vec2(tidev, .5));
	vec4 planetColor = texture2D(tex, texCoordv);
	gl_FragColor = mix(planetColor, hsvColor, heatAlpha);
}
*/}),
		uniforms : {
			tex : 0,
			hsvTex : 1
		}
	});

/*
per-pixel heat shader:
//will be used by all per-pixel shaders
attribute vec3 vertex;
attribute vec2 texCoord;
uniform mat4 mvMat;
uniform mat4 objMat;
uniform mat4 projMat;
varying vec2 texCoordv;
varying vec4 worldv;
void main() {
	tidev = tide;
	texCoordv = texCoord;
	worldv = objMat * vec4(vertex, 1.);
	gl_Position = projMat * mvMat * worldv;
	gl_Position.z = depthfunction(gl_Position);
}

tangent tidal shader:
varying vec2 texCoordv;
varying vec4 worldv;
uniform float heatAlpha;
uniform sampler2D tex;
uniform sampler2D hsvTex;
void main() {
	float tidev = something(worldv);	//do calculations with this.  fix it if the precision isn't good enough.
	vec4 hsvColor = texture2D(hsvTex, vec2(tidev, .5));
	vec4 planetColor = texture2D(tex, texCoordv);
	gl_FragColor = mix(planetColor, hsvColor, heatAlpha);
}
*/

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

	pointObj = new glutil.SceneObject({
		mode : gl.POINTS,
		shader : colorShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				count : 1,
				usage : gl.DYNAMIC_DRAW
			}),
		},
		uniforms : {
			pointSize : 4
		},
		parent : null
	});

	//init stars now that shaders are made 
	initStars();

	var planetsDone = 0;

	for (var planetIndex_ = 0; planetIndex_ < planets.length; ++planetIndex_) { (function(){
		var planetIndex = planetIndex_;
		var planet = planets[planetIndex];
		var planetClassPrototype = planet.init.prototype;
		initPlanetColorSchRadiusAngle(planet);
		
		var checkDone = function() {
			
			initPlanetSceneLatLonLineObjs(planet);
			
			++planetsDone;
			if (planetsDone == planets.length) {
				initOrbitPaths();
			}
		};

		// load texture
		if (primaryPlanetHorizonIDs.indexOf(planet.id) !== -1) {
			var img = new Image();
			img.onload = function() {
				planetClassPrototype.tex = new glutil.Texture2D({
					flipY : true,
					data : img,
					minFilter : gl.LINEAR_MIPMAP_LINEAR,
					magFilter : gl.LINEAR,
					generateMipmap : true
				});
				//checkDone();
			};
			img.onerror = function() {
				console.log('failed to find texture for planet '+planet.name);
				//checkDone();
			};
			img.src = 'textures/'+planet.name.toLowerCase()+'.png';
		} else {
			//checkDone();
		}
		
		//or... don't wait for textures 
		checkDone();
	})(); }

	//while we're here, load the rings
	{
		var ringShader = new ModifiedDepthShaderProgram({
			context : gl,
			vertexCode : mlstr(function(){/*
#define M_PI 3.1415926535897931 
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float minRadius;
uniform float maxRadius;
varying vec2 texCoordv;
void main() {
	texCoordv = vertex;
	//vertex is the uv of the texcoord wrapped around the annulus 
	//u coord is radial, v coord is angular 
	float rad = mix(minRadius, maxRadius, vertex.x);
	float theta = 2. * M_PI * vertex.y;
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);
	vec4 vtx4 = mvMat * vec4(cosTheta * rad, sinTheta * rad, 0., 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
varying vec2 texCoordv;
uniform sampler2D tex;
void main() {
	gl_FragColor = texture2D(tex, texCoordv);
}
*/})
		});
		
		var img = new Image();
		img.onload = function() {
			var planetClassPrototype = planets[planets.indexes.Saturn].init.prototype;
			planetClassPrototype.ringTex = new glutil.Texture2D({
				data : img,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				generateMipmap : true
			});

			var vertexes = [];
			var res = 200;
			for (var i = 0; i < res; ++i) {
				var f = i / (res - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}
		
			//and a ring object
			planetClassPrototype.ringObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				shader : ringShader,
				attrs : {
					vertex : new glutil.ArrayBuffer({
						dim : 2,
						data : vertexes
					})
				},
				uniforms : {
					tex : 0,
					minRadius : 7.4e+7,	//meters
					maxRadius : 1.4025e+8
				},
				texs : [planetClassPrototype.ringTex],
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null
			});
		};
		img.src = 'textures/saturn-rings.png';
	}
}

var orbitShader;
function initPlanetOrbitPathObj(planet) {
	var planetClassPrototype = planet.init.prototype;
	var planetIndex = planet.index;

	// based on planet position and velocity, find plane of orbit
	if (planet.parent === undefined) {
		console.log(planet.name+' has no orbit parent');
		planetClassPrototype.orbitAxis = [0,0,1];
		planetClassPrototype.orbitBasis = [[1,0,0],[0,1,0],[0,0,1]];
		return;
	}

	var parentPlanet = planets[planet.parent];
	if (parentPlanet === undefined) {
		console.log(planet.name+' has an invalid orbit planet');
		planetClassPrototype.orbitAxis = [0,0,1];
		planetClassPrototype.orbitBasis = [[1,0,0],[0,1,0],[0,0,1]];
		return;
	}

	//if it's a comet then 
	//calculate pos and vel and mass by parameters ... ?
	// or just put an orbit mesh there?
	if (planet.isComet || planet.isAsteroid) {
	
		var eccentricity = assert(planet.orbitData.eccentricity);
		if (eccentricity >= 1) {
			console.log("WARNING: omitting hyperbolic orbit of comet "+planet.name);
			return;
		}

		var pericenterDistance, semiMajorAxis;
		if (planet.isComet) {
			pericenterDistance = assert(planet.orbitData.perihelionDistance);
			semiMajorAxis = pericenterDistance / (1 - eccentricity);
		} else if (planet.isAsteroid) {
			semiMajorAxis = assert(planet.orbitData.semiMajorAxis);
			pericenterDistance = semiMajorAxis * (1 - eccentricity);
		}

		var gravitationalParameter = gravitationalConstant * parentPlanet.mass;	//assuming the comet mass is negligible, since the comet mass is not provided
		var semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;
		var orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed  / gravitationalParameter) / (60*60*24);	//julian day
		
		var longitudeOfAscendingNode = assert(planet.orbitData.longitudeOfAscendingNode);
		var cosAscending = Math.cos(longitudeOfAscendingNode);
		var sinAscending = Math.sin(longitudeOfAscendingNode);

		var argumentOfPericenter = assert(planet.orbitData.argumentOfPerihelion);
		var cosPericenter = Math.cos(argumentOfPericenter);
		var sinPericenter = Math.sin(argumentOfPericenter);

		var inclination = assert(planet.orbitData.inclination);
		var cosInclination = Math.cos(inclination);
		var sinInclination = Math.sin(inclination);
	
		var A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
				 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
				 semiMajorAxis * sinPericenter * sinInclination];
		var B = [-semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
				 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
				 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * cosPericenter * sinInclination];

		var timeOfPeriapsisCrossing, eccentricAnomaly;
		if (planet.isComet) {
			//how long until it crosses the periapsis
			// solve for eccentric anomaly...
			timeOfPeriapsisCrossing = planet.orbitData.timeOfPerihelionPassage;	//julian day
			var tau = (julianDate - timeOfPeriapsisCrossing) / orbitalPeriod;	//unitless
			var pathEccentricAnomaly = 2 * Math.PI * tau;
			var pathCosEccentricAnomaly = Math.cos(pathEccentricAnomaly);
			var pathSinEccentricAnomaly = Math.sin(pathEccentricAnomaly);
			var posX = A[0] * (pathCosEccentricAnomaly - eccentricity) + B[0] * pathSinEccentricAnomaly;
			var posY = A[1] * (pathCosEccentricAnomaly - eccentricity) + B[1] * pathSinEccentricAnomaly;
			var posZ = A[2] * (pathCosEccentricAnomaly - eccentricity) + B[2] * pathSinEccentricAnomaly;
			var velX = (A[0] * -pathSinEccentricAnomaly + B[0] * pathCosEccentricAnomaly) * 2 * Math.PI / orbitalPeriod;	//m/day
			var velY = (A[1] * -pathSinEccentricAnomaly + B[1] * pathCosEccentricAnomaly) * 2 * Math.PI / orbitalPeriod;
			var velZ = (A[2] * -pathSinEccentricAnomaly + B[2] * pathCosEccentricAnomaly) * 2 * Math.PI / orbitalPeriod;
			planet.pos[0] = posX + parentPlanet.pos[0];
			planet.pos[1] = posY + parentPlanet.pos[1];
			planet.pos[2] = posZ + parentPlanet.pos[2];
			planet.vel[0] = velX + parentPlanet.vel[0];
			planet.vel[1] = velY + parentPlanet.vel[1];
			planet.vel[2] = velZ + parentPlanet.vel[2];
			vec3.copy(initPlanets[planetIndex].pos, planet.pos);
			vec3.copy(initPlanets[planetIndex].vel, planet.vel);
			velX = velX / (60*60*24);	//m/s
			velY = velY / (60*60*24);
			velZ = velZ / (60*60*24);
			var posDotVel = posX * velX + posY * velY + posZ * velZ;	//m^2/s

			//rather than assume the comet is at its closest approach, calculate it using the 'timeOfPerihelionPassage' / 'epoch' info
			var distanceToParent = vec3.length(planet.pos);
			var cosEccentricAnomaly = (1 - distanceToParent / semiMajorAxis) / eccentricity;						//unitless
			var sinEccentricAnomaly = posDotVel / (eccentricity * Math.sqrt(gravitationalParameter * semiMajorAxis));	//m^2/s / sqrt(m^3/s^2 * m) = m^2/s / sqrt(m^4/s^2) = m^2/s / (m^2/s) = unitless
			eccentricAnomaly = Math.atan2(sinEccentricAnomaly, cosEccentricAnomaly);	//radians (unitless)
		} else if (planet.isAsteroid) {
			//solve Newton Rhapson
			//f(E) = M - E + e sin E = 0
			var meanAnomaly = planet.orbitData.meanAnomaly;
			eccentricAnomaly = meanAnomaly;
			for (var i = 0; i < 10; ++i) {
				var func = meanAnomaly - eccentricAnomaly + eccentricity * Math.sin(eccentricAnomaly);
				var deriv = -1 + eccentricity * Math.cos(eccentricAnomaly);	//has zeroes ...
				var delta = func / deriv;
				if (Math.abs(delta) < 1e-15) break;
				eccentricAnomaly -= delta;
			}
			var sinEccentricAnomaly = Math.sin(eccentricAnomaly);	
			var cosEccentricAnomaly = Math.cos(eccentricAnomaly);
			timeOfPeriapsisCrossing = -(eccentricAnomaly - eccentricity * sinEccentricAnomaly) / Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (60*60*24);	//julian day
		
			var posX = A[0] * (cosEccentricAnomaly - eccentricity) + B[0] * sinEccentricAnomaly;
			var posY = A[1] * (cosEccentricAnomaly - eccentricity) + B[1] * sinEccentricAnomaly;
			var posZ = A[2] * (cosEccentricAnomaly - eccentricity) + B[2] * sinEccentricAnomaly;
			var velX = (A[0] * -sinEccentricAnomaly + B[0] * cosEccentricAnomaly) * 2 * Math.PI / orbitalPeriod;	//m/day
			var velY = (A[1] * -sinEccentricAnomaly + B[1] * cosEccentricAnomaly) * 2 * Math.PI / orbitalPeriod;
			var velZ = (A[2] * -sinEccentricAnomaly + B[2] * cosEccentricAnomaly) * 2 * Math.PI / orbitalPeriod;
			planet.pos[0] = posX + parentPlanet.pos[0];
			planet.pos[1] = posY + parentPlanet.pos[1];
			planet.pos[2] = posZ + parentPlanet.pos[2];
			planet.vel[0] = velX + parentPlanet.vel[0];
			planet.vel[1] = velY + parentPlanet.vel[1];
			planet.vel[2] = velZ + parentPlanet.vel[2];
			vec3.copy(initPlanets[planetIndex].pos, planet.pos);
			vec3.copy(initPlanets[planetIndex].vel, planet.vel);	
		}

		planetClassPrototype.keplerianOrbitalElements = {
			semiMajorAxis : semiMajorAxis,
			eccentricity : eccentricity,
			eccentricAnomaly : eccentricAnomaly,
			longitudeOfAscendingNode : longitudeOfAscendingNode,
			argumentOfPericenter : argumentOfPericenter,
			inclination : inclination,
			timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
			orbitalPeriod : orbitalPeriod
		};
		
		//iterate around the eccentric anomaly to reconstruct the path
		var vertexes = [];
		var res = 250;
		for (var i = 0; i < res; ++i) {
			var theta = i / (res - 1) * 2 * Math.PI;
			var pathEccentricAnomaly = eccentricAnomaly + theta;
			var pathCosEccentricAnomaly = Math.cos(pathEccentricAnomaly);
			var pathSinEccentricAnomaly = Math.sin(pathEccentricAnomaly);
			
			var vtxPosX = A[0] * (pathCosEccentricAnomaly - eccentricity) + B[0] * pathSinEccentricAnomaly;
			var vtxPosY = A[1] * (pathCosEccentricAnomaly - eccentricity) + B[1] * pathSinEccentricAnomaly;
			var vtxPosZ = A[2] * (pathCosEccentricAnomaly - eccentricity) + B[2] * pathSinEccentricAnomaly;

			//offset by center (parent planet), calculated error (which is pretty small), and then from our own position (for subsequent offsetting/rendering)
			vtxPosX += parentPlanet.pos[0] /*- checkPosToPosX*/ - planet.pos[0];
			vtxPosY += parentPlanet.pos[1] /*- checkPosToPosY*/ - planet.pos[1];
			vtxPosZ += parentPlanet.pos[2] /*- checkPosToPosZ*/ - planet.pos[2];
	
			//add to buffer
			vertexes.push(vtxPosX);
			vertexes.push(vtxPosY);
			vertexes.push(vtxPosZ);
		
			//pack transparency info into the vertex
			//comets (and asteroids?) aren't evaluated perfectly, so I'm getting a bit of an offest ... TODO fixme
			var alphaAngleOffset = 238/250;
			var alpha = ((i / (res-1) + alphaAngleOffset) % 1) * .75 + .25;
			vertexes.push(alpha);
		}

		planetClassPrototype.orbitPathObj = new glutil.SceneObject({
			mode : gl.LINE_STRIP,
			shader : orbitShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({
					dim : 4,
					data : vertexes
				})
			},
			uniforms : {
				color : planetClassPrototype.color
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			pos : [0,0,0],
			angle : [0,0,0,1],
			parent : null
		});	

		return;
	}

	//consider position relative to orbitting parent
	// should I be doing the same thing with the velocity?  probably...
	var posX = planet.pos[0] - parentPlanet.pos[0];
	var posY = planet.pos[1] - parentPlanet.pos[1];
	var posZ = planet.pos[2] - parentPlanet.pos[2];

	//convert from m/day to m/s to coincide with the units of our gravitational constant
	var velX = (planet.vel[0] - parentPlanet.vel[0]) / (60 * 60 * 24);
	var velY = (planet.vel[1] - parentPlanet.vel[1]) / (60 * 60 * 24);
	var velZ = (planet.vel[2] - parentPlanet.vel[2]) / (60 * 60 * 24);

	var posDotVel = posX * velX + posY * velY + posZ * velZ;	//m^2/s

	var angularMomentumX = posY * velZ - posZ * velY; //m^2/s
	var angularMomentumY = posZ * velX - posX * velZ; 
	var angularMomentumZ = posX * velY - posY * velX; 
	var angularMomentumMagSq = angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY + angularMomentumZ * angularMomentumZ;		//m^4/s^2
	var angularMomentumMag = Math.sqrt(angularMomentumMagSq);
	if (angularMomentumMag < 1e-9) {
		planetClassPrototype.orbitAxis = [0,0,1];
	} else {
		var axisX = angularMomentumX / angularMomentumMag;
		var axisY = angularMomentumY / angularMomentumMag;
		var axisZ = angularMomentumZ / angularMomentumMag;
		planetClassPrototype.orbitAxis = [axisX, axisY, axisZ];
	}

	var basisX = [0,0,0];
	var basisY = [0,0,0];
	var basisZ = planetClassPrototype.orbitAxis;
	//TODO use the inclination and longitudeOfAscendingNode
	calcBasis(basisX, basisY, basisZ);
	//a[j][i] = a_ij, so our indexing is backwards, but our storage is column-major
	planetClassPrototype.orbitBasis = [basisX, basisY, basisZ];

	//now decompose the relative position in the coordinates of the orbit basis
	//i've eliminated all but one of the rotation degrees of freedom ...

	//http://www.mathworks.com/matlabcentral/fileexchange/31333-orbital-elements-from-positionvelocity-vectors/content/vec2orbElem.m

	var velSq = velX * velX + velY * velY + velZ * velZ;		//(m/s)^2
	var distanceToParent = Math.sqrt(posX * posX + posY * posY + posZ * posZ);		//m
	var gravitationalParameter = gravitationalConstant * ((planet.mass || 0) + parentPlanet.mass);	//m^3 / (kg s^2) * kg = m^3 / s^2
	var specificOrbitalEnergy  = .5 * velSq - gravitationalParameter / distanceToParent;		//m^2 / s^2 - m^3 / s^2 / m = m^2/s^2, supposed to be negative for elliptical orbits
	var semiMajorAxis = -.5 * gravitationalParameter / specificOrbitalEnergy;		//m^3/s^2 / (m^2/s^2) = m
	var semiLatusRectum = angularMomentumMagSq / gravitationalParameter;			//m^4/s^2 / (m^3/s^2) = m
	var eccentricity = Math.sqrt(1 - semiLatusRectum / semiMajorAxis);				//unitless (assuming elliptical orbit)
	
	var cosEccentricAnomaly = (1 - distanceToParent / semiMajorAxis) / eccentricity;						//unitless
	var sinEccentricAnomaly = posDotVel / (eccentricity * Math.sqrt(gravitationalParameter * semiMajorAxis));	//m^2/s / sqrt(m^3/s^2 * m) = m^2/s / sqrt(m^4/s^2) = m^2/s / (m^2/s) = unitless
	var eccentricAnomaly = Math.atan2(sinEccentricAnomaly, cosEccentricAnomaly);	//radians (unitless)
	
	var sinInclination = Math.sqrt(angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY) / angularMomentumMag;	//unitless
	var cosInclination = angularMomentumZ / angularMomentumMag;	//unitless
	var inclination = Math.atan2(sinInclination, cosInclination);

	var sinPericenter = ((velX * angularMomentumY - velY * angularMomentumX) / gravitationalParameter - posZ / distanceToParent) / (eccentricity * sinInclination);
	var cosPericenter = (angularMomentumMag * velZ / gravitationalParameter - (angularMomentumX * posY - angularMomentumY * posX) / (angularMomentumMag * distanceToParent)) / (eccentricity * sinInclination);
	var argumentOfPericenter = Math.atan(sinPericenter, cosPericenter);

	var cosAscending = -angularMomentumY / (angularMomentumMag * sinInclination);
	var sinAscending = angularMomentumX / (angularMomentumMag * sinInclination);
	var longitudeOfAscendingNode = Math.atan2(sinAscending, cosAscending);

	var semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;	//m^3
	var orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed  / gravitationalParameter) / (60*60*24);	//julian day
	var timeOfPeriapsisCrossing = -(eccentricAnomaly - eccentricity * sinEccentricAnomaly) / Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (60*60*24);	//julian day

	var A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
			 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
			 semiMajorAxis * sinPericenter * sinInclination];
	var B = [-semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * cosPericenter * sinInclination];

	/*
	to convert back:
	pos[i] = A[i] * (cosEccentricAnomaly - eccentricity) + B[i] * sinEccentricAnomaly
	rDot[i] = (-A[i] * sinEccentricAnomaly + B[i] * cosEccentricAnomaly) * Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (1 - eccentricity * cosEccentricAnomaly) 
	*/
	var checkPosX = A[0] * (cosEccentricAnomaly - eccentricity) + B[0] * sinEccentricAnomaly;
	var checkPosY = A[1] * (cosEccentricAnomaly - eccentricity) + B[1] * sinEccentricAnomaly;
	var checkPosZ = A[2] * (cosEccentricAnomaly - eccentricity) + B[2] * sinEccentricAnomaly;

	var checkPosToPosX = checkPosX - posX;
	var checkPosToPosY = checkPosY - posY;
	var checkPosToPosZ = checkPosZ - posZ;
	var checkPosToPosDist = Math.sqrt(checkPosToPosX * checkPosToPosX + checkPosToPosY * checkPosToPosY + checkPosToPosZ * checkPosToPosZ);
	var checkPosError = checkPosToPosDist / distanceToParent;
	if (checkPosError === checkPosError) {
		if (checkPosError > 1e-5) {	//only report significant error
			console.log(planet.name+' error of reconstructed position '+ checkPosError);
		}
	} else {	//NaN? debug!
		
	/*		
		for (k in planetClassPrototype) {
			var v = planetClassPrototype[k];
			if (k != 'name' && typeof(v) != 'function') {
				console.log(planetClassPrototype.name, k, v);
			}
		}
	*/
		console.log(planet.name+' has no orbit info.  mass: '+planet.mass+' radius: '+planet.radius);
	}

	planetClassPrototype.keplerianOrbitalElements = {
		relVelSq : velSq,
		gravitationalParameter : gravitationalParameter,
		specificOrbitalEnergy : specificOrbitalEnergy,
		distanceToParent : distanceToParent,
		semiMajorAxis : semiMajorAxis,
		semiLatusRectum : semiLatusRectum,
		eccentricity : eccentricity,
		eccentricAnomaly : eccentricAnomaly,
		inclination : inclination,
		argumentOfPericenter : argumentOfPericenter,
		longitudeOfAscendingNode : longitudeOfAscendingNode,
		timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
		orbitalPeriod : orbitalPeriod
	};	

	//not NaN, we successfully reconstructed the position
	if (checkPosError === checkPosError) {
		//iterate around the eccentric anomaly to reconstruct the path
		var vertexes = [];
		var res = 250;
		
		for (var i = 0; i < res; ++i) {
			var theta = i / (res - 1) * 2 * Math.PI;
			var pathEccentricAnomaly = eccentricAnomaly + theta;
			var pathCosEccentricAnomaly = Math.cos(pathEccentricAnomaly);
			var pathSinEccentricAnomaly = Math.sin(pathEccentricAnomaly);
			
			var vtxPosX = A[0] * (pathCosEccentricAnomaly - eccentricity) + B[0] * pathSinEccentricAnomaly;
			var vtxPosY = A[1] * (pathCosEccentricAnomaly - eccentricity) + B[1] * pathSinEccentricAnomaly;
			var vtxPosZ = A[2] * (pathCosEccentricAnomaly - eccentricity) + B[2] * pathSinEccentricAnomaly;

			//offset by center (parent planet), calculated error (which is pretty small), and then from our own position (for subsequent offsetting/rendering)
			var sun = planets
			//if centered around the orbitting planet ...
			vtxPosX += parentPlanet.pos[0] - checkPosToPosX - planet.pos[0];
			vtxPosY += parentPlanet.pos[1] - checkPosToPosY - planet.pos[1];
			vtxPosZ += parentPlanet.pos[2] - checkPosToPosZ - planet.pos[2];
			//global positions...
			//vtxPosX += parentPlanet.pos[0];
			//vtxPosY += parentPlanet.pos[1];
			//vtxPosZ += parentPlanet.pos[2];
	
			//add to buffer
			vertexes.push(vtxPosX);
			vertexes.push(vtxPosY);
			vertexes.push(vtxPosZ);
		
			//pack transparency info into the vertex
			var alpha = (i / (res-1)) * .75 + .25;
			vertexes.push(alpha);
		}
		
		planetClassPrototype.orbitPathObj = new glutil.SceneObject({
			mode : gl.LINE_STRIP,
			shader : orbitShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({
					dim : 4,
					data : vertexes
				})
			},
			uniforms : {
				color : planetClassPrototype.color
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			pos : [0,0,0],
			angle : [0,0,0,1],
			parent : null
		});
	}

}

function initOrbitPaths() {
	var gravWellShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec4 vertex;
varying float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex.xyz, 1.);
	alpha = vertex.w;
	gl_Position = projMat * vtx4;
	gl_PointSize = 4.;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
varying float alpha;
void main() {
	gl_FragColor = vec4(color.xyz, color.w * alpha);
}
*/})
	});


	orbitShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec4 vertex;
varying float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex.xyz, 1.);
	alpha = vertex.w;
	gl_Position = projMat * vtx4;
	gl_PointSize = 4.;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
varying float alpha;
void main() {
	gl_FragColor = vec4(color.xyz, color.w * alpha);
}
*/})
	});

	//TODO update these when we integrate!
	var calcOrbitPathStartTime = Date.now();
	for (var i = 0; i < planets.length; ++i) {
		initPlanetOrbitPathObj(planets[i]);
	}
	var calcOrbitPathEndTime = Date.now();

	var setPlanetTiltAngleToMoonOrbitPlane = function(planetName, moonName) {
		var planetClassPrototype = Planets.prototype.planetClasses[Planets.prototype.indexes[planetName]].prototype;
		var moonClassPrototype = Planets.prototype.planetClasses[Planets.prototype.indexes[moonName]].prototype;
		assert(planetClassPrototype.tiltAngle === undefined);
		planetClassPrototype.tiltAngle = [0,0,0,1];
		quat.rotateZ(planetClassPrototype.tiltAngle, planetClassPrototype.tiltAngle, moonClassPrototype.keplerianOrbitalElements.longitudeOfAscendingNode);
		quat.rotateX(planetClassPrototype.tiltAngle, planetClassPrototype.tiltAngle, moonClassPrototype.keplerianOrbitalElements.inclination);
	};

	//accepts degrees
	//TODO start at orbit axis plane rather than earth's (ie J2000) orbit axis plane
	var setPlanetTiltAngleToFixedValue = function(planetName, inclination, tiltDirection) {
		if (tiltDirection === undefined) tiltDirection = 0;
		var planetClassPrototype = Planets.prototype.planetClasses[Planets.prototype.indexes[planetName]].prototype;
		planetClassPrototype.tiltAngle = [0,0,0,1];
		quat.rotateZ(planetClassPrototype.tiltAngle, planetClassPrototype.tiltAngle, Math.rad(tiltDirection));
		quat.rotateX(planetClassPrototype.tiltAngle, planetClassPrototype.tiltAngle, Math.rad(inclination));
	};

	setPlanetTiltAngleToFixedValue('Mercury', 2.11/60);		//TODO tilt from mercury orbit plane.  until then it's off
	setPlanetTiltAngleToFixedValue('Venus', 177.3);			//TODO tilt from venus orbit plane.  until then, this measures 175 degrees.
	setPlanetTiltAngleToFixedValue('Earth', 23.44, 180);
	setPlanetTiltAngleToMoonOrbitPlane('Mars', 'Phobos');		//ours: 25.79, exact: 25.19
	setPlanetTiltAngleToMoonOrbitPlane('Jupiter', 'Metis');		//ours: 3.12, exact: 3.13
	setPlanetTiltAngleToMoonOrbitPlane('Saturn', 'Atlas');		//ours: 26.75, exact: 26.73
	setPlanetTiltAngleToMoonOrbitPlane('Uranus', 'Cordelia');	//ours: 97.71, exact: 97.77
	setPlanetTiltAngleToMoonOrbitPlane('Neptune', 'Galatea');	//ours: 28.365, exact: 28.32
	setPlanetTiltAngleToMoonOrbitPlane('Pluto', 'Charon');		//ours: 119, exact: 123

	//looks like low grav wells run into fp accuracy issues
	//how about extracting the depth and storing normalized values?
	var calcGravWellStartTime = Date.now();
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		var planet = planets[planetIndex];
		var planetClassPrototype = planet.init.prototype;

		/*
		gravity wells
	
		calculate spacetime embedding radius
		  for radial distance r, radius R, Schwarzschild radius Rs 
		 inside the planet  (r <= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R (r/R)^2 ))
		 outside the planet (r >= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R)) + sqrt(4Rs(r - Rs)) - sqrt(4Rs(R - Rs))
		*/
		var R = planetClassPrototype.radius;
		var Rs = planetClassPrototype.schwarzschildRadius;
		var R_Rs = R / Rs;
		var Rs_R = Rs / R;
		var R_sqrt_R_Rs = R * Math.sqrt(R_Rs);

		//extract the normalized scalar to post-multiply after transformation (so small fields will not lose accuracy)
		planetClassPrototype.gravityWellScalar = Math.sqrt(R / Rs) - Math.sqrt(R / Rs - 1);

		//populate first vertex
		var gravWellVtxs = [0,0,0,1];
		
		var gravWellIndexes = [];
		
		var rimax = 60;		//was 200 r and 60 th, but I added a lot of planets.  need to occlude these based on distance/angle ...
		var thimax = 60;
		
		for (var ri = 1; ri < rimax; ++ri) {
			var r = R * Math.pow(100, ri / rimax * (gravityWellRadialMaxLog100 - gravityWellRadialMinLog100) + gravityWellRadialMinLog100);
			//max radial dist is R * Math.pow(100, gravityWellRadialMaxLog100)

			var z;
			if (r <= R) {
				var r_R = r / R;
				z = R_sqrt_R_Rs * (1 - Math.sqrt(1 - Rs_R * r_R * r_R));
			} else {
				z = R_sqrt_R_Rs * (1 - Math.sqrt(1 - Rs_R))
					+ Math.sqrt(4 * Rs * (r - Rs))
					- Math.sqrt(4 * Rs * (R - Rs));
			}
			z /= planetClassPrototype.gravityWellScalar;
			
			for (var thi = 0; thi < thimax; ++thi) {
				var th = 2 * Math.PI * thi / thimax;
		
				var x = r * Math.cos(th);
				var y = r * Math.sin(th);
				
				//it would be great to include influences of other planets' gravity wells ...
				// but that would require us performing our calculation at some point in 3D space -- for the sake of the radial coordinate
				// I could approximate that if I knew the approximate plane of orbit of most all the planets ...
				// from there, take samples based on radial distance within that plane between planets ...
				// I could remap planes on a per-planet basis depending on which you are orbitting ...
				// or I could always try for something 3D ... exxhagerated pinch lattice vectors ...
				// to do this it might be best to get the chebyshev interval calculations working, or somehow calculate the ellipses of rotation of each planet
				//TODO also: recalculate the gravity well mesh when the planets change, to watch it in realtime	

				gravWellVtxs.push(x * planetClassPrototype.orbitBasis[0][0] + y * planetClassPrototype.orbitBasis[1][0] + z * planetClassPrototype.orbitBasis[2][0]);
				gravWellVtxs.push(x * planetClassPrototype.orbitBasis[0][1] + y * planetClassPrototype.orbitBasis[1][1] + z * planetClassPrototype.orbitBasis[2][1]);
				gravWellVtxs.push(x * planetClassPrototype.orbitBasis[0][2] + y * planetClassPrototype.orbitBasis[1][2] + z * planetClassPrototype.orbitBasis[2][2]);
				var tau = ri/(rimax-1);
				gravWellVtxs.push(1 - tau * tau * tau * tau * tau * tau * tau * tau);

				if (ri == 1) {
					gravWellIndexes.push(0);
					gravWellIndexes.push(1 + thi + thimax * (ri-1));
				}
				if (ri < rimax-1) {
					gravWellIndexes.push(1 + thi + thimax * (ri-1));
					gravWellIndexes.push(1 + thi + thimax * ri);	//plus one radial
				}
				gravWellIndexes.push(1 + thi + thimax * (ri-1));
				gravWellIndexes.push(1 + ((thi+1)%thimax) + thimax * (ri-1));	//plus one tangential
			}
		}
		planetClassPrototype.gravWellObj = new glutil.SceneObject({
			mode : gl.LINES,
			indexes : new glutil.ElementArrayBuffer({
				data : gravWellIndexes
			}),
			shader : gravWellShader,
			attrs : {
				vertex : new glutil.ArrayBuffer({
					dim : 4,
					data : gravWellVtxs
				})
			},
			uniforms : {
				color : [1,1,1,.2]
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			parent : null
		});
		//scenegraph is a mess
		//static objects have mvMat pointed to glutil.scene.mvMat
		//but dynamic objects only give control over position and angle
		//!static implies creating the object's matrices, but !static also implies overwriting them every draw call...
		planetClassPrototype.localMat = mat4.create();
		planetClassPrototype.gravWellObj.mvMat = mat4.create();
		planetClassPrototype.gravWellObj.uniforms.mvMat = planetClassPrototype.gravWellObj.mvMat;
	}
	var calcGravWellEndTime = Date.now();

	//on my machine this was 769ms
	console.log('calc orbit path time ',calcOrbitPathEndTime-calcOrbitPathStartTime,'ms');
	
	//on my machine this was 386ms
	console.log('calc grav well time ',calcGravWellEndTime-calcGravWellStartTime,'ms');

	//looks like doing these in realtime will mean toning the detail down a bit ...

	new glutil.TextureCube({
		flipY : true,
		generateMipmap : true,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR_MIPMAP_LINEAR,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		},
		urls : skyTexFilenames,
		onload : function(side,url,image) {
			if (image.width > glMaxCubeMapTextureSize || image.height > glMaxCubeMapTextureSize) {
				throw "cube map size "+image.width+"x"+image.height+" cannot exceed "+glMaxCubeMapTextureSize;
			}
		},
		done : function() {
			initSkyCube(this);
		}
	});

	initScene();
}

//init the "sky" cubemap (the galaxy background) once the texture for it loads
function initSkyCube(skyTex) {
	var cubeShader = new ModifiedDepthShaderProgram({
		context : gl,
		vertexCode : mlstr(function(){/*
attribute vec3 vertex;
varying vec3 vertexv;
uniform mat4 projMat;

void main() {
	vertexv = vertex;
	gl_Position = projMat * vec4(vertex, 1.);
}
*/}),
		fragmentCode : mlstr(function(){/*
precision mediump float;
varying vec3 vertexv;
uniform samplerCube skyTex;

uniform vec4 angle;
uniform vec4 viewAngle;
vec3 quatRotate( vec4 q, vec3 v ){ 
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
void main() {
	vec3 dir = vertexv;
	dir = quatRotate(viewAngle, dir);
	dir = quatRotate(vec4(angle.xyz, -angle.w), dir);
	gl_FragColor = .3 * textureCube(skyTex, dir);
	gl_FragColor.w = 1.; 
}
*/}),
		uniforms : {
			skyTex : 0
		}
	});

	var cubeVtxArray = new Float32Array(3*8);
	for (var i = 0; i < 8; i++) {
		cubeVtxArray[0+3*i] = glutil.view.zNear*10*(2*(i&1)-1);
		cubeVtxArray[1+3*i] = glutil.view.zNear*10*(2*((i>>1)&1)-1);
		cubeVtxArray[2+3*i] = glutil.view.zNear*10*(2*((i>>2)&1)-1);
	}

	var cubeIndexBuf = new glutil.ElementArrayBuffer({
		data : [
			5,7,3,3,1,5,		// <- each value has the x,y,z in the 0,1,2 bits (off = 0, on = 1)
			6,4,0,0,2,6,
			2,3,7,7,6,2,
			4,5,1,1,0,4,
			6,7,5,5,4,6,
			0,1,3,3,2,0
		]
	});

	/*
	right ascention is "to the right of the vernal equinox" ...
	but we're at the Autumnal equinox and NASA's J2000 coordinates say we're at x+ ... 
	shouldn't NASA's J2000 coordinate system be rotated 180' about the z axis, so x+ is the vernal equinox and x- is the Autumnal equinox?

	my galaxy texture is centered at x+ and lies in the xy plane
	*/
	//"Reconsidering the galactic coordinate system", Jia-Cheng Liu, Zi Zhu, and Hong Zhang, Oct 20, 2010
	//eqn 10
	var alpha = 2*Math.PI/24 * (12 + 1/60 * (51 + 1/60 * 26.27549)) + Math.PI;
	var delta = -2*Math.PI/180 * (27 + 1/60 * (7 + 1/60 * 41.7043));
	var theta = 2*Math.PI/180 * 122.93191857;
	var rz = [0, 0, Math.sin(.5*alpha), Math.cos(.5*alpha)];
	var ry = [0, Math.sin(.5*delta), 0, Math.cos(.5*delta)];
	var rx = [Math.sin(.5*theta), 0, 0, Math.cos(.5*theta)];
	var j2000ToGalactic = [0,0,0,1];
	quat.mul(j2000ToGalactic, ry, rx);
	quat.mul(j2000ToGalactic, rz, j2000ToGalactic);
	skyCubeObj = new glutil.SceneObject({
		mode : gl.TRIANGLES,
		indexes : cubeIndexBuf,
		shader : cubeShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				data : cubeVtxArray
			})
		},
		uniforms : {
			angle : j2000ToGalactic,
			viewAngle : glutil.view.angle
		},
		texs : [skyTex],
		parent : null
	});
}

function initScene() {
	//gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL);
	gl.clearColor(0,0,0,0);

	var trackPlanet = planets[planets.indexes.Earth];
	orbitPlanetIndex = trackPlanet.index;
	orbitTargetDistance = 2. * trackPlanet.radius;
	refreshOrbitTargetDistanceText();
	orbitDistance = orbitTargetDistance;
	$('#orbitPlanetText').text(trackPlanet.name);
	$('#hoverPlanetText').text(trackPlanet.name);

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
			quat.multiply(glutil.view.angle, glutil.view.angle, tmpQ);
			quat.normalize(glutil.view.angle, glutil.view.angle);
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
			quat.mul(glutil.view.angle, glutil.view.angle, r);
			quat.normalize(glutil.view.angle, glutil.view.angle);
		}
	}

	var orbitPlanet = planets[orbitPlanetIndex];

	/* converage angle on target planet * /

	var delta = vec3.create();
	vec3.scale(delta, glutil.view.pos, -1);
	var fwd = vec3.create();
	vec3.quatZAxis(fwd, glutil.view.angle);
	vec3.scale(fwd, fwd, -1);
	var axis = vec3.create();
	vec3.cross(axis, fwd, delta);
	vec3.scale(axis, axis, 1/vec3.length(delta));	//divide out length of delta so axis length is sin(theta)
	var sinTheta = vec3.length(axis);
	var theta = 0;
	if (sinTheta > 1e-3) {
		vec3.scale(axis, axis, 1/sinTheta);	//normalize axis
		theta = Math.asin(Math.clamp(sinTheta, -1,1));
		var q = quat.create();
		var convergeAngleCoeff = .5;
		quat.setAxisAngle(q, axis, theta * convergeAngleCoeff);
		quat.mul(glutil.view.angle, glutil.view.angle, q);
	}
	/**/

	// track ball orbit
	//assumes z points away from the planet

	var orbitCenter;
	if (orbitGeodeticLocation !== undefined) {
		planetGeodeticToSolarSystemBarycentric(
			orbitCenter,
			orbitPlanet, 
			orbitGeodeticLocation.lat, 
			orbitGeodeticLocation.lon, 
			orbitGeodeticLocation.height);
	} else {
		orbitCenter = orbitPlanet.pos;
	}
	var viewAngleZAxis = vec3.create();
	vec3.quatZAxis(viewAngleZAxis, glutil.view.angle);
	vec3.scale(glutil.view.pos, viewAngleZAxis, orbitDistance);
	vec3.add(glutil.view.pos, glutil.view.pos, orbitOffset);
	var orbitConvergeCoeff = .9;
	vec3.scale(orbitOffset, orbitOffset, orbitConvergeCoeff);
	{
		var logDist = Math.log(orbitDistance);
		var logTarget = Math.log(orbitTargetDistance);
		var coeff = .05;
		var newLogDist = (1 - coeff) * logDist + coeff * logTarget;
		orbitDistance = Math.exp(newLogDist);
	}

/*
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
	if (!integrationPaused) {
		// if we're integrating ...
		planets = integrate.run(julianDate, planets, integrateTimeStep, integrateFunction, integrationMethod, integrationArgs);
		julianDate += integrateTimeStep;
		refreshCurrentTimeText();
	}
	
	// set the angle for the time
	{
		//Nowhere in any of this do I seem to be taking tiltAngle into account ...
		var qdst = planets[planets.indexes.Earth].angle;
		quat.identity(qdst);
		quat.rotateZ(qdst, qdst, (julianDate % 1) * 2 * Math.PI);
	}
	
	//if we are close enough to the planet then rotate with it
	if (false &&
		lastJulianDate !== julianDate && 
		orbitPlanetIndex == planets.indexes.Earth &&	//only for earth at the moment ...
		orbitDistance < orbitPlanet.radius * 10) 
	{
		var deltaJulianDate = julianDate - lastJulianDate;
		var deltaAngle = quat.create();
		quat.rotateZ(deltaAngle, deltaAngle, deltaJulianDate * 2 * Math.PI);
		vec3TransformQuat(glutil.view.pos, glutil.view.pos, deltaAngle);
		quat.multiply(glutil.view.angle, deltaAngle, glutil.view.angle); 
	}
	lastJulianDate = julianDate;

	glutil.scene.setupMatrices();
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	drawScene();
	glutil.clearAlpha();

	requestAnimFrame(update);
}

