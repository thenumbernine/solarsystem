//all the unused init code for various other system goes here

import {cfg} from './globals.js';

// with a little help from
// http://www.cv.nrao.edu/~rfisher/Ephemerides/ephem_descr.html
// ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/jpl_eph.cc
// http://www.astro-phys.com/js/astro/api.js
class EphemerisData {

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

	interp(coeff, numCoeffs, intervalLength, time) {
		let pc = new Float64Array(18);//ffi.new('double[18]')
		let vc = new Float64Array(18);//ffi.new('double[18]')
	
		if (time < 0 || time > 1) throw 'out of bounds';
		// tc is the normalized chebyshev time (-1 <= tc <= 1)
		let tc = 2 * time - 1;
		
		pc[0] = 1;
		pc[1] = tc;

		let twot = tc + tc
		for (let i = 2; i < numCoeffs; ++i) {
			pc[i] = twot * pc[i-1] - pc[i-2];
		}
		
		let pos = 0;
		for (let i = numCoeffs-1; i >= 0; --i) {
			pos += pc[i] * coeff[i];
		}
		
		vc[0] = 0;
		vc[1] = 1;
		vc[2] = twot + twot;
		for (let i=3; i < numCoeffs; ++i) {
			vc[i] = twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
		}
		let vel = 0;
		for (let i = numCoeffs-1; i >= 1; --i) {
			vel += vc[i] * coeff[i];
		}
		vel *= 2 / intervalLength;
		return [pos, vel];
	}

	posVel(planetIndex, timeOrigin, timeOffset) {
		timeOffset = +timeOffset;
		let coeffBuffer = this.getCoeffBuffer(timeOrigin, timeOffset);
		let coeff, numCoeffs, subintervalLength, subintervalFrac = getCoeffSubinterval(planetIndex, coeffBuffer, timeOrigin, timeOffset);
		let pos = [];
		let vel = [];
		let tmp;
		for (let i = 0; i < 3; ++i) {
			tmp = interp(coeff, numCoeffs, subintervalLength, subintervalFrac);
			pos[i] = tmp[0];
			vel[i] = tmp[1];
			coeff = coeff + numCoeffs
		}
		return [pos, vel];
	}

	mercury(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Mercury, timeOrigin, timeOffset); }
	venus(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Venus, timeOrigin, timeOffset); }
	mars(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Mars, timeOrigin, timeOffset); }
	jupiter(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Jupiter, timeOrigin, timeOffset); }
	saturn(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Saturn, timeOrigin, timeOffset); }
	uranus(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Uranus, timeOrigin, timeOffset); }
	neptune(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Neptune, timeOrigin, timeOffset); }
	pluto(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Pluto, timeOrigin, timeOffset); }
	sun(timeOrigin, timeOffset) { return this.posVel(this.objIndexForName.Sun, timeOrigin, timeOffset); }

	earth(timeOrigin, timeOffset) {
		let tmp = this.posVel(this.objIndexForName.EM_Bary, timeOrigin, timeOffset);
		let earthMoonPos = tmp[0];
		let earthMoonVel = tmp[1] 
		let tmp = this.posVel(this.objIndexForName.GeoCMoon, timeOrigin, timeOffset);
		let geoMoonPos = tmp[0];
		let geoMoonVel = tmp[1];
		let scale = 1 / (1 + this.hdr.emrat);
		
		let earthPos = [
			earthMoonPos[0] - geoMoonPos[0] * scale,
			earthMoonPos[1] - geoMoonPos[1] * scale,
			earthMoonPos[2] - geoMoonPos[2] * scale
		];
		let earthVel = [
			earthMoonVel[0] - geoMoonVel[0] * scale,
			earthMoonVel[1] - geoMoonVel[1] * scale,
			earthMoonVel[2] - geoMoonVel[2] * scale
		];

		return [earthPos, earthVel];
	}

	moon(timeOrigin, timeOffset) {
		let tmp = this.posVel(this.objIndexForName.EM_Bary, timeOrigin, timeOffset);
		let earthMoonPos = tmp[0];
		let earthMoonVel = tmp[1];
		let tmp = this.posVel(this.objIndexForName.GeoCMoon, timeOrigin, timeOffset);
		let geoMoonPos = tmp[0];
		let geoMoonVel = tmp[1];
		let scale = 1 / (1 + this.hdr.emrat);
		
		let earthPos = [
			earthMoonPos[0] - geoMoonPos[0] * scale,
			earthMoonPos[1] - geoMoonPos[1] * scale,
			earthMoonPos[2] - geoMoonPos[2] * scale
		];
		let earthVel = [
			earthMoonVel[0] - geoMoonVel[0] * scale,
			earthMoonVel[1] - geoMoonVel[1] * scale,
			earthMoonVel[2] - geoMoonVel[2] * scale
		];
		let moonPos = [
			geoMoonPos[0] + earthPos[0],
			geoMoonPos[1] + earthPos[1],
			geoMoonPos[2] + earthPos[2]
		];
		let moonVel = [
			geoMoonVel[0] + earthVel[0],
			geoMoonVel[1] + earthVel[1],
			geoMoonVel[2] + earthVel[2]
		];
		
		return [moonPos, moonVel];
	}

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
}
EphemerisData.prototype.objIndexForName = {};
EphemerisData.prototype.objNames.forEach((objName,objIndex) => {
	EphemerisData.prototype.objIndexForName[objName] = objIndex;
});
let ephemerisData = new EphemerisData();

function initAstroPhysPlanetsByRecords() {
	//getting records ...
	let url = 'http://www.astro-phys.com/api/de406/records';
	//let url = 'astro-phys-records.json';
	console.log('reading from '+url);
	fetch(url)
	.then(response => {
		if (!response.ok) return Promise.reject('not ok');
		response.json()
		.then(data => {
			//copied from astro-api.  i need to change the apicall but keep the function the same
			let results = data.results;
			let start = data.start;
			let end = data.end;
			for(let body in results) {
				let nchunks = results[body].length
				let c_step = (end - start) / nchunks;
				for(let i = 0; i < nchunks; ++i) {
					let c_start = start + i * c_step;
					let c_end = c_start + c_step;
					let coeffs = results[body][i];
					results[body][i] = new astro.CoeffSet(coeffs, c_start, c_end);
				}
			}
			
			let record = new astro.Record(results, start, end);
			cfg.julianDate = data.date;
			initPlanets = planets.clone();
			initJulianDate = cfg.julianDate;
			cfg.refreshCurrentTimeText();
			
			//record.getStates is giving back the whole coefficients ... ?	
			let positions = record.getPositions(cfg.julianDate);
			planets = new Planets();
			for (let i = 0; i < planets.length; ++i) {
				let planet = planets[i];
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
	}).catch(console.log);
}

//not used at the moment
// convert from json query object from http://www.astro-phys.com
function makePlanetsFromAstroPhysState(astroPhysPlanets) {
	let planets = new Planets();
	for (let i = 0; i < planets.length; ++i) {
		let planet = planets[i];
		let astroPhysPlanet = astroPhysPlanets[planet.name.toLowerCase()];
	
		//convert km to m
		planet.pos[0] = astroPhysPlanet[0][0] * 1000;
		planet.pos[1] = astroPhysPlanet[0][1] * 1000;
		planet.pos[2] = astroPhysPlanet[0][2] * 1000;
		planet.vel[0] = astroPhysPlanet[1][0] * 1000;
		planet.vel[1] = astroPhysPlanet[1][1] * 1000;
		planet.vel[2] = astroPhysPlanet[1][2] * 1000;
	}
	return planets;
}

function initAstroPhysPlanetsByStates() {
	//getting single sets of states ...
	//returning a single state
	let planetNames = [];
	primaryPlanetHorizonIDs.forEach(horizonID => {
		let planetClassPrototype = Planets.prototype.planetClassForHorizonID[horizonID].prototype;
		planetNames.push(planetClassPrototype.name.toLowerCase());
	});

	//works, but not for the extra horizon data
	if (false) {
		let url = 'http://www.astro-phys.com/api/de406/states?bodies=' + planetNames.join(',');
		//let url = 'astro-phys-state.json';
		console.log('reading from '+url);
		fetch(url)
		.then(response => {
			if (!response.ok) return Promise.reject('not ok');
			response.json()
			.then(d => {
				planets = makePlanetsFromAstroPhysState(d.results);
				cfg.julianDate = d.date;
				
				initPlanets = planets.clone();
				initJulianDate = cfg.julianDate;
				cfg.refreshCurrentTimeText();
				
				init2();
			});
		}).catch(console.log);
	}

	if (false) { // astro.getStates not working as well as my own manual call ...
		console.log('planetNames '+JSON.stringify(planetNames));
		console.log('julianDate '+cfg.julianDate);
		astro.getStates(planetNames, cfg.julianDate, function(results, date) {
			console.log('results '+JSON.stringify(results));
			console.log('date '+date);
			planets = makePlanetsFromAstroPhysState(results);
			cfg.julianDate = date;
			
			initPlanets = planets.clone();
			initJulianDate = cfg.julianDate;
			cfg.refreshCurrentTimeText();
			
			init2();
		});
	}
}

//TODO if this gets reinstated then go from clone to copy
//
//ephemeris data is a few hundred megs
//so i'm not putting it on my serverspace
//so don't use makePlanetsFromEphemeris
//
// convert an ephemeris dataset
// date is julian date
//	dir is only used if solarsys.ephemerisData hasn't been initialized
function makePlanetsFromEphemeris(date, denum, dir) {
	if (!ephemerisData.hasInitialized()) {
		ephemerisData.initialize(denum, dir);
	}
	let planets = new Planets();
	for (let i=0; i < planets.length; ++i) {
		let planet = planets[i];
		let tmp = ephemerisData[planet.name.toLowerCase()](date);	// returns pos & vel in terms of km and km/julian day
		let pos = tmp[0];
		let vel = tmp[1];
		
		//convert km to m
		planet.pos[0] = pos[0] * 1000;
		planet.pos[1] = pos[1] * 1000;
		planet.pos[2] = pos[2] * 1000;
		planet.vel[0] = vel[0] * 1000;
		planet.vel[1] = vel[1] * 1000;
		planet.vel[2] = vel[2] * 1000;
	}
	return planets;
}

//TODO if this gets reinstated then go from clone to copy
function updateTimeForEphemeris() {
	if (targetJulianDate !== undefined) {
		let logDate = Math.log(cfg.julianDate);
		let logTarget = Math.log(targetJulianDate);
		let coeff = .5;
		let deltaLog = logTarget - logDate;
		if (deltaLog < 0) {
			let newLogDate = logDate + coeff * deltaLog;
			cfg.julianDate = Math.exp(newLogDate);
		} else {
			cfg.julianDate = targetJulianDate;
			targetJulianDate = undefined; 
		}
		planets = makePlanetsFromEphemeris(cfg.julianDate);
	} else
}


// routines for calculating the earth's angle ...
// ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/astrtime.cc
function getLeapSec(mjd, utc) {
	// oh great there is a leap-second file...
	return 0;
}
function tai2utc(taimjd, tai) {
	let mjd = taimjd;
	let ls = getLeapSec(taimjd, tai);
	let utc = tai - ls / 86400;
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
	let tai = tt - 32.184 / 86400;
	let taimjd = ttmjd;
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
	let ut1mjd = mjd;
	let ut1 = utc + getUt1Offset(mjd, utc) / 86400;
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
	let tmp = utc2ut1(mjd, utc);
	let ut1mjd = tmp[0];
	let ut1 = tmp[1];
	return iauGmst82(2400000.5, ut1mjd + ut1) / (2 * Math.pi);
}
function utc2tai(mjd, utc) {
	let taimjd = mjd;
	let tai = utc;
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
	let gmst = utc2gmst(mjd, utc);
	let tmp = utc2tt(mjd, utc);
	let ttmjd = tmp[0];
	let tt = tmp[1];
	return gmst + iauEqeq94(2400000.5, ttmjd + tt) / (2 * Math.pi);
}
function utc2last(mjd,utc,lon) {
	let lon = 0;	// longitude of our observatory ... from space!
	let last = utc2gast(mjd, utc) + lon / (2 * Math.pi);
	if (last < 0) ++last;	//modulo?
	if (last > 1) --last;
	return last;
}
function getEarthAngle(jd) {
	let jdofs = cfg.julianDate - 2400000.5;
	let ttmjd = Math.floor(jdofs);
	let tt = jdofs - ttmjd;
	let tmp = tt2utc(ttmjd, tt);
	let mjd = tmp[0];
	let utc = tmp[1];
	let arg = utc2last(mjd,utc) * Math.pi * 2;
	return deg;
}
