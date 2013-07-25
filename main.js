var Planet = makeClass({

	init : function(args) {
		if (args !== undefined) {
			if (args.pos !== undefined) {
				this.pos = vec3.clone(args.pos);
			}
			if (args.vel !== undefined) {
				this.vel = vec3.clone(args.vel);
			}
		}
		if (this.pos === undefined) {
			this.pos = vec3.create();
		}
		if (this.vel === undefined) {
			this.vel = vec3.create();
		}
	},

	add : function(other) {
		var a = this;
		var b = other;
		assert(getmetatable(a) == getmetatable(b))
		var p = new a.init();	//init is the ctor is the class
		vec3.add(p.pos, a.pos, b.pos);
		vec3.add(p.vel, a.vel, b.vel);
		return p;
	},

	//scalar multiply
	mul : function(other) {
		var a = this;
		var b = other;
		var p = new a.init();
		vec3.scale(p.pos, a.pos, b);
		vec3.scale(p.vel, a.vel, b);
		return p;
	},

	toString : function() {
		return 'Planet('+JSON.stringify(this)+')';
	},
	
	geodeticPosition : function(lat, lon, height) {
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
			var x = NPlusH * cosPhi * Math.cos(lambda);
			var y = NPlusH * cosPhi * Math.sin(lambda);
			var z = (N * (1 - eccentricitySquared) + height) * Math.sin(phi);
			return [x, y, z];
		} else {
			var NPlusH = equatorialRadius + height;
			var x = NPlusH * cosPhi * Math.cos(lambda);
			var y = NPlusH * cosPhi * Math.sin(lambda);
			var z = NPlusH * Math.sin(phi);
			return [x, y, z];
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
		makeClass({super:Planet, name:'jupiter', mass:1.89813e+27, radius:69911e+3, equatorialRadius:714192e+3, inverseFlattening:1/0.06487}),
		makeClass({super:Planet, name:'saturn', mass:5.68319e+26, radius:58232e+3, equatorialRadius:60268e+3, inverseFlattening:1/0.09796}),
		makeClass({super:Planet, name:'uranus', mass:8.68103e+25, radius:25362e+3, equatorialRadius:25559e+3, inverseFlattening:1/0.02293}),
		makeClass({super:Planet, name:'neptune', mass:1.0241e+26, radius:24624e+3, equatorialRadius:24766e+3, inverseFlattening:1/0.0171}),
		makeClass({super:Planet, name:'pluto', mass:1.314e+22, radius:1151e+3}),
	],

	init : function() {
		for (var i=0; i < this.planetClasses.length; ++i) {
			this[i] = new (this.planetClasses[i])();
			this[i].index = i;	// does anyone use this anymore?
		}
		this.length = this.planetClasses.length;
	},

	// convert from json query object from http://www.astro-phys.com
	fromAstroPhys : function(astroPhysPlanets) {
		var planets = new Planets();	//why isn't this new'ing?
		for (var i = 0; i < planets.length; ++i) {
			var planet = planets[i];
			var name = planet.name;
			var astroPhysPlanet = astroPhysPlanets[name];
			planet.pos = vec3.create();
			planet.vel = vec3.create();
			vec3.scale(planet.pos, astroPhysPlanet[0], 1000);	// in m
			vec3.scale(planet.vel, astroPhysPlanet[1], 1000);	// in m/s
		}
		return planets;
	},

//ephemeris data is a few hundred megs
//so i'm not putting it on my serverspace
//so don't use fromEphemeris
	// convert an ephemeris dataset
	// date is julian date
	//	dir is only used if solarsys.eph hasn't been initialized
	fromEphemeris : function(date, denum, dir) {
		throw "not implemented";
		if (!eph.hasInitialized()) {
			eph.initialize(denum, dir);
		}
		var planets = new Planets();
		for (var i=0; i < planets.length; ++i) {
			var planet = planets[i];
			var tmp = eph[planet.name](date);	// returns pos & vel in terms of km and km/julian day
			var pos = tmp[0];
			var vel = tmp[1];
			planet.pos = vec3.fromValues(pos[0] * 1000, pos[1] * 1000, pos[2] * 1000);
	 		planet.vel = vec3.fromValues(vel[0] * 1000, pos[1] * 1000, pos[2] * 1000);
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
			l += vec3.length(p.pos) + vec3.length(p.vel);
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
var orbitZoomFactor = .9;	// upon mousewheel

var mouse;

var gravitationConstant = 6.6738480e-11;	// m^3 / (kg * s^20

var tidalMin = 0;
var tidalMax = 1;
function calcTidalForce(srcPlanet, pos) {
	var accel = vec3.create();
	var srcPlanetToPos = vec3.create();
	vec3.sub(srcPlanetToPos, pos, srcPlanet.pos);
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		var planet = planets[planetIndex];
		if (planet != srcPlanet) {
			var x = vec3.create();
			vec3.sub(x, pos, planet.pos);
			var xLength = vec3.length(x);
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
	}
	return accel;
}

function calcGravitationForce(pos) {
	var accel = vec3.create();
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		var planet = planets[planetIndex];
		var x = vec3.create();
		vec3.sub(x, pos, planet.pos);
		var xLength = vec3.length(x);
		var xToTheSecond = xLength * xLength;
		var xToTheThird = xLength * xToTheSecond;
		for (var i = 0; i < 3; ++i) {
			accel[i] -= x[i] * (gravitationConstant * planet.mass / xToTheThird);
		}
	}
	return accel;
}

function planetCartesianToSolarSystemBarycentric(planet, x) {
	x = vec3.clone(x);
	vec3.transformQuat(x, planet.angle, x);		// right now planet angles aren't stored in the planet state (for adaptive integration's sake -- how to weight and combine time + space measurements)
	
	// now rotate by axial tilt (along
	var tiltAngle = planet.tiltAngle;
	if (tiltAngle !== undefined) {
		vec3.transformQuat(x, tiltAngle, x);
	}
	
	vec3.add(x, x, planet.pos);

	return x;
}

function planetGeodeticToSolarSystemBarycentric(planet, lat, lon, height) {
	var x = planet.geodeticPosition(lat, lon, height);		// position relative to the planet center
	x = planetCartesianToSolarSystemBarycentric(planet, x);
	return x;
}

function mouseRay() {
	var w = canvas.width;
	var h = canvas.height;
	var aspectRatio = w / h;
	// ray intersect
	var fx = mouse.xf * 2 - 1
	var fy = mouse.yf * 2 - 1
	var tanFovY = Math.tan(GL.view.fovY * Math.PI / 360);
	var v = vec3.fromValues(fx * aspectRatio * tanFovX, fy * tanFovY, -1);
	vec3.transformQuat(v, GL.view.angle, v);
	vec3.normalize(v, v);
	return v;
}

function chooseNewPlanet(mouseDir,doChoose) {
	var bestDot = 0;
	var bestPlanet = undefined; 
	for (var i = 0; i < planets.length; ++i) {
		var planet = planets[i];
		var delta = vec3.create();
		vec3.sub(delta, planet.pos, viewPos);
		vec3.normalize(delta, delta);
		var dot = vec3.dot(delta, mouseDir);
		if (dot > bestDot) {
			bestDot = dot;
			bestPlanet = planet;
		}
	}
	if (bestPlanet !== undefined) {
		planetText.text(bestPlanet.name);
		if (bestPlanet.index !== orbitPlanetIndex && doChoose) {
			orbitPlanetIndex = bestPlanet.index;
			orbitDistance = vec3.distance(viewPos, bestPlanet.pos);
			orbitTargetDistance = 2 * bestPlanet.radius;
		}
	}
}

// routines for calculating the earth's angle ...
// ftp://ftp.cv.nrao.edu/NRAO-staff/rfisher/SSEphem/astrtime.cc
function getLeapSec(mjd, utc) {
	// oh great there is a leap-second file...
	return 0
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
var panel;
var hsvTex;
var colorShader;
var texShader;
var hsvShader;

var showTide = false;
var colorBarHSVRange = 2/3;	// how much of the rainbow to use

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
	GL.resize();
}


$(document).ready(init1);

function init1() {
	panel = $('#panel');

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

	$(window).resize(resize);
	resize();

	var planetNames = [];
	$.each(Planets.prototype.planetClasses, function(_,planetClass) {
		planetNames.push(planetClass.prototype.name);
	});
	//var url = 'http://www.astro-phys.com/api/de406/states?date=' + encodeURIComponent(dateStr) + '&bodies=' + planetNames.join(',');
	var url = 'astro-phys-state.json';
	console.log('reading from '+url);
	$.ajax({
		url : url
	}).done(function(d) {
		//console.log('results');
		//console.log(JSON.stringify(d.results));
		planets = Planets.prototype.fromAstroPhys(d.results);
		
		julianDate = d.date;
		init2();
	});
}

function init2() {
	colorShader = new GL.ShaderProgram({
		vertexCodeID : 'color-vsh',
		fragmentCodeID : 'color-fsh',
		uniforms : {
			color : [1,1,1]
		}
	});

	texShader : new GL.ShaderProgram({
		vertexCodeID : 'tex-vsh',
		fragmentCodeID : 'tex-fsh',
		uniforms : {
			tex : 0
		}
	});

	hsvShader : new GL.ShaderProgram({
		vertexCodeID : 'heat-vsh',
		fragmentCodeID : 'heat-fsh',
		uniforms : {
			hsvTex : 0
		}
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
	var checkDone = function() {
		++planetsDone;
		if (planetsDone == planets.length) {
			init3();
		}
	};
	$.each(planets, function(i,planet) {
		var planetClassPrototype = planet.init.prototype;
		planetClassPrototype.color = colors[planet.name];
		// load texture
		var img = new Image();
		img.onload = function() {
			console.log('loading tex for planet '+planet.name);
			checkDone();
			planetClassPrototype.tex = new GL.Texture2D({
				flipY : true,
				data : img,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				generateMipmap : true
			});
		};
		img.onerror = function() {
			console.log('failed to find texture for planet '+planet.name);
			checkDone();
		};
		img.src = 'textures/'+planet.name+'.png';
		
		planetClassPrototype.angle = quat.create();			// rotation ... only used for earth at the moment

		planetClassPrototype.elementArray = [];
		planetClassPrototype.vertexArray = [];
		planetClassPrototype.texCoordArray = [];
		planetClassPrototype.tideArray = [];

		var latdiv = Math.floor((latMax-latMin)/latStep);
		var londiv = Math.floor((lonMax-lonMin)/lonStep);
		for (var loni=0; loni <= londiv; ++loni) {
			var lon = lonMin + loni * lonStep;
			for (var lati=0; lati <= latdiv; ++lati) {
				var lat = latMin + lati * latStep;
				
				var vtx = planetClassPrototype.geodeticPosition(lat, lon, 0);
				for (var j = 0; j < 3; ++j) {
					planetClassPrototype.vertexArray.push(vtx[j]);
				}
				
				planetClassPrototype.texCoordArray.push(lon / 360 + .5);
				planetClassPrototype.texCoordArray.push(-lat / 180 + .5);

				planetClassPrototype.tideArray.push(0);

				if (loni < londiv && lati < latdiv) {
					for (var j = 0; j < quad.length; ++j) {
						var ofs = quad[j];
						var index = (lati + ofs[0]) + (latdiv + 1) * (loni + ofs[1]);
						planetClassPrototype.elementArray.push(index);
					}
				}
			}
		}

		planetClassPrototype.indexBuffer = new GL.ElementArrayBuffer({data:planetClassPrototype.elementArray});
		planetClassPrototype.vertexBuffer = new GL.ArrayBuffer({data:planetClassPrototype.vertexArray});
		planetClassPrototype.texCoordBUffer = new GL.ArrayBuffer({dim:2, data:planetClassPrototype.texCoordArray});
		planetClassPrototype.tideBuffer = new GL.ArrayBuffer({dim:1, data:planetClassPrototype.tideArray, keep:true, usage:gl.DYNAMIC_DRAW});

		planetClassPrototype.sceneObj = new GL.SceneObject({
			mode : gl.TRIANGLES,
			indexes : planetClassPrototype.indexBuffer,
			attrs : {
				vertex : planetClassPrototype.vertexBuffer,
				texCoords : planetClassPrototype.texCoordBUffer,
				tide : planetClassPrototype.tideBuffer
			},
			shader : planetClassPrototype.tex !== undefined ? texShader : colorShader,
			texs : [planetClassPrototype.tex],
			uniforms : {
				color : planetClassPrototype.color
			},
			pos : vec3.create(),
			angle : quat.create()
		});
	});
}

function init3() {
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL);
	gl.clearColor(0,0,0,0);

	hsvTex = new GL.HSVTexture(256);

	var earth = planets[planets.indexes.earth];
	orbitPlanetIndex = earth.index;
	orbitTargetDistance = 2. * earth.radius;
	orbitDistance = orbitTargetDistance;

	mouse = new Mouse3D({
		pressObj : canvas,
		mousedown : function() {
		},
		move : function() {
			var mouseDir = mouseRay();
			chooseNewPlanet(mouseDir, false);
		},
		zoom : function(zoomChange) {
			var scale = Math.exp(-orbitZoomFactor * zoomChange);
			orbitTargetDistance *= scale;
		},
		click : function() {
			var mouseDir = mouseRay();
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
			quat.mul(viewAngle, viewAngle, r);
			quat.normalize(viewAngle, viewAngle);
		}
	}

	// track ball orbit
	
	var earth = planets[planets.indexes.earth];
	
	var orbitCenter;
	if (orbitGeodeticLocation !== undefined) {
		orbitCenter = planetGeodeticToSolarSystemBarycentric(
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
	}
*/

	//drawScene:
	for (var planetIndex = 0; planetIndex < planets.length; ++planetIndex) {
		var planet = planets[planetIndex];
		var planetClassPrototype = planet.init.prototype;
	
		planet.visRatio = planet.radius / vec3.distance(planet.pos, GL.view.pos);
		//if (planet.visRatio >= .005) draw the whole planet, otherwise just a dot
		
		//update scene object
		vec3.copy(planetClassPrototype.sceneObj.pos, planet.pos);
		quat.copy(planetClassPrototype.sceneObj.angle, planet.angle);
	
		//update tide attributes
		if (!showTide) {
			planetClassPrototype.sceneObj.shader = planetClassPrototype.tex !== undefined ? texShader : colorShader;
			planetClassPrototype.sceneObj.texs[0] = planetClassPrototype.tex;
		} else {
			planetClassPrototype.sceneObj.shader = hsvShader;
			planetClassPrototype.sceneObj.texs[0] = hsvTex;
			if (planetClassPrototype.lastTideCalcDate !== julianDate) {
				planetClassPrototype.lastTideCalcDate = julianDate;
				var tidalMin = undefined;
				var tidalMax = undefined;
				var vertexIndex = 0;
				for (var tideIndex = 0; tideIndex < planet.tideArray.length; ++tideIndex) {
					var planetX = vec3.fromValues(
						planet.vertexArray[vertexIndex++],
						planet.vertexArray[vertexIndex++],
						planet.vertexArray[vertexIndex++]);
					var x = planetCartesianToSolarSystemBarycentric(planet, planetX);
					var accel = calcTidalForce(planet, x);
					//var accel = calcGravitationalForce(x);
					var norm = vec3.create();
					vec3.sub(norm, x, planet.pos);
					vec3.normalize(norm, norm);
					//var toTheMoon = (planets[planets.indexes.moon].pos - x):normalize()
					//var normCrossMoon = norm:cross(toTheMoon)	//points upwards, tangent, right angle to norm and moon
					//var tangentTowardsMoon = normCrossMoon:cross(norm)
					//var tidalAccel = accel:dot(tangentTowardsMoon)
					//var tidalAccel = accel:dot(norm)	// normal component
					//var tidalAccel = accel:dot(toTheMoon)	// moonward component
					var dot = vec3.dot(accel, norm);
					var proj = vec3.create();
					vec3.scale(proj, norm, dot);
					var tidalAccel = vec3.distance(accel, proj);	// tangent component
					//var tidalAccel = accel:length()	// magnitude
					var t = tidalAccel;
					if (tidalMin === undefined || t < tidalMin) tidalMin = t;
					if (tidalMax === undefined || t > tidalMax) tidalMax = t;
					planet.tideArray[tideIndex] = t
				}
				for (var i = 0; i < planet.tideArray.length; ++i) {
					planet.tideArray[i] = (255/256 - (planet.tideArray[i] - tidalMin) / (tidalMax - tidalMin) * 254/256) * colorBarHSVRange;
				}
				planet.tideBuffer.updateData();
				gl.bindBuffer(gl.ARRAY_BUFFER, null);
			}
		}
	}
	//should draw these manually while debugging ...
	GL.draw();
//	requestAnimFrame(update);
}

