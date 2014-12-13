//getting rid of the old way incrementally
var CALCULATE_TIDES_WITH_GPU = true;

//TODO put this in js/gl-util-kernel.js


GLUtil.prototype.oninit.push(function() {
	var glutil = this;
	glutil.KernelShader = makeClass({
		super : glutil.ShaderProgram,
		init : function(args) {
			
			var varyingCodePrefix = 'varying vec2 pos;\n';

			var fragmentCodePrefix = '';
			var uniforms = {};
			if (args.uniforms !== undefined) {
				$.each(args.uniforms, function(uniformName, uniformType) {
					if ($.isArray(uniformType)) {
						//save initial value
						uniforms[uniformName] = uniformType[1];
						uniformType = uniformType[0];
					}
					fragmentCodePrefix += 'uniform '+uniformType+' '+uniformName+';\n';
				});
			}
			if (args.texs !== undefined) {
				for (var i = 0; i < args.texs.length; ++i) {
					var v = args.texs[i];
					var name, vartype;
					if (typeof(v) == 'string') {
						name = v;
						vartype = 'sampler2D';
					} else {
						name = v[0];
						vartype = v[1];
					}
					fragmentCodePrefix += 'uniform '+vartype+' '+name+';\n';
					uniforms[name] = i;
				}
			}


			if (!glutil.KernelShader.prototype.kernelVertexShader) {
				glutil.KernelShader.prototype.kernelVertexShader = new glutil.VertexShader({
					code : 
						glutil.vertexPrecision + 
						varyingCodePrefix +
						mlstr(function(){/*
attribute vec2 vertex;
attribute vec2 texCoord;
void main() {
	pos = texCoord; 
	gl_Position = vec4(vertex, 0., 1.);
}
*/})
				});	
			}

			args.vertexShader = glutil.KernelShader.prototype.kernelVertexShader;
			args.fragmentCode = glutil.fragmentPrecision + varyingCodePrefix + fragmentCodePrefix + args.code;
			delete args.code;
			args.uniforms = uniforms;	
			glutil.KernelShader.super.call(this, args);
		}
	});
});


var toGLSLFloat = function(x) {
	x = ''+x;
	if (x.indexOf('.') == -1) x = x + '.';
	return x;
};

function setCSSRotation(obj, degrees) {
	obj.css('-webkit-transform','rotate('+degrees+'deg)');
	obj.css('-moz-transform','rotate('+degrees+'deg)');
	obj.css('transform','rotate('+degrees+'deg)');
}

var glutil;

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

var Planet = makeClass({

	init : function(args) {
		if (args !== undefined) {
			if (args.pos !== undefined) {
				this.pos = [args.pos[0], args.pos[1], args.pos[2]];
			}
			if (args.vel !== undefined) {
				this.vel = [args.vel[0], args.vel[1], args.vel[2]];
			}
			if (args.angle !== undefined) {
				this.angle = [args.angle[0], args.angle[1], args.angle[2], args.angle[3]];
			}
		}
		if (this.pos === undefined) {
			this.pos = [0,0,0];
		}
		if (this.vel === undefined) {
			this.vel = [0,0,0];
		}
		if (this.angle === undefined) {
			this.angle = [0,0,0,1];
		}
		this.tiltAngle = [0,0,0,1];
	},

	clone : function() {
		var p = new this.init();
		vec3.copy(p.pos, this.pos);
		vec3.copy(p.vel, this.vel);
		quat.copy(p.angle, this.angle);
		quat.copy(p.tiltAngle, this.tiltAngle);
		return p;
	},

	copy : function(src) {
		this.pos[0] = src.pos[0];
		this.pos[1] = src.pos[1];
		this.pos[2] = src.pos[2];
		this.vel[0] = src.vel[0];
		this.vel[1] = src.vel[1];
		this.vel[2] = src.vel[2];
	},

	//no longer used for mesh construction -- that all goes on in the shader
	//this is only used for getting positions for updating the tidal array calculations
	// and for (geosynchronously) orbitting geodetic locations (which is currently disabled)
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
			destX[2] = (N * (1 - eccentricitySquared) + height) * sinPhi;
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


var StarSystem = makeClass({
	init : function() {
		this.pos = [0,0,0];
		this.vel = [0,0,0];
		this.angle = [0,0,0,1];
		this.planets = [];
		this.stars = [];
	},

	//this calls 'createPlanetsFBOTex' which shouldn't be called until after WebGL init
	doneBuildingPlanets : function() {
		this.buildIndexes();
		this.mapParents();
		this.createPlanetsFBOTex();
	},

	/*
	need the following texture per-planet:
		position (and maybe velocity) x2 for front and back buffers
		mass
		keplerian orbital elements
			A
			B
			eccentricity
			eccentricAnomaly
	
		for now just do position and mass 
		-- and manually update them
		-- and use them for the planetSurfaceCalcuationShader
	
	TODO what to do about comets ...
	for now just don't add them.
	they don't get any tidal/gravitational force calculations anywyas.
		
	eventually set up a mask tex to hold what planets the user flags on/off
	*/
	createPlanetsFBOTex : function() {
		
		this.planetStateTex = new glutil.Texture2D({
			internalFormat : gl.RGBA,		//xyz = pos, w = mass.  double this up if you need more precision
			type : gl.FLOAT,
			width : 1,	//might double this if we need more accuracy
			height : this.planets.length,	//npo2 ...
			magFilter : gl.NEAREST,
			minFilter : gl.NEAREST,
			wrap : {
				s : gl.CLAMP_TO_EDGE,
				t : gl.CLAMP_TO_EDGE
			}
		});
	},

	//builds solarSystem.indexes[planetName]
	//and remaps solarSystem.planets[i].parent from a name to an index (why not a pointer?)
	buildIndexes : function() {
		this.indexes = {};
		for (var i = 0; i < this.planets.length; ++i) {
			this.planets[i].index = i;
			var planet = this.planets[i];
			this.indexes[planet.name] = i;
			//while we're here...
			planet.starSystem = this;
		}
	},

	//map parent field from name to index (or should it be to object?)
	mapParents : function() {
		//convert parent from name to class (or undefined if no such name exists)
		for (var i = 0; i < this.planets.length; ++i) {
			var planet = this.planets[i];
			if (planet.parent !== undefined) {
				assert(typeof(planet.parent) === 'string');
				var index = assertExists(this.indexes, planet.parent);
				planet.parent = assertExists(this.planets, index);
			}
		}
	},

	//this is the old Planets behavior
	// which I might recreate
	// combination of Array and integration functions
	clonePlanets : function() {
		var planets = [];
		for (var i = 0; i < this.planets.length; ++i) {
			planets[i] = this.planets[i].clone();
		}
		return planets;
	},

	//static
	copyPlanets : function(dest, src) {
		assert(dest.length == src.length);
		for (var i = 0; i < src.length; ++i) {
			dest[i].copy(src[i]);
		}
	}

});

//our solar system
var SolarSystem = makeClass({
	super : StarSystem,
	name : 'Solar System',
	init : function() {
		SolarSystem.super.apply(this, arguments);

		//add our initial planets ...
		this.planets.push(mergeInto(new Planet(), {id:10, name:'Sun', mass:1.9891e+30, radius:6.960e+8, type:'star'}));
		this.planets.push(mergeInto(new Planet(), {id:199, name:'Mercury', parent:'Sun', mass:3.302e+23, radius:2.440e+6, equatorialRadius:2440e+3, rotationPeriod:58.6462, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:299, name:'Venus', parent:'Sun', mass:4.8685e+24, radius:6.0518e+6, equatorialRadius:6051.893e+3, rotationPeriod:-243.0185, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:399, name:'Earth', parent:'Sun', mass:5.9736e+24, radius:6.37101e+6, equatorialRadius:6378.136e+3, inverseFlattening:298.257223563, rotationPeriod:1, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:301, name:'Moon', parent:'Earth', mass:7.349e+22, radius:1.73753e+6, rotationPeriod:30.25, rotationOffset:2.3247785636564, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:499, name:'Mars', parent:'Sun', mass:6.4185e+23, radius:3.3899e+6, equatorialRadius:3397e+3, inverseFlattening:154.409, rotationPeriod:24.622962/24, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:599, name:'Jupiter', parent:'Sun', mass:1.89813e+27, radius:6.9911e+7, equatorialRadius:71492e+3, inverseFlattening:1/0.06487, ringRadiusRange:[102200000,227000000], type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:699, name:'Saturn', parent:'Sun', mass:5.68319e+26, radius:5.8232e+7, equatorialRadius:60268e+3, inverseFlattening:1/0.09796, ringRadiusRange:[74510000,140390000], type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:799, name:'Uranus', parent:'Sun', mass:8.68103e+25, radius:2.5362e+7, equatorialRadius:25559e+3, inverseFlattening:1/0.02293, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:899, name:'Neptune', parent:'Sun', mass:1.0241e+26, radius:2.4624e+7, equatorialRadius:24766e+3, inverseFlattening:1/0.0171, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:999, name:'Pluto', parent:'Sun', mass:1.314e+22, radius:1.151e+6, type:'planet'}));

		//sun is our only star
		this.stars.push(this.planets[0]);

		//start out with the current planets
		var currentIDs = {};
		for (var i = 0; i < this.planets.length; ++i) {
			currentIDs[this.planets[i].id] = true;
		}

		//I need to fix up my export script ...
		if (horizonsDynamicData.coords.length != horizonsStaticData.length) throw 'static to dynamic data lengths differ: dynamic has '+horizonsDynamicData.coords.length+' while static has '+horizonsStaticData.length;
		for (var i = 0; i < horizonsDynamicData.coords.length; ++i) {
			var dynamicData = horizonsDynamicData.coords[i];
			var staticData = horizonsStaticData[i];
			assert(dynamicData.id == staticData.id, "got some bad data");
			if (dynamicData.name.match(/^L[1-5]/)) {
				//391-395 are Lagrangian points - skip
			} else if (dynamicData.name.match(/Barycenter$/)) {
				//1-9 are Barycenter points - skip
			} else {
				if (currentIDs[dynamicData.id]) {
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
					this.planets.push(mergeInto(new Planet(), {
						id:dynamicData.id,
						name:dynamicData.name,
						mass:staticData.mass,
						radius:staticData.radius,
						equatorialRadius:staticData.equatorialRadius,
						inverseFlattening:staticData.inverseFlattening,
						parent:staticData.parent
					}));
					currentIDs[dynamicData.id] = true;
				}
			}
		}


		//used for getting dynamic data, and for texture loading
		this.planetForHorizonID = {};
		for (var i = 0; i < this.planets.length; ++i) {
			var planet = this.planets[i];
			this.planetForHorizonID[planet.id] = planet;
		}


		//extra horizon data
		//has to be telnetted out of jpl and that takes about 5 mins for everything
		//some dude must have typed in every one of the 200 info cards by hand
		// because there's nothing consistent about formatting or variable names
		//so I'm thinking cron job to update and then integrate to extrapolate for later time.

		for (var i = 0; i < horizonsDynamicData.coords.length; ++i) {
			var dynamicData = horizonsDynamicData.coords[i];
			var planet = this.planetForHorizonID[dynamicData.id];
			if (planet) {	//excluding the BCC and Ln points
				for (var j = 0; j < 3; ++j) {
					//convert km to m
					planet.pos[j] = dynamicData.pos[j] * 1000;
					planet.vel[j] = dynamicData.vel[j] * 1000;
				}
			}
		}

		//once we're done making planets, make a copy of the init
		this.initPlanets = this.clonePlanets();
	}
});

var starSystems = [];
var starSystemForNames = {};
var solarSystem;	//the one and only.  don't construct until after WebGL init so we can populate our float tex for the planets


//TODO merge with starSystems[] ... keep the StarField for point rendering of all StarSystems (or make it a Galaxy object, honestly, that's where thignsn are going)
// and remove StarInField ... make that just StarSystem (even for zero-planet systems)
//
//only instanciate these for the named stars.  87 in all.
var StarInField = makeClass({
	init : function(args) {
	}
});

var StarField = makeClass({});
var starfield = undefined;


var orbitPathResolution = 500;
var ringResolution = 200;

var julianDate = 0;
var lastJulianDate = 0;
var initJulianDate = 0;
var targetJulianDate = undefined;
var dateTime = Date.now();

// track ball motion variables
var mouseOverTarget;
var orbitStarSystem;	//only do surface calculations for what star system we are in
var orbitTarget;
var orbitGeodeticLocation;
var orbitDistance;
var orbitOffset = [0,0,0];
var orbitTargetDistance;
var orbitZoomFactor = .0003;	// upon mousewheel

var mouse;
var mouseDir;
var skyCubeObj;

var skyTexFilenamePrefixes = [
	'textures/sky-visible-cube-xp-',
	'textures/sky-visible-cube-xn-',
	'textures/sky-visible-cube-yp-',
	'textures/sky-visible-cube-yn-',
	'textures/sky-visible-cube-zp-',
	'textures/sky-visible-cube-zn-'
];

var glMaxCubeMapTextureSize;

var speedOfLight = 299792458;	// m/s
var gravitationalConstant = 6.6738480e-11;	// m^3 / (kg * s^2)

//1 = c m/s <-> c = s/m
var secondsPerMeter = speedOfLight;

//1 = G m^3 / (kg s^2) <=> G m^3 / (c m)^2 = kg <=> G/c^2 = kg/m
var kilogramsPerMeter = gravitationalConstant / (speedOfLight * speedOfLight);

var Metric = makeClass();

//low gravity newtonian approximation
var NewtonApproximateMetric = makeClass({
	super : Metric,

	/*
	geodesic calculation:
	x''^u = -Conn^u_ab x'^a x'^b

	Newtonian approximation:
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
			var r3 = r * r2;
			accel[0] -= x * (gravitationalConstant * planet.mass / r3);
			accel[1] -= y * (gravitationalConstant * planet.mass / r3);
			accel[2] -= z * (gravitationalConstant * planet.mass / r3);
		}
	},

	/*
	geodesic deviation:
	x''^i = R^i_jkl n^j dx^k n^l

	Newtonian approximation:
	R^i_tjt = R^i_ttj = phi_,ij
	But what if phi changes wrt time? then phi_,tt is nonzero, right? How does our Riemann metric change?
	*/
	calcTidal : function(accel, pos, srcPlanet) {
		var x = [];
		var n = [];

		n[0] = pos[0] - srcPlanet.pos[0];
		n[1] = pos[1] - srcPlanet.pos[1];
		n[2] = pos[2] - srcPlanet.pos[2];

		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.index === srcPlanet.index) continue;
			if (planet.mass === undefined) continue;

			x[0] = pos[0] - planet.pos[0];
			x[1] = pos[1] - planet.pos[1];
			x[2] = pos[2] - planet.pos[2];
			var r = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
			var r2 = r * r;
			var r3 = r * r2;
			var r4 = r * r3;
			var r5 = r * r4;

			var xDotN = x[0] * n[0] + x[1] * n[1] + x[2] * n[2];

			for (i = 0; i < 3; ++i) {
				accel[i] += gravitationalConstant * planet.mass * (3 * xDotN * x[i] / r5 - n[i] / r3);
			}
		}
	}
});

//rotation-less spherical body
var SchwarzschildMetric = makeClass({
	super : Metric,

	/*
	geodesic calculation:
	x''^u = -Conn^u_ab x'^a x'^b
	*/
	calcGravity : function(accel, pos) {
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			if (!planetInfluences[planetIndex]) continue;
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.mass === undefined) continue;
			var M = planet.mass;
			var x = pos[0] - planet.pos[0];
			var y = pos[1] - planet.pos[1];
			var z = pos[2] - planet.pos[2];
			var r = Math.sqrt(x*x + y*y + z*z);

			var M2 = M * M;
			var M3 = M2 * M;
			var x2 = x * x;
			var y2 = y * y;
			var z2 = z * z;
			var x3 = x2 * x;
			var y3 = y2 * y;
			var z3 = z2 * z;
			var r2 = r * r;
			var r3 = r * r2;
			var r5 = r3 * r2;
			var r6 = r5 * r;
			var r8 = r6 * r2;
			var r9 = r6 * r3;

			//TODO simplify more!
			//something needs to factor out, and until it does you'll get spikes in the evaluation where x^2 = y^2 = z^2
			accel[0] -= speedOfLight * (4 * x * y2 * z2 * M3 + (-2 * r3 * x * z2-2 * r3 * x * y2) * M2 + r6 * x * M)/(40 * x2 * y2 * z2 * M3 + ((-12 * r3 * y2-12 * r3 * x2) * z2-12 * r3 * x2 * y2) * M2 + 2 * r8 * M + r9);
			accel[1] -= speedOfLight * (4 * x2 * y * z2 * M3 + (-2 * r3 * y * z2-2 * r3 * x2 * y) * M2 + r6 * y * M)/(40 * x2 * y2 * z2 * M3 + ((-12 * r3 * y2-12 * r3 * x2) * z2-12 * r3 * x2 * y2) * M2 + 2 * r8 * M + r9);
			accel[2] -= speedOfLight * (4 * x2 * y2 * z * M3 + (-2 * r3 * z * y2-2 * r3 * z * x2) * M2 + r6 * z * M)/(40 * x2 * y2 * z2 * M3 + ((-12 * r3 * y2-12 * r3 * x2) * z2-12 * r3 * x2 * y2) * M2 + 2 * r8 * M + r9);
		}
	}
});

//uncharged, rotating spherical body
var KerrMetric = makeClass({
	super : Metric
});

var metric = new NewtonApproximateMetric();
//var metric = new SchwarzschildMetric();


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
	vec3TransformQuat(destX, srcX, planet.angle);
	destX[0] += planet.pos[0];
	destX[1] += planet.pos[1];
	destX[2] += planet.pos[2];
}

function planetGeodeticToSolarSystemBarycentric(destX, planet, lat, lon, height) {
	planet.geodeticPosition(destX, lat, lon, height);		// position relative to the planet center
	planetCartesianToSolarSystemBarycentric(destX, destX, planet);
}

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

function chooseNewOrbitObject(mouseDir, doChoose) {
	var bestDot = Math.cos(Math.deg(10));
	var bestTarget = undefined;
	/*
	list contains:
	.length
	[i]
		.hide
		.pos[j]
	*/
	var processList = function(list) {
		for (var i = 0; i < list.length; ++i) {
			var target = list[i];

			//no need to select in starfield the system we're already orbitting
			if (list === starfield && target == orbitStarSystem) continue;

			if (target.hide) continue;
			var deltaX = target.pos[0] - glutil.view.pos[0] - orbitTarget.pos[0];
			var deltaY = target.pos[1] - glutil.view.pos[1] - orbitTarget.pos[1];
			var deltaZ = target.pos[2] - glutil.view.pos[2] - orbitTarget.pos[2];
			var deltaLength = Math.sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
			deltaX /= deltaLength;
			deltaY /= deltaLength;
			deltaZ /= deltaLength;
			var dot = deltaX * mouseDir[0] + deltaY * mouseDir[1] + deltaZ * mouseDir[2];
			if (dot > bestDot) {
				bestDot = dot;
				bestTarget = target;
			}
		}
	};

	//only check the star system we're in for planet clicks
	processList(orbitStarSystem.planets);

	//for larger-scale clicks, use the starfield
/* gets annoying, trying to click on a planet and catching a star.  just use the side menu for selecting new stars. 
	if (showStars) {
		processList(starSystems);
	}
*/
	
	mouseOverTarget = undefined;
	if (bestTarget !== undefined) {
		$('#hoverTargetText').text(bestTarget.name);
		mouseOverTarget = bestTarget;
		if (bestTarget !== orbitTarget && doChoose) {
			setOrbitTarget(bestTarget);
			refreshMeasureText();
		}
	}
}

function refreshMeasureText() {
	$('#measureMin').text(orbitTarget.measureMin === undefined ? '' : (orbitTarget.measureMin.toExponential() + ' m/s^2'));
	$('#measureMax').text(orbitTarget.measureMax === undefined ? '' : (orbitTarget.measureMax.toExponential() + ' m/s^2'));
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

var fbo;
var hsvTex;
var colorShader;
var latLonShader;
var planetHeatMapAttrShader;
var planetHeatMapTexShader;
var orbitPathShader;

var pointObj;
var planetSceneObj;
var planetLatLonObj;

var colorIndexMin = 2.;
var colorIndexMax = -.4;
var colorIndexTex;
var colorIndexShader;	//used for rendering point stars

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
var showPlanetsAsDistantPoints = true;
var showOrbits = true;
var showStars = true;
var starsVisibleMagnitudeBias = 4;
var gravityWellScaleNormalized = true;
var gravityWellScaleFixed = false;
var gravityWellScaleFixedValue = 2000;
var gravityWellRadialMinLog100 = -1;
var gravityWellRadialMaxLog100 = 2;

var heatAlpha = .5;
var colorBarHSVRange = 2/3;	// how much of the rainbow to use

var depthConstant = 1e-6;//2 / Math.log(1e+7 + 1);

var integrationPaused = true;
var defaultIntegrateTimeStep = 1/(24*60);
var integrateTimeStep = defaultIntegrateTimeStep;

var updatePlanetClassSceneObj;
(function(){
	var x = [];
	var accel = [];
	var norm = [];
	var proj = [];
	updatePlanetClassSceneObj = function(planet) {
		var planetShaders = getPlanetShadersForNumberOfStars(planet.starSystem.stars.length);

		//update tide attributes
		var useOverlay = displayMethod != 'None';
		if (!useOverlay) {
			if (planet.tex === undefined) {
				planet.sceneObj.shader = planetShaders.colorShader;
				planet.sceneObj.texs.length = 0;
			} else {
				if (planet.ringObj) {
					planet.sceneObj.shader = planetShaders.ringShadowShader;
					planet.sceneObj.texs.length = 2;
					if (planet.ringTransparencyTex !== undefined) {
						planet.sceneObj.texs[1] = planet.ringTransparencyTex;
					} else if (planet.ringColorTex !== undefined) {
						planet.sceneObj.texs[1] = planet.ringColorTex;
					} else {
						//...or instead of throwing, should I just assume full-alpha?
						throw 'planet has a ringObj but no ring texture';
					}
				} else {
					planet.sceneObj.shader = planetShaders.texShader;
					planet.sceneObj.texs.length = 1;
				}
				planet.sceneObj.texs[0] = planet.tex;
			}
		} else {
if (!CALCULATE_TIDES_WITH_GPU) {
			planet.sceneObj.shader = planetHeatMapAttrShader;
			planet.sceneObj.texs.length = 2;
			planet.sceneObj.texs[0] = planet.tex;
			planet.sceneObj.texs[1] = hsvTex;
} else {
			planet.sceneObj.shader = planetHeatMapTexShader;
			planet.sceneObj.texs.length = 3;
			planet.sceneObj.texs[0] = planet.tex;
			planet.sceneObj.texs[1] = planet.tideTex;
			planet.sceneObj.texs[2] = hsvTex;
}
			//and update calculated variable if it is out of date ...
			if (planet.lastMeasureCalcDate !== julianDate) {
				planet.lastMeasureCalcDate = julianDate;
				
				// old way -- calculate on CPU, upload to vertex buffer 
if (!CALCULATE_TIDES_WITH_GPU) {
				
				var measureMin = undefined;
				var measureMax = undefined;
				var vertexIndex = 0;
				for (var tideIndex = 0; tideIndex < planet.tideBuffer.data.length; ++tideIndex) {
					var lat = planet.sceneObj.attrs.vertex.data[vertexIndex++];
					var lon = planet.sceneObj.attrs.vertex.data[vertexIndex++];
					planetGeodeticToSolarSystemBarycentric(x, planet, lat, lon, 0);

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

					if (measureMin === undefined || t < measureMin) measureMin = t;
					if (measureMax === undefined || t > measureMax) measureMax = t;
					planet.tideBuffer.data[tideIndex] = t;
				}
				planet.measureMin = measureMin;
				planet.measureMax = measureMax;
				for (var i = 0; i < planet.tideBuffer.data.length; ++i) {
					planet.tideBuffer.data[i] =
						(255/256 - (planet.tideBuffer.data[i] - measureMin)
							/ (measureMax - measureMin) * 254/256) * colorBarHSVRange;
				}
				planet.tideBuffer.updateData();
				//if it updated...
				if (planet == orbitTarget) {
					refreshMeasureText();
				}

} else {		//new way -- update planet state buffer to reflect position & mass
				// (could store tide values in a texture then reduce to find min/max, but this means a texture per planet ... that's lots of textures ...)
				// (could compute this on the fly, but then there's no easy way to find the min/max ... )
				// (could store as texture, and just keep the texture small ... around the size that the tideBuffer already is ... 36x72 )

				//fbo render to the tide tex to calcuate the float values 
		
				//look in the shader's code for which does what
				var flags;
				switch (displayMethod) {
				case 'Tangent Tidal':
					flags = [true, false, false, false];
					break;
				case 'Normal Tidal':
					flags = [false, true, false, false];
					break;
				case 'Total Tidal':
					flags = [true, true, false, false];
					break;
				case 'Tangent Gravitational':
					flags = [false, false, true, false];
					break;
				case 'Normal Gravitational':
					flags = [false, false, false, true];
					flags[3] = true;
					break;
				case 'Total Gravitational':
					flags = [false, false, true, true];
					break;
				case 'Tangent Total':
					flags = [true, false, true, false];
					break;
				case 'Normal Total':
					flags = [false, true, false, true];
					break;
				case 'Total':
					flags = [true, true, true, true];
					break;
				}

				fbo.setColorAttachmentTex2D(0, planet.tideTex);
				gl.viewport(0, 0, tideTexWidth, tideTexHeight);
				gl.disable(gl.DEPTH_TEST);
				gl.disable(gl.CULL_FACE);
				fbo.draw({
					callback : function() {
						//clear is working, but the render is not ...
						//gl.clearColor(.5,0,0,1);
						//gl.clear(gl.COLOR_BUFFER_BIT);
						quadObj.draw({
							shader : planetSurfaceCalculationShader,
							uniforms : {
								pos : planet.pos,
								angle : planet.angle,
								equatorialRadius : planet.equatorialRadius !== undefined ? planet.equatorialRadius : planet.radius,
								inverseFlattening : planet.inverseFlattening !== undefined ? planet.inverseFlattening : .5,
								sourcePlanetIndex : planet.index,
								planetStateTexHeight : orbitStarSystem.planetStateTex.height,
								flags : flags
							},
							texs : [orbitStarSystem.planetStateTex],
						});
					}
				});

				//...then min/max reduce

				var reduce = function(kernelShader, src) {
					var width = tideTexWidth;
					var height = tideTexHeight;
					var dstIndex = 0; 
					var current = src;

					while (width > 1 && height > 1) {
						width >>= 1;
						height >>= 1;
						if (!width || !height) throw 'got a non-square size... TODO add support for this'; 
						gl.viewport(0, 0, width, height);
						
						fbo.bind();
						gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tideReduceTexs[dstIndex].obj, 0);
						fbo.check();
						gl.clear(gl.COLOR_BUFFER_BIT);
						quadObj.draw({
							shader : kernelShader,
							uniforms : {
								texsize : [tideTexWidth, tideTexHeight], 
								viewsize : [width, height]
							},
							texs : [current]
						});
						fbo.unbind();

						current = tideReduceTexs[dstIndex];
						dstIndex = (dstIndex + 1) & 1;
					}
			
					//'current' has our texture

					//now that the viewport is 1x1, run the encode shader on it
					fbo.bind();
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, encodeTempTex.obj, 0);
					fbo.check();
					gl.viewport(0, 0, encodeTempTex.width, encodeTempTex.height);
					quadObj.draw({
						shader : encodeShader[0],
						texs : [current]
					});

					var cflUint8Result = new Uint8Array(4);
					gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, cflUint8Result);
					fbo.unbind();
					
					var cflFloat32Result = new Float32Array(cflUint8Result.buffer);
					var result = cflFloat32Result[0];
					return result;
				};

				planet.measureMin = reduce(minReduceShader, planet.tideTex);
				planet.measureMax = reduce(maxReduceShader, planet.tideTex);
//console.log('measure min', planet.measureMin, 'max', planet.measureMax);
				if (planet == orbitTarget) {
					refreshMeasureText();
				}

				//planet.sceneObj.uniforms.forceMin = planet.measureMax;
				//planet.sceneObj.uniforms.forceMax = (planet.measureMin - planet.measureMax) / colorBarHSVRange + planet.measureMax;
planet.forceMin = planet.measureMax;
planet.forceMax = (planet.measureMin - planet.measureMax) / colorBarHSVRange + planet.measureMax;

				gl.viewport(0, 0, glutil.canvas.width, glutil.canvas.height);
				gl.enable(gl.DEPTH_TEST);
				gl.enable(gl.CULL_FACE);
				
				//...then use the float buffer, min, and max, to do the hsv overlay rendering
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
	var viewfwd = vec3.create();

	//used for new gpu update of tide tex
	var updatePlanetStateBuffer = new Float32Array(1);

	drawScene = function() {

		//update the planet state texture
		// right now it's only used for tide calcs so
		//1) only update it if tide calcs are on
		//2) only include planets that are enabled for tide calcs
		// more discernment later when i make it general purpose
		if (displayMethod != 'None') {
		
			var targetSize = orbitStarSystem.planetStateTex.width * orbitStarSystem.planetStateTex.height * 4;
			
			if (updatePlanetStateBuffer.length != targetSize) {
				updatePlanetStateBuffer = new Float32Array(targetSize);
			}

			//gather pos and mass
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
				updatePlanetStateBuffer[0 + 4 * planetIndex] = planet.pos[0];
				updatePlanetStateBuffer[1 + 4 * planetIndex] = planet.pos[1];
				updatePlanetStateBuffer[2 + 4 * planetIndex] = planet.pos[2];
				
				if (!planetInfluences[planetIndex]) {
					//if we're not using this planet then set the mass to zero
					// works the same as not using it
					updatePlanetStateBuffer[3 + 4 * planetIndex] = 0;
				} else {
					updatePlanetStateBuffer[3 + 4 * planetIndex] = planet.mass === undefined ? 0 : planet.mass;
				}
			}

			//update orbit star system's planet state tex
			orbitStarSystem.planetStateTex.bind();
			
			gl.texSubImage2D(
				gl.TEXTURE_2D,
				0,
				0,
				0,
				orbitStarSystem.planetStateTex.width,
				orbitStarSystem.planetStateTex.height,
				gl.RGBA,
				gl.FLOAT,
				updatePlanetStateBuffer);
			
			orbitStarSystem.planetStateTex.unbind();
		}
		
		mat4.identity(glutil.scene.mvMat);

		quat.conjugate(viewAngleInv, glutil.view.angle);
		mat4.fromQuat(invRotMat, viewAngleInv);
		mat4.multiply(glutil.scene.mvMat, glutil.scene.mvMat, invRotMat);

		//TODO pull from matrix
		vec3.quatZAxis(viewfwd, glutil.view.angle);
		vec3.scale(viewfwd, viewfwd, -1);

		if (skyCubeObj) {
			gl.disable(gl.DEPTH_TEST);
			skyCubeObj.draw();
			gl.enable(gl.DEPTH_TEST);
		}

		vec3.scale(viewPosInv, glutil.view.pos, -1);
		mat4.translate(glutil.scene.mvMat, glutil.scene.mvMat, viewPosInv);


		if (showStars) {
			if (starfield !== undefined && starfield.sceneObj !== undefined) {
				gl.disable(gl.DEPTH_TEST);
				starfield.sceneObj.uniforms.visibleMagnitudeBias = starsVisibleMagnitudeBias;
				starfield.sceneObj.draw();
				gl.enable(gl.DEPTH_TEST);
			}
		}

		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;
			if (planet.pos === undefined) continue;	//comets don't have pos yet, but I'm working on that
			if (orbitTarget.pos === undefined) continue;
		
			/*
			if (planet.isComet) continue;
			if (planet.isAsteroid) continue;
			if (planet.parent !== solarSystem.planets[solarSystem.indexes.Sun] &&
				planet !== solarSystem.planets[solarSystem.indexes.Sun] &&
				planet !== solarSystem.planets[solarSystem.indexes.Moon]) continue;
			*/

			if (showLinesToOtherPlanets && orbitStarSystem == solarSystem) {
				if (orbitTarget !== planet) {
					//while here, update lines

					vec3.sub(delta, planet.pos, orbitTarget.pos);

					var dist = vec3.length(delta);
					vec3.scale(delta, delta, 1/dist);

					lineObj.attrs.vertex.data[0] = delta[0] * orbitTarget.radius;
					lineObj.attrs.vertex.data[1] = delta[1] * orbitTarget.radius;
					lineObj.attrs.vertex.data[2] = delta[2] * orbitTarget.radius;
					lineObj.attrs.vertex.data[3] = delta[0] * orbitDistance;
					lineObj.attrs.vertex.data[4] = delta[1] * orbitDistance;
					lineObj.attrs.vertex.data[5] = delta[2] * orbitDistance;

					lineObj.attrs.vertex.updateData();
					lineObj.draw({uniforms : { color : planet.color }});
				}
			}

			if (showVelocityVectors) {
				vec3.sub(delta, planet.pos, orbitTarget.pos);
				lineObj.attrs.vertex.data[0] = delta[0];
				lineObj.attrs.vertex.data[1] = delta[1];
				lineObj.attrs.vertex.data[2] = delta[2];
				lineObj.attrs.vertex.data[3] = delta[0] + planet.vel[0] * velocityVectorScale;
				lineObj.attrs.vertex.data[4] = delta[1] + planet.vel[1] * velocityVectorScale;
				lineObj.attrs.vertex.data[5] = delta[2] + planet.vel[2] * velocityVectorScale;
				lineObj.attrs.vertex.updateData();
				lineObj.draw({uniforms : { color : planet.color }});
			}

			if (showRotationAxis) {
				vec3.sub(delta, planet.pos, orbitTarget.pos);
				var axis = [0,0,1];
				vec3.quatZAxis(axis, planet.angle);
				lineObj.attrs.vertex.data[0] = delta[0] + axis[0] * 2 * planet.radius;
				lineObj.attrs.vertex.data[1] = delta[1] + axis[1] * 2 * planet.radius;
				lineObj.attrs.vertex.data[2] = delta[2] + axis[2] * 2 * planet.radius;
				lineObj.attrs.vertex.data[3] = delta[0] + axis[0] * -2 * planet.radius;
				lineObj.attrs.vertex.data[4] = delta[1] + axis[1] * -2 * planet.radius;
				lineObj.attrs.vertex.data[5] = delta[2] + axis[2] * -2 * planet.radius;
				lineObj.attrs.vertex.updateData();
				lineObj.draw({uniforms : { color : planet.color }});
			}

			if (showOrbitAxis && planet.orbitAxis) {
				vec3.sub(delta, planet.pos, orbitTarget.pos);
				lineObj.attrs.vertex.data[0] = delta[0] + planet.orbitAxis[0] * 2 * planet.radius;
				lineObj.attrs.vertex.data[1] = delta[1] + planet.orbitAxis[1] * 2 * planet.radius;
				lineObj.attrs.vertex.data[2] = delta[2] + planet.orbitAxis[2] * 2 * planet.radius;
				lineObj.attrs.vertex.data[3] = delta[0] + planet.orbitAxis[0] * -2 * planet.radius;
				lineObj.attrs.vertex.data[4] = delta[1] + planet.orbitAxis[1] * -2 * planet.radius;
				lineObj.attrs.vertex.data[5] = delta[2] + planet.orbitAxis[2] * -2 * planet.radius;
				lineObj.attrs.vertex.updateData();
				lineObj.draw({uniforms : { color : planet.color }});
			}
		}

		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			//update vis ratio
			var dx = planet.pos[0] - glutil.view.pos[0] - orbitTarget.pos[0];
			var dy = planet.pos[1] - glutil.view.pos[1] - orbitTarget.pos[1];
			var dz = planet.pos[2] - glutil.view.pos[2] - orbitTarget.pos[2];
			planet.visRatio = planet.radius / Math.sqrt(dx * dx + dy * dy + dz * dz);
		}

		//draw sphere planets
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			if (planet.sceneObj && (planet.visRatio >= planetPointVisRatio)) {
				updatePlanetClassSceneObj(planet);

				//update scene object
				//don't forget one is shared among all planets
				vec3.sub(planet.sceneObj.uniforms.pos, planet.pos, orbitTarget.pos);
				quat.copy(planet.sceneObj.uniforms.angle, planet.angle);

				//only allow ring obj if there's a radii and a sphere
				//... this way i can just copy over the pos and angle
				if (planet.ringObj !== undefined) {
					vec3.sub(planet.ringObj.pos, planet.pos, orbitTarget.pos);
					quat.copy(planet.ringObj.angle, planet.sceneObj.uniforms.angle);
				}

				//calculate sun position for lighting
				for (var starIndex = 0; starIndex < orbitStarSystem.stars.length; ++starIndex) {
					var star = orbitStarSystem.stars[starIndex];
					var tmp = [];
					vec3.sub(tmp, star.pos, planet.pos);
					planet.sceneObj.uniforms.sunDir[0+3*starIndex] = tmp[0];
					planet.sceneObj.uniforms.sunDir[1+3*starIndex] = tmp[1];
					planet.sceneObj.uniforms.sunDir[2+3*starIndex] = tmp[2];
				}

				//webkit bug
				planet.sceneObj.shader.use();
				gl.uniform3fv(
					gl.getUniformLocation(
						planet.sceneObj.shader.obj, 'sunDir[0]'
					), planet.sceneObj.uniforms.sunDir);

				//update ellipsoid parameters
				planet.sceneObj.uniforms.equatorialRadius = planet.equatorialRadius !== undefined ? planet.equatorialRadius : planet.radius;
				planet.sceneObj.uniforms.inverseFlattening = planet.inverseFlattening !== undefined ? planet.inverseFlattening : .5;
				if (planet.ringObj !== undefined) {
					planet.sceneObj.uniforms.ringMinRadius = planet.ringRadiusRange[0];
					planet.sceneObj.uniforms.ringMaxRadius = planet.ringRadiusRange[1];
				}
				planet.sceneObj.uniforms.color = planet.color;

if (!CALCULATE_TIDES_WITH_GPU) {
				planet.sceneObj.attrs.tide = planet.tideBuffer;
}

//something is overriding this between the update calc and here 
// making me need to re-assign it ...
planet.sceneObj.uniforms.forceMin = planet.forceMin;
planet.sceneObj.uniforms.forceMax = planet.forceMax;

				//TODO - ambient for each star
				// TODO even more - absolute magnitude, radiance, HDR, etc
				planet.sceneObj.uniforms.ambient = planet.type == 'star' ? 1 : .1;

				planet.sceneObj.draw();

				if (showLatAndLonLines) {
					vec3.copy(planet.latLonObj.pos, planet.sceneObj.uniforms.pos);
					quat.copy(planet.latLonObj.angle, planet.sceneObj.uniforms.angle);
					planet.latLonObj.uniforms.equatorialRadius = planet.equatorialRadius !== undefined ? planet.equatorialRadius : planet.radius;
					planet.latLonObj.uniforms.inverseFlattening = planet.inverseFlattening !== undefined ? planet.inverseFlattening : 1.;
					planet.latLonObj.draw();
				}
			}
		}

		//draw rings and transparent objects last
		//disable depth writing so orbits are drawn in front and behind them
		//do this last so (without depth writing) other planets in the background don't show up in front of the rings
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			if (planet.sceneObj && (planet.visRatio >= planetPointVisRatio)) {
				if (planet.ringObj !== undefined) {
					gl.disable(gl.CULL_FACE);
					gl.depthMask(false);
					//to provide to the ring shader:
					//1) vector from planet to the sun (for fwd-vs-back scattering)
					//2) radius and eccentricity of planet (for self-shadows)
					//3) axis of planet rotation / normal of ring plane (for unlit side) ... can be derived from the angle, which is copied into the ring object

					//TODO cache rotation axis -- for here and for showing axis (like we already do for orbit axis)
					var axis = [];
					var sunDir = [];
					vec3.sub(sunDir, planet.pos, orbitStarSystem.planets[orbitStarSystem.indexes.Sun].pos);
					vec3.normalize(sunDir, sunDir);
					vec3.quatZAxis(axis, planet.ringObj.angle);
					var axisDotSun = vec3.dot(axis, sunDir);
					var lookingAtLitSide = vec3.dot(viewfwd, axis) * axisDotSun;
					lookingAtLitSide = (lookingAtLitSide < 0 ? -1 : 1) * Math.pow(Math.abs(lookingAtLitSide), 1/4);	//preserve sign, so raise to odd power
					planet.ringObj.uniforms.lookingAtLitSide = Math.clamp(.5 + .5 * lookingAtLitSide, 0, 1);
					var viewDotSun = vec3.dot(viewfwd, sunDir);
					planet.ringObj.uniforms.backToFrontLitBlend = .5 - .5 * viewDotSun;	//clamp(sqrt(sqrt(dot( normalize(sun.pos - planet.pos), axis ) * dot( viewfwd, axis ))), 0., 1.) * -.5 + .5


					//have to recalculate these because the uniforms are shared, so they've been overwritten
					//...and the ring objs have to be drawn last for transparency reasons...

					//calculate sun position for lighting
					for (var starIndex = 0; starIndex < orbitStarSystem.stars.length; ++starIndex) {
						var star = orbitStarSystem.stars[starIndex];
						var tmp = [];
						vec3.sub(tmp, star.pos, planet.pos);
						planet.ringObj.uniforms.sunDir[0+3*starIndex] = tmp[0];
						planet.ringObj.uniforms.sunDir[1+3*starIndex] = tmp[1];
						planet.ringObj.uniforms.sunDir[2+3*starIndex] = tmp[2];
					}
					vec3.sub(planet.ringObj.uniforms.pos, planet.pos, orbitTarget.pos);
					quat.copy(planet.ringObj.uniforms.angle, planet.angle);

					//webkit bug
					planet.ringObj.shader.use();
					gl.uniform3fv(
						gl.getUniformLocation(
							planet.ringObj.shader.obj, 'sunDir[0]'
						), planet.ringObj.uniforms.sunDir);


					planet.ringObj.draw();

					gl.depthMask(true);
					gl.enable(gl.CULL_FACE);
				}
			}
		}

		//draw point planets
		// goes slow when comets are included
		//TODO make a buffer with a vertex per-planet rather than changing a single vertex
		for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
			var planet = orbitStarSystem.planets[planetIndex];
			if (planet.hide) continue;

			if (!planet.sceneObj || planet.visRatio < planetPointVisRatio) {
				if (showPlanetsAsDistantPoints) {
					vec3.sub(pointObj.attrs.vertex.data, planet.pos, orbitTarget.pos);
					pointObj.attrs.vertex.updateData();
					pointObj.draw({
						uniforms : {
							color : planet.color,
							pointSize : 4
						}
					});
				} else {	//if (showStars) {
					//draw as a pixel with apparent magnitude
					//use the starfield shader
					vec3.sub(pointObj.attrs.vertex.data, planet.pos, orbitTarget.pos);
					pointObj.attrs.vertex.data[3] = planet.absoluteMagnitude || 10;
					pointObj.attrs.vertex.updateData();
					pointObj.attrs.colorIndex.data[0] = planet.colorIndex || 0;
					pointObj.attrs.colorIndex.updateData();
					pointObj.draw({
						shader : colorIndexShader,
						//disable depth test too?
						blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
						uniforms : {
							pointSize : 1
						},
						texs : [colorIndexTex]
					});
				}
			}
		}

		if (mouseOverTarget !== undefined) {
			var planet = mouseOverTarget;
			if (planet !== undefined) {
				gl.disable(gl.DEPTH_TEST);
				pointObj.attrs.vertex.data[0] = planet.pos[0] - orbitTarget.pos[0];
				pointObj.attrs.vertex.data[1] = planet.pos[1] - orbitTarget.pos[1];
				pointObj.attrs.vertex.data[2] = planet.pos[2] - orbitTarget.pos[2];
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
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
				if (planet.hide) continue;

				if (planet.orbitPathObj) {

					var semiMajorAxis = planet.keplerianOrbitalElements.semiMajorAxis;
					var eccentricity = planet.keplerianOrbitalElements.eccentricity;
					var distPeriapsis = semiMajorAxis * (1 + eccentricity);	//largest distance from the parent planet

					//vector from view to parent planet
					var parentPlanet = planet.parent;
					vec3.sub(delta, parentPlanet.pos, orbitTarget.pos);
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
								vec3.sub(planet.orbitPathObj.pos, planet.parent.pos, orbitTarget.pos);
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
			for (var planetIndex = 0; planetIndex < orbitStarSystem.planets.length; ++planetIndex) {
				var planet = orbitStarSystem.planets[planetIndex];
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
					orbitBasis[i] = planet.orbitBasis[0][i];
					orbitBasis[4+i] = planet.orbitBasis[1][i];
					orbitBasis[8+i] = planet.orbitBasis[2][i];
				}
				//gravWellObjBasis = orbitBasis * gravityWellZScale * zOffsetByRadius * orbitBasis^-1 * planetPosBasis

				//calc this for non-sceneObj planets.  maybe I should store it as a member variable?
				var relPos = [
					planet.pos[0] - orbitTarget.pos[0],
					planet.pos[1] - orbitTarget.pos[1],
					planet.pos[2] - orbitTarget.pos[2]
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
					gravityWellTargetZScale = 1 / orbitTarget.gravityWellScalar;
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
	};
})();

//quad of tris
var quad = [[0,0],[0,1],[1,1],[1,1],[1,0],[0,0]];
var latitudeMin = -90;
var latitudeMax = 90;
var latitudeStep = 5;
var longitudeMin = -180;
var longitudeMax = 180;
var longitudeStep = 5;
var latitudeDivisions = Math.floor((latitudeMax-latitudeMin)/latitudeStep);
var longitudeDivisions = Math.floor((longitudeMax-longitudeMin)/longitudeStep);
var tideTexWidth = 128;
var tideTexHeight = 128;

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;

	//fix side panel heights
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('height', window.innerHeight);
	});

	//fix info panel height
	if (showBodyInfo) {
		var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
		$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);
	}

	glutil.resize();

	//TODO fix mouseline to work with this
	//glutil.view.fovY = Math.min(1, canvas.height / canvas.width) * 90;

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
	for (var starSystemIndex = 0; starSystemIndex < starSystems.length; ++starSystemIndex) {
		var starSystem = starSystems[starSystemIndex];
		for (var planetIndex = 0; planetIndex < starSystem.planets.length; ++planetIndex) {
			var planet = starSystem.planets[planetIndex];
			planet.lastMeasureCalcDate = undefined;
		}
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
//local debugging
if (astro) {
	$('#currentTimeText').text(astro.julianToCalendar(julianDate));
}
}

var primaryPlanetHorizonIDs = [10, 199, 299, 301, 399, 499, 599, 699, 799, 899, 999];

var ModifiedDepthShaderProgram;

$(document).ready(init1);

var slideDuration = 500;
var slideWidth = 300;
var currentOpenSidePanelID = undefined;
var showBodyInfo = false;
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
	$('#menu').animate(
		{
			left : slideWidth
		}, {
			duration : slideDuration,
			step : function(now, fx) {
				var degrees = now / slideWidth * 180;
				setCSSRotation($(this), degrees);
			}
		}
	);
	sidePanel.animate({left:0}, {duration:slideDuration});
	currentOpenSidePanelID = sidePanelID;
}

function hideSidePanel(sidePanelID, dontMoveOpenButton) {
	if (sidePanelID === currentOpenSidePanelID) currentOpenSidePanelID = undefined;
	var sidePanel = $('#'+sidePanelID);
	if (!dontMoveOpenButton) {
		$('#menu').animate(
			{
				left : 0
			},
			{
				duration : slideDuration,
				step : function(now, fx) {
					var degrees = (now / slideWidth) * 180;
					setCSSRotation($(this), degrees);
				}
			}
		);
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
		'displayOptionsSidePanel',
		'overlaySidePanel',
		'solarSystemSidePanel',
		'smallBodiesSidePanel',
		'starSystemsSidePanel',
		'controlsSidePanel'
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
		{buttonID:'mainButtonDisplayOptions', divID:'displayOptionsSidePanel'},
		{buttonID:'mainButtonOverlay', divID:'overlaySidePanel'},
		{buttonID:'mainButtonSolarSystem', divID:'solarSystemSidePanel'},
		{buttonID:'mainButtonSmallBodies', divID:'smallBodiesSidePanel'},
		{buttonID:'mainButtonStarSystems', divID:'starSystemsSidePanel'},
		{buttonID:'mainButtonControls', divID:'controlsSidePanel'},
	], function(i, info) {
		$('#'+info.buttonID).click(function() {
			showSidePanel(info.divID);
		});
	});

	$('#reset').click(function() {
		integrationPaused = true;
		integrateTimeStep = defaultIntegrateTimeStep;
		for (var i = 0; i < starSystems.length; ++i) {
			var starSystem = starSystems[i];
			for (var j = 0; j < 1/*starSystem.planets.length*/; ++j) {
				var planet = starSystem.planets[j];
				var initPlanet = starSystem.initPlanets[j];
				vec3.copy(planet.pos, initPlanet.pos);
				vec3.copy(planet.vel, initPlanet.vel);
				quat.copy(planet.angle, initPlanet.angle);
			}
		}
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
			//  but, true to the original design, I actually want more detail up front.  Maybe I'll map it to get more detail up front than it already has?
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
console.log('glMaxCubeMapTextureSize', glMaxCubeMapTextureSize);

	/*glutil.view.angle[0] = -0.4693271591372717;
	glutil.view.angle[1] = 0.7157221264895661;
	glutil.view.angle[2] = -0.4298784661116332;
	glutil.view.angle[3] = 0.28753844912098436;*/
	glutil.view.zNear = 1e+4;
	glutil.view.zFar = 1e+25;

	fbo = new glutil.Framebuffer();
	gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.obj);
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);

	//now that GL is initialized, and before we are referencing planets, 
	// build the solar sytsem

	solarSystem = new SolarSystem();
	starSystemForNames[solarSystem.name] = solarSystem;
	solarSystem.index = starSystems.length;
	starSystems.push(solarSystem);
	solarSystem.doneBuildingPlanets();


	// overlay side panel


	var overlaySidePanelContents = $('#overlaySidePanelContents');
	$('<span>', {text:'Overlay:'}).appendTo(overlaySidePanelContents);
	$('<br>').appendTo(overlaySidePanelContents);
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
			.appendTo(overlaySidePanelContents);
		if (thisDisplayMethod == displayMethod) radio.attr('checked', 'checked');
		$('<span>', {text:thisDisplayMethod}).appendTo(overlaySidePanelContents);
		$('<br>').appendTo(overlaySidePanelContents);
	});
	$('<br>').appendTo(overlaySidePanelContents);
	$('<span>', {text:'Influencing Planets:'}).appendTo(overlaySidePanelContents);
	$('<br>').appendTo(overlaySidePanelContents);

	//add radio buttons hierarchically ...
	var overlayControlsForPlanets = {};

	var HierarchicalCheckboxControl = makeClass({
		/*
		args:
			title		<- used to identify this checkbox
			change		<- callback upon changing the checkbox value
			clickTitle	<- callback upon clicking the title
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
				text : args.title,
				css : {
					textDecoration : 'underline',
					cursor : 'pointer',
				},
				click : function(e) {
					args.clickTitle.call(thiz);
				}
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

	//TODO add star systems
	$.each(solarSystem.planets, function(planetIndex, planet) {
		//if any other planet doesn't have recorded mass then skip it
		if (planet.mass === undefined) return;

		var parentPlanet = planet.parent;
		if (parentPlanet !== undefined) {
			if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
		}

		var controls = new HierarchicalCheckboxControl({
			title : planet.name,
			isChecked : true,
			change : function() {
				planetInfluences[this.args.planetIndex] = this.checkbox.is(':checked');
				invalidateForces();
			},
			clickTitle : function() {
				setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
			},
			planetIndex : planetIndex
		});

		//add to parent or to the side panel
		if (parentPlanet === undefined) {
			controls.div.appendTo(overlaySidePanelContents);
		} else {
			overlayControlsForPlanets[parentPlanet.index].addChild(controls);
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
		if (planetIndex !== solarSystem.indexes.Sun &&
			planetIndex !== solarSystem.indexes.Earth)
		{
			controls.childDiv.hide();
		}
	}

	$('<br>').appendTo(overlaySidePanelContents);


	// display options side panel


	$.each([
		'showLinesToOtherPlanets',
		'showVelocityVectors',
		'showRotationAxis',
		'showOrbitAxis',
		'showLatAndLonLines',
		'showGravityWell',
		'showOrbits',
		'showStars',
		'showPlanetsAsDistantPoints',
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
		'starsVisibleMagnitudeBias'
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


	var solarSystemSidePanel = $('#solarSystemSidePanel');

	//add radio buttons hierarchically ...
	var celestialBodiesControlsForPlanets = {};

	var cometParent;

	$.each(solarSystem.planets, function(planetIndex,planet) {

		var parentPlanet = planet.parent;
		if (parentPlanet !== undefined) {
			if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
		}

		var controls = new HierarchicalCheckboxControl({
			title : planet.name,
			isChecked : !planet.hide,
			change : function() {
				solarSystem.planets[this.args.planetIndex].hide = !this.checkbox.is(':checked');
			},
			clickTitle : function() {
				setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
			},
			planetIndex : planetIndex
		});

		if (parentPlanet === undefined) {
			controls.div.appendTo($('#celestialBodiesVisiblePlanets'));
		} else {
			celestialBodiesControlsForPlanets[parentPlanet.index].addChild(controls);
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
		if (planetIndex !== solarSystem.indexes.Sun &&
			planetIndex !== solarSystem.indexes.Earth)
		{
			controls.childDiv.hide();
		}
	}

	$('<br>').appendTo($('#celestialBodiesVisiblePlanets'));

	//these are added to the end of the result
	//they should get greyed upon new query (search, prev, next click)
	//they should be regenerated upon new content
	//they should be re-enabled upon error
	var nextButton = undefined;
	var prevButton = undefined;

	var searchResults = [];

	var searchLastID = 0;	//uid to inc so if we search twice, the first can be invalidated

	//dataSource = 'remote' for remote queries, 'local' for the currently-selected planets
	var processSearch = function(pageIndex, dataSource) {
		var button = $('#celestialBodiesSearch');
		var searchText = $('#celestialBodiesSearchText');
		var searchStr = searchText.val();
		button.prop('disabled', 1);
		searchText.val('searching...');
		searchText.prop('disabled', 1);

		if (prevButton) prevButton.prop('disabled', 1);
		if (nextButton) nextButton.prop('disabled', 1);

		var searchID = ++searchLastID;

		var processResults = function(results) {
			if (searchID < searchLastID-1) return;	//someone else has searched

			searchText.val(searchStr);
			searchText.prop('disabled', 0);
			button.prop('disabled', 0);

			$('#celestialBodiesSearchToggleAll').prop('checked', 0);	//disable too?

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

				var titleSpan;
				var checkbox = $('<input>', {
					type : 'checkbox',
					change : function() {

						if (!$(this).is(':checked')) {	//uncheck checkbox => remove planet

							//only remove if it's already there
							if (solarSystem.indexes[name] !== undefined) {
								var index = solarSystem.indexes[name];
								var planet = solarSystem.planets[index];
								solarSystem.initPlanets.splice(index, 1);
								solarSystem.planets.splice(index, 1);
								//now remap indexes
								for (var i = index; i < solarSystem.planets.length; ++i) {
									solarSystem.planets[index].index = i;
								}
								if (orbitTarget === planet) {
									setOrbitTarget(solarSystem.planets[solarSystem.indexes.Sun]);
								}
								//TODO destruct WebGL geometry?  or is it gc'd automatically?
								//now rebuild indexes
								solarSystem.indexes = {};
								for (var i = 0; i < solarSystem.planets.length; ++i) {
									solarSystem.indexes[solarSystem.planets[i].name] = i;
								}
							}

							titleSpan.css({textDecoration:'', cursor:''});

						} else {	//check checkbox => add planet

							//only add if it's not there
							if (solarSystem.indexes[name] === undefined) {

								//add the row to the bodies

								var index = solarSystem.planets.length;
								var planet = mergeInto(new Planet(), {
									name : name,
									isComet : row.bodyType == 'comet',
									isAsteroid : row.bodyType == 'numbered asteroid' || row.bodyType == 'unnumbered asteroid',
									sourceData : row,
									parent : solarSystem.planets[solarSystem.indexes.Sun],
									starSystem : solarSystem,
									index : index
								});

								solarSystem.planets.push(planet);

								//add copy to initPlanets for when we get reset() working
								solarSystem.initPlanets[index] = planet.clone();

								solarSystem.indexes[planet.name] = index;

								initPlanetColorSchRadiusAngle(planet);
								initPlanetSceneLatLonLineObjs(planet);
								initPlanetOrbitPathObj(planet, false);
							}

							titleSpan.css({textDecoration:'underline', cursor:'pointer'});
						}
					}
				})
					.prop('checked', solarSystem.indexes[name] !== undefined)
					.appendTo(rowDiv);

				titleSpan = $('<span>', {
					text : name,
					click : function() {
						var targetPlanet = solarSystem.planets[solarSystem.indexes[name]];
						if (targetPlanet !== undefined) setOrbitTarget(targetPlanet);
					}
				}).appendTo(rowDiv);
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
					processSearch(pageIndex+dir, dataSource);
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
		};

		if (dataSource == 'remote') {
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
				$('#celestialBodiesSearchWarning').after(warning);
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
			}).done(processResults);
		} else if (dataSource == 'local') {
			var rows = [];
			var searchingComets = $('#celestialBodiesSearchComets').prop('checked');
			var searchingNumbered = $('#celestialBodiesSearchNumbered').prop('checked');
			var searchingUnnumbered = $('#celestialBodiesSearchUnnumbered').prop('checked');
			for (var i = 0; i < solarSystem.planets.length; ++i) {
				var planet = solarSystem.planets[i];
				var row = planet.sourceData;
				if (row) {
					if ((row.bodyType == 'comet' && searchingComets) ||
						(row.bodyType == 'numbered asteroid' && searchingNumbered) ||
						(row.bodyType == 'unnumbered asteroid' && searchingUnnumbered))
					{
						rows.push(row);
					}
				}
			}
			
/*hack for debugging*/
rows = [
	   //test case for hyperbolic
		{
			"perihelionDistance" : 209232954020.56,
			"inclination" : 2.2519519080902,
			"timeOfPerihelionPassage" : 2456955.82823,
			"argumentOfPerihelion" : 0.042602963442406,
			"bodyType" : "comet",
			"orbitSolutionReference" : "JPL 101",
			"epoch" : 56691,
			"eccentricity" : 1.00074241,
			"idNumber" : "C",
			"name" : "2013 A1 (Siding Spring)",
			"longitudeOfAscendingNode" : 5.2530270564044
		},
		/*
target C1/2013 Siding Spring data: on 2014-11-15 10:30:00
position m: 	64031628815.774, -196629902235.24, 57392197050.865
velocity m/day: -1916173862.297, -455287182.34414, 2315832701.4279
		*/
		
		//test case for elliptical
		{
			"perihelionDistance" : 87661077532.81,
			"inclination" : 2.8320181936429,
			"timeOfPerihelionPassage" : 2446467.39532,
			"argumentOfPerihelion" : 1.9431185149437,
			"bodyType" : "comet",
			"orbitSolutionReference" : "JPL J863/77",
			"epoch" : 49400,
			"eccentricity" : 0.96714291,
			"idNumber" : "1P",
			"name" : "Halley",
			"longitudeOfAscendingNode" : 1.0196227452785
		}
];
/**/
			
			processResults({rows:rows, count:rows.length});
		} else {
			throw "got an unknown data source request " + dataSource;
		}
	};

	$('#celestialBodiesSearchText').keydown(function(e){
		if (e.keyCode == 13) {
			$('#celestialBodiesSearch').trigger('click');
		}
	});

	//change a check box, immediately update search results
	$.each([
		$('#celestialBodiesSearchComets'),
		$('#celestialBodiesSearchNumbered'),
		$('#celestialBodiesSearchUnnumbered'),
		$('#celestialBodiesSearchVisible')
	], function(_,checkbox) {
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
		processSearch(0,
			$('#celestialBodiesSearchVisible').prop('checked') ? 'local' : 'remote'
		);
	});
	$('#celestialBodiesSearch').trigger('click');	//fire one off


	// rest of the init


	$(window).resize(resize);
	resize();

	//julianDate = horizonsDynamicData.julianDate;	//get most recent recorded date from planets
	//TODO get current julian date from client
	// integrate forward by the small timestep between the two
	julianDate = calendarToJulian(new Date());


	initJulianDate = julianDate;
	refreshCurrentTimeText();

	solarSystem.copyPlanets(solarSystem.initPlanets, solarSystem.planets);

	init2();
}

//call this once we get back our star data
function initStarsControls() {
	//TODO maybe some pages or something for this, like the asteroid/smallbody search?
	for (var i = 0; i < starSystems.length; ++i) {
		(function(){
			var starSystem = starSystems[i];
			$('<div>', {
				css : {
					textDecoration : 'underline',
					cursor : 'pointer',
					paddingLeft : '10px'
				},
				click : function() {
					setOrbitTarget(starSystem);
				},
				text : starSystem.name
			}).appendTo($('#starSystemContents'));
		})();
	}
}

function init2() {

	var imgs = [];
	$.each(primaryPlanetHorizonIDs, function(_,horizonID) {
		var planet = solarSystem.planetForHorizonID[horizonID];
		imgs.push('textures/'+planet.name.toLowerCase()+'.png');
	});
	imgs.push('textures/jupiter-rings-color.png');
	imgs.push('textures/saturn-rings-color.png');
	imgs.push('textures/saturn-rings-transparency.png');
	imgs.push('textures/saturn-rings-back-scattered.png');
	imgs.push('textures/saturn-rings-forward-scattered.png');
	imgs.push('textures/saturn-rings-unlit-side.png');
	//pick the cubemap textures once we know our max cubemap texture size -- 512 on firefox =( but 1024 on chrome =)
	for (var i = 0; i < skyTexFilenamePrefixes.length; ++i) {
		imgs.push(skyTexFilenamePrefixes[i] + Math.min(1024, glMaxCubeMapTextureSize) + '.png');
	}
	console.log('loading '+imgs.join(', '));
	$(imgs).preload(function(){
		$('#loadingDiv').hide();
	}, function(percent){
		$('#loading').attr('value', parseInt(100*percent));
	});

	$('#infoPanel').bind('touchmove', function(e) {
		if (!showBodyInfo) {
			if (!e.target.classList.contains('scrollable')) {
				//only touch devices will scroll the div even when browsers know not to
				//so for touch scroll, explicitly tell it not to scroll
				e.preventDefault();
			}
		}
	});

	setCSSRotation($('#toggleBodyInfo'), -90);
	//$('#infoPanel').css('height', $('#infoDiv').offset().top - $('#infoPanel').offset().top);
	//$('#infoPanel').css('bottom', -15);
	$('#toggleBodyInfo').click(function() {
		if (!showBodyInfo) {
			showBodyInfo = true;
			var infoDivTop = $('#infoPanel').offset().top;
			var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
			$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);
			console.log('fixing body info top at', infoDivTop);
			console.log('interpolating to dest top', infoDivDestTop);

			$('#infoPanel')
				//go from bottom-aligned to top-algined
				.css('bottom', 'auto')
				.css('top', infoDivTop)
				.css('overflow', 'scroll')
				//interpolate up to the time controls ... or some fixed distance to it
				.stop(true)
				.animate(
					{
						top : infoDivDestTop+'px'
					},
					{
						duration : slideDuration,
						step : function(now, fx) {
							var frac = 1 - (now - infoDivDestTop) / (infoDivTop - infoDivDestTop);
							var degrees = 180 * frac - 90;
							setCSSRotation($('#toggleBodyInfo'), degrees);
						}
					}
				);
		} else {
			showBodyInfo = false;
			$('#infoPanel').css('height', '104px');
			var infoDivBottom = window.innerHeight - ($('#infoPanel').offset().top + $('#infoPanel').height());
			var infoDivDestBottom = -15;
			console.log('fixing body info bottom at', infoDivBottom);
			console.log('interpolating to dest bottom ', infoDivDestBottom);

			$('#infoPanel')
				//go from bottom-aligned to top-algined
				.css('top', 'auto')
				.css('bottom', infoDivBottom)
				.css('overflow', 'visible')
				//interpolate up to the time controls ... or some fixed distance to it
				.stop(true)
				.animate(
					{
						bottom : infoDivDestBottom+'px'
					},
					{
						duration : slideDuration,
						step : function(now, fx) {
							var frac = (now - infoDivDestBottom) / (infoDivBottom - infoDivDestBottom);
							var degrees = 180 * frac - 90;
							setCSSRotation($('#toggleBodyInfo'), degrees);
						}
					}
				);

		}
	});

	$('#menu').show();
	$('#timeDiv').show();
	//$('#infoPanel').show();	//don't show here -- instead show upon first setOrbitTarget, when the css positioning gets fixed

	init3();
}

function initStars() {
	var xhr = new XMLHttpRequest();
	xhr.open('GET', 'hyg/stardata.f32', true);
	xhr.responseType = 'arraybuffer';
	/* if we want a progress bar ...
	xhr.onprogress = function(e) {
		if (e.total) {
			progress.attr('value', parseInt(e.loaded / e.total * 100));
		}
	};
	*/
	var numElem = 5;
	xhr.onload = function(e) {
		var arrayBuffer = this.response;
		var data = new DataView(arrayBuffer);

		var floatBuffer = new Float32Array(data.byteLength / Float32Array.BYTES_PER_ELEMENT);
		var len = floatBuffer.length;
		for (var j = 0; j < len; ++j) {
			var x = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);

			if (j % numElem < 3) {
				//convert xyz from parsec coordinates to meters ... max float is 10^38, so let's hope (/warn) if an incoming value is close to 10^22
				if (Math.abs(x) > 1e+20) {
					console.log('star '+Math.floor(j/numElem)+' has coordinate that position exceeds fp resolution');
				}
				x *= 3.08567758e+16;
			}

			floatBuffer[j] = x;
		}

		//now that we have the float buffer ...
		StarField.prototype.buffer = new glutil.ArrayBuffer({data : floatBuffer, dim : numElem});
		StarField.prototype.sceneObj = new glutil.SceneObject({
			mode : gl.POINTS,
			shader : colorIndexShader,
			texs : [colorIndexTex],
			attrs : {
				vertex : new glutil.Attribute({buffer : StarField.prototype.buffer, size : 4, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 0}),	//xyz abs-mag
				colorIndex : new glutil.Attribute({buffer : StarField.prototype.buffer, size : 1, stride : numElem * Float32Array.BYTES_PER_ELEMENT, offset : 4 * Float32Array.BYTES_PER_ELEMENT})	//color-index
			},
			uniforms : {
				pointSize : 1,
				color : [1,1,1,1],
				visibleMagnitudeBias : starsVisibleMagnitudeBias,
				colorIndexTex : 0
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			parent : null
		});

		//assign after all prototype buffer stuff is written, so StarField can call Star can use it during ctor
		starfield = new StarField();

		//now that we've built all our star system data ... add it to the star field
		if (starSystems.length > 1) addStarSystemsToStarField();

	};
	xhr.send();
}

function initExoplanets() {
	var processResults = function(results) {
		//process results
		$.each(results.systems, function(i,systemInfo) {
			var systemName = assertExists(systemInfo, 'name');

			var rightAscension = systemInfo.rightAscension;
			if (rightAscension === undefined) {
				console.log('failed to find right ascension for system '+systemName);
				return;
			}
			var declination = systemInfo.declination;
			if (declination === undefined) {
				console.log('failed to find declination for system '+systemName);
				return;
			}
			var cosRA = Math.cos(rightAscension);
			var sinRA = Math.sin(rightAscension);
			var cosDec = Math.cos(declination);
			var sinDec = Math.sin(declination);
			//convert to coordinates
			var pos = [];
			pos[0] = cosRA * cosDec;
			pos[1] = sinRA * cosDec;
			pos[2] = sinDec;
			//rotate for earth's tilt
			var epsilon = Math.rad(23 + 1/60*(26 + 1/60*(21.4119)));
			var cosEps = Math.cos(epsilon);
			var sinEps = Math.sin(epsilon);
			var yn = cosEps * pos[1] + sinEps * pos[2];
			pos[2] = -sinEps * pos[1] + cosEps * pos[2];
			pos[1] = yn;
			//distance
			var distance = systemInfo.distance;
			if (distance === undefined) {
				console.log('failed to find distance for system '+systemName);
				return;
			}
			pos[0] *= distance;
			pos[1] *= distance;
			pos[2] *= distance;


			var starSystem = new StarSystem();
			starSystem.name = systemName;
			starSystem.sourceData = systemInfo;
			vec3.copy(starSystem.pos, pos);

			//TODO visual magnitude, which is undocumented but some have


			$.each(assertExists(systemInfo, 'bodies'), function(j,bodyInfo) {

				var name = assertExists(bodyInfo, 'name');
				var radius = bodyInfo.radius;
				if (radius === undefined) {
					if (bodyInfo.type !== 'barycenter') {
						//console.log('no radius for body '+name);
						//if planets don't have radii they can be estimated by the planet's mass or distance from sun or both?
						radius = 7e+7;	// use a jupiter radius
					}
				}

				var mass = bodyInfo.mass;
				if (mass === undefined) {
					if (bodyInfo.density && bodyInfo.radius) {
						mass = 4/3 * Math.PI * bodyInfo.density * radius * radius * radius;
					}
				}
				if (mass === undefined) {
					//this prevents us from an orbit ..
					//console.log('no mass for body '+name);
					mass = 2e+27;	//use a jupiter mass
				}

				var body = mergeInto(new Planet(), {
					type : assertExists(bodyInfo, 'type'),
					name : name,
					mass : mass,
					radius : radius,	//for planets/suns this is the radius of the planet.  does it exist for barycenters?
					//visualMagnitude : bodyInfo.visualMagnitude,	// TODO convert to absolute magnitude for star shader?
					//spectralType : bodyInfo.spectralType,			// or this one?
					//temperature : bodyInfo.temperature,			// or this one?
					sourceData : bodyInfo,
					parent : bodyInfo.parent,
					starSystem : starSystem,
					hide : bodyInfo.type === 'barycenter',
					isExoplanet : true
				});

				//hacks for orbit
				bodyInfo.longitudeOfAscendingNode = bodyInfo.longitudeOfAscendingNode || 0;
				bodyInfo.inclination = bodyInfo.inclination || 0;
				bodyInfo.eccentricity = bodyInfo.eccentricity || 0;

				starSystem.planets.push(body);
				if (body.type == 'star') starSystem.stars.push(body);
			});

			starSystem.doneBuildingPlanets();

			for (var i = 0; i < starSystem.planets.length; ++i) {
				var body = starSystem.planets[i];

				//further hacks for orbit, now that the parent pointer has been established
				if (body.sourceData.semiMajorAxis === undefined) {
					if (body.parent &&
						body.parent.type === 'barycenter' &&
						//this is what I don't like about the open exoplanet database:
						//sometimes 'separation' or 'semimajoraxis' is the distance to the parent
						// sometimes it's the distance to the children
						(body.parent.sourceData.separation !== undefined ||
						body.parent.sourceData.semiMajorAxis !== undefined))
					{
						body.sourceData.semiMajorAxis = body.parent.sourceData.semiMajorAxis ||
							body.parent.sourceData.separation || 0;
						//TODO also remove semiMajorAxis from barycenter so it doesn't get offset from the planet's center?
						// also -- do this all in exoplanet file preprocessing?
					}
				}
				//longitude of periapsis = longitude of ascending node + argument of periapsis
				//argument of periapsis = longitude of periapsis - longitude of ascending node
				if (body.sourceData.longitudeOfPeriapsis === undefined) {
					body.sourceData.longitudeOfPeriapsis = 0;
				}
				body.sourceData.argumentOfPerihelion = body.sourceData.longitudeOfPeriapsis - body.sourceData.longitudeOfAscendingNode;
				if (body.sourceData.argumentOfPerihelion == 0) {
					body.sourceData.argumentOfPerihelion = Math.PI * 2 * i / starSystem.planets.length;
				}

				vec3.add(body.pos, body.pos, starSystem.pos);

				initPlanetColorSchRadiusAngle(body);
				initPlanetSceneLatLonLineObjs(body);

				if (body.pos[0] !== body.pos[0]) {
					console.log('system '+starSystem.name+' planet '+body.name+' has bad pos');
				}

				initPlanetOrbitPathObj(body, false);
			}

			starSystem.initPlanets = starSystem.clonePlanets();
			starSystemForNames[starSystem.name] = starSystem;
			starSystem.index = starSystems.length;
			starSystems.push(starSystem);
		});

		//now that we've built all our star system data ... add it to the star field
		if (starfield !== undefined) addStarSystemsToStarField();
	};

	//just over a meg, so I might as well ajax it
	$.ajax({
		url : 'exoplanet/openExoplanetCatalog.json',
		dataType : 'json',
		timeout : 5000
	}).error(function() {
		console.log('failed to get exoplanets, trying again...');
		setTimeout(function() {
			initExoplanets();
		}, 5000);
	}).done(processResults);
}

function addStarSystemsToStarField() {
console.log('adding star systems to star fields and vice versa');	
	assert(starfield !== undefined);
	assert(starSystems.length > 1);

	//add buffer points

	var array = starfield.buffer.data;	//5 fields: x y z absmag colorIndex
	array = Array.prototype.slice.call(array);	//to js array

	for (var i = 0; i < starSystems.length; ++i) {
		var starSystem = starSystems[i];
		starSystem.starfieldIndex = array.length / 5;
		array.push(starSystem.pos[0]);
		array.push(starSystem.pos[1]);
		array.push(starSystem.pos[2]);
		array.push(5);	//abs mag
		array.push(0);	//color index
	}

	starfield.buffer.setData(array, gl.STATIC_DRAW, true);

	//then add named stars to the starSystem array
	// if they're not there ... and I'm pretty sure they're all not there ...

	for (var i = 0; i < namedStars.length; ++i) {
		var args = namedStars[i];

		if (({
			Sun:1,
			"Kapteyn's Star":1,
			"Rigel Kentaurus A":1,
			"Rigel Kentaurus B":1
		})[args.name]) continue;	//hmm... there are some double-ups ...
		
		var starSystem = new StarSystem();

		starSystem.name = args.name;

		//index in the stars.buffer, stride of stars.buffer.dim
		//preserved since the original 1000 or however many from the HYG database are not moved
		starSystem.starfieldIndex = args.index;

		//note stars.buffer holds x y z abs-mag color-index
		starSystem.pos = [
			starfield.buffer.data[0 + starSystem.starfieldIndex * starfield.buffer.dim],
			starfield.buffer.data[1 + starSystem.starfieldIndex * starfield.buffer.dim],
			starfield.buffer.data[2 + starSystem.starfieldIndex * starfield.buffer.dim]
		];

		//starSystem.sourceData = ...

		//now add a single star to the new star system ...

		var body = mergeInto(new Planet(), {
			type : 'star',
			name : args.name,
			//mass : undefined,
			//radius : undefined,
			//sourceData : {},
			//parent : undefined,
			starSytsem : starSystem,
			hide : false,
			isExoplanet : true
		});
		vec3.add(body.pos, body.pos, starSystem.pos);
		initPlanetColorSchRadiusAngle(body);
		initPlanetSceneLatLonLineObjs(body);
		initPlanetOrbitPathObj(body, false);
		
		starSystem.planets.push(body);
		starSystem.stars.push(body);
		
		starSystem.doneBuildingPlanets();
	
		starSystem.initPlanets = starSystem.clonePlanets();
		starSystemForNames[starSystem.name] = starSystem;
		starSystem.index = starSystems.length;
		starSystems.push(starSystem);

	}

	initStarsControls();
}

function initPlanetColorSchRadiusAngle(planet) {
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
		planet.color = [Math.random(), Math.random(), Math.random(), 1];
		vec3.normalize(planet.color, planet.color);
	} else {
		planet.color = [color[0], color[1], color[2], 1];
	}
	planet.schwarzschildRadius = 2 * planet.mass * kilogramsPerMeter;
}

//TODO no more multiple copies of tide array per-planet
function initPlanetSceneLatLonLineObjs(planet) {
	if (planet.radius === undefined) {
		//only/always use a point/basis/etc?
	} else {
		var tideArray = [];
		tideArray.length = planetSceneObj.attrs.vertex.count;
		for (var i = 0; i < tideArray.length; ++i) tideArray[i] = 0;

		planet.sceneObj = planetSceneObj;
		planet.latLonObj = planetLatLonObj;

if (!CALCULATE_TIDES_WITH_GPU) {
		//old way, per-vertex storage, updated by CPU
		planet.tideBuffer = new glutil.ArrayBuffer({dim : 1, data : tideArray, usage : gl.DYNAMIC_DRAW});
} else {
		//new way, per-texel storage, updated by GPU FBO kernel
		//TODO pull these by request rather than allocating them per-planet (since we only ever see one or two or maybe 10 or 20 at a time .. never all 180+ local and even more)
		planet.tideTex = new glutil.Texture2D({
			internalFormat : gl.RGBA,
			type : gl.FLOAT,
			width : tideTexWidth,
			height : tideTexHeight,
			magFilter : gl.LINEAR,
			minFilter : gl.NEAREST,
			wrap : {
				s : gl.CLAMP_TO_EDGE,
				t : gl.CLAMP_TO_EDGE
			}
		});
}
	}
}


//GLSL code generator
function unravelForLoop(varname, start, end, code) {
	var lines = [];
	for (var i = start; i <= end; ++i) {
		lines.push('#define '+varname+' '+i);
		lines.push(code);
		lines.push('#undef '+varname);
	}
	return lines.join('\n')+'\n';
};

var geodeticPositionCode = mlstr(function(){/*

//takes in a vec2 of lat/lon
//returns the ellipsoid coordinates in the planet's frame of reference
//I could move the uniforms here, but then the function would be imposing on the shader
#define M_PI 3.141592653589793115997963468544185161590576171875
vec3 geodeticPosition(vec2 latLon) {
	float phi = latLon.x * M_PI / 180.;
	float lambda = latLon.y * M_PI / 180.;
	float cosPhi = cos(phi);
	float sinPhi = sin(phi);
	float eccentricitySquared = (2. * inverseFlattening - 1.) / (inverseFlattening * inverseFlattening);
	float sinPhiSquared = sinPhi * sinPhi;
	float N = equatorialRadius / sqrt(1. - eccentricitySquared * sinPhiSquared);
	const float height = 0.;
	float NPlusH = N + height;	//plus height, if you want?  no one is using height at the moment.  heightmaps someday...
	return vec3(
		NPlusH * cosPhi * cos(lambda),
		NPlusH * cosPhi * sin(lambda),
		(N * (1. - eccentricitySquared) + height) * sinPhi);
}
*/});

//request this per solar system.  rebuild if we need, return from cache if we don't.
var planetShadersPerNumberOfStarsCache = {};
function getPlanetShadersForNumberOfStars(numberOfStars) {
	if (numberOfStars <= 0) numberOfStars = 1;	//huh, I guess I have a star system with no stars ... "CFBDSIR2149 / CFBDSIR J214947.2-040308.9 / CFBDS J214947-040308"
	var shaders = planetShadersPerNumberOfStarsCache[numberOfStars];
	if (shaders !== undefined) return shaders;

	shaders = {};
	shaders.colorShader = new ModifiedDepthShaderProgram({
		vertexCode :
'#define NUM_STARS '+numberOfStars+'\n'+
		mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs:
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
//to fragment shader:
varying vec3 lightDir[NUM_STARS];		//light position
varying vec3 normal;		//surface normal

*/}) + geodeticPositionCode + mlstr(function(){/*

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
*/}) + unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'+
		mlstr(function(){/*
uniform vec4 color;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
uniform float ambient;
void main() {
	float litLum = 0.;
*/}) + unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
+ mlstr(function(){/*
	float luminance = min(1., litLum);
	gl_FragColor = color * max(ambient, sqrt(luminance));
}
*/}),
		uniforms : {
			color : [1,1,1,1],
			pointSize : 4
		}
	});

	shaders.texShader = new ModifiedDepthShaderProgram({
		vertexCode :
'#define NUM_STARS '+numberOfStars+'\n'+
		mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
//to fragment shader:
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;

*/}) + geodeticPositionCode + mlstr(function(){/*

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex);
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
*/}) + unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'+
		mlstr(function(){/*
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
uniform sampler2D tex;
uniform float ambient;
void main() {
	float litLum = 0.;
*/}) + unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
+ mlstr(function(){/*
	float luminance = min(1., litLum);
	gl_FragColor = texture2D(tex, texCoordv) * max(ambient, sqrt(luminance));
}
*/}),
		uniforms : {
			tex : 0
		}
	});

	shaders.ringShadowShader = new ModifiedDepthShaderProgram({
		vertexCode :
'#define NUM_STARS '+numberOfStars+'\n'+
		mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
//to fragment shader:
varying vec3 modelVertexv;
varying vec3 normal;
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
*/}) + geodeticPositionCode + mlstr(function(){/*

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

void main() {
	//vertex is really the lat/lon in degrees
	modelVertexv = geodeticPosition(vertex);
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 worldVertex = quatRotate(angle, modelVertexv) + pos;
	normal = quatRotate(angle, normalize(modelVertexv));
*/}) + unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - worldVertex);')
+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(worldVertex, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'+
		mlstr(function(){/*
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
varying vec3 modelVertexv;
uniform sampler2D tex;
uniform sampler2D ringTransparencyTex;
uniform float ambient;
uniform float ringMinRadius;
uniform float ringMaxRadius;
uniform vec4 angle;

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

float ringIntersect(vec3 startPos, vec3 dir) {
	if (dot(startPos, dir) < 0.) return -1.;	//occluded by planet
	float t = -startPos.z / dir.z;
	if (t < 0.) return -1.;	//trace intersects backwards
	vec2 intersect = startPos.xy + t * dir.xy;
	float r = length(intersect);
	return (r - ringMinRadius) / (ringMaxRadius - ringMinRadius);
}

void main() {
	float luminance = 1.;

	//inverse rotate lightDir[0]
	//TODO how to combine this per-light-source ... depenends on the intensity of each ...
	vec3 lightDirInModelSpace = quatRotate(vec4(-angle.xyz, angle.w), lightDir[0]);
	float intersectPos = ringIntersect(modelVertexv, lightDirInModelSpace);
	if (intersectPos >= 0. && intersectPos <= 1.) {
		luminance *= texture2D(ringTransparencyTex, vec2(intersectPos, .5)).r;
	}

	float litLum = 0.;
*/}) + unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
+ mlstr(function(){/*
	luminance *= min(1., litLum);

	gl_FragColor.rgb = texture2D(tex, texCoordv).rgb * max(ambient, sqrt(luminance));
	gl_FragColor.a = 1.;
}
*/}),
		uniforms : {
			tex : 0,
			ringTransparencyTex : 1
		}
	});

	planetShadersPerNumberOfStarsCache[numberOfStars] = shaders;
	return shaders;
}

function init3() {
	hsvTex = new glutil.HSVTexture(256);

	colorShader = new ModifiedDepthShaderProgram({
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

	latLonShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;	//lat/lon pairs
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;
uniform float inverseFlattening;

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

		*/}) + geodeticPositionCode + mlstr(function(){/*

void main() {
	vec3 modelVertex = geodeticPosition(vertex);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
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
			color : [1,1,1,1]
		}
	});

if (!CALCULATE_TIDES_WITH_GPU) {
	
	//renders a heat map from the float values of the 'tide' attribute
	planetHeatMapAttrShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
attribute float tide;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
varying float tidev;
varying vec2 texCoordv;

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

*/}) + geodeticPositionCode + mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex);
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	tidev = tide;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
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

} else { //if (CALCULATE_TIDES_WITH_GPU)

	//renders a heat map from the float values of the 'tide' texture 
	planetHeatMapTexShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
varying vec2 texCoordv;

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

*/}) + geodeticPositionCode + mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex);
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
precision highp sampler2D;

varying vec2 texCoordv;
uniform sampler2D tex;
uniform sampler2D tideTex;
uniform sampler2D hsvTex;
uniform float heatAlpha;
uniform float forceMin, forceMax;	//for clamping range
void main() {
	float gradientValue = (texture2D(tideTex, texCoordv).r - forceMin) / (forceMax - forceMin);
	vec4 hsvColor = texture2D(hsvTex, vec2(gradientValue, .5));
	vec4 planetColor = texture2D(tex, texCoordv);
	gl_FragColor = mix(planetColor, hsvColor, heatAlpha);
}
*/}),
		uniforms : {
			tex : 0,
			tideTex : 1,
			hsvTex : 2
		}
	});



	//this is the FBO kernel update shader
	//currently tidal and gravitational calculations are done by cpu and mapped to hsv tex
	// next step would be to calculate them on GPU according to world position
	// (possibly need double precision in GLSL ... )
	// ( use https://thasler.com/blog/?p=93 for emulated double precision in GLSL code)
	planetSurfaceCalculationShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec2 vertex;
attribute vec2 texCoord;
varying vec2 texCoordv;
void main() {
	texCoordv = texCoord;
	gl_Position = vec4(vertex, 0., 1.);
}
		*/}),
		fragmentCode : mlstr(function(){/*
precision highp sampler2D;

varying vec2 texCoordv;

uniform vec3 pos;					//planet position, relative to the sun
uniform vec4 angle;					//planet angle

uniform float equatorialRadius;
uniform float inverseFlattening;

uniform int sourcePlanetIndex;		//index of the source planet.

uniform int planetStateTexHeight;
uniform sampler2D planetStateTex;	//texture of planet states.  currently [x y z mass]

// x = tangent tidal
// y = normal tidal
// z = tangent gravitational
// w = normal gravitational
uniform bvec4 flags; 

*/}) 
+ 'const float gravitationalConstant = '+gravitationalConstant+';	// m^3 / (kg * s^2)\n'
+ geodeticPositionCode + mlstr(function(){/*

vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}


//while the double precision hack works great for fractions, it doesn't add too much to the extent of the range of the number
//1e+19 is still an upper bound for our floating point numbers
//soo ... scale everything down here, and up later
const float precisionBias = 1e-4;


//https://thasler.com/blog/?p=93

struct double1 {
	vec2 v;	//v.x == hi, v.y == lo
};

double1 double1_set(float a) {
	double1 d;
	d.v = vec2(a, 0.);
	return d;
}

double1 double1_set2(float hi, float lo) {
	double1 d;
	d.v = vec2(hi, lo);
	return d;
}

double1 double1_add(double1 dsa, double1 dsb) {
	double1 dsc;
	float t1 = dsa.v.x + dsb.v.x;
	float e = t1 - dsa.v.x;
	float t2 = ((dsb.v.x - e) + (dsa.v.x - (t1 - e))) + dsa.v.y + dsb.v.y;
	dsc.v.x = t1 + t2;
	dsc.v.y = t2 - (dsc.v.x - t1);
	return dsc;
}

double1 double1_sub(double1 dsa, double1 dsb) {
	double1 dsc;
	float t1 = dsa.v.x - dsb.v.x;
	float e = t1 - dsa.v.x;
	float t2 = ((-dsb.v.x - e) + (dsa.v.x - (t1 - e))) + dsa.v.y - dsb.v.y;
	dsc.v.x = t1 + t2;
	dsc.v.y = t2 - (dsc.v.x - t1);
	return dsc;
}

double1 double1_mul(double1 dsa, double1 dsb) {
	double1 dsc;
	const float split = 8193.;

	float cona = dsa.v.x * split;
	float conb = dsb.v.x * split;
	float a1 = cona - (cona - dsa.v.x);
	float b1 = conb - (conb - dsb.v.x);
	float a2 = dsa.v.x - a1;
	float b2 = dsb.v.x - b1;

	float c11 = dsa.v.x * dsb.v.x;
	float c21 = a2 * b2 + (a2 * b1 + (a1 * b2 + (a1 * b1 - c11)));

	float c2 = dsa.v.x * dsb.v.y + dsa.v.y * dsb.v.x;

	float t1 = c11 + c2;
	float e = t1 - c11;
	float t2 = dsa.v.y * dsb.v.y + ((c2 - e) + (c11 - (t1 - e))) + c21;

	dsc.v.x = t1 + t2;
	dsc.v.y = t2 - (dsc.v.x - t1);

	return dsc;
}

double1 double1_div(double1 dsa, double1 dsb) {
	const float split = 8193.;
	float s1 = dsa.v.x / dsb.v.x;
	float cona = s1 * split;
	float conb = dsb.v.x * split;
	float a1 = cona - (cona - s1);
	float b1 = conb - (conb - dsb.v.x);
	float a2 = s1 - a1;
	float b2 = dsb.v.x - b1;
	float c11 = s1 * dsb.v.x;
	float c21 = (((a1 * b1 - c11) + a1 * b2) + a2 * b1) + a2 * b2;
	float c2 = s1 * dsb.v.y;
	float t1 = c11 + c2;
	float e = t1 - c11;
	float t2 = ((c2 - e) + (c11 - (t1 - e))) + c21;
	float t12 = t1 + t2;
	float t22 = t2 - (t12 - t1);
	float t11 = dsa.v.x - t12;
	e = t11 - dsa.v.x;
	float t21 = ((-t12 - e) + (dsa.v.x - (t11 - e))) + dsa.v.y - t22;
	float s2 = (t11 + t21) / dsb.v.x;
	float dsc_hi = s1 + s2;
	return double1_set2(dsc_hi, s2 - (dsc_hi - s1));
}

double1 double1_sqrt(double1 a) {
	if (a.v.x == 0.) return double1_set(0.);
	float t1 = 1. / sqrt(a.v.x);
	float t2 = a.v.x * t1;
	double1 s0 = double1_mul(double1_set(t2), double1_set(t2));
	double1 s1 = double1_sub(a, s0);
	float t3 = 0.5 * s1.v.x * t1;
	s0 = double1_set(t2);
	s1 = double1_set(t3);
	return double1_add(s0, s1);
}

//same as double1 but trying to take advantage of parallelization
//upper bound: 1.3e+19

struct double3 {
	vec3 hi, lo;
};

double3 double3_zero() {
	double3 d;
	d.hi = d.lo = vec3(0., 0., 0.);
	return d;
}

double3 double3_set(vec3 v) {
	double3 d;
	d.hi = v;
	d.lo = vec3(0., 0., 0.);
	return d;
}

double3 double3_set2(vec3 hi, vec3 lo) {
	double3 d;
	d.hi = hi;
	d.lo = lo;
	return d;
}

double3 double3_add(double3 dsa, double3 dsb) {
	double3 dsc;
	vec3 t1 = dsa.hi + dsb.hi;
	vec3 e = t1 - dsa.hi;
	vec3 t2 = ((dsb.hi - e) + (dsa.hi - (t1 - e))) + dsa.lo + dsb.lo;
	dsc.hi = t1 + t2;
	dsc.lo = t2 - (dsc.hi - t1);
	return dsc;
}

double3 double3_sub(double3 dsa, double3 dsb) {
	double3 dsc;
	vec3 t1 = dsa.hi - dsb.hi;
	vec3 e = t1 - dsa.hi;
	vec3 t2 = ((-dsb.hi - e) + (dsa.hi - (t1 - e))) + dsa.lo - dsb.lo;
	dsc.hi = t1 + t2;
	dsc.lo = t2 - (dsc.hi - t1);
	return dsc;
}

double3 double3_mul(double3 dsa, double3 dsb) {
	double3 dsc;
	const vec3 split = vec3(8193.);

	vec3 cona = dsa.hi * split;
	vec3 conb = dsb.hi * split;
	vec3 a1 = cona - (cona - dsa.hi);
	vec3 b1 = conb - (conb - dsb.hi);
	vec3 a2 = dsa.hi - a1;
	vec3 b2 = dsb.hi - b1;

	vec3 c11 = dsa.hi * dsb.hi;
	vec3 c21 = a2 * b2 + (a2 * b1 + (a1 * b2 + (a1 * b1 - c11)));

	vec3 c2 = dsa.hi * dsb.lo + dsa.lo * dsb.hi;

	vec3 t1 = c11 + c2;
	vec3 e = t1 - c11;
	vec3 t2 = dsa.lo * dsb.lo + ((c2 - e) + (c11 - (t1 - e))) + c21;

	dsc.hi = t1 + t2;
	dsc.lo = t2 - (dsc.hi - t1);

	return dsc;
}

double3 double3_scale(double3 dsa, double1 dsb) {
	double1 x = double1_mul(double1_set2(dsa.hi.x, dsa.lo.x), dsb);
	double1 y = double1_mul(double1_set2(dsa.hi.y, dsa.lo.y), dsb);
	double1 z = double1_mul(double1_set2(dsa.hi.z, dsa.lo.z), dsb);
	return double3_set2(vec3(x.v.x, y.v.x, z.v.x), vec3(x.v.y, y.v.y, z.v.y));
}

double1 double3_dot(double3 dsa, double3 dsb) {
	double3 dsc = double3_mul(dsa, dsb);
	return double1_add(
		double1_set2(dsc.hi.x, dsc.lo.x),
		double1_add(
			double1_set2(dsc.hi.y, dsc.lo.y),
			double1_set2(dsc.hi.z, dsc.lo.z)
		)
	);
}

double1 double3_length(double3 dsa) {
	return double1_sqrt(double3_dot(dsa, dsa));
}

double3 calcGravityAccel(double3 solarSystemVertex, double3 planetPos, double1 planetMass) {
	double3 x = double3_sub(solarSystemVertex, planetPos);
	double1 r = double3_length(x);
	double1 r2 = double1_mul(r, r);	//zero -- too big
	double1 r3 = double1_mul(r, r2);
	return double3_scale(
			x, 
			double1_div(
				double1_mul(
					double1_set(-gravitationalConstant * precisionBias), 
					planetMass
				), 
				r3
			)
		);
}

double3 calcTidalAccel(double3 solarSystemVertex, double3 solarSystemNormal, double3 planetPos, double1 planetMass) {
	double3 x = double3_sub(solarSystemVertex, planetPos);
	double1 r = double3_length(x);
	double1 r2 = double1_mul(r, r);
	double1 r3 = double1_mul(r, r2);
	double1 xDotN = double3_dot(x, solarSystemNormal);
	return double3_scale(
			double3_sub(
				double3_scale(
					x, 
					double1_div(
						double1_mul(
							double1_set(3.), 
							xDotN
						), 
						r2
					)
				), 
				solarSystemNormal
			),
			double1_div(
				double1_mul(
					double1_set(gravitationalConstant * precisionBias),
					planetMass
				),
				r3
			)
		);
}

float calcForceValue(vec3 solarSystemVertex_s, vec3 solarSystemNormal_s) {
	//here's where we have to cycle through every planet and perform our calculation on each of them
	//then sum the results ...
	// I can guarantee already that the # of planets will exceed the GLSL max uniforms ...
	//that means we'll have to upload the positions and masses via textures ...
	//this would fit well with my plans to eventually do the keplerian orbital element pos calcs on the GPU

	double3 solarSystemVertex = double3_set(solarSystemVertex_s * precisionBias);
	double3 solarSystemNormal = double3_set(solarSystemNormal_s);

	double3 tideAccel = double3_zero();
	double3 gravAccel = double3_zero();
	
	float delta_i = 1. / float(planetStateTexHeight);
	//only works with constant upper bounds ...
	const int NEEDLESSUPPERBOUND = 1024;
	for (int i = 0; i < NEEDLESSUPPERBOUND; ++i) {
		//...so I just have to put another upper bound here ...
		if (i >= planetStateTexHeight) break;
		vec4 planetState = texture2D(planetStateTex, vec2(.5, (float(i)+.5)*delta_i));
		if (planetState.w == 0.) continue;

		double3 planetPos = double3_set(planetState.xyz * precisionBias);
		double1 planetMass = double1_set(planetState.w * precisionBias);
		if ((flags[0] || flags[1]) && (i != sourcePlanetIndex)) {
			tideAccel = double3_add(tideAccel, calcTidalAccel(solarSystemVertex, solarSystemNormal, planetPos, planetMass));
		}
		if (flags[2] || flags[3]) {
			gravAccel = double3_add(gravAccel, calcGravityAccel(solarSystemVertex, planetPos, planetMass));
		}
	}

	double3 tideNormal = double3_scale(solarSystemNormal, double3_dot(tideAccel, solarSystemNormal));
	double3 tideTangent = double3_sub(tideAccel, tideNormal);
	double3 gravNormal = double3_scale(solarSystemNormal, double3_dot(gravAccel, solarSystemNormal));
	double3 gravTangent = double3_sub(gravAccel, gravNormal);

	double3 accel = double3_zero();
	if (flags[0]) accel = double3_add(accel, tideNormal);
	if (flags[1]) accel = double3_add(accel, tideTangent);
	if (flags[2]) accel = double3_add(accel, gravNormal);
	if (flags[3]) accel = double3_add(accel, gravTangent);
	return double3_length(accel).v.x;
}

void main() {
	vec2 latLonCoord = ((texCoordv - vec2(.5, .5)) * vec2(360., 180.)).yx;
	vec3 planetVertex = geodeticPosition(latLonCoord);	//planet's local frame
	vec3 solarSystemRotatedVertex = quatRotate(angle, planetVertex);	//...rotated into the solar system (but still centered at earth origin)
	vec3 solarSystemVertex = solarSystemRotatedVertex + pos;	///...translated to be relative to the SSB 
	vec3 solarSystemNormal = normalize(solarSystemRotatedVertex);	//normalize() doesn't correctly handle ellipses.  you're supposed to scale the normal by [1/sx, 1/sy, 1/sz] or something to correct for the flattening distortion of the normals.
	float forceValue = calcForceValue(solarSystemVertex, solarSystemNormal);
	gl_FragColor = vec4(forceValue);
}
*/}),
		uniforms : {
			planetStateTex : 0
		}
	});
}



	//tex used for performing min/max across tide tex
	tideReduceTexs = [];
	$.each([0,1],function() {
		tideReduceTexs.push(new glutil.Texture2D({
			width : tideTexWidth,
			height : tideTexHeight,
			internalFormat : gl.RGBA,
			format : gl.RGBA,
			type : gl.FLOAT,
			minFilter : gl.NEAREST,
			magFilter : gl.NEAREST,
			wrap : {
				s : gl.CLAMP_TO_EDGE,
				t : gl.CLAMP_TO_EDGE
			}//,
			//data : initialDataF32	//why is this needed?
		}));
	});

	encodeTempTex = new glutil.Texture2D({
		internalFormat : gl.RGBA,
		format : gl.RGBA,
		type : gl.UNSIGNED_BYTE,
		width : tideTexWidth,
		height : tideTexHeight,
		minFilter : gl.NEAREST,
		magFilter : gl.NEAREST,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		}
	});

	minReduceShader = new glutil.KernelShader({
		code : mlstr(function(){/*
void main() {
vec2 intPos = pos * viewsize - .5;

float a = texture2D(srcTex, (intPos * 2. + .5) / texsize).x;
float b = texture2D(srcTex, (intPos * 2. + vec2(1., 0.) + .5) / texsize).x;
float c = texture2D(srcTex, (intPos * 2. + vec2(0., 1.) + .5) / texsize).x;
float d = texture2D(srcTex, (intPos * 2. + vec2(1., 1.) + .5) / texsize).x;
float e = min(a,b);
float f = min(c,d);
float g = min(e,f);
gl_FragColor = vec4(g, 0., 0., 0.);
}
*/}),
		uniforms : {
			texsize : 'vec2',
			viewsize : 'vec2'
		},
		texs : ['srcTex']
	});

	//should I just double this up in the one reduciton pass?
	//it would mean twice as many operations anyways
	maxReduceShader = new glutil.KernelShader({
		code : mlstr(function(){/*
void main() {
vec2 intPos = pos * viewsize - .5;

float a = texture2D(srcTex, (intPos * 2. + .5) / texsize).x;
float b = texture2D(srcTex, (intPos * 2. + vec2(1., 0.) + .5) / texsize).x;
float c = texture2D(srcTex, (intPos * 2. + vec2(0., 1.) + .5) / texsize).x;
float d = texture2D(srcTex, (intPos * 2. + vec2(1., 1.) + .5) / texsize).x;
float e = max(a,b);
float f = max(c,d);
float g = max(e,f);
gl_FragColor = vec4(g, 0., 0., 0.);
}
*/}),
		uniforms : {
			texsize : 'vec2',
			viewsize : 'vec2'
		},
		texs : ['srcTex']
	});



	//http://lab.concord.org/experiments/webgl-gpgpu/webgl.html
	encodeShader = [];
	for (var channel = 0; channel < 4; ++channel) {
		encodeShader[channel] = new glutil.KernelShader({
			code : mlstr(function(){/*
float shift_right(float v, float amt) {
	v = floor(v) + 0.5;
	return floor(v / exp2(amt));
}

float shift_left(float v, float amt) {
	return floor(v * exp2(amt) + 0.5);
}

float mask_last(float v, float bits) {
	return mod(v, shift_left(1.0, bits));
}

float extract_bits(float num, float from, float to) {
	from = floor(from + 0.5);
	to = floor(to + 0.5);
	return mask_last(shift_right(num, from), to - from);
}

vec4 encode_float(float val) {
	if (val == 0.0)
		return vec4(0, 0, 0, 0);
	float sign = val > 0.0 ? 0.0 : 1.0;
	val = abs(val);
	float exponent = floor(log2(val));
	float biased_exponent = exponent + 127.0;
	float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0;
	
	float t = biased_exponent / 2.0;
	float last_bit_of_biased_exponent = fract(t) * 2.0;
	float remaining_bits_of_biased_exponent = floor(t);
	
	float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0;
	float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0;
	float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0;
	float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0;
	return vec4(byte4, byte3, byte2, byte1);
}

void main() {
	vec4 data = texture2D(tex, pos);
	gl_FragColor = encode_float(data[$channel]);
}
*/}).replace(/\$channel/g, channel),
			texs : ['tex']
		});
	}

	//going by http://stackoverflow.com/questions/21977786/star-b-v-color-index-to-apparent-rgb-color
	//though this will be helpful too: http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color.html
	var colorIndexTexWidth = 1024;
	colorIndexTex = new glutil.Texture2D({
		width : colorIndexTexWidth,
		height : 1,
		internalFormat : gl.RGBA,
		format : gl.RGBA,	//TODO function callback support for format == gl.RGB
		type : gl.UNSIGNED_BYTE,
		magFilter : gl.LINEAR,
		minFilter : gl.NEAREST,
		wrap : gl.CLAMP_TO_EDGE,
		data : function(texX, texY) {
//console.log('texX',texX);
			var frac = (texX + .5) / colorIndexTexWidth;
//console.log('frac',frac);
			var bv = frac * (colorIndexMax - colorIndexMin) + colorIndexMin;
//console.log('bv',bv);

//from here on out is the OP of the above article
if (false) {

			var t = 4600 * ((1 / ((.92 * bv) + 1.7)) + (1 / ((.92 * bv) + .62)));
//console.log('t',t);

			var x, y = 0;

			if (t>=1667 && t<=4000) {
				x = ((-0.2661239 * Math.pow(10,9)) / Math.pow(t,3)) + ((-0.2343580 * Math.pow(10,6)) / Math.pow(t,2)) + ((0.8776956 * Math.pow(10,3)) / t) + 0.179910;
			} else if (t > 4000 && t <= 25000) {
				x = ((-3.0258469 * Math.pow(10,9)) / Math.pow(t,3)) + ((2.1070379 * Math.pow(10,6)) / Math.pow(t,2)) + ((0.2226347 * Math.pow(10,3)) / t) + 0.240390;
			}
//console.log('x',x);

			if (t >= 1667 && t <= 2222) {
				y = -1.1063814 * Math.pow(x,3) - 1.34811020 * Math.pow(x,2) + 2.18555832 * x - 0.20219683;
			} else if (t > 2222 && t <= 4000) {
				y = -0.9549476 * Math.pow(x,3) - 1.37418593 * Math.pow(x,2) + 2.09137015 * x - 0.16748867;
			} else if (t > 4000 && t <= 25000) {
				y = 3.0817580 * Math.pow(x,3) - 5.87338670 * Math.pow(x,2) + 3.75112997 * x - 0.37001483;
			}
//console.log('y',y);

			//the rest is found at http://en.wikipedia.org/wiki/SRGB

			// xyY to XYZ, Y = 1
			var Y = (y == 0)? 0 : 1;
			var X = (y == 0)? 0 : (x * Y) / y;
			var Z = (y == 0)? 0 : ((1 - x - y) * Y) / y;
//console.log('X',X);
//console.log('Y',Y);
//console.log('Z',Z);

			var R_linear = 3.2406 * X - 1.5372 * Y - .4986 * Z;
			var G_linear = -.9689 * X + 1.8758 * Y + .0415 * Z;
			var B_linear = .0557 * X - .2040 * Y + 1.0570 * Z;
//console.log('R_linear',R_linear);
//console.log('G_linear',G_linear);
//console.log('B_linear',B_linear);

			var srgbGammaAdjust = function(C_linear) {
				var a = .055;
				if (C_linear <= .0031308) return 12.92 * C_linear;
				return (1 + a) * Math.pow(C_linear, 1/0.5) - a;
			};

			var R_srgb = srgbGammaAdjust(R_linear);
			var G_srgb = srgbGammaAdjust(G_linear);
			var B_srgb = srgbGammaAdjust(B_linear);
//console.log('R_srgb',R_srgb);
//console.log('G_srgb',G_srgb);
//console.log('B_srgb',B_srgb);

			var result = [
				Math.clamp(R_srgb, 0, 1),
				Math.clamp(G_srgb, 0, 1),
				Math.clamp(B_srgb, 0, 1),
				1];

			return result;

}	//...and this is the reply ...
if (true) {

			var t;  r=0.0; g=0.0; b=0.0; if (t<-0.4) t=-0.4; if (t> 2.0) t= 2.0;
			if ((bv>=-0.40)&&(bv<0.00)) { t=(bv+0.40)/(0.00+0.40); r=0.61+(0.11*t)+(0.1*t*t); }
			else if ((bv>= 0.00)&&(bv<0.40)) { t=(bv-0.00)/(0.40-0.00); r=0.83+(0.17*t)		  ; }
			else if ((bv>= 0.40)&&(bv<2.10)) { t=(bv-0.40)/(2.10-0.40); r=1.00				   ; }
			if ((bv>=-0.40)&&(bv<0.00)) { t=(bv+0.40)/(0.00+0.40); g=0.70+(0.07*t)+(0.1*t*t); }
			else if ((bv>= 0.00)&&(bv<0.40)) { t=(bv-0.00)/(0.40-0.00); g=0.87+(0.11*t)		  ; }
			else if ((bv>= 0.40)&&(bv<1.60)) { t=(bv-0.40)/(1.60-0.40); g=0.98-(0.16*t)		  ; }
			else if ((bv>= 1.60)&&(bv<2.00)) { t=(bv-1.60)/(2.00-1.60); g=0.82		 -(0.5*t*t); }
			if ((bv>=-0.40)&&(bv<0.40)) { t=(bv+0.40)/(0.40+0.40); b=1.00				   ; }
			else if ((bv>= 0.40)&&(bv<1.50)) { t=(bv-0.40)/(1.50-0.40); b=1.00-(0.47*t)+(0.1*t*t); }
			else if ((bv>= 1.50)&&(bv<1.94)) { t=(bv-1.50)/(1.94-1.50); b=0.63		 -(0.6*t*t); }

			return [r,g,b,1];
}
		}
	});

	colorIndexShader = new ModifiedDepthShaderProgram({
		vertexCode :
'#define COLOR_INDEX_MIN '+toGLSLFloat(colorIndexMin)+'\n'+
'#define COLOR_INDEX_MAX '+toGLSLFloat(colorIndexMax)+'\n'+
		mlstr(function(){/*
attribute vec4 vertex;
attribute float colorIndex;
varying float alpha;
varying vec3 color;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float visibleMagnitudeBias;
uniform sampler2D colorIndexTex;
#define M_LOG_10			2.3025850929940459010936137929093092679977416992188
#define PARSECS_PER_M		3.2407792910106957712004544136882149907718250416875e-17
void main() {
	vec4 vtx4 = mvMat * vec4(vertex.xyz, 1.);

	//calculate apparent magnitude, convert to alpha
	float distanceInParsecs = length(vtx4.xyz) * PARSECS_PER_M;	//in parsecs
	float absoluteMagnitude = vertex.w;
	float apparentMagnitude = absoluteMagnitude - 5. * (1. - log(distanceInParsecs) / M_LOG_10);

	//something to note: these power functions give us a sudden falloff,
	// but providing pre-scaled vertices and scaling them afterwards gives us an even more sudden falloff
	// one last fix would be to never scale the vertices and render different scaled objects at separate coordinate systems
	alpha = pow(100., -.2*(apparentMagnitude - visibleMagnitudeBias));

	//calculate color
	color = texture2D(colorIndexTex, vec2((colorIndex - COLOR_INDEX_MIN) / (COLOR_INDEX_MAX - COLOR_INDEX_MIN), .5)).rgb;

	//TODO point sprite / point spread function?
	gl_PointSize = 1.;

	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
	fragmentCode : mlstr(function(){/*
varying vec3 color;
varying float alpha;
void main() {
	gl_FragColor = vec4(color, alpha);
}
*/})
	});


	$('#overlaySlider').slider({
		range : 'max',
		width : '200px',
		min : 0,
		max : 100,
		value : 100 * heatAlpha,
		slide : function(event, ui) {
			heatAlpha = ui.value / 100;
if (!CALCULATE_TIDES_WITH_GPU) {
			gl.useProgram(planetHeatMapAttrShader.obj);
			gl.uniform1f(planetHeatMapAttrShader.uniforms.heatAlpha.loc, heatAlpha);
} else {
			gl.useProgram(planetHeatMapTexShader.obj);
			gl.uniform1f(planetHeatMapTexShader.uniforms.heatAlpha.loc, heatAlpha);
}
			gl.useProgram(null);
		}
	});
if (!CALCULATE_TIDES_WITH_GPU) {
	gl.useProgram(planetHeatMapAttrShader.obj);
	gl.uniform1f(planetHeatMapAttrShader.uniforms.heatAlpha.loc, heatAlpha);
} else {
	gl.useProgram(planetHeatMapTexShader.obj);
	gl.uniform1f(planetHeatMapTexShader.uniforms.heatAlpha.loc, heatAlpha);
}
	gl.useProgram(null);

	pointObj = new glutil.SceneObject({
		mode : gl.POINTS,
		shader : colorShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 4,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			}),
			colorIndex : new glutil.ArrayBuffer({
				dim : 1,
				count : 1,
				usage : gl.DYNAMIC_DRAW
			})
		},
		uniforms : {
			pointSize : 4
		},
		parent : null
	});

	//init our planet shaders

	(function(){
		var triIndexArray = [];
		var latLonIndexArray = [];
		var vertexArray = [];

		for (var loni=0; loni <= longitudeDivisions; ++loni) {
			var lon = longitudeMin + loni * longitudeStep;
			for (var lati=0; lati <= latitudeDivisions; ++lati) {
				var lat = latitudeMin + lati * latitudeStep;

				vertexArray.push(lat);
				vertexArray.push(lon);

				if (loni < longitudeDivisions && lati < latitudeDivisions) {
					for (var j = 0; j < quad.length; ++j) {
						var ofs = quad[j];
						var index = (lati + ofs[0]) + (latitudeDivisions + 1) * (loni + ofs[1]);
						triIndexArray.push(index);
					}
					//if we're using 5 div step then every 6 will be 30 degrees
					if ((loni * longitudeStep) % 30 == 0) {
						latLonIndexArray.push(lati + (latitudeDivisions+1) * loni);
						latLonIndexArray.push(lati+1 + (latitudeDivisions+1) * loni);
					}
					if ((lati * latitudeStep) % 30 == 0) {
						latLonIndexArray.push(lati + (latitudeDivisions+1) * loni);
						latLonIndexArray.push(lati + (latitudeDivisions+1) * (loni+1));
					}
				}
			}
		}

		//these could be merged if you don't mind evaluating the cos() and sin() in-shader
		var vertexBuffer = new glutil.ArrayBuffer({dim : 2, data : vertexArray});

		planetSceneObj = new glutil.SceneObject({
			mode : gl.TRIANGLES,
			indexes : new glutil.ElementArrayBuffer({data : triIndexArray}),
			attrs : {
				vertex : vertexBuffer
			},
			uniforms : {
				color : [1,1,1,1],
				pos : [0,0,0],
				angle : [0,0,0,1],
				sunDir : [0,0,0, 0,0,0, 0,0,0, 0,0,0],
				ambient : .3
			},
			texs : [],
			parent : null,
			static : true
		});

		planetLatLonObj = new glutil.SceneObject({
			mode : gl.LINES,
			indexes : new glutil.ElementArrayBuffer({
				data : latLonIndexArray
			}),
			shader : latLonShader,
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
	})();

	//init stars now that shaders are made
	initStars();

	lineObj = new glutil.SceneObject({
		mode : gl.LINES,
		shader : colorShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				count : 2,
				usage : gl.DYNAMIC_DRAW
			})
		},
		uniforms : {
			color : [1,1,1,1]
		},
		parent : null,
		static : true
	});

	quadObj = new glutil.SceneObject({
		mode : gl.TRIANGLE_STRIP,
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 2,
				data : [-1,-1, 1,-1, -1,1, 1,1]
			}),
			texCoord : new glutil.ArrayBuffer({
				dim : 2,
				data : [0,0, 1,0, 0,1, 1,1]
			})
		},
		parent : null,
		static : true
	});	
	var planetsDone = 0;

	for (var planetIndex_ = 0; planetIndex_ < solarSystem.planets.length; ++planetIndex_) { (function(){
		var planetIndex = planetIndex_;
		var planet = solarSystem.planets[planetIndex];
		initPlanetColorSchRadiusAngle(planet);

		var checkDone = function() {

			initPlanetSceneLatLonLineObjs(planet);

			++planetsDone;
			if (planetsDone == solarSystem.planets.length) {
				initOrbitPaths();
			}
		};

		// load texture
		if (primaryPlanetHorizonIDs.indexOf(planet.id) !== -1) {
			var img = new Image();
			img.onload = function() {
				planet.tex = new glutil.Texture2D({
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

	//ring texture for Jupiter
	//http://www.celestiamotherlode.net/catalog/jupiter.php
	(function(){
		var jupiterRingShader = new ModifiedDepthShaderProgram({
			//vertex code matches Saturn
			vertexCode : mlstr(function(){/*
#define M_PI 3.1415926535897931
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float ringMinRadius;
uniform float ringMaxRadius;

uniform vec3 pos;
uniform vec4 angle;

varying vec2 texCoordv;
varying vec3 worldPosv;	//varying in world coordinates

vec3 quatRotate( vec4 q, vec3 v ){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

void main() {
	texCoordv = vertex;

	//vertex is the uv of the texcoord wrapped around the annulus
	//u coord is radial, v coord is angular
	float rad = mix(ringMinRadius, ringMaxRadius, vertex.x);
	float theta = 2. * M_PI * vertex.y;
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);

	//rotate the planet-local-frame vector into modelview space
	vec4 modelPos = vec4(cosTheta * rad, sinTheta * rad, 0., 1.);

	//rotate the planet-local-frame vector into modelview space
	vec4 eyePos = mvMat * modelPos;
	worldPosv = quatRotate(angle, modelPos.xyz) + pos;

	gl_Position = projMat * eyePos;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
#define NUM_STARS 1
varying vec2 texCoordv;
varying vec3 worldPosv;

uniform sampler2D colorTex;

//vector from planet to sun
uniform vec3 sunDir[NUM_STARS];
uniform vec3 pos;

uniform float planetRadius;

float sphereIntersect(vec3 startPos, vec3 dir, vec3 spherePos, float sphereRadius) {
	vec3 sphereToStart = startPos - spherePos;
	float a = dot(dir, dir);
	float b = 2. * dot(sphereToStart, dir);
	float c = dot(sphereToStart, sphereToStart) - sphereRadius * sphereRadius;
	float discr = b*b - 4.*a*c;
	if (discr < 0.) return 1.;	//full trace, no intersection
	float sqrtDiscr = sqrt(discr);
	float invDenom = .5 / a;
	float t1 = (-b - sqrtDiscr) * invDenom;
	float t2 = (-b + sqrtDiscr) * invDenom;
	if (t1 >= 0. && t1 <= 1.) return t1;
	if (t2 >= 0. && t2 <= 1.) return t2;
	return 1.;
}

void main() {
	float luminance = 1.;
	//notice that I have to scale down the meters here for shader accuracy to work
	luminance = min(luminance, step(1., sphereIntersect(worldPosv * 1e-8, sunDir[0] * 1e-8, pos * 1e-8, planetRadius * 1e-8)));

	gl_FragColor = texture2D(colorTex, texCoordv);
	gl_FragColor.rgb *= sqrt(luminance);
}
*/})
		});

		var planet = solarSystem.planets[solarSystem.indexes.Jupiter];

		var texSrcInfo = [
			{field:'ringColorTex', url:'textures/jupiter-rings-color.png', format:gl.RGBA},
		];
		var numLoaded = 0;
		var onLoadRingTex = function(url) {
			++numLoaded;
			if (numLoaded < texSrcInfo.length) return;
			if (numLoaded > texSrcInfo.length) throw "already created the rings!";

			//done! create the ring object
			var vertexes = [];
			for (var i = 0; i < ringResolution; ++i) {
				var f = i / (ringResolution - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}

			var texs = [];
			for (var i = 0; i < texSrcInfo.length; ++i) {
				texs[i] = planet[texSrcInfo[i].field];
			}

			//and a ring object
			planet.ringObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				shader : jupiterRingShader,
				attrs : {
					vertex : new glutil.ArrayBuffer({
						dim : 2,
						data : vertexes
					})
				},
				uniforms : {
					colorTex : 0,
					ringMinRadius : planet.ringRadiusRange[0],
					ringMaxRadius : planet.ringRadiusRange[1],
					planetRadius : planet.radius,
					sunDir : [0,0,0, 0,0,0, 0,0,0, 0,0,0],
					pos : [0,0,0],
					angle : [0,0,0,1]
				},
				blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
				texs : texs,
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null
			});
		};
		$.each(texSrcInfo, function(i,info) {
			planet[info.field] = new glutil.Texture2D({
				url : info.url,
				format : info.format,
				internalFormat : info.format,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				wrap : {
					s :  gl.CLAMP_TO_EDGE,
					t : gl.REPEAT
				},
				generateMipmap : true,
				onload : onLoadRingTex
			});
		});
	})();


	//Saturn's rings
	(function(){
		var saturnRingShader = new ModifiedDepthShaderProgram({
			vertexCode : mlstr(function(){/*
#define M_PI 3.1415926535897931
attribute vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float ringMinRadius;
uniform float ringMaxRadius;

//these are actually baked into the mvMat, but I need them separate for converting coordinates to world coordinates (rather than eye coordinates)
//pos is in the scene, which is offset by the oribitting planet
//techinically the offset doesn't matter since i'm using it to offset vertices from the planet, then sphere-intersect-testing rays from those vertices against the planet
uniform vec3 pos;
uniform vec4 angle;

varying vec2 texCoordv;
varying vec3 worldPosv;	//varying in world coordinates

vec3 quatRotate( vec4 q, vec3 v ){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

void main() {
	texCoordv = vertex;

	//vertex is the uv of the texcoord wrapped around the annulus
	//u coord is radial, v coord is angular
	float rad = mix(ringMinRadius, ringMaxRadius, vertex.x);
	float theta = 2. * M_PI * vertex.y;
	float cosTheta = cos(theta);
	float sinTheta = sin(theta);

	vec4 modelPos = vec4(cosTheta * rad, sinTheta * rad, 0., 1.);

	//rotate the planet-local-frame vector into modelview space
	vec4 eyePos = mvMat * modelPos;
	worldPosv = quatRotate(angle, modelPos.xyz) + pos;

	gl_Position = projMat * eyePos;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
#define NUM_STARS 1
//described here:
//http://www.mmedia.is/~bjj/data/s_rings/index.html

varying vec2 texCoordv;
varying vec3 worldPosv;

uniform sampler2D colorTex;
uniform sampler2D backScatteredTex;	//from the sun looking at the rings
uniform sampler2D forwardScatteredTex;	//looking at the rings towards the sun
uniform sampler2D transparencyTex;
uniform sampler2D unlitSideTex;

uniform float lookingAtLitSide;
uniform float backToFrontLitBlend;

//vector from planet to sun
uniform vec3 sunDir[NUM_STARS];
uniform vec3 pos;

uniform float planetRadius;

float sphereIntersect(vec3 startPos, vec3 dir, vec3 spherePos, float sphereRadius) {
	vec3 sphereToStart = startPos - spherePos;
	float a = dot(dir, dir);
	float b = 2. * dot(sphereToStart, dir);
	float c = dot(sphereToStart, sphereToStart) - sphereRadius * sphereRadius;
	float discr = b*b - 4.*a*c;
	if (discr < 0.) return 1.;	//full trace, no intersection
	float sqrtDiscr = sqrt(discr);
	float invDenom = .5 / a;
	float t1 = (-b - sqrtDiscr) * invDenom;
	float t2 = (-b + sqrtDiscr) * invDenom;
	if (t1 >= 0. && t1 <= 1.) return t1;
	if (t2 >= 0. && t2 <= 1.) return t2;
	return 1.;
}

void main() {
	vec3 color = texture2D(colorTex, texCoordv).rgb;
	gl_FragColor.rgb = color;

	float backScattered = texture2D(backScatteredTex, texCoordv).r;
	float forwardScattered = texture2D(forwardScatteredTex, texCoordv).r;
	float litSide = mix(backScattered, forwardScattered, backToFrontLitBlend);

	float unlitSide = texture2D(unlitSideTex, texCoordv).r;

	//shadow from the sun
	//raytrace from worldPosv along sunDir (non-normalized)
	//if it intersects with a sphere at pos with radius planetRadius then we're in shadow
	//if you want to do a soft shadow based on the sun's radius and distance (length of sunDir) then you can

	float luminance = 1.;
	//notice that I have to scale down the meters here for shader accuracy to work
	luminance = min(luminance, step(1., sphereIntersect(worldPosv * 1e-8, sunDir[0] * 1e-8, pos * 1e-8, planetRadius * 1e-8)));
	luminance *= mix(unlitSide, litSide, lookingAtLitSide);

	gl_FragColor.rgb *= sqrt(luminance);

	float transparency = texture2D(transparencyTex, texCoordv).r;
	gl_FragColor.a = transparency;
}
*/})
		});

		var planet = solarSystem.planets[solarSystem.indexes.Saturn];

		var texSrcInfo = [
			{field:'ringColorTex', url:'textures/saturn-rings-color.png', format:gl.RGBA},
			{field:'ringBackScatteredTex', url:'textures/saturn-rings-back-scattered.png', format:gl.LUMINANCE},
			{field:'ringForwardScatteredTex', url:'textures/saturn-rings-forward-scattered.png', format:gl.LUMINANCE},
			{field:'ringTransparencyTex', url:'textures/saturn-rings-transparency.png', format:gl.LUMINANCE},
			{field:'ringUnlitSideTex', url:'textures/saturn-rings-unlit-side.png', format:gl.LUMINANCE}
		];
		var numLoaded = 0;
		var onLoadRingTex = function(url) {
			++numLoaded;
			if (numLoaded < texSrcInfo.length) return;
			if (numLoaded > texSrcInfo.length) throw "already created the rings!";

			//done! create the ring object
			var vertexes = [];
			for (var i = 0; i < ringResolution; ++i) {
				var f = i / (ringResolution - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}

			var texs = [];
			for (var i = 0; i < texSrcInfo.length; ++i) {
				texs[i] = planet[texSrcInfo[i].field];
			}

			//and a ring object
			planet.ringObj = new glutil.SceneObject({
				mode : gl.TRIANGLE_STRIP,
				shader : saturnRingShader,
				attrs : {
					vertex : new glutil.ArrayBuffer({
						dim : 2,
						data : vertexes
					})
				},
				uniforms : {
					colorTex : 0,
					backScatteredTex : 1,
					forwardScatteredTex : 2,
					transparencyTex : 3,
					unlitSideTex : 4,
					ringMinRadius : planet.ringRadiusRange[0],
					ringMaxRadius : planet.ringRadiusRange[1],
					planetRadius : planet.radius,
					lookingAtLitSide : 1,
					backToFrontLitBlend : 0,
					sunDir : [0,0,0, 0,0,0, 0,0,0, 0,0,0],
					pos : [0,0,0],
					angle : [0,0,0,1]
				},
				blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
				texs : texs,
				pos : [0,0,0],
				angle : [0,0,0,1],
				parent : null
			});
		};
		$.each(texSrcInfo, function(i,info) {
			planet[info.field] = new glutil.Texture2D({
				url : info.url,
				format : info.format,
				internalFormat : info.format,
				minFilter : gl.LINEAR_MIPMAP_LINEAR,
				magFilter : gl.LINEAR,
				wrap : {
					s :  gl.CLAMP_TO_EDGE,
					t : gl.REPEAT
				},
				generateMipmap : true,
				onload : onLoadRingTex
			});
		});
	})();
}

function initPlanetOrbitPathObj(planet, useVectorState) {
	var planetIndex = planet.index;

	// based on planet position and velocity, find plane of orbit
	if (planet.parent === undefined) {
		//console.log(planet.name+' has no orbit parent');
		planet.orbitAxis = [0,0,1];
		planet.orbitBasis = [[1,0,0],[0,1,0],[0,0,1]];
		return;
	}

	var parentPlanet = planet.parent;
	if (parentPlanet === undefined) {
		console.log(planet.name+' has an invalid orbit planet');
		planet.orbitAxis = [0,0,1];
		planet.orbitBasis = [[1,0,0],[0,1,0],[0,0,1]];
		return;
	}

	//calculated variables:
	var eccentricity = undefined;
	var eccentricAnomaly = undefined;
	var orbitType = undefined;
	var meanAnomaly = undefined;
	var A, B;

	//used by planets to offset reconstructed orbit coordinates to exact position of planet
	var checkPosToPosX = 0;
	var checkPosToPosY = 0;
	var checkPosToPosZ = 0;

	//if we're not using vector state then we're using keplerian orbital elements
	if (!useVectorState) {
		/*
		we have:
			eccentricity
			semiMajorAxis
			pericenterDistance
		we compute:
		*/

		eccentricity = assertExists(planet.sourceData, 'eccentricity');

		var parabolicEccentricityEpsilon = 1e-7;
		if (Math.abs(eccentricity - 1) < parabolicEccentricityEpsilon) {
			orbitType = 'parabolic';
		} else if (eccentricity > 1) {
			orbitType = 'hyperbolic';
		} else {
			orbitType = 'elliptic';
		}

		var pericenterDistance;
		var semiMajorAxis;
		if (planet.isComet) {
			pericenterDistance = assert(planet.sourceData.perihelionDistance);
			if (orbitType !== 'parabolic') {
				semiMajorAxis = pericenterDistance / (1 - eccentricity);
			}	//otherwise, if it is parabolic, we don't get the semi-major axis ...
		} else if (planet.isAsteroid) {
			semiMajorAxis = assert(planet.sourceData.semiMajorAxis);
			pericenterDistance = semiMajorAxis * (1 - eccentricity);
		} else if (planet.isExoplanet) {
			semiMajorAxis = planet.sourceData.semiMajorAxis || 0;
		}

		var gravitationalParameter = gravitationalConstant * parentPlanet.mass;	//assuming the comet mass is negligible, since the comet mass is not provided
		var semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;
		
		//orbital period is only defined for circular and elliptical orbits (not parabolic or hyperbolic)
		var orbitalPeriod = undefined;
		if (orbitType === 'elliptic') {
			orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed / gravitationalParameter) / (60*60*24);	//julian day
		}

		var longitudeOfAscendingNode = assertExists(planet.sourceData, 'longitudeOfAscendingNode');
		var cosAscending = Math.cos(longitudeOfAscendingNode);
		var sinAscending = Math.sin(longitudeOfAscendingNode);

		var argumentOfPericenter = assertExists(planet.sourceData, 'argumentOfPerihelion');
		var cosPericenter = Math.cos(argumentOfPericenter);
		var sinPericenter = Math.sin(argumentOfPericenter);

		var inclination = assertExists(planet.sourceData, 'inclination');
		var cosInclination = Math.cos(inclination);
		var sinInclination = Math.sin(inclination);

		var oneMinusEccentricitySquared = 1 - eccentricity * eccentricity;
		//magnitude of A is a 
		A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
			 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
			 semiMajorAxis * sinPericenter * sinInclination];
		//magnitude of B is a * sqrt(|1 - e^2|)
		B = [semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * -(cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * cosPericenter * sinInclination];
		//inner product: A dot B = 0

		var timeOfPeriapsisCrossing;
		if (planet.isComet) {
			timeOfPeriapsisCrossing = planet.sourceData.timeOfPerihelionPassage;	//julian day
		}

		if (orbitType === 'parabolic') {
			eccentricAnomaly = Math.tan(argumentOfPericenter / 2);
			meanAnomaly = eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3; 
		} else if (orbitType === 'hyperbolic') {
			assert(timeOfPeriapsisCrossing !== undefined);	//only comets are hyperbolic, and all comets have timeOfPeriapsisCrossing defined
			meanAnomaly = Math.sqrt(-gravitationalParameter / semiMajorAxisCubed) * timeOfPeriapsisCrossing * 60*60*24;	//in seconds
		} else if (orbitType === 'elliptic') {
			//in theory I can say 
			//eccentricAnomaly = Math.acos((eccentricity + Math.cos(argumentOfPericenter)) / (1 + eccentricity * Math.cos(argumentOfPericenter)));
			// ... but Math.acos has a limited range ...

			if (planet.isComet) {
				timeOfPeriapsisCrossing = planet.sourceData.timeOfPerihelionPassage;	//julian day
				var timeSinceLastPeriapsisCrossing = julianDate - timeOfPeriapsisCrossing;
				meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * Math.PI / orbitalPeriod;
			} else if (planet.isAsteroid) {
				meanAnomaly = planet.sourceData.meanAnomaly;
			//} else {
			//	throw 'here';
			}
		} else {
			throw 'here';
		}
		
		//solve Newton Rhapson
		//for elliptical orbits:
		//	f(E) = M - E + e sin E = 0
		//	f'(E) = -1 + e cos E
		//for parabolic oribts:
		//	f(E) = M - E - E^3 / 3
		//	f'(E) = -1 - E^2
		//for hyperbolic orbits:
		//	f(E) = M - e sinh(E) - E
		//	f'(E) = -1 - e cosh(E)
		eccentricAnomaly = meanAnomaly;
		for (var i = 0; i < 10; ++i) {
			var func, deriv;
			if (orbitType === 'parabolic') {	//parabolic
				func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3;
				deriv = -1 - eccentricAnomaly * eccentricAnomaly;
			} else if (orbitType === 'elliptic') { 	//elliptical
				func = meanAnomaly - eccentricAnomaly + eccentricity * Math.sin(eccentricAnomaly);
				deriv = -1 + eccentricity * Math.cos(eccentricAnomaly);	//has zeroes ...
			} else if (orbitType === 'hyperbolic') {	//hyperbolic
				func = meanAnomaly + eccentricAnomaly - eccentricity  * Math.sinh(eccentricAnomaly);
				deriv = 1 - eccentricity * Math.cosh(eccentricAnomaly);
			} else {
				throw 'here';
			}

			var delta = func / deriv;
			if (Math.abs(delta) < 1e-15) break;
			eccentricAnomaly -= delta;
		}
	
		var sinEccentricAnomaly = Math.sin(eccentricAnomaly);
		var cosEccentricAnomaly = Math.cos(eccentricAnomaly);
		
		//parabolas and hyperbolas don't define orbitalPeriod
		// so no need to recalculate it
		if (orbitalPeriod !== undefined && meanAnomaly !== undefined) {
			timeOfPeriapsisCrossing = meanAnomaly * orbitalPeriod / (2 * Math.PI); //if it is a comet then we're just reversing the calculation above ...
		}

/*
r = radial distance
a = semi major axis
E = eccentric anomaly
nu = true anomaly

r^2 = a^2 (cos(E) - e)^2 + a^2 |1 - e^2| sin(E)^2
r^2 = a^2 (cos(E)^2 - 2 e cos(E) + e^2 + |1 - e^2| sin(E)^2)

=== for eccentric orbits ... e < 1 <=> e^2 < 1 <=> 1 - e^2 > 0
r^2 = a^2 (cos(E)^2 - 2 e cos(E) + e^2 + (1 - e^2) sin(E)^2)
r^2 = a^2 (cos(E)^2 - 2 e cos(E) + e^2 + sin(E)^2 - e^2 sin(E)^2)
r^2 = a^2 (1 - 2 e cos(E) + e^2 - e^2 sin(E)^2)
r^2 = a^2 (1 - 2 e cos(E) + e^2 - e^2 (1 - cos(E)^2))
r^2 = a^2 (1 - 2 e cos(E) + e^2 cos(E)^2)
r^2 = a^2 (1 - e cos(E))^2
r = a (1 - e cos(E))		true according to http://en.wikipedia.org/wiki/Eccentric_anomaly

r = a (1 - e^2) / (1 + e cos(nu)) is stated by http://www.bogan.ca/orbits/kepler/orbteqtn.html
therefore should be true:
(1 - e cos(E)) = (1 - e^2) / (1 + e cos(nu))
1 + e cos(nu) = (1 - e^2) / (1 - e cos(E))
e cos(nu) = (1 - e^2 - 1 + e cos(E)) / (1 - e cos(E))
cos(nu) = (cos(E) - e) / (1 - e cos(E))
...which checks out with Wikipedia: http://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly

r = a (1 - e^2) / (1 + e cos(nu))
r/a = (1 - e^2) / (1 + e cos(nu))
a/r = (1 + e cos(nu)) / (1 - e^2)
2 a/r - 1 = (1 + e cos(nu)) / (1 - e^2) - (1 - e^2)/(1 - e^2)
2 a/r - r/r = e (e + cos(nu)) / (1 - e^2)
(2 a - r)/(r a) = [e (e + cos(nu))] / [a (1 - e^2)]
v^2 = mu (2/r - 1/a) = mu (e^2 + e cos(nu)) / (a (1 - e^2))
the first of these is true: v^2 = mu (2/r - 1/a) according to both bogan.ca and wikipedia 

the second: v^2 = mu e (e + cos(nu)) / (a (1 - e^2)) should match v^2 = mu / (a (1 - e^2)) (1 - 2 cos(nu) + e^2)
is true only for e^2 + e cos(nu) = 1 - 2 cos(nu) + e^2
	... e cos(nu) = 1 - 2 cos(nu)
	... (2 + e) cos(nu) = 1
	... cos(nu) = 1/(e + 2)	...which doesn't look true ... so maybe that bogan.ca second velocity equation v^2 = (mu/p)(1 - 2 cos(nu) + e^2) is wrong?

=== for hyperbolic orbits ... e > 1 <=> e^2 > 1 <=> 1 - e^2 < 0
r^2 = a^2 (cos(E) 
	... we should get to 

*/
		var dt_dE;
		if (orbitType == 'parabolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 + eccentricAnomaly * eccentricAnomaly);
		} else if (orbitType == 'elliptic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 - eccentricity * Math.cos(eccentricAnomaly));
		} else if (orbitType == 'hyperbolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (eccentricity * Math.cosh(eccentricAnomaly) - 1);
		}
		var dE_dt = 1/dt_dE;
		//finally using http://en.wikipedia.org/wiki/Kepler_orbit like I should've in the first place ...
		var coeffA, coeffB;
		var coeffDerivA, coeffDerivB;
		if (orbitType == 'parabolic') {
			//...?
		} else if (orbitType == 'elliptic') { 
			coeffA = Math.cos(eccentricAnomaly) - eccentricity;
			coeffB = Math.sin(eccentricAnomaly);
			coeffDerivA = -Math.sin(eccentricAnomaly) * dE_dt;
			coeffDerivB = Math.cos(eccentricAnomaly) * dE_dt;
		} else if (orbitType == 'hyperbolic') {
			coeffA = eccentricity - Math.cosh(eccentricAnomaly);
			coeffB = Math.sinh(eccentricAnomaly);
			coeffDerivA = -Math.sinh(eccentricAnomaly) * dE_dt;
			coeffDerivB = Math.cosh(eccentricAnomaly) * dE_dt;
		}
		
		var posX = A[0] * coeffA + B[0] * coeffB;
		var posY = A[1] * coeffA + B[1] * coeffB;
		var posZ = A[2] * coeffA + B[2] * coeffB;
		//v^2 = (a^2 sin^2(E) + a^2 |1 - e^2| cos^2(E)) * mu/a^3 / (1 - e cos(E))
		//v^2 = (sin^2(E) + |1 - e^2| cos^2(E)) * mu/(a (1 - e cos(E)))
		
		//v^2 should be = mu/(a(1-e^2)) * (1 + e^2 - 2 cos(nu))	<- for nu = true anomaly
		//v^2 should also be = mu (2/r - 1/a) = mu (2a - r) / (r a)
		//... then (2 a - r) / (r a) should = (1 - 2 cos(nu) + e^2) / (a (1 - e^2))
		//...	2 a/r - 1 should = (1 - 2 cos(nu) + e^2) / ((1 + e) (1 - e))
		//...	2 a/r should = (1-e^2)/(1-e^2) + (1 - 2 cos(nu) + e^2) / (1-e^2)
		//...	2 a/r should = (2 - 2 cos(nu)) / (1-e^2)
		//...	a/r should = (1 - cos(nu)) / (1 - e^2)
		//...	r/a should = (1 - e^2) / (1 - cos(nu))
		//...	r should = a (1 - e^2) / (1 - cos(nu))
		var velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
		var velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
		var velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;
		planet.pos[0] = posX + parentPlanet.pos[0];
		planet.pos[1] = posY + parentPlanet.pos[1];
		planet.pos[2] = posZ + parentPlanet.pos[2];
		planet.vel[0] = velX + parentPlanet.vel[0];
		planet.vel[1] = velY + parentPlanet.vel[1];
		planet.vel[2] = velZ + parentPlanet.vel[2];
		vec3.copy(solarSystem.initPlanets[planetIndex].pos, planet.pos);
		vec3.copy(solarSystem.initPlanets[planetIndex].vel, planet.vel);

		planet.keplerianOrbitalElements = {
			semiMajorAxis : semiMajorAxis,
			eccentricity : eccentricity,
			eccentricAnomaly : eccentricAnomaly,
			longitudeOfAscendingNode : longitudeOfAscendingNode,
			argumentOfPericenter : argumentOfPericenter,
			inclination : inclination,
			timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
			meanAnomaly : meanAnomaly,
			orbitType : orbitType,
			orbitalPeriod : orbitalPeriod,	//only exists for elliptical orbits
			A : A,
			B : B
		};

	//using vector state (pos & vel) as provided by Horizons
	} else {

		orbitType = 'elliptic';	//better be ...

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
			planet.orbitAxis = [0,0,1];
		} else {
			var axisX = angularMomentumX / angularMomentumMag;
			var axisY = angularMomentumY / angularMomentumMag;
			var axisZ = angularMomentumZ / angularMomentumMag;
			planet.orbitAxis = [axisX, axisY, axisZ];
		}

		var basisX = [0,0,0];
		var basisY = [0,0,0];
		var basisZ = planet.orbitAxis;
		//TODO use the inclination and longitudeOfAscendingNode
		calcBasis(basisX, basisY, basisZ);
		//a[j][i] = a_ij, so our indexing is backwards, but our storage is column-major
		planet.orbitBasis = [basisX, basisY, basisZ];

		//now decompose the relative position in the coordinates of the orbit basis
		//i've eliminated all but one of the rotation degrees of freedom ...

		//http://www.mathworks.com/matlabcentral/fileexchange/31333-orbital-elements-from-positionvelocity-vectors/content/vec2orbElem.m
		//http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

		var velSq = velX * velX + velY * velY + velZ * velZ;		//(m/s)^2
		var distanceToParent = Math.sqrt(posX * posX + posY * posY + posZ * posZ);		//m
		var gravitationalParameter = gravitationalConstant * ((planet.mass || 0) + parentPlanet.mass);	//m^3 / (kg s^2) * kg = m^3 / s^2
		var specificOrbitalEnergy  = .5 * velSq - gravitationalParameter / distanceToParent;		//m^2 / s^2 - m^3 / s^2 / m = m^2/s^2, supposed to be negative for elliptical orbits
		var semiMajorAxis = -.5 * gravitationalParameter / specificOrbitalEnergy;		//m^3/s^2 / (m^2/s^2) = m
		var semiLatusRectum = angularMomentumMagSq / gravitationalParameter;			//m^4/s^2 / (m^3/s^2) = m
		eccentricity = Math.sqrt(1 - semiLatusRectum / semiMajorAxis);				//unitless (assuming elliptical orbit)

		var cosEccentricAnomaly = (1 - distanceToParent / semiMajorAxis) / eccentricity;						//unitless
		var sinEccentricAnomaly = posDotVel / (eccentricity * Math.sqrt(gravitationalParameter * semiMajorAxis));	//m^2/s / sqrt(m^3/s^2 * m) = m^2/s / sqrt(m^4/s^2) = m^2/s / (m^2/s) = unitless
		eccentricAnomaly = Math.atan2(sinEccentricAnomaly, cosEccentricAnomaly);	//radians (unitless)

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

		A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
				 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
				 semiMajorAxis * sinPericenter * sinInclination];
		B = [-semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
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

		checkPosToPosX = checkPosX - posX;
		checkPosToPosY = checkPosY - posY;
		checkPosToPosZ = checkPosZ - posZ;
		var checkPosToPosDist = Math.sqrt(checkPosToPosX * checkPosToPosX + checkPosToPosY * checkPosToPosY + checkPosToPosZ * checkPosToPosZ);
		var checkPosError = checkPosToPosDist / distanceToParent;
		if (checkPosError === checkPosError) {
			if (checkPosError > 1e-5) {	//only report significant error
				console.log(planet.name+' error of reconstructed position '+ checkPosError);
			}
		} else {	//NaN? debug!

		/*
			for (k in planet) {
				var v = planet[k];
				if (k != 'name' && typeof(v) != 'function') {
					console.log(planet.name, k, v);
				}
			}
		*/
			console.log(planet.name+' has no orbit info.  mass: '+planet.mass+' radius: '+planet.radius);
		}
				
		meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * Math.PI / orbitalPeriod;

		planet.keplerianOrbitalElements = {
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
			meanAnomaly : meanAnomaly,
			orbitType : orbitType,
			orbitalPeriod : orbitalPeriod,
			A : A,
			B : B
		};

		//not NaN, we successfully reconstructed the position
		if (checkPosError !== checkPosError) return;
	}

	//iterate around the eccentric anomaly to reconstruct the path
	var vertexes = [];
	for (var i = 0; i < orbitPathResolution; ++i) {
		var frac = i / (orbitPathResolution - 1);
		var theta = frac * 2 * Math.PI;
		var pathEccentricAnomaly = eccentricAnomaly + theta;

		//matches above
		var coeffA, coeffB;
		if (orbitType == 'parabolic') {
			//...?
		} else if (orbitType == 'elliptic') { 
			coeffA = Math.cos(pathEccentricAnomaly) - eccentricity;
			coeffB = Math.sin(pathEccentricAnomaly);
			//t = a sqrt(a/mu) (E - e sin(E))
		} else if (orbitType == 'hyperbolic') {
			coeffA = eccentricity - Math.cosh(pathEccentricAnomaly);
			coeffB = Math.sinh(pathEccentricAnomaly);
			//t = a sqrt(a/mu) (e sinh(E) - E)
		}

		var vtxPosX = A[0] * coeffA + B[0] * coeffB;
		var vtxPosY = A[1] * coeffA + B[1] * coeffB;
		var vtxPosZ = A[2] * coeffA + B[2] * coeffB;

		//add to buffer
		vertexes.push(vtxPosX);
		vertexes.push(vtxPosY);
		vertexes.push(vtxPosZ);

		//pack transparency info into the vertex
		var alpha = frac;
		vertexes.push(alpha);
	}

	planet.orbitPathObj = new glutil.SceneObject({
		mode : gl.LINE_STRIP,
		shader : assert(orbitPathShader, "failed to find orbitPathShader for planet "+planet.name),
		attrs : {
			vertex : new glutil.ArrayBuffer({
				dim : 4,
				data : vertexes
			})
		},
		uniforms : {
			color : planet.color,
			fractionOffset : 0,
		},
		blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
		pos : [0,0,0],
		angle : [0,0,0,1],
		parent : null
	});

	/*
	integration can be simulated along Keplerian orbits using the A and B vectors ...
		time advanced = simulated date - data snapshot date
		orbit fraction advanced = time advanced / orbital period
		planet position and velocity can then be computed using A and B's evaluation and derivatives
		then the alpha of the orbit needs to be updated ...
	*/
}

function recomputePlanetsAlongOrbit() {
	var timeAdvanced = julianDate - initJulianDate;
	var starSystem = orbitStarSystem;
	for (var i = 0; i < starSystem.planets.length; ++i) {	//or do it for all systems?
		var planet = starSystem.planets[i];
		if (planet.parent) {
			var ke = planet.keplerianOrbitalElements;
			var orbitType = ke.orbitType;

			var meanMotion = undefined;
			if (ke.orbitalPeriod !== undefined) {	//elliptical
				meanMotion = 2 * Math.PI / ke.orbitalPeriod;
			} else if (ke.timeOfPeriapsisCrossing !== undefined) {	//hyperbolic ... shouldn't be using mean motion?
				meanMotion = ke.meanAnomaly / (julianDate - ke.timeOfPeriapsisCrossing);
			} else {
				throw 'here';
			}
		
			//TODO don't use meanMotion for hyperbolic orbits
			var fractionOffset = timeAdvanced * meanMotion / (2 * Math.PI); 
			var theta = timeAdvanced * meanMotion;
			var pathEccentricAnomaly = ke.eccentricAnomaly + theta;
			var A = ke.A;
			var B = ke.B;
		
			//matches above
			var dt_dE;
			var semiMajorAxisCubed = ke.semiMajorAxis * ke.semiMajorAxis * ke.semiMajorAxis;
			if (orbitType == 'parabolic') {
				dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (1 + pathEccentricAnomaly * pathEccentricAnomaly);
			} else if (orbitType == 'elliptic') {
				dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (1 - ke.eccentricity * Math.cos(pathEccentricAnomaly));
			} else if (orbitType == 'hyperbolic') {
				dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (ke.eccentricity * Math.cosh(pathEccentricAnomaly) - 1);
			}
			var dE_dt = 1/dt_dE;
			var coeffA, coeffB;
			var coeffDerivA, coeffDerivB;
			if (orbitType == 'parabolic') {
				//...?
			} else if (orbitType == 'elliptic') { 
				coeffA = Math.cos(pathEccentricAnomaly) - ke.eccentricity;
				coeffB = Math.sin(pathEccentricAnomaly);
				coeffDerivA = -Math.sin(pathEccentricAnomaly) * dE_dt;
				coeffDerivB = Math.cos(pathEccentricAnomaly) * dE_dt;
			} else if (orbitType == 'hyperbolic') {
				coeffA = ke.eccentricity - Math.cosh(pathEccentricAnomaly);
				coeffB = Math.sinh(pathEccentricAnomaly);
				coeffDerivA = -Math.sinh(pathEccentricAnomaly) * dE_dt;
				coeffDerivB = Math.cosh(pathEccentricAnomaly) * dE_dt;
			}
			var posX = A[0] * coeffA + B[0] * coeffB;
			var posY = A[1] * coeffA + B[1] * coeffB;
			var posZ = A[2] * coeffA + B[2] * coeffB;
			var velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
			var velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
			var velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;
			
			planet.pos[0] = posX + planet.parent.pos[0];
			planet.pos[1] = posY + planet.parent.pos[1];
			planet.pos[2] = posZ + planet.parent.pos[2];
			planet.vel[0] = velX + planet.parent.vel[0];
			planet.vel[1] = velY + planet.parent.vel[1];
			planet.vel[2] = velZ + planet.parent.vel[2];
			planet.orbitPathObj.uniforms.fractionOffset = fractionOffset;
		}
	}
}

function downloadOrbitPaths() {
	var s = '';
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		var lines = [];

		var planet = solarSystem.planets[i];
		if (planet.orbitPathObj) {
			var name = '';
			if (planet.id) name += planet.id + ' ';
			name += planet.name;
			if (planet.parent) name += ', moon of ' + planet.parent.name;

			lines.push('#####');
			lines.push('# ' + name);

			//json.stringify chokes on float32arrays
			var alpha = [];
			var vsrc = planet.orbitPathObj.attrs.vertex.data;
			if (vsrc.length % 4 != 0) throw 'this shouldnt be';
			var vtx = [];
			for (var j = 0; j < vsrc.length; j += 4) {
				for (var k = 0; k < 3; ++k) {
					vtx[k] = vsrc[j+k] + planet.pos[k];	//position .. offset it by the planet's position
				}
				lines.push('v ' + vtx[0] + ' ' + vtx[1] + ' ' + vtx[2]);
				alpha[j/4] = vsrc[j+3];	//alpha
			}
			for (var j = 0; j < alpha.length; ++j) {
				lines.push('vt ' + alpha[j] + ' .5');
			}
			for (var j = 1; j <= alpha.length; ++j) {
				var jn = (j % alpha.length) + 1;
				lines.push('l v'+j+'/vt'+j+' v'+jn+'/vt'+jn);
			}
		}

		s += lines.join('\n') + '\n';
	}

	console.log('size of orbit path data is '+s.length);

	//TODO open in new tab?  can chrome handle opening each piece in a new tab at the same time?
	//chrome can't handle the file, it's too big ...
	var chopped = 3;	//TODO calculate by data size and by browser max download size ... which varies on the browser
	var i = 2;{//for (var i = 0; i < chopped; ++i) {	//too many of these at once, even chopped down correctly, will crash Chrome
		var start = Math.floor(i*s.length/chopped);
		var end = Math.floor((i+1)*s.length/chopped);
		var sub = s.substring(start, end);
		var a = document.createElement('a');
		console.log('size of chopped piece is '+sub.length);
		a.href = "data:text/plain,"+encodeURIComponent(sub);
		a.target = '_blank';
		a.click();
	}
}

function initOrbitPaths() {
	var gravWellShader = new ModifiedDepthShaderProgram({
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

	//for rendering the orbital path
	orbitPathShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec4 vertex;
varying float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex.xyz, 1.);
	alpha = vertex.w;
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
		fragmentCode : mlstr(function(){/*
uniform vec4 color;
uniform float fractionOffset;
varying float alpha;
void main() {
	float a = mod(alpha - fractionOffset, 1.);
	a = fract(a + 1.);
	a = a * .75 + .25;
	gl_FragColor = vec4(color.xyz, color.w * a);
}
*/})
	});

	//shader for recomputing planet positions
	//associated with texture per-planet that stores position (and maybe velocity)

	//TODO update these when we integrate!
	var calcOrbitPathStartTime = Date.now();
	//init solar system planets based on Horizons data
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		initPlanetOrbitPathObj(solarSystem.planets[i], true);
	}
	var calcOrbitPathEndTime = Date.now();

	var setPlanetTiltAngleToMoonOrbitPlane = function(planetName, moonName) {
		var planet = solarSystem.planets[solarSystem.indexes[planetName]];
		var moon = solarSystem.planets[solarSystem.indexes[moonName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.longitudeOfAscendingNode);
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.inclination);
		quat.copy(solarSystem.initPlanets[solarSystem.indexes[planetName]].tiltAngle, planet.tiltAngle);
	};

	//accepts degrees
	//TODO start at orbit axis plane rather than earth's (ie J2000) orbit axis plane
	var setPlanetTiltAngleToFixedValue = function(planetName, inclination, tiltDirection) {
		if (tiltDirection === undefined) tiltDirection = 0;
		var planet = solarSystem.planets[solarSystem.indexes[planetName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, Math.rad(tiltDirection));
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, Math.rad(inclination));
		quat.copy(solarSystem.initPlanets[solarSystem.indexes[planetName]].tiltAngle, planet.tiltAngle);
	};

	setPlanetTiltAngleToFixedValue('Mercury', 2.11/60);		//TODO tilt from mercury orbit plane.  until then it's off
	setPlanetTiltAngleToFixedValue('Venus', 177.3);			//TODO tilt from venus orbit plane.  until then, this measures 175 degrees.
	setPlanetTiltAngleToFixedValue('Earth', 23 + 1/60*(26 + 1/60*(21.4119)), 180);
	setPlanetTiltAngleToMoonOrbitPlane('Mars', 'Phobos');		//ours: 25.79, exact: 25.19
	setPlanetTiltAngleToMoonOrbitPlane('Jupiter', 'Metis');		//ours: 3.12, exact: 3.13
	setPlanetTiltAngleToMoonOrbitPlane('Saturn', 'Atlas');		//ours: 26.75, exact: 26.73
	setPlanetTiltAngleToMoonOrbitPlane('Uranus', 'Cordelia');	//ours: 97.71, exact: 97.77
	setPlanetTiltAngleToMoonOrbitPlane('Neptune', 'Galatea');	//ours: 28.365, exact: 28.32
	setPlanetTiltAngleToMoonOrbitPlane('Pluto', 'Charon');		//ours: 119, exact: 123

	//looks like low grav wells run into fp accuracy issues
	//how about extracting the depth and storing normalized values?
	var calcGravWellStartTime = Date.now();
	for (var planetIndex = 0; planetIndex < solarSystem.planets.length; ++planetIndex) {
		var planet = solarSystem.planets[planetIndex];

		/*
		gravity wells

		calculate spacetime embedding radius
		  for radial distance r, radius R, Schwarzschild radius Rs
		 inside the planet  (r <= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R (r/R)^2 ))
		 outside the planet (r >= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R)) + sqrt(4Rs(r - Rs)) - sqrt(4Rs(R - Rs))
		*/
		var R = planet.radius;
		var Rs = planet.schwarzschildRadius;
		var R_Rs = R / Rs;
		var Rs_R = Rs / R;
		var R_sqrt_R_Rs = R * Math.sqrt(R_Rs);

		//extract the normalized scalar to post-multiply after transformation (so small fields will not lose accuracy)
		planet.gravityWellScalar = Math.sqrt(R / Rs) - Math.sqrt(R / Rs - 1);

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
			z /= planet.gravityWellScalar;

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

				gravWellVtxs.push(x * planet.orbitBasis[0][0] + y * planet.orbitBasis[1][0] + z * planet.orbitBasis[2][0]);
				gravWellVtxs.push(x * planet.orbitBasis[0][1] + y * planet.orbitBasis[1][1] + z * planet.orbitBasis[2][1]);
				gravWellVtxs.push(x * planet.orbitBasis[0][2] + y * planet.orbitBasis[1][2] + z * planet.orbitBasis[2][2]);
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
		planet.gravWellObj = new glutil.SceneObject({
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
		planet.localMat = mat4.create();
		planet.gravWellObj.mvMat = mat4.create();
		planet.gravWellObj.uniforms.mvMat = planet.gravWellObj.mvMat;
	}
	var calcGravWellEndTime = Date.now();

	//on my machine this was 769ms
	console.log('calc orbit path time ',calcOrbitPathEndTime-calcOrbitPathStartTime,'ms');

	//on my machine this was 386ms
	console.log('calc grav well time ',calcGravWellEndTime-calcGravWellStartTime,'ms');

	//looks like doing these in realtime will mean toning the detail down a bit ...
	//if the image is too big, how do we downsample the skymap without lagging the whole browser?  1) browsers start using LuaJIT (not going to happen, stupid JS) 2) provide pre-computed sampled down versions.

	(function(){
		new glutil.TextureCube({
			flipY : true,
			generateMipmap : true,
			magFilter : gl.LINEAR,
			minFilter : gl.LINEAR_MIPMAP_LINEAR,
			wrap : {
				s : gl.CLAMP_TO_EDGE,
				t : gl.CLAMP_TO_EDGE
			},
			urls : skyTexFilenamePrefixes.map(function(filename) {
				return filename + Math.min(1024, glMaxCubeMapTextureSize) +'.png';
			}),
			done : function() {
				//firefox says something is binding this cubemap before we even load it
				// which can only mean i'm leaving a binding attached while initializing it (or that Firefox is buggy)
				//Chrome doesn't complain.
				console.log('loaded sky cubemap!');
				initSkyCube(this);
			}
		});
	})();
	//only do this after the orbit shader is made
	initExoplanets();

	initScene();
}

//init the "sky" cubemap (the galaxy background) once the texture for it loads
function initSkyCube(skyTex) {
	var cubeShader = new ModifiedDepthShaderProgram({
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

//TODO scale this from low value in solar system (.3 or so) to a large value at the range of the star field (1) before fading into the universe model
const float brightness = .3;

//uniform vec4 angle;
uniform vec4 viewAngle;
vec3 quatRotate( vec4 q, vec3 v ){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}

//"Reconsidering the galactic coordinate system", Jia-Cheng Liu, Zi Zhu, and Hong Zhang, Oct 20, 2010 eqn 9
//..but this is in equatorial coordinates, so rotate to get to ecliptic
const mat3 nj = mat3(	//represented transposed
	-0.054875539390, 0.494109453633, -0.867666135681,	//x-axis column
	-0.873437104725, -0.444829594298, -0.198076389622,	//y-axis column
	-0.483834991775, 0.746982248696, 0.455983794523);	//z-axis column
#define M_COS_EPSILON	0.9177546256839811400496387250314000993967056274414
#define M_SIN_EPSILON	0.3971478906347805648557880431326339021325111389160
const mat3 equatorialToEcliptical = mat3(	//represented transposed
	1., 0., 0.,
	0., M_COS_EPSILON, M_SIN_EPSILON,
	0., -M_SIN_EPSILON, M_COS_EPSILON);
void main() {
	vec3 dir = vertexv;
	dir = quatRotate(viewAngle, dir);
	dir = nj * equatorialToEcliptical * dir;
	gl_FragColor.rgb = brightness * textureCube(skyTex, dir).rgb;
	gl_FragColor.a = 1.;
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

	//my galaxy texture is centered at x+ and lies in the xy plane
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
			viewAngle : glutil.view.angle
		},
		texs : [skyTex],
		parent : null
	});
}

function setOrbitTarget(newTarget) {
	var selectingNewSystem = false;
	if (newTarget === undefined) {
		newTarget = orbitStarSystem.stars[0];
	}
	if (newTarget.isa(StarSystem)) {
		var targetSystem = newTarget;
		var i = 0;
		for (; i < targetSystem.planets.length; ++i) {
			if (targetSystem.planets[i].type !== 'barycenter') {
				newTarget = targetSystem.planets[i];
				selectingNewSystem = true;
				break;
			}
		}
		if (i == targetSystem.planets.length) {
			//assume we get something?
			//will this get reached for the named planets? if so, give them a star!
			console.log("failed to find a planet in the system to select!");
		}
	}
	if (newTarget.isa(Planet)) {
		//if we're changing star systems...
		if (newTarget.starSystem !== orbitStarSystem) {
			selectingNewSystem = true;
			orbitStarSystem = newTarget.starSystem;
		}
	}

	if (selectingNewSystem) {
		orbitTargetDistance = Math.max(100000, newTarget.radius || newTarget.equatorialRadius || 0);
		for (var i = 0; i < orbitStarSystem.planets.length; ++i) {
			var planet = orbitStarSystem.planets[i];
			orbitTargetDistance = Math.max(orbitTargetDistance,
				(planet.sourceData || {}).semiMajorAxis ||
				(planet.keplerianOrbitalElements || {}).semiMajorAxis ||
				0);
			refreshOrbitTargetDistanceText();
		}
		recomputePlanetsAlongOrbit();
	}

	if (orbitTarget !== undefined) vec3.sub(orbitOffset, orbitTarget.pos, newTarget.pos);
	orbitTarget = newTarget

	//refresh info div

	$('#orbitTargetText').text(orbitTarget.name);
	$('#infoDiv').empty();
	$('#infoDiv').append($('<hr>'));
	$('#infoDiv').append($('<div>', {text:'Type: '+orbitTarget.type}));
	if (orbitTarget.mass !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Mass: '+orbitTarget.mass+' kg'}));
	}
	if (orbitTarget.equatorialRadius !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Equatorial Radius: '+orbitTarget.equatorialRadius+' m'}));
	} else if (orbitTarget.radius !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Radius: '+orbitTarget.radius+' m'}));
	}
	if (orbitTarget.inverseFlattening !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Inverse Flattening: '+orbitTarget.inverseFlattening}));
	}
	if (orbitTarget.rotationPeriod !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Rotation Period: '+orbitTarget.rotationPeriod+' days'}));
	}
	if (orbitTarget.ringRadiusRange !== undefined) {
		$('#infoDiv').append($('<div>', {text:'Ring Min Radius: '+orbitTarget.ringRadiusRange[0]+' m'}));
		$('#infoDiv').append($('<div>', {text:'Ring Max Radius: '+orbitTarget.ringRadiusRange[1]+' m'}));
	}
	if (orbitTarget.keplerianOrbitalElements) {
		if (orbitTarget.keplerianOrbitalElements.semiMajorAxis) {
			$('#infoDiv').append($('<div>', {text:'Semi-Major Axis: '+orbitTarget.keplerianOrbitalElements.semiMajorAxis+' m'}));
		}
		if (orbitTarget.keplerianOrbitalElements.orbitType) {
			$('#infoDiv').append($('<div>', {text:'Orbit Type: '+orbitTarget.keplerianOrbitalElements.orbitType}));
		}
		if (orbitTarget.keplerianOrbitalElements.eccentricity) {
			$('#infoDiv').append($('<div>', {text:'Eccentricity: '+orbitTarget.keplerianOrbitalElements.eccentricity}));
		}
		if (orbitTarget.keplerianOrbitalElements.eccentricAnomaly) {
			$('#infoDiv').append($('<div>', {html:'Eccentric Anomaly: '+Math.deg(orbitTarget.keplerianOrbitalElements.eccentricAnomaly)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.longitudeOfAscendingNode) {
			$('#infoDiv').append($('<div>', {html:'Longitude of Ascending Node: '+Math.deg(orbitTarget.keplerianOrbitalElements.longitudeOfAscendingNode)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.argumentOfPericenter) {
			$('#infoDiv').append($('<div>', {html:'Argument of Pericenter: '+Math.deg(orbitTarget.keplerianOrbitalElements.argumentOfPericenter)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.inclination) {
			$('#infoDiv').append($('<div>', {html:'Inclination: '+Math.deg(orbitTarget.keplerianOrbitalElements.inclination)+'&deg;'}));
		}
		if (orbitTarget.keplerianOrbitalElements.timeOfPeriapsisCrossing) {
			$('#infoDiv').append($('<div>', {html:'Time of Periapsis Crossing: '+orbitTarget.keplerianOrbitalElements.timeOfPeriapsisCrossing+' days'}));
		}
		if (orbitTarget.keplerianOrbitalElements.orbitalPeriod) {
			$('#infoDiv').append($('<div>', {text:'Orbital Period: '+orbitTarget.keplerianOrbitalElements.orbitalPeriod+' days'}));
		}
	}
	$('#infoDiv').append($('<br>'));
	$('#infoDiv').append($('<br>'));

	setTimeout(function() {
		//refresh offset if the infoPanel is hidden
		$('#infoPanel').css('bottom', -15);
		if (!showBodyInfo) {
			$('#infoPanel').css('height', '104px');
		} else {
			var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
			$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);
		}
		setTimeout(function() {
			$('#infoPanel').show();
		}, 1);
	}, 1);
}

function initScene() {
	//gl.blendFunc(gl.SRC_ALPHA, gl.ONE);
	gl.enable(gl.DEPTH_TEST);
	gl.enable(gl.CULL_FACE);
	gl.depthFunc(gl.LEQUAL);
	gl.clearColor(0,0,0,0);

	//assign our initial orbitting solar system
	orbitStarSystem = solarSystem;

	var trackPlanetName = 'Earth';
	if ($.url().param('target') !== undefined) {
		trackPlanetName = $.url().param('target');
	}

	var trackPlanet = solarSystem.planets[solarSystem.indexes[trackPlanetName]];
	setOrbitTarget(trackPlanet);
	orbitTargetDistance = 2. * orbitTarget.radius;
	refreshOrbitTargetDistanceText();
	orbitDistance = orbitTargetDistance;
	$('#hoverTargetText').text(orbitTarget.name);

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
			chooseNewOrbitObject(mouseDir, false);
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
			chooseNewOrbitObject(mouseDir, true);
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
			orbitTarget,
			orbitGeodeticLocation.lat,
			orbitGeodeticLocation.lon,
			orbitGeodeticLocation.height);
	} else {
		orbitCenter = orbitTarget.pos;
	}
	var viewAngleZAxis = vec3.create();
	vec3.quatZAxis(viewAngleZAxis, glutil.view.angle);
	vec3.scale(glutil.view.pos, viewAngleZAxis, orbitDistance + (orbitTarget.equatorialRadius || orbitTarget.radius || 0));
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

	if (!integrationPaused) {
		julianDate += integrateTimeStep;
		refreshCurrentTimeText();
	}
	if (julianDate !== lastJulianDate) {
		//recompute position by delta since init planets
		recomputePlanetsAlongOrbit();

		//recompute angle based on sidereal rotation period
		for (var i = 0; i < orbitStarSystem.planets.length; ++i) {
			var planet = orbitStarSystem.planets[i];
			if (planet.rotationPeriod) {
				var zrot = [0,0,0,1];
				quat.rotateZ(zrot, zrot, ((julianDate / planet.rotationPeriod) % 1) * 2 * Math.PI + (planet.rotationOffset || 0));
				quat.multiply(planet.angle, planet.tiltAngle, zrot);
			} else {
				quat.copy(planet.angle, planet.tiltAngle);
			}
		}
	}

	//if we are close enough to the planet then rotate with it
	/*if (lastJulianDate !== julianDate &&
		orbitTarget == solarSystem.planets[solarSystem.indexes.Earth] &&	//only for earth at the moment ...
		orbitDistance < orbitTarget.radius * 10)
	{
		var deltaJulianDate = julianDate - lastJulianDate;
		var deltaAngle = quat.create();
		quat.rotateZ(deltaAngle, deltaAngle, deltaJulianDate * 2 * Math.PI);
		vec3TransformQuat(glutil.view.pos, glutil.view.pos, deltaAngle);
		quat.multiply(glutil.view.angle, deltaAngle, glutil.view.angle);
	}*/

	lastJulianDate = julianDate;

	glutil.scene.setupMatrices();
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	drawScene();
	glutil.clearAlpha();

	requestAnimFrame(update);
}

