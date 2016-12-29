//tried an experiment of doing surface calculations on the GPU
//it ran a lot faster than doing them in CPU for JS ... but the floating point accuracy was too low to get any good results back, even with double precision functions
//I might try worker threads later...
var CALCULATE_TIDES_WITH_GPU = false;

var colorBarHSVRange = 2/3;	// how much of the rainbow to use
var tideTexWidth = 128;
var tideTexHeight = 128;

var CalcTides = makeClass({
	init3 : function() {
		this.hsvTex = new glutil.HSVTexture(256);
	},

	invalidateForces : function() {
		//invalidate all
		for (var starSystemIndex = 0; starSystemIndex < starSystems.length; ++starSystemIndex) {
			var starSystem = starSystems[starSystemIndex];
			for (var planetIndex = 0; planetIndex < starSystem.planets.length; ++planetIndex) {
				var planet = starSystem.planets[planetIndex];
				planet.lastMeasureCalcDate = undefined;
			}
		}
	}
});

var CalcTidesCPU = makeClass({
	super : CalcTides,
	init3 : function() {
		CalcTidesCPU.superProto.init3.apply(this, arguments);

		//renders a heat map from the float values of the 'tide' attribute
		this.planetHeatMapAttrShader = new ModifiedDepthShaderProgram({
			vertexCode : mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
attribute float tide;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;
varying float tidev;
varying vec2 texCoordv;

*/}) + geodeticPositionCode + quatRotateCode + mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
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
	},

	initPlanetSceneObj : function(planet) {
		//old way, per-vertex storage, updated by CPU
		
		var tideArray = [];
		tideArray.length = planetSceneObj.attrs.vertex.count;
		for (var i = 0; i < tideArray.length; ++i) tideArray[i] = 0;
	
		planet.tideBuffer = new glutil.ArrayBuffer({dim : 1, data : tideArray, usage : gl.DYNAMIC_DRAW});
	},

	updateHeatMapAlpha : function(heatAlpha) {
		gl.useProgram(this.planetHeatMapAttrShader.obj);
		gl.uniform1f(this.planetHeatMapAttrShader.uniforms.heatAlpha.loc, heatAlpha);
		gl.useProgram(null);
	},

	// old way -- calculate on CPU, upload to vertex buffer 
	updatePlanetSceneObj : function(planet) {
		//TODO what if planet.tex === undefined?
		planet.sceneObj.shader = this.planetHeatMapAttrShader;
		planet.sceneObj.texs.length = 2;
		planet.sceneObj.texs[0] = planet.tex;
		planet.sceneObj.texs[1] = this.hsvTex;
			
		//and update calculated variable if it is out of date ...
		if (planet.lastMeasureCalcDate !== julianDate) {
			planet.lastMeasureCalcDate = julianDate;

			var measureMin = undefined;
			var measureMax = undefined;
			var vertexIndex = 0;
			for (var tideIndex = 0; tideIndex < planet.tideBuffer.data.length; ++tideIndex) {
				var lat = planet.sceneObj.attrs.vertex.data[vertexIndex++];
				var lon = planet.sceneObj.attrs.vertex.data[vertexIndex++];
var x = [];
				planetGeodeticToSolarSystemBarycentric(x, planet, lat, lon, 0);

				var t = calcMetricForce(x, planet);

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
		}
	},

	drawPlanet : function(planet) {
		planet.sceneObj.attrs.tide = planet.tideBuffer;
	}
});

var CalcTidesGPU = makeClass({
	super : CalcTides,
	init3 : function() {
		CalcTidesCPU.superProto.init3.apply(this, arguments);
		
		this.fbo = new glutil.Framebuffer();
		gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo.obj);
		gl.bindFramebuffer(gl.FRAMEBUFFER, null);
		
		//renders a heat map from the float values of the 'tide' texture 
		this.planetHeatMapTexShader = new ModifiedDepthShaderProgram({
			vertexCode : mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 pos;
uniform vec4 angle;
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;
varying vec2 texCoordv;

*/}) + geodeticPositionCode + quatRotateCode + mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
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
uniform float scaleExaggeration;

uniform int sourcePlanetIndex;		//index of the source planet.

uniform int planetStateTexHeight;
uniform sampler2D planetStateTex;	//texture of planet states.  currently [x y z mass]

// x = tangent tidal
// y = normal tidal
// z = tangent gravitational
// w = normal gravitational
uniform bvec4 flags; 

*/}) 
			+ 'const float gravitationalConstant = '+(gravitationalConstant*1e+11)+';	// m^3 / (kg * s^2)\n'
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

//while the double precision hack works great for fractions, it doesn't add too much to the extent of the range of the number
//1e+19 is still an upper bound for our floating point numbers
//soo ... scale everything down here, and up later

//pluto sits at 4.9e+12 m from the sun, but distances are cubed (and cannot exceed 1e+19) so ...
const float distancePrecisionBias = 1e-6;

//planet masses are all less than the sun at 2e+30
// so scale mass down by 1e-18 or so to put it in the same range as distance
const float massPrecisionBias = 1.;	//baking this directly into the uploaded mass: 1e-18

//...and gravity is pretty out there ... 6e-11
const float gravityPrecisionBias = 1.;//baking this directly in the constant: 1e+11

//also note, m only shows up next to g, so maybe merge them into the numbers as well?
//maybe upload M*G=mu gravitational parameters of planets instead of masses?

//finally the force rescaling
// equal to 2*d-g-m to renormalize everything
const float forcePrecisionBias = 1e-5;

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

//returns magnitude g+m-2d
double3 calcGravityAccel(
	double3 solarSystemVertex, 	//d
	double3 planetPos, 			//d
	double1 planetMass			//m
) {
	double3 x = double3_sub(solarSystemVertex, planetPos);	//d
	double1 r = double3_length(x);							//d
	double1 r2 = double1_mul(r, r);							//2d
	double1 r3 = double1_mul(r, r2);						//3d
	return double3_scale(									//g+m-2d
			x,												//d
			double1_div(									//g+m-3d
				double1_mul(								//g+m
					double1_set(-gravitationalConstant * gravityPrecisionBias), //g
					planetMass								 //m
				), 
				r3											//3d
			)
		);
}

//returns magnitude g+m-2d
double3 calcTidalAccel(
	double3 solarSystemVertex,	//d
	double3 solarSystemNormal,	//0
	double3 planetPos,			//d
	double1 planetMass			//m
) {
	double3 x = double3_sub(solarSystemVertex, planetPos);	//d
	double1 r = double3_length(x);							//d
	double1 r2 = double1_mul(r, r);							//2d
	double1 r3 = double1_mul(r, r2);						//3d
	double1 xDotN = double3_dot(x, solarSystemNormal);		//d
	return double3_scale(						//g+m-2d
			double3_sub(						//0
				double3_scale(					//0
					x, 							//d
					double1_div(				//-d
						double1_mul(			//d
							double1_set(3.), 	//0
							xDotN				//d
						), 
						r2						//2d
					)
				), 
				solarSystemNormal				//0
			),
			double1_mul(							//g+m-2d
				double1_set(distancePrecisionBias),	//d
				double1_div(						//g+m-3d
					double1_mul(					//g+m
						double1_set(gravitationalConstant * gravityPrecisionBias),	//g
						planetMass					//m
					),
					r3								//3d
				)
			)
		);
}

float calcForceValue(vec3 solarSystemVertex_s, vec3 solarSystemNormal_s) {
	//here's where we have to cycle through every planet and perform our calculation on each of them
	//then sum the results ...
	// I can guarantee already that the # of planets will exceed the GLSL max uniforms ...
	//that means we'll have to upload the positions and masses via textures ...
	//this would fit well with my plans to eventually do the keplerian orbital element pos calcs on the GPU

	double3 solarSystemVertex = double3_set(solarSystemVertex_s * distancePrecisionBias);
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

		double3 planetPos = double3_set(planetState.xyz * distancePrecisionBias);
		double1 planetMass = double1_set(planetState.w * massPrecisionBias);	//might have to pull the rescale outside the shader ...
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

	double1 force = double3_length(accel);	//units of g + m - 2 * d
	force = double1_mul(		//0
		force, 					//g+m-2d
		double1_set(forcePrecisionBias));	//2d-g-m
	return force.v.x;
}

void main() {
	vec2 latLonCoord = ((texCoordv - vec2(.5, .5)) * vec2(360., 180.)).yx;
	vec3 planetVertex = geodeticPosition(latLonCoord) * scaleExaggeration;	//planet's local frame
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
	
		//tex used for performing min/max across tide tex
		console.log('init tide reduce texs...');
		this.tideReduceTexs = [];
		var thiz = this;
		$.each([0,1],function() {
			thiz.tideReduceTexs.push(new glutil.Texture2D({
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
	
		console.log('init encode temp tex...');
		this.encodeTempTex = new glutil.Texture2D({
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

		console.log('init min reduce shader...');
		this.minReduceShader = new glutil.KernelShader({
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
		this.maxReduceShader = new glutil.KernelShader({
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
		console.log('init encode shader...');
		this.encodeShader = [];
		for (var channel = 0; channel < 4; ++channel) {
			this.encodeShader[channel] = new glutil.KernelShader({
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
	},
	
	initPlanetSceneObj : function(planet) {
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
	},
	
	updateHeatMapAlpha : function(heatAlpha) {
		gl.useProgram(this.planetHeatMapTexShader.obj);
		gl.uniform1f(this.planetHeatMapTexShader.uniforms.heatAlpha.loc, heatAlpha);
		gl.useProgram(null);
	},

	//new way -- update planet state buffer to reflect position & mass
	updatePlanetSceneObj : function(planet) {
		planet.sceneObj.shader = this.planetHeatMapTexShader;
		planet.sceneObj.texs.length = 3;
		planet.sceneObj.texs[0] = planet.tex;
		planet.sceneObj.texs[1] = planet.tideTex;
		planet.sceneObj.texs[2] = this.hsvTex;
		
		//and update calculated variable if it is out of date ...
		if (planet.lastMeasureCalcDate !== julianDate) {
			planet.lastMeasureCalcDate = julianDate;


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
	
			var viewport = gl.getParameter(gl.VIEWPORT);

			this.fbo.setColorAttachmentTex2D(0, planet.tideTex);
			gl.viewport(0, 0, tideTexWidth, tideTexHeight);
			gl.disable(gl.DEPTH_TEST);
			gl.disable(gl.CULL_FACE);
			this.fbo.draw({
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
							scaleExaggeration : planetScaleExaggeration,
							flags : flags
						},
						texs : [orbitStarSystem.planetStateTex],
					});
				}
			});

			//...then min/max reduce

			var thiz = this;
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
					
					thiz.fbo.bind();
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, thiz.tideReduceTexs[dstIndex].obj, 0);
					thiz.fbo.check();
					gl.clear(gl.COLOR_BUFFER_BIT);
					quadObj.draw({
						shader : kernelShader,
						uniforms : {
							texsize : [tideTexWidth, tideTexHeight], 
							viewsize : [width, height]
						},
						texs : [current]
					});
					thiz.fbo.unbind();

					current = thiz.tideReduceTexs[dstIndex];
					dstIndex = (dstIndex + 1) & 1;
				}
		
				//'current' has our texture

				//now that the viewport is 1x1, run the encode shader on it
				thiz.fbo.bind();
				gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, thiz.encodeTempTex.obj, 0);
				thiz.fbo.check();
				gl.viewport(0, 0, thiz.encodeTempTex.width, thiz.encodeTempTex.height);
				quadObj.draw({
					shader : thiz.encodeShader[0],
					texs : [current]
				});

				var cflUint8Result = new Uint8Array(4);
				gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, cflUint8Result);
				thiz.fbo.unbind();
				
				var cflFloat32Result = new Float32Array(cflUint8Result.buffer);
				var result = cflFloat32Result[0];
				return result;
			};

			planet.measureMin = reduce(this.minReduceShader, planet.tideTex);
			planet.measureMax = reduce(this.maxReduceShader, planet.tideTex);
//console.log('measure min', planet.measureMin, 'max', planet.measureMax);
			if (planet == orbitTarget) {
				refreshMeasureText();
			}

			//planet.sceneObj.uniforms.forceMin = planet.measureMax;
			//planet.sceneObj.uniforms.forceMax = (planet.measureMin - planet.measureMax) / colorBarHSVRange + planet.measureMax;
planet.forceMin = planet.measureMax;
planet.forceMax = (planet.measureMin - planet.measureMax) / colorBarHSVRange + planet.measureMax;

			gl.viewport.apply(gl, viewport);
			gl.enable(gl.DEPTH_TEST);
			gl.enable(gl.CULL_FACE);
			
			//...then use the float buffer, min, and max, to do the hsv overlay rendering

		}
	},
	
	drawPlanet : function(planet) {}
});

var calcTides = CALCULATE_TIDES_WITH_GPU
	? new CalcTidesGPU()
	: new CalcTidesCPU();
