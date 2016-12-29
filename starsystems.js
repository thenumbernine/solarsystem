/*
this hold static data used by all the star systems
kind of like starfield
but that's just a point cloud
*/

var starSystems = [];
var starSystemForNames = {};
var solarSystem;	//the one and only.  don't construct until after WebGL init so we can populate our float tex for the planets

var planetSceneObj;
var planetLatLonObj;

var latitudeMin = -90;
var latitudeMax = 90;
var latitudeStep = 5;
var longitudeMin = -180;
var longitudeMax = 180;
var longitudeStep = 5;
var latitudeDivisions = Math.floor((latitudeMax-latitudeMin)/latitudeStep);
var longitudeDivisions = Math.floor((longitudeMax-longitudeMin)/longitudeStep);

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
	var orbitPathShader = new ModifiedDepthShaderProgram({
		vertexCode : mlstr(function(){/*
attribute vec4 vertex;
varying float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform vec3 A, B;
uniform float eccentricity;
uniform float eccentricAnomaly;
uniform int orbitType;

float cosh(float x) { return .5 * (exp(x) + exp(-x)); }
float sinh(float x) { return .5 * (exp(x) - exp(-x)); }

#define TWO_PI 6.283185307179586231995926937088
void main() {
	float frac = vertex.x;
	float theta = frac * TWO_PI;
	float pathEccentricAnomaly = eccentricAnomaly + theta;

	float coeffA, coeffB;
	if (orbitType == 0) {	//elliptic
		coeffA = cos(pathEccentricAnomaly) - eccentricity;
		coeffB = sin(pathEccentricAnomaly);
	} else {	//if (orbitType == 1) {	//hyperbolic
		coeffA = eccentricity - cosh(pathEccentricAnomaly);
		coeffB = sinh(pathEccentricAnomaly);
	}
	vec3 pos = A * coeffA + B * coeffB;	
	
	alpha = frac; 

	gl_Position = projMat * (mvMat * vec4(pos, 1.));
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

	//iterate around the eccentric anomaly to reconstruct the path
	var vertexes = [];
	for (var i = 0; i < orbitPathResolution; ++i) {
		vertexes.push(i / (orbitPathResolution - 1));
	}
	
	orbitPathSceneObj = new glutil.SceneObject({
		mode : gl.LINE_STRIP,
		shader : orbitPathShader,
		attrs : {
			vertex : new glutil.ArrayBuffer({dim : 1, data : vertexes})
		},
		blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
		pos : [0,0,0],
		angle : [0,0,0,1],
		parent : null
	});


	//shader for recomputing planet positions
	//associated with texture per-planet that stores position (and maybe velocity)

	//TODO update these when we integrate!
	var calcOrbitPathStartTime = Date.now();
	//init solar system planets based on Horizons data
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		calcKeplerianOrbitalElements(solarSystem.planets[i], true);
	}
	var calcOrbitPathEndTime = Date.now();

	var setPlanetTiltAngleToMoonOrbitPlane = function(planetName, moonName) {
		var planet = solarSystem.planets[solarSystem.indexes[planetName]];
		var moon = solarSystem.planets[solarSystem.indexes[moonName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.longitudeOfAscendingNode);
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.inclination);
		quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
	};

	//accepts degrees
	//TODO start at orbit axis plane rather than earth's (ie J2000) orbit axis plane
	var setPlanetTiltAngleToFixedValue = function(planetName, inclination, tiltDirection) {
		if (tiltDirection === undefined) tiltDirection = 0;
		var planet = solarSystem.planets[solarSystem.indexes[planetName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, Math.rad(tiltDirection));
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, Math.rad(inclination));
		quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
	};

	setPlanetTiltAngleToFixedValue('Mercury', 2.11/60);		//TODO tilt from mercury orbit plane.  until then it's off
	setPlanetTiltAngleToFixedValue('Venus', 177.3);			//TODO tilt from venus orbit plane.  until then, this measures 175 degrees.
	setPlanetTiltAngleToFixedValue('Earth', 23 + 1/60*(26 + 1/60*(21.4119)), 180);
	setPlanetTiltAngleToMoonOrbitPlane('Mars', 'Phobos');		//ours: 25.79, exact: 25.19
	setPlanetTiltAngleToMoonOrbitPlane('Jupiter', 'Metis');		//ours: 3.12, exact: 3.13
	setPlanetTiltAngleToMoonOrbitPlane('Saturn', 'Atlas');		//ours: 26.75, exact: 26.73
	setPlanetTiltAngleToMoonOrbitPlane('Uranus', 'Cordelia');	//ours: 97.71, exact: 97.77
	setPlanetTiltAngleToMoonOrbitPlane('Neptune', 'Galatea');	//ours: 28.365, exact: 28.32
	setPlanetTiltAngleToMoonOrbitPlane('Pluto', 'Charon');		//ours: 119, exact: 119

	//now set all moon tilt angles to the orbit axis 
	var Sun = solarSystem.planets[0];
	for (var i = 0; i < solarSystem.planets.length; ++i) {
		var planet = solarSystem.planets[i];
		if (planet.parent !== undefined &&
			planet.parent.parent !== undefined &&
			planet.parent.parent == Sun)
		{
			quat.rotateZ(planet.tiltAngle, planet.tiltAngle, planet.keplerianOrbitalElements.longitudeOfAscendingNode);
			quat.rotateX(planet.tiltAngle, planet.tiltAngle, planet.keplerianOrbitalElements.inclination);
			quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
		}
	}

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
				// I could remap planes on a per-planet basis depending on which you are orbiting ...
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
}



//collection of all star systems
//kind of like starfield, but that's a point cloud
var StarSystems = makeClass({

	initSolarSystem : function() {
		solarSystem = new SolarSystem();
		starSystemForNames[solarSystem.name] = solarSystem;
		solarSystem.index = starSystems.length;
		starSystems.push(solarSystem);
		solarSystem.doneBuildingPlanets();
		solarSystem.magnitude = solarSystem.planets[solarSystem.indexes.Sun].magnitude;
	},

	initPlanetSceneObjs : function() {
		console.log('init planet shaders...');
		var quad = [[0,0],[0,1],[1,1],[1,1],[1,0],[0,0]];
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
		
		for (var i = 0; i < solarSystem.planets.length; ++i) { 
			var planet = solarSystem.planets[i];
			planet.initColorSchRadiusAngle();
			planet.initSceneLatLonLineObjs();
		}
	
		initOrbitPaths();
	
	//ring texture for Jupiter
	//http://www.celestiamotherlode.net/catalog/jupiter.php
	(function(){
		var jupiterRingShader = new ModifiedDepthShaderProgram({
			//vertex code matches Saturn
			vertexCode :
quatRotateCode
+ mlstr(function(){/*
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
			vertexCode : 
quatRotateCode
+ mlstr(function(){/*
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

	
	},

	planetShadersForNumStars : {},

	//request this per solar system.  rebuild if we need, return from cache if we don't.
	getPlanetShadersForNumberOfStars : function(numberOfStars) {
		if (numberOfStars <= 0) numberOfStars = 1;	//huh, I guess I have a star system with no stars ... "CFBDSIR2149 / CFBDSIR J214947.2-040308.9 / CFBDS J214947-040308"
		var shaders = this.planetShadersForNumStars[numberOfStars];
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
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
varying vec3 lightDir[NUM_STARS];		//light position
varying vec3 normal;		//surface normal

*/}) 
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
			+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
uniform vec4 color;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
uniform float ambient;
void main() {
	float litLum = 0.;
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
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
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;

*/}) 
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
			+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
uniform sampler2D tex;
uniform float ambient;
void main() {
	float litLum = 0.;
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
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
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
varying vec3 modelVertexv;
varying vec3 normal;
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];

*/}) 	
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	modelVertexv = geodeticPosition(vertex) * scaleExaggeration;
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 worldVertex = quatRotate(angle, modelVertexv) + pos;
	normal = quatRotate(angle, normalize(modelVertexv));
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - worldVertex);')
			+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(worldVertex, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ quatRotateCode
			+ mlstr(function(){/*
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
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
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

		this.planetShadersForNumStars[numberOfStars] = shaders;
		return shaders;
	}
});

var starSystemsExtra = new StarSystems(); 
