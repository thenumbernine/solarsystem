/*
this hold static data used by all the star systems
kind of like starfield
but that's just a point cloud
*/

let starSystems = [];
let starSystemForNames = {};
let solarSystem;	//the one and only.  don't construct until after WebGL init so we can populate our float tex for the planets

//ugly ugly singletons
let starSystemsHasGotResults = false;

let planetSceneObj;
let planetLatLonObj;

let latitudeMin = -90;
let latitudeMax = 90;
let latitudeStep = 5;
let longitudeMin = -180;
let longitudeMax = 180;
let longitudeStep = 5;
let latitudeDivisions = Math.floor((latitudeMax-latitudeMin)/latitudeStep);
let longitudeDivisions = Math.floor((longitudeMax-longitudeMin)/longitudeStep);

function initOrbitPaths() {
	let gravWellShader = new ModifiedDepthShaderProgram({
		vertexCode : `
in vec4 vertex;
out float alpha;
uniform mat4 mvMat;
uniform mat4 projMat;
void main() {
	vec4 vtx4 = mvMat * vec4(vertex.xyz, 1.);
	alpha = vertex.w;
	vtx4 = flatEarthXForm(vtx4);
	gl_Position = projMat * vtx4;
	gl_PointSize = 4.;
	gl_Position.z = depthfunction(gl_Position);
}
`,
		fragmentCode : `
uniform vec4 color;
in float alpha;
out vec4 fragColor;
void main() {
	fragColor = vec4(color.xyz, color.w * alpha);
}
`
	});

	//for rendering the orbital path
	let orbitPathShader = new ModifiedDepthShaderProgram({
		vertexCode : `
in vec4 vertex;
out float alpha;
uniform mat4 viewMatInv;
uniform mat4 localMat;
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

	//I have to separate the modelTransform and viewTransform for any object that has a non-identity modelTransform that wants to use the flatEarthXForm
	gl_Position = projMat * (viewMatInv * flatEarthXForm(localMat * vec4(pos, 1.)));
	gl_Position.z = depthfunction(gl_Position);
}
`,
		fragmentCode : `
uniform vec4 color;
uniform float fractionOffset;
in float alpha;
out vec4 fragColor;
void main() {
	float a = mod(alpha - fractionOffset, 1.);
	a = fract(a + 1.);
	a = a * .75 + .25;
	fragColor = vec4(color.xyz, color.w * a);
}
`
	});

	//iterate around the eccentric anomaly to reconstruct the path
	let vertexes = [];
	for (let i = 0; i < orbitPathResolution; ++i) {
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
	let calcOrbitPathStartTime = Date.now();
	//init solar system planets based on Horizons data
	for (let i = 0; i < solarSystem.planets.length; ++i) {
		solarSystem.planets[i].calcKOEFromPosVel();
	}
	let calcOrbitPathEndTime = Date.now();

	let setPlanetTiltAngleToMoonOrbitPlane = function(planet, moonName) {
		let moon = solarSystem.planets[solarSystem.indexes[moonName]];
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.longitudeOfAscendingNode);
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, moon.keplerianOrbitalElements.inclination);
		quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
	};

	//accepts degrees
	//TODO start at orbit axis plane rather than earth's (ie J2000) orbit axis plane
	let setPlanetTiltAngleToFixedValue = function(planet, inclination, tiltDirection) {
		if (tiltDirection === undefined) tiltDirection = 0;
		quat.rotateZ(planet.tiltAngle, planet.tiltAngle, Math.rad(tiltDirection));
		quat.rotateX(planet.tiltAngle, planet.tiltAngle, Math.rad(inclination));
		quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
	};

	let setBy = {
		Mercury:{inclination:2.11/60},//TODO tilt from mercury orbit plane.  until then it's off
		Venus:{inclination:177.3},//TODO tilt from venus orbit plane.  until then, this measures 175 degrees.
		Earth:{inclination:23 + 1/60*(26 + 1/60*(21.4119)), tiltDirection:180},
		Mars:{toMoon:'Phobos'},//ours: 25.79, exact: 25.19
		Jupiter:{toMoon:'Metis'},//ours: 3.12, exact: 3.13
		Saturn:{toMoon:'Atlas'},//ours: 26.75, exact: 26.73
		Uranus:{toMoon:'Cordelia'},//ours: 97.71, exact: 97.77
		Neptune:{toMoon:'Galatea'},//ours: 28.365, exact: 28.32
		Pluto:{toMoon:'Charon'},//ours: 119, exact: 119
	};

	//now set all moon tilt angles to the orbit axis 
	let Sun = solarSystem.planets[0];
	for (let i = 0; i < solarSystem.planets.length; ++i) {
		let planet = solarSystem.planets[i];
		let sb = setBy[planet.name];
		if (sb && sb.inclination !== undefined) {
			setPlanetTiltAngleToFixedValue(planet, sb.inclination, sb.tiltDirection);
		} else if (sb && sb.toMoon) {
			setPlanetTiltAngleToMoonOrbitPlane(planet, sb.toMoon);
		} else if (planet.keplerianOrbitalElements !== undefined
				&& planet.keplerianOrbitalElements.longitudeOfAscendingNode !== undefined
				&& planet.keplerianOrbitalElements.inclination !== undefined
		) {
			quat.rotateZ(planet.tiltAngle, planet.tiltAngle, planet.keplerianOrbitalElements.longitudeOfAscendingNode);
			quat.rotateX(planet.tiltAngle, planet.tiltAngle, planet.keplerianOrbitalElements.inclination);
			quat.copy(solarSystem.initPlanets[planet.index].tiltAngle, planet.tiltAngle);
		}
	}

	//looks like low grav wells run into fp accuracy issues
	//how about extracting the depth and storing normalized values?
	let calcGravWellStartTime = Date.now();
	for (let planetIndex = 0; planetIndex < solarSystem.planets.length; ++planetIndex) {
		let planet = solarSystem.planets[planetIndex];

		/*
		gravity wells

		calculate spacetime embedding radius
		  for radial distance r, radius R, Schwarzschild radius Rs
		 inside the planet  (r <= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R (r/R)^2 ))
		 outside the planet (r >= R): z(r) = R sqrt(R/Rs) (1 - sqrt(1 - Rs/R)) + sqrt(4Rs(r - Rs)) - sqrt(4Rs(R - Rs))
		*/
		let R = planet.radius;
		let Rs = planet.schwarzschildRadius;
		let R_Rs = R / Rs;
		let Rs_R = Rs / R;
		let R_sqrt_R_Rs = R * Math.sqrt(R_Rs);

		//extract the normalized scalar to post-multiply after transformation (so small fields will not lose accuracy)
		planet.gravityWellScalar = Math.sqrt(R / Rs) - Math.sqrt(R / Rs - 1);

		//populate first vertex
		let gravWellVtxs = [0,0,0,1];

		let gravWellIndexes = [];

		let rimax = 60;		//was 200 r and 60 th, but I added a lot of planets.  need to occlude these based on distance/angle ...
		let thimax = 60;

		for (let ri = 1; ri < rimax; ++ri) {
			let r = R * Math.pow(100, ri / rimax * (gravityWellRadialMaxLog100 - gravityWellRadialMinLog100) + gravityWellRadialMinLog100);
			//max radial dist is R * Math.pow(100, gravityWellRadialMaxLog100)

			let z;
			if (r <= R) {
				let r_R = r / R;
				z = R_sqrt_R_Rs * (1 - Math.sqrt(1 - Rs_R * r_R * r_R));
			} else {
				z = R_sqrt_R_Rs * (1 - Math.sqrt(1 - Rs_R))
					+ Math.sqrt(4 * Rs * (r - Rs))
					- Math.sqrt(4 * Rs * (R - Rs));
			}
			z /= planet.gravityWellScalar;

			for (let thi = 0; thi < thimax; ++thi) {
				let th = 2 * Math.PI * thi / thimax;

				let x = r * Math.cos(th);
				let y = r * Math.sin(th);

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
				let tau = ri/(rimax-1);
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
	let calcGravWellEndTime = Date.now();

	//on my machine this was 769ms
	console.log('calc orbit path time ',calcOrbitPathEndTime-calcOrbitPathStartTime,'ms');

	//on my machine this was 386ms
	console.log('calc grav well time ',calcGravWellEndTime-calcGravWellStartTime,'ms');
}



//collection of all star systems
//kind of like starfield, but that's a point cloud
let StarSystems = makeClass({

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
		let quad = [[0,0],[0,1],[1,1],[1,1],[1,0],[0,0]];
		let triIndexArray = [];
		let latLonIndexArray = [];
		let vertexArray = [];

		for (let loni=0; loni <= longitudeDivisions; ++loni) {
			let lon = longitudeMin + loni * longitudeStep;
			for (let lati=0; lati <= latitudeDivisions; ++lati) {
				let lat = latitudeMin + lati * latitudeStep;

				vertexArray.push(lat);
				vertexArray.push(lon);

				if (loni < longitudeDivisions && lati < latitudeDivisions) {
					for (let j = 0; j < quad.length; ++j) {
						let ofs = quad[j];
						let index = (lati + ofs[0]) + (latitudeDivisions + 1) * (loni + ofs[1]);
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
		let vertexBuffer = new glutil.ArrayBuffer({dim : 2, data : vertexArray});

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
		
		for (let i = 0; i < solarSystem.planets.length; ++i) { 
			let planet = solarSystem.planets[i];
			planet.initColorSchRadiusAngle();
			planet.initSceneLatLonLineObjs();
		}
	
		initOrbitPaths();
	
	//ring texture for Jupiter
	//http://www.celestiamotherlode.net/catalog/jupiter.php
	(function(){
		let jupiterRingShader = new ModifiedDepthShaderProgram({
			//vertex code matches Saturn
			vertexCode :
quatRotateCode
+ `
#define M_PI 3.1415926535897931
in vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float ringMinRadius;
uniform float ringMaxRadius;
uniform vec3 pos;
uniform vec4 angle;
out vec2 texCoordv;
out vec3 worldPosv;	//in world coordinates

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
	vec4 eyePos = mvMat * flatEarthXForm(modelPos);
	worldPosv = quatRotate(angle, modelPos.xyz) + pos;
	
	gl_Position = projMat * eyePos;
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode : `
#define NUM_STARS 1
in vec2 texCoordv;
in vec3 worldPosv;

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

out vec4 fragColor;
void main() {
	float luminance = 1.;
	//notice that I have to scale down the meters here for shader accuracy to work
	luminance = min(luminance, step(1., sphereIntersect(worldPosv * 1e-8, sunDir[0] * 1e-8, pos * 1e-8, planetRadius * 1e-8)));

	fragColor = texture(colorTex, texCoordv);
	fragColor.rgb *= sqrt(luminance);
}
`
		});

		let planet = solarSystem.planets[solarSystem.indexes.Jupiter];

		let texSrcInfo = [
			{field:'ringColorTex', url:'textures/jupiter-rings-color.png', format:gl.RGBA},
		];
		let numLoaded = 0;
		let onLoadRingTex = function(url) {
			++numLoaded;
			if (numLoaded < texSrcInfo.length) return;
			if (numLoaded > texSrcInfo.length) throw "already created the rings!";

			//done! create the ring object
			let vertexes = [];
			for (let i = 0; i < ringResolution; ++i) {
				let f = i / (ringResolution - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}

			let texs = [];
			for (let i = 0; i < texSrcInfo.length; ++i) {
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
		let saturnRingShader = new ModifiedDepthShaderProgram({
			vertexCode : 
quatRotateCode
+ `
#define M_PI 3.1415926535897931
in vec2 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float ringMinRadius;
uniform float ringMaxRadius;

//these are actually baked into the mvMat, but I need them separate for converting coordinates to world coordinates (rather than eye coordinates)
//pos is in the scene, which is offset by the oribitting planet
//techinically the offset doesn't matter since i'm using it to offset vertices from the planet, then sphere-intersect-testing rays from those vertices against the planet
uniform vec3 pos;
uniform vec4 angle;

out vec2 texCoordv;
out vec3 worldPosv;	//in world coordinates

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
	vec4 eyePos = mvMat * flatEarthXForm(modelPos);
	worldPosv = quatRotate(angle, modelPos.xyz) + pos;

	gl_Position = projMat * eyePos;
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode : `
#define NUM_STARS 1
//described here:
//http://www.mmedia.is/~bjj/data/s_rings/index.html

in vec2 texCoordv;
in vec3 worldPosv;

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

out vec4 fragColor;
void main() {
	vec3 color = texture(colorTex, texCoordv).rgb;
	fragColor.rgb = color;

	float backScattered = texture(backScatteredTex, texCoordv).r;
	float forwardScattered = texture(forwardScatteredTex, texCoordv).r;
	float litSide = mix(backScattered, forwardScattered, backToFrontLitBlend);

	float unlitSide = texture(unlitSideTex, texCoordv).r;

	//shadow from the sun
	//raytrace from worldPosv along sunDir (non-normalized)
	//if it intersects with a sphere at pos with radius planetRadius then we're in shadow
	//if you want to do a soft shadow based on the sun's radius and distance (length of sunDir) then you can

	float luminance = 1.;
	//notice that I have to scale down the meters here for shader accuracy to work
	luminance = min(luminance, step(1., sphereIntersect(worldPosv * 1e-8, sunDir[0] * 1e-8, pos * 1e-8, planetRadius * 1e-8)));
	luminance *= mix(unlitSide, litSide, lookingAtLitSide);

	fragColor.rgb *= sqrt(luminance);

	float transparency = texture(transparencyTex, texCoordv).r;
	fragColor.a = transparency;
}
`
		});

		let planet = solarSystem.planets[solarSystem.indexes.Saturn];

		let texSrcInfo = [
			{field:'ringColorTex', url:'textures/saturn-rings-color.png', format:gl.RGBA},
			{field:'ringBackScatteredTex', url:'textures/saturn-rings-back-scattered.png', format:gl.LUMINANCE},
			{field:'ringForwardScatteredTex', url:'textures/saturn-rings-forward-scattered.png', format:gl.LUMINANCE},
			{field:'ringTransparencyTex', url:'textures/saturn-rings-transparency.png', format:gl.LUMINANCE},
			{field:'ringUnlitSideTex', url:'textures/saturn-rings-unlit-side.png', format:gl.LUMINANCE}
		];
		let numLoaded = 0;
		let onLoadRingTex = function(url) {
			++numLoaded;
			if (numLoaded < texSrcInfo.length) return;
			if (numLoaded > texSrcInfo.length) throw "already created the rings!";

			//done! create the ring object
			let vertexes = [];
			for (let i = 0; i < ringResolution; ++i) {
				let f = i / (ringResolution - 1);
				vertexes.push(1);
				vertexes.push(f);
				vertexes.push(0);
				vertexes.push(f);
			}

			let texs = [];
			for (let i = 0; i < texSrcInfo.length; ++i) {
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
		let shaders = this.planetShadersForNumStars[numberOfStars];
		if (shaders !== undefined) return shaders;

		shaders = {};
		shaders.colorShader = new ModifiedDepthShaderProgram({
			vertexCode :`
#define NUM_STARS `+numberOfStars+`

in vec2 vertex;		//lat/lon pairs:
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
out vec3 lightDir[NUM_STARS];		//light position
out vec3 normal;		//surface normal

`
+ geodeticPositionCode 
+ quatRotateCode 
+ `

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
`
+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
+ `
	vec4 vtx4 = mvMat * flatEarthXForm(vec4(vtx3, 1.));
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode :`
#define NUM_STARS `+numberOfStars+`

uniform vec4 color;
in vec3 lightDir[NUM_STARS];
in vec3 normal;
uniform float ambient;
out vec4 fragColor;
void main() {
	float litLum = 0.;
` + unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
+ `
	float luminance = min(1., litLum);
	fragColor = color * max(ambient, sqrt(luminance));
}
`,
			uniforms : {
				color : [1,1,1,1],
				pointSize : 4
			}
		});

		shaders.texShader = new ModifiedDepthShaderProgram({
			vertexCode :
`#define NUM_STARS `+numberOfStars+`

in vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
out vec2 texCoordv;
out vec3 lightDir[NUM_STARS];
out vec3 normal;

`
+ geodeticPositionCode 
+ quatRotateCode 
+ `

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
` + unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
+ `
	vec4 vtx4 = mvMat * flatEarthXForm(vec4(vtx3, 1.));
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode :
`#define NUM_STARS `+numberOfStars+`

in vec2 texCoordv;
in vec3 lightDir[NUM_STARS];
in vec3 normal;
uniform sampler2D tex;
uniform float ambient;
out vec4 fragColor;
void main() {
	float litLum = 0.;
` + unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
+ `
	float luminance = min(1., litLum);
	fragColor = texture(tex, texCoordv) * max(ambient, sqrt(luminance));
}
`,
			uniforms : {
				tex : 0
			}
		});

		shaders.ringShadowShader = new ModifiedDepthShaderProgram({
			vertexCode :
`#define NUM_STARS `+numberOfStars+`

in vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
out vec3 modelVertexv;
out vec3 normal;
out vec2 texCoordv;
out vec3 lightDir[NUM_STARS];

`
+ geodeticPositionCode 
+ quatRotateCode 
+ `

void main() {
	//vertex is really the lat/lon in degrees
	modelVertexv = geodeticPosition(vertex) * scaleExaggeration;
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 worldVertex = quatRotate(angle, modelVertexv) + pos;
	normal = quatRotate(angle, normalize(modelVertexv));
`
+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - worldVertex);')
+ `
	vec4 vtx4 = mvMat * flatEarthXForm(vec4(worldVertex, 1.));
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode :
`#define NUM_STARS `+numberOfStars + `
` + quatRotateCode
+ `
in vec2 texCoordv;
in vec3 lightDir[NUM_STARS];
in vec3 normal;
in vec3 modelVertexv;
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

out vec4 fragColor;
void main() {
	float luminance = 1.;

	//inverse rotate lightDir[0]
	//TODO how to combine this per-light-source ... depenends on the intensity of each ...
	vec3 lightDirInModelSpace = quatRotate(vec4(-angle.xyz, angle.w), lightDir[0]);
	float intersectPos = ringIntersect(modelVertexv, lightDirInModelSpace);
	if (intersectPos >= 0. && intersectPos <= 1.) {
		luminance *= texture(ringTransparencyTex, vec2(intersectPos, .5)).r;
	}

	float litLum = 0.;
` + unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
+ `	luminance *= min(1., litLum);

	fragColor.rgb = texture(tex, texCoordv).rgb * max(ambient, sqrt(luminance));
	fragColor.a = 1.;
}
`,
			uniforms : {
				tex : 0,
				ringTransparencyTex : 1
			}
		});

		this.planetShadersForNumStars[numberOfStars] = shaders;
		return shaders;
	},

	initExoplanets : function() {
		//just over a meg, so I might as well ajax it
		let exoplanetURL = 'exoplanet/openExoplanetCatalog.json';
		let thiz = this;
		$.ajax({
			url : exoplanetURL,
			dataType : 'json'
		}).error(function() {
			console.log('failed to get exoplanets from '+exoplanetURL+' , trying again...');
			setTimeout(function() {
				thiz.initExoplanets();
			}, 5000);
		}).done(function(data) {
			thiz.processResults(data);
		});
	},

	processResults : function(results) {
		starSystemsHasGotResults = true;

		//process results
		$.each(results.systems, function(i,systemInfo) {
			let systemName = assertExists(systemInfo, 'name');

			let rightAscension = systemInfo.rightAscension;
			if (rightAscension === undefined) {
				console.log('failed to find right ascension for system '+systemName);
				return;
			}
			let declination = systemInfo.declination;
			if (declination === undefined) {
				console.log('failed to find declination for system '+systemName);
				return;
			}
			let cosRA = Math.cos(rightAscension);
			let sinRA = Math.sin(rightAscension);
			let cosDec = Math.cos(declination);
			let sinDec = Math.sin(declination);
			//convert to coordinates
			let pos = [];
			pos[0] = cosRA * cosDec;
			pos[1] = sinRA * cosDec;
			pos[2] = sinDec;
			//rotate for earth's tilt
			let epsilon = Math.rad(23 + 1/60*(26 + 1/60*(21.4119)));
			let cosEps = Math.cos(epsilon);
			let sinEps = Math.sin(epsilon);
			let yn = cosEps * pos[1] + sinEps * pos[2];
			pos[2] = -sinEps * pos[1] + cosEps * pos[2];
			pos[1] = yn;
			//distance
			let distance = systemInfo.distance;
			if (distance === undefined) {
//				console.log('failed to find distance for system '+systemName);
				return;
			}
			pos[0] *= distance;
			pos[1] *= distance;
			pos[2] *= distance;


			let starSystem = new StarSystem();
			starSystem.name = systemName;
			starSystem.sourceData = systemInfo;
			vec3.copy(starSystem.pos, pos);

			//TODO absoluate magnitude of the collective system (sum of all parts?)

			let minAbsMag = undefined;

			$.each(assertExists(systemInfo, 'bodies'), function(j, bodyInfo) {

				let name = assertExists(bodyInfo, 'name');
				let radius = bodyInfo.radius;
				if (radius === undefined) {
					if (bodyInfo.type !== 'barycenter') {
						//console.log('no radius for body '+name);
						//if planets don't have radii they can be estimated by the planet's mass or distance from sun or both?
						radius = 7e+7;	// use a jupiter radius
					}
				}

				let mass = bodyInfo.mass;
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

				let body = mergeInto(new Planet(), {
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
				bodyInfo.meanAnomaly = Math.PI * 2 * (j + Math.random()) / systemInfo.bodies.length;

				starSystem.planets.push(body);
				if (body.type == 'star') starSystem.stars.push(body);
				
				if (bodyInfo.visualMagnitude !== undefined) {
					//absolute magnitude based on visual magnitude
					body.magnitude = bodyInfo.visualMagnitude - 5 * (Math.log10(distance / metersPerUnits.pc) - 1);
					//...???
					minAbsMag = minAbsMag === undefined ? body.magnitude : Math.min(minAbsMag, body.magnitude);
				}
			});
			if (minAbsMag !== undefined) {
				starSystem.magnitude = minAbsMag;
			}

			starSystem.doneBuildingPlanets();

			for (let i = 0; i < starSystem.planets.length; ++i) {
				let body = starSystem.planets[i];

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
				body.sourceData.argumentOfPeriapsis = body.sourceData.longitudeOfPeriapsis - body.sourceData.longitudeOfAscendingNode;
				if (body.sourceData.argumentOfPeriapsis == 0) {
					body.sourceData.argumentOfPeriapsis = Math.PI * 2 * i / starSystem.planets.length;
				}

				vec3.add(body.pos, body.pos, starSystem.pos);

				body.initColorSchRadiusAngle();
				body.initSceneLatLonLineObjs();

				if (body.pos[0] !== body.pos[0]) {
					console.log('system '+starSystem.name+' planet '+body.name+' has bad pos');
				}

				body.getKOEFromSourceData();
			}

			starSystem.initPlanets = starSystem.clonePlanets();
			starSystemForNames[starSystem.name] = starSystem;
			starSystem.index = starSystems.length;
			starSystems.push(starSystem);
		});

		//now that we've built all our star system data ... add it to the star field
		starfield.addStarSystems();
	},

	//call this once we get back our star data
	initStarsControls : function() {
		//TODO maybe some pages or something for this, like the asteroid/smallbody search?
		for (let i = 0; i < starSystems.length; ++i) {
			(function(){
				let starSystem = starSystems[i];
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
});

let starSystemsExtra = new StarSystems(); 
