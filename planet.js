/*
http://en.wikipedia.org/wiki/Eccentric_anomaly
http://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
http://en.wikipedia.org/wiki/Kepler_orbit
http://en.wikipedia.org/wiki/Longitude_of_the_periapsis
http://en.wikipedia.org/wiki/Mean_longitude
http://en.wikipedia.org/wiki/Orbital_elements
http://en.wikipedia.org/wiki/True_anomaly#From_the_eccentric_anomaly
http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
http://spaceweather.com/
http://www.bogan.ca/orbits/kepler/orbteqtn.html
http://www.mathworks.com/matlabcentral/fileexchange/31333-orbital-elements-from-positionvelocity-vectors/content/vec2orbElem.m
https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
https://space.stackexchange.com/questions/8911/determining-orbital-position-at-a-future-point-in-time

TODO theres some vel errors in here
in calcKOEFromPosVel and calcOrbitBasis the vel is converted from m/day to m/s (to work with gravitationalConstant)
 but later it looks like the Vinti A, B expect m/day
 and later the parent-rel vel is computed by adding the KOE-computed vel to the parent vel
 ... pretty sure one of these is in m/s and the other in m/day
though atm i don't think I use vel but TODO fix this
*/
import {vec3, quat} from '/js/gl-matrix-3.4.1/index.js';
import {mathRad, mathDeg, assertExists} from '/js/util.js';
import {cfg} from './globals.js';
import {kilogramsPerMeter, gravitationalConstant} from './units.js';
import {starSystemsExtra} from './starsystems.js';
import {calcTides} from './calctides.js';
import {calcBasis} from './vec.js';

class Planet {
	constructor(args) {
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

		this.maxSemiMajorAxisOfEllipticalSatellites = 0;	//used for spherical bounding volumes ... excluding parabolic/hyperbolic satellites ...
	}

	clone() {
		const p = new Planet();
		vec3.copy(p.pos, this.pos);
		vec3.copy(p.vel, this.vel);
		quat.copy(p.angle, this.angle);
		quat.copy(p.tiltAngle, this.tiltAngle);
		return p;
	}

	copy(src) {
		this.pos[0] = src.pos[0];
		this.pos[1] = src.pos[1];
		this.pos[2] = src.pos[2];
		this.vel[0] = src.vel[0];
		this.vel[1] = src.vel[1];
		this.vel[2] = src.vel[2];
	}

	//TODO the math ... don't be lazy
	fromGeodeticPosition(pos) {
		const x = pos[0];
		const y = pos[1];
		const z = pos[2];
		//if (this.inverseFlattening === undefined) {
		const r2 = Math.sqrt(x*x + y*y);
		const r = Math.sqrt(r2*r2 + z*z);
		const phi = Math.atan2(z, r2);
		const lambda = Math.atan2(y, x);

		const equatorialRadius = this.equatorialRadius !== undefined ? this.equatorialRadius : this.radius;
		if (equatorialRadius === undefined) {
			return {
				lat : 0,
				lon : 0,
				height : 0
			};
		}
		return {
			lat : mathDeg(phi),
			lon : mathDeg(lambda),
			height : r - equatorialRadius
		}
		//} else it is a nonlinear problem and I will think about it later 
	}
	
	//no longer used for mesh construction -- that all goes on in the shader
	//this is only used for getting positions for updating the tidal array calculations
	// and for (geosynchronously) orbiting geodetic locations (which is currently disabled)
	geodeticPosition(destX, lat, lon, height) {
		const phi = mathRad(lat);
		const lambda = mathRad(lon);
		const cosPhi = Math.cos(phi);
		const sinPhi = Math.sin(phi);

		const equatorialRadius = this.equatorialRadius !== undefined ? this.equatorialRadius : this.radius;
		if (equatorialRadius === undefined) {
			destX[0] = 0;
			destX[1] = 0;
			destX[2] = 0;
		} else if (this.inverseFlattening !== undefined) {
			const eccentricitySquared = (2 * this.inverseFlattening - 1) / (this.inverseFlattening * this.inverseFlattening);
			const sinPhiSquared = sinPhi * sinPhi;
			const N = equatorialRadius / Math.sqrt(1 - eccentricitySquared * sinPhiSquared);
			const NPlusH = N + height;
			destX[0] = NPlusH * cosPhi * Math.cos(lambda);
			destX[1] = NPlusH * cosPhi * Math.sin(lambda);
			destX[2] = (N * (1 - eccentricitySquared) + height) * sinPhi;
		} else {
			const NPlusH = equatorialRadius + height;
			destX[0] = NPlusH * cosPhi * Math.cos(lambda);
			destX[1] = NPlusH * cosPhi * Math.sin(lambda);
			destX[2] = NPlusH * Math.sin(phi);
		}
	}

	geodeticNormal(lat, lon) {
		const phi = mathRad(lat)
		const lambda = mathRad(lon)
		const cosPhi = Math.cos(phi)
		const sinPhi = Math.sin(phi)

		const equatorialRadius = this.equatorialRadius !== undefined
			? this.equatorialRadius
			: this.radius;
		if (equatorialRadius === undefined) {
			throw "don't know how to calculate this planet's surface "+this;
		}
		if (this.inverseFlattening !== undefined) {
			//mind you this is the planet shape eccentricity, not to be confused with the orbit path eccentricity
			const eccentricitySquared = (2 * this.inverseFlattening - 1) / (this.inverseFlattening * this.inverseFlattening);
			const sinPhiSquared = sinPhi * sinPhi;
			const N = equatorialRadius / Math.sqrt(1 - eccentricitySquared * sinPhiSquared);
			const oneMinusEccSq = 1 - eccentricitySquared;
			const NPlusH = N + height;
			const NPrime = sinPhi * eccentricitySquared / (equatorialRadius * equatorialRadius) * N * N * N;
			const sinPhi = Math.sin(phi);
			const cosPhi = Math.cos(phi);
			const sinLambda = Math.sin(lambda);
			const cosLambda = Math.cos(lambda);
			const nx = oneMinusEccSq * NPlusH * cosPhi * cosLambda * (NPrime * sinPhi + NPlusH * cosPhi);
			const ny = oneMinusEccSq * NPlusH * cosPhi * sinLambda * (NPrime * sinPhi + NPlusH * cosPhi);
			const nz = -NPlusH * cosPhi * (NPrime * cosPhi - NPlusH * sinPhi);
			const l = Math.sqrt(nx*nx + ny*ny + nz*nz);
			return [nx/l, ny/l, nz/l];
		} else {
			const x = cosPhi * Math.cos(lambda);
			const y = cosPhi * Math.sin(lambda);
			const z = Math.sin(phi);
			return [x, y, z]
		}
	}

	initColorSchRadiusAngle() {
		const colors = {
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
		const color = colors[this.name];
		if (!color) {
			//console.log("failed to find color for "+this.name);
			this.color = [Math.random(), Math.random(), Math.random(), 1];
			vec3.normalize(this.color, this.color);
		} else {
			this.color = [color[0], color[1], color[2], 1];
		}
		this.schwarzschildRadius = 2 * this.mass * kilogramsPerMeter;
	}

	//TODO no more multiple copies of tide array per-planet
	initSceneLatLonLineObjs() {
		if (this.radius === undefined) {
			//only/always use a point/basis/etc?
		} else {
			this.sceneObj = starSystemsExtra.planetSceneObj;
			calcTides.initPlanetSceneObj(this);
		}
	}

	updateSceneObj() {
		const planetShaders = starSystemsExtra.getPlanetShadersForNumberOfStars(this.starSystem.stars.length);

		//update tide attributes
		const useOverlay = cfg.overlayShowOrbitTarget && cfg.displayMethod != 'None';
		if (!useOverlay) {
			if (this.tex === undefined) {
				this.sceneObj.shader = planetShaders.colorShader;
				this.sceneObj.texs.length = 0;
			} else {
				if (this.ringObj) {
					this.sceneObj.shader = planetShaders.ringShadowShader;
					this.sceneObj.texs.length = 2;
					if (this.ringTransparencyTex !== undefined) {
						this.sceneObj.texs[1] = this.ringTransparencyTex;
					} else if (this.ringColorTex !== undefined) {
						this.sceneObj.texs[1] = this.ringColorTex;
					} else {
						//...or instead of throwing, should I just assume full-alpha?
						throw 'planet has a ringObj but no ring texture';
					}
				} else {
					this.sceneObj.shader = planetShaders.texShader;
					this.sceneObj.texs.length = 1;
				}
				this.sceneObj.texs[0] = this.tex;
			}
		} else {
			calcTides.updatePlanetSceneObj(this);
		}
	}

	calcOrbitBasis() {
		// based on this position and velocity, find plane of orbit
		let parentBody = this.parent;
		if (parentBody === undefined) {
			//console.log(this.name+' has no orbit parent');
			this.orbitAxis = [0,0,1];
			this.orbitBasis = [[1,0,0],[0,1,0],[0,0,1]];
			return;
		} else if (parentBody.name.match(/Barycenter$/g)) {
			parentBody = parentBody.parent;
		}

		//consider position relative to orbiting parent
		const posX = this.pos[0] - parentBody.pos[0];
		const posY = this.pos[1] - parentBody.pos[1];
		const posZ = this.pos[2] - parentBody.pos[2];

		//convert from m/day to m/s to coincide with the units of our gravitational constant
		const velX = (this.vel[0] - parentBody.vel[0]) / (60 * 60 * 24);
		const velY = (this.vel[1] - parentBody.vel[1]) / (60 * 60 * 24);
		const velZ = (this.vel[2] - parentBody.vel[2]) / (60 * 60 * 24);

		const angularMomentumX = posY * velZ - posZ * velY; //m^2/s
		const angularMomentumY = posZ * velX - posX * velZ;
		const angularMomentumZ = posX * velY - posY * velX;
		const angularMomentumMagSq = angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY + angularMomentumZ * angularMomentumZ;		//m^4/s^2
		const angularMomentumMag = Math.sqrt(angularMomentumMagSq);

		if (angularMomentumMag < 1e-9) {
			this.orbitAxis = [0,0,1];
		} else {
			const axisX = angularMomentumX / angularMomentumMag;
			const axisY = angularMomentumY / angularMomentumMag;
			const axisZ = angularMomentumZ / angularMomentumMag;
			this.orbitAxis = [axisX, axisY, axisZ];
		}

		const basisX = [0,0,0];
		const basisY = [0,0,0];
		const basisZ = this.orbitAxis;
		calcBasis(basisX, basisY, basisZ);
		//a[j][i] = a_ij, so our indexing is backwards, but our storage is column-major
		this.orbitBasis = [basisX, basisY, basisZ];
	}

	calcKOEFromPosVel() {
//console.log("this.mass", this.mass);		
		// based on this position and velocity, find plane of orbit
		let parentBody = this.parent;
		if (parentBody === undefined) {
			console.log('planet',this.name,'has no parent');
			this.calcOrbitBasis();
			return;
		}
//console.log("parentBody.mass", parentBody.mass);		

		//consider position relative to orbiting parent
		// should I be doing the same thing with the velocity?  probably...
		const posX = this.pos[0] - parentBody.pos[0];
		const posY = this.pos[1] - parentBody.pos[1];
		const posZ = this.pos[2] - parentBody.pos[2];

		//convert from m/day to m/s to coincide with the units of our gravitational constant
		const velX = (this.vel[0] - parentBody.vel[0]) / (60 * 60 * 24);
		const velY = (this.vel[1] - parentBody.vel[1]) / (60 * 60 * 24);
		const velZ = (this.vel[2] - parentBody.vel[2]) / (60 * 60 * 24);

		const posDotVel = posX * velX + posY * velY + posZ * velZ;	//m^2/s

		const angularMomentumX = posY * velZ - posZ * velY; //m^2/s
		const angularMomentumY = posZ * velX - posX * velZ;
		const angularMomentumZ = posX * velY - posY * velX;
		const angularMomentumMagSq = angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY + angularMomentumZ * angularMomentumZ;		//m^4/s^2
		const angularMomentumMag = Math.sqrt(angularMomentumMagSq);
		
		//now decompose the relative position in the coordinates of the orbit basis
		//i've eliminated all but one of the rotation degrees of freedom ...

		//http://www.mathworks.com/matlabcentral/fileexchange/31333-orbital-elements-from-positionvelocity-vectors/content/vec2orbElem.m
		//http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto

		const velSq = velX * velX + velY * velY + velZ * velZ;		//(m/s)^2
		const distanceToParent = Math.sqrt(posX * posX + posY * posY + posZ * posZ);		//m
//console.log("gravitationalConstant", gravitationalConstant);
		const gravitationalParameter = gravitationalConstant * ((this.mass || 0) + parentBody.mass);	//m^3 / (kg s^2) * kg = m^3 / s^2
//console.log("gravitationalParameter", gravitationalParameter);
		const specificOrbitalEnergy  = .5 * velSq - gravitationalParameter / distanceToParent;		//m^2 / s^2 - m^3 / s^2 / m = m^2/s^2, supposed to be negative for elliptical orbits
		const semiMajorAxis = -.5 * gravitationalParameter / specificOrbitalEnergy;		//m^3/s^2 / (m^2/s^2) = m
//console.log("semiMajorAxis", semiMajorAxis);
		const semiLatusRectum = angularMomentumMagSq / gravitationalParameter;			//m^4/s^2 / (m^3/s^2) = m
		const eccentricity = Math.sqrt(1 - semiLatusRectum / semiMajorAxis);		//e, unitless (assuming elliptical orbit)

		let orbitType = undefined;
		const parabolicEccentricityEpsilon = 1e-7;
//		if (Math.abs(eccentricity - 1) < parabolicEccentricityEpsilon) {
//console.log("danger! got a parabolic orbit for ",this," based on eccentricity epsilon",Math.abs(eccentricity-1));
//			orbitType = 'parabolic';
//		} else 
		if (eccentricity > 1) {
			orbitType = 'hyperbolic';
		} else {
			orbitType = 'elliptic';
		}

		const cosEccentricAnomaly = (1 - distanceToParent / semiMajorAxis) / eccentricity;						//unitless
		const sinEccentricAnomaly = posDotVel / (eccentricity * Math.sqrt(gravitationalParameter * semiMajorAxis));	//m^2/s / sqrt(m^3/s^2 * m) = m^2/s / sqrt(m^4/s^2) = m^2/s / (m^2/s) = unitless
		const eccentricAnomaly = Math.atan2(sinEccentricAnomaly, cosEccentricAnomaly);	//E, in radians (unitless)

		const sinInclination = Math.sqrt(angularMomentumX * angularMomentumX + angularMomentumY * angularMomentumY) / angularMomentumMag;	//unitless
		const cosInclination = angularMomentumZ / angularMomentumMag;	//unitless
		const inclination = Math.atan2(sinInclination, cosInclination);	//i

		const sinPericenter = ((velX * angularMomentumY - velY * angularMomentumX) / gravitationalParameter - posZ / distanceToParent) / (eccentricity * sinInclination);
		const cosPericenter = (angularMomentumMag * velZ / gravitationalParameter - (angularMomentumX * posY - angularMomentumY * posX) / (angularMomentumMag * distanceToParent)) / (eccentricity * sinInclination);
		const argumentOfPeriapsis = Math.atan2(sinPericenter, cosPericenter);	//omega

		const cosAscending = -angularMomentumY / (angularMomentumMag * sinInclination);
		const sinAscending = angularMomentumX / (angularMomentumMag * sinInclination);
		const longitudeOfAscendingNode = Math.atan2(sinAscending, cosAscending);	//Omega

		const semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;	//m^3
		let orbitalPeriod;
		if (orbitType == 'elliptic') {
			orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed  / gravitationalParameter) / (60*60*24);	//julian day
//override for known planets?  don't forget to factor in the barycenter 
if (this.orbitalPeriod !== undefined) {
	console.log('for planet',this.name,'orbitalPeriod calculated',orbitalPeriod,'provided',this.orbitalPeriod);
	orbitalPeriod = this.orbitalPeriod;
}
		}

		const longitudeOfPeriapsis = longitudeOfAscendingNode + argumentOfPeriapsis;	//omega-bar
		
		let meanAnomaly;
		if (orbitType == 'parabolic') {
			meanAnomaly = eccentricAnomaly + eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3;
		} else if (orbitType == 'hyperbolic') {
			meanAnomaly = eccentricity * Math.sinh(eccentricAnomaly) - eccentricAnomaly;
		} else if (orbitType == 'elliptic') {
			meanAnomaly = eccentricAnomaly - eccentricity * sinEccentricAnomaly;
		}

		const meanLongitude = meanAnomaly + longitudeOfPeriapsis;
		
		//can I do this?
		const epoch = cfg.initJulianDate;
		const meanAnomalyAtEpoch = meanAnomaly;

		const timeOfPeriapsisCrossing = -meanAnomaly / Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (60*60*24);	//julian day

		const A = [semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
				 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
				 semiMajorAxis * sinPericenter * sinInclination];
		const B = [-semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
				 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
				 semiMajorAxis * Math.sqrt(1 - eccentricity * eccentricity) * cosPericenter * sinInclination];
//console.log("A", A);
//console.log("B", B);

		/*
		to convert back:
		pos[i] = A[i] * (cosEccentricAnomaly - eccentricity) + B[i] * sinEccentricAnomaly
		rDot[i] = (-A[i] * sinEccentricAnomaly + B[i] * cosEccentricAnomaly) * Math.sqrt(gravitationalParameter / semiMajorAxisCubed) / (1 - eccentricity * cosEccentricAnomaly)
		*/
//console.log("pos", posX, posY, posZ);   
		const checkPosX = A[0] * (cosEccentricAnomaly - eccentricity) + B[0] * sinEccentricAnomaly;
		const checkPosY = A[1] * (cosEccentricAnomaly - eccentricity) + B[1] * sinEccentricAnomaly;
		const checkPosZ = A[2] * (cosEccentricAnomaly - eccentricity) + B[2] * sinEccentricAnomaly;
//console.log("checkPos", checkPosX, checkPosY, checkPosZ);   
		//used by planets to offset reconstructed orbit coordinates to exact position of this
		const checkPosToPosX = checkPosX - posX;
		const checkPosToPosY = checkPosY - posY;
		const checkPosToPosZ = checkPosZ - posZ;
//console.log("checkPosToPos", checkPosToPosX, checkPosToPosY, checkPosToPosZ);   
		const checkPosToPosDist = Math.sqrt(checkPosToPosX * checkPosToPosX + checkPosToPosY * checkPosToPosY + checkPosToPosZ * checkPosToPosZ);
//console.log("checkPosToPosDist", checkPosToPosDist);  	
		const checkPosError = checkPosToPosDist / distanceToParent;
		if (checkPosError === checkPosError) {
			if (checkPosError > 1e-5) {	//only report significant error
				console.log(this.name+' error of reconstructed position '+ checkPosError);
			}
		} else {	//NaN? debug!
			console.log(this.name+' has no orbit info.  mass: '+this.mass+' radius: '+this.radius);
		}

		this.keplerianOrbitalElements = {
			relVelSq : velSq,
			gravitationalParameter : gravitationalParameter,
			specificOrbitalEnergy : specificOrbitalEnergy,
			distanceToParent : distanceToParent,
			semiMajorAxis : semiMajorAxis,
			semiLatusRectum : semiLatusRectum,
			inclination : inclination,
			argumentOfPeriapsis : argumentOfPeriapsis,
			longitudeOfAscendingNode : longitudeOfAscendingNode,
			timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
			meanAnomaly : meanAnomaly,
			orbitalPeriod : orbitalPeriod,
			//should I do this?
			epoch : epoch,
			meanAnomalyAtEpoch : meanAnomalyAtEpoch,
			//the following are used for the orbit path shader:
			orbitType : orbitType,
			eccentricity : eccentricity,
			eccentricAnomaly : eccentricAnomaly,
			A : A,
			B : B,
			//this is used for drawing but not in the shader
			fractionOffset : 0
		};

		//not NaN, we successfully reconstructed the position
		if (checkPosError !== checkPosError) {
			console.log('check-position error was NaN for planet ',this);
			this.calcOrbitBasis();
			return;
		}
		
		//for elliptic orbits,
		// we can accumulate & store the largest semi major axis of all children
		//but for parabolic/hyperbolic ...
		// there is no bound ... so maybe those objects should be stored/rendered separately?
		if (parentBody !== undefined && orbitType == 'elliptic') {
			parentBody.maxSemiMajorAxisOfEllipticalSatellites = Math.max(semiMajorAxis, parentBody.maxSemiMajorAxisOfEllipticalSatellites);
		}

		this.renderOrbit = true;

		this.calcOrbitBasis();
	}

	getKOEFromSourceData() {

		// based on this position and velocity, find plane of orbit
		const parentBody = this.parent;
		if (parentBody === undefined) {
			this.calcOrbitBasis();
			return;
		}

		//used by planets to offset reconstructed orbit coordinates to exact position of this
		const checkPosToPosX = 0;
		const checkPosToPosY = 0;
		const checkPosToPosZ = 0;

		/*
		we have:
			eccentricity
			semiMajorAxis
			pericenterDistance
		we compute:
		*/

		const eccentricity = assertExists(this.sourceData, 'eccentricity');

		let orbitType = undefined;
		const parabolicEccentricityEpsilon = 1e-7;
		if (Math.abs(eccentricity - 1) < parabolicEccentricityEpsilon) {
			orbitType = 'parabolic';
		} else if (eccentricity > 1) {
			orbitType = 'hyperbolic';
		} else {
			orbitType = 'elliptic';
		}

		let pericenterDistance;
		let semiMajorAxis = undefined;
		if (this.isComet) {
			pericenterDistance = assert(this.sourceData.perihelionDistance);
			if (orbitType !== 'parabolic') {
				semiMajorAxis = pericenterDistance / (1 - eccentricity);
			}	//otherwise, if it is parabolic, we don't get the semi-major axis ...
		} else if (this.isAsteroid) {
			semiMajorAxis = assert(this.sourceData.semiMajorAxis);
			pericenterDistance = semiMajorAxis * (1 - eccentricity);
		} else if (this.isExoplanet) {
			semiMajorAxis = this.sourceData.semiMajorAxis || 0;
		}

		const gravitationalParameter = gravitationalConstant * parentBody.mass;	//assuming the comet mass is negligible, since the comet mass is not provided
		const semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis;
		
		//orbital period is only defined for circular and elliptical orbits (not parabolic or hyperbolic)
		let orbitalPeriod = undefined;
		if (orbitType === 'elliptic') {
			orbitalPeriod = 2 * Math.PI * Math.sqrt(semiMajorAxisCubed / gravitationalParameter) / (60*60*24);	//julian day
		}

		const longitudeOfAscendingNode = assertExists(this.sourceData, 'longitudeOfAscendingNode');
		const cosAscending = Math.cos(longitudeOfAscendingNode);
		const sinAscending = Math.sin(longitudeOfAscendingNode);

		const argumentOfPeriapsis = assertExists(this.sourceData, 'argumentOfPeriapsis');
		const cosPericenter = Math.cos(argumentOfPeriapsis);
		const sinPericenter = Math.sin(argumentOfPeriapsis);

		const inclination = assertExists(this.sourceData, 'inclination');
		const cosInclination = Math.cos(inclination);
		const sinInclination = Math.sin(inclination);

		const oneMinusEccentricitySquared = 1 - eccentricity * eccentricity;
		//magnitude of A is a 
		const A = [
			semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
			 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
			 semiMajorAxis * sinPericenter * sinInclination
		];
		//magnitude of B is a * sqrt(|1 - e^2|)
		const B = [
			semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * -(cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
			 semiMajorAxis * Math.sqrt(Math.abs(oneMinusEccentricitySquared)) * cosPericenter * sinInclination
		];
		//inner product: A dot B = 0

		let timeOfPeriapsisCrossing;
		if (this.isComet) {
			timeOfPeriapsisCrossing = this.sourceData.timeOfPerihelionPassage;	//julian day
		}

		let meanAnomaly = this.sourceData.meanAnomaly;
		let eccentricAnomaly = this.sourceData.eccentricAnomaly;
		const epoch = this.sourceData.epoch;
		if (orbitType === 'parabolic') {
			if (eccentricAnomaly === undefined || eccentricAnomaly === null) {
				eccentricAnomaly = Math.tan(argumentOfPeriapsis / 2);
			}
			if (meanAnomaly === undefined || meanAnomaly === null) {
				meanAnomaly = eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3; 
			}
		} else if (orbitType === 'hyperbolic') {
			assert(timeOfPeriapsisCrossing !== undefined);	//only comets are hyperbolic, and all comets have timeOfPeriapsisCrossing defined
			if (meanAnomaly === undefined || meanAnomaly === null) {
				meanAnomaly = Math.sqrt(-gravitationalParameter / semiMajorAxisCubed) * timeOfPeriapsisCrossing * 60*60*24;	//in seconds
			}
		} else if (orbitType === 'elliptic') {
			//in theory I can say 
			//eccentricAnomaly = Math.acos((eccentricity + Math.cos(argumentOfPeriapsis)) / (1 + eccentricity * Math.cos(argumentOfPeriapsis)));
			// ... but Math.acos has a limited range ...

			if (meanAnomaly === undefined || meanAnomaly === null) {
				if (this.isComet) {
					timeOfPeriapsisCrossing = this.sourceData.timeOfPerihelionPassage;	//julian day
					const timeSinceLastPeriapsisCrossing = cfg.initJulianDate - timeOfPeriapsisCrossing;
					meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * Math.PI / orbitalPeriod;
				} else if (this.isAsteroid) {
					meanAnomaly = this.sourceData.meanAnomalyAtEpoch + 2 * Math.PI / orbitalPeriod * (cfg.initJulianDate - epoch);
				} else if (this.isExoplanet) {
					meanAnomaly = this.sourceData.meanAnomaly !== undefined ? this.sourceData.meanAnomaly : 0;
				} else {
				//	throw 'here';
				}
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
		if (eccentricAnomaly === undefined || eccentricAnomaly === null) {
			eccentricAnomaly = meanAnomaly;
			for (let i = 0; i < 10; ++i) {
				let func, deriv;
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

				const delta = func / deriv;
				if (Math.abs(delta) < 1e-15) break;
				eccentricAnomaly -= delta;
			}
		}

		const sinEccentricAnomaly = Math.sin(eccentricAnomaly);
		const cosEccentricAnomaly = Math.cos(eccentricAnomaly);
		
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
		let dt_dE;
		if (orbitType == 'parabolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 + eccentricAnomaly * eccentricAnomaly);
		} else if (orbitType == 'elliptic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 - eccentricity * Math.cos(eccentricAnomaly));
		} else if (orbitType == 'hyperbolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (eccentricity * Math.cosh(eccentricAnomaly) - 1);
		}
		const dE_dt = 1/dt_dE;
		//finally using http://en.wikipedia.org/wiki/Kepler_orbit like I should've in the first place ...
		let coeffA, coeffB;
		let coeffDerivA, coeffDerivB;
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
		
		const posX = A[0] * coeffA + B[0] * coeffB;
		const posY = A[1] * coeffA + B[1] * coeffB;
		const posZ = A[2] * coeffA + B[2] * coeffB;
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
		const velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
		const velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
		const velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;
		this.pos[0] = posX + parentBody.pos[0];
		this.pos[1] = posY + parentBody.pos[1];
		this.pos[2] = posZ + parentBody.pos[2];
		this.vel[0] = velX + parentBody.vel[0];
		this.vel[1] = velY + parentBody.vel[1];
		this.vel[2] = velZ + parentBody.vel[2];
		const solarSystem = starSystemsExtra.solarSystem;
		vec3.copy(solarSystem.initPlanets[this.index].pos, this.pos);
		vec3.copy(solarSystem.initPlanets[this.index].vel, this.vel);

		this.keplerianOrbitalElements = {
			semiMajorAxis : semiMajorAxis,
			eccentricity : eccentricity,
			eccentricAnomaly : eccentricAnomaly,
			longitudeOfAscendingNode : longitudeOfAscendingNode,
			argumentOfPeriapsis : argumentOfPeriapsis,
			inclination : inclination,
			timeOfPeriapsisCrossing : timeOfPeriapsisCrossing,
			meanAnomaly : meanAnomaly,
			meanAnomalyAtEpoch : this.sourceData.meanAnomalyAtEpoch,
			epoch : epoch,
			orbitType : orbitType,
			orbitalPeriod : orbitalPeriod,	//only exists for elliptical orbits
			A : A,
			B : B
		};
		
		//for elliptic orbits,
		// we can accumulate & store the largest semi major axis of all children
		//but for parabolic/hyperbolic ...
		// there is no bound ... so maybe those objects should be stored/rendered separately?
		if (parentBody !== undefined && orbitType == 'elliptic') {
			parentBody.maxSemiMajorAxisOfEllipticalSatellites = Math.max(semiMajorAxis, parentBody.maxSemiMajorAxisOfEllipticalSatellites);
		}

		this.renderOrbit = true;

		this.calcOrbitBasis();
	}

	updatePosVel() {
		const timeAdvanced = cfg.julianDate - cfg.initJulianDate;
		if (!this.parent) return;
		const ke = this.keplerianOrbitalElements;
		const orbitType = ke.orbitType;

		//https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time

		let meanAnomaly, meanMotion;
		if (orbitType == 'elliptic') {
			assertExists(ke, 'orbitalPeriod');
			assertExists(ke, 'meanAnomalyAtEpoch');	//not always there.  esp not in planets.
			meanMotion = 2 * Math.PI / ke.orbitalPeriod;
			meanAnomaly = ke.meanAnomalyAtEpoch + meanMotion * (cfg.julianDate - ke.epoch);
		} else if (orbitType == 'hyperbolic') {
			assertExists(ke, 'timeOfPeriapsisCrossing');
			meanAnomaly = ke.meanAnomaly;
			meanMotion = ke.meanAnomaly / (cfg.julianDate - ke.timeOfPeriapsisCrossing);
		} else if (orbitType == 'parabolic') {
console.log('parabolic orbit for planet',this);			
			throw 'got a parabolic orbit';
		} else {
			throw 'here';
		}

		const eccentricity = ke.eccentricity;

		//solve eccentricAnomaly from meanAnomaly via Newton Rhapson
		//for elliptical orbits:
		//	f(E) = M - E + e sin E = 0
		//	f'(E) = -1 + e cos E
		//for parabolic oribts:
		//	f(E) = M - E - E^3 / 3
		//	f'(E) = -1 - E^2
		//for hyperbolic orbits:
		//	f(E) = M - e sinh(E) - E
		//	f'(E) = -1 - e cosh(E)
		let eccentricAnomaly = meanAnomaly;
		for (let i = 0; i < 10; ++i) {
			let func, deriv;
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

			const delta = func / deriv;
			if (Math.abs(delta) < 1e-15) break;
			eccentricAnomaly -= delta;
		}

		//TODO don't use meanMotion for hyperbolic orbits
		const fractionOffset = timeAdvanced * meanMotion / (2 * Math.PI); 
		// TODO where did I get this from?  I think it's more like 'theta == eccentricAnomaly' cuz I'm doubling the value by adding this ...
		//const theta = timeAdvanced * meanMotion;
		const pathEccentricAnomaly = eccentricAnomaly;	// + theta;
		const A = ke.A;
		const B = ke.B;

		//matches above
		let dt_dE;
		const semiMajorAxisCubed = ke.semiMajorAxis * ke.semiMajorAxis * ke.semiMajorAxis;
		if (orbitType == 'parabolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (1 + pathEccentricAnomaly * pathEccentricAnomaly);
		} else if (orbitType == 'elliptic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (1 - ke.eccentricity * Math.cos(pathEccentricAnomaly));
		} else if (orbitType == 'hyperbolic') {
			dt_dE = Math.sqrt(semiMajorAxisCubed / ke.gravitationalParameter) * (ke.eccentricity * Math.cosh(pathEccentricAnomaly) - 1);
		}
		const dE_dt = 1/dt_dE;
		let coeffA, coeffB;
		let coeffDerivA, coeffDerivB;
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
		const posX = A[0] * coeffA + B[0] * coeffB;
		const posY = A[1] * coeffA + B[1] * coeffB;
		const posZ = A[2] * coeffA + B[2] * coeffB;
		const velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
		const velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
		const velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;
		
		this.pos[0] = posX + this.parent.pos[0];
		this.pos[1] = posY + this.parent.pos[1];
		this.pos[2] = posZ + this.parent.pos[2];
		this.vel[0] = velX + this.parent.vel[0];
		this.vel[1] = velY + this.parent.vel[1];
		this.vel[2] = velZ + this.parent.vel[2];

		//ke == this.keplerianOrbitalElements
		ke.meanAnomaly = meanAnomaly;
		ke.eccentricAnomaly = eccentricAnomaly;
		ke.fractionOffset = fractionOffset;
	}
}

/*
TODO Apsidal precession:
https://en.wikipedia.org/wiki/Apsidal_precession

epsilon = 24 pi^3 a^2 / (T^2 c^2 (1 + e^2))
a = semi-major axis,		[a] = m
T = period of revolution,	[T] = s/revolution ... [revolution] = 1
c = speed of light,			[c] = m/s
e = eccentricity of orbit,	[e] = 1 = m/m
epsilon = precession during one period of revolution, in radians
[epsilon] = [a]^2 / ([T]^2 [c]^2 (1 + [e]^2))
[epsilon] = m^2 / (s^2/rev^2 m^2/s^2) = rev^2 = 1
*/

export {Planet};
