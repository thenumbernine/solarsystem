import {getIDs} from '/js/util.js';
import {
	eclipticalToGalacticTransform,
	equatorialToEclipticalTransform
} from './vec.js';

const ids = getIDs();
window.ids = ids;
const urlparams = new URLSearchParams(window.location.search);

const cfg = {
	showLinesToOtherPlanets : false,
	showVelocityVectors : false,
	velocityVectorScale : 30,
	showRotationAxis : false,
	showOrbitAxis : false,
	showEllipseAxis : false,
	showLatAndLonLines : false,
	showGravityWell : false,
	showNames : true,
	showOrbits : true,

	showSmallBodies : true,
	allowSelectSmallBodies : true,

	showStars : true,
	starPointSizeScale : 3,
	starPointSizeBias : -3,
	starPointAlpha : 1,
	allowSelectStars : true,
	bubbleStartFadeDistInLyr : .25,
	bubbleStopFadeDistInLyr : 1.25,
	drawConstellationColorScalar : 20,
	drawConstellationPointSizeMax : 10,

	showPlanetsAsDistantPoints : true,

	showGalaxies : true,
	allowSelectGalaxies : true,
	//too small (1e+3) and the depth buffer precision gets destroyed - which is important for color-picking
	//too big (1e+20) and the scene gets zfar clipped away
	//used by milkyway.js and galaxies.js
	interGalacticRenderScale : 1e+15,

	planetScaleExaggeration : 1,

	gravityWellScaleNormalized : true,
	gravityWellScaleFixed : false,
	gravityWellScaleFixedValue : 2000,
	gravityWellRadialMinLog100 : -1,
	gravityWellRadialMaxLog100 : 2,

	overlayShowOrbitTarget : true,
	overlayShowCurrentPosition : false,

	heatAlpha : .5,

	integrationPaused : true,
	defaultIntegrateTimeStep : 1/(24*60),

	// in case you want to see things in flat earth mode ( just a half theta mapping of spherical coordinates, then flatten planet z)
	targetFlatEarthCoeff : 0,
	flatEarthConvCoeff : .03,
	flatEarthCoeff : 0,
	flatEarthRelativeEarthPos : [0,0,0],
	flatEarthRelativeEarthNorthDir : [0,0,1],

	orbitPathResolution : 500,
	ringResolution : 200,

	julianDate : 0,
	lastJulianDate : 0,
	initJulianDate : 0,

	// track ball motion variables
	mouseOverTarget : undefined,
	orbitStarSystem : undefined,	//only do surface calculations for what star system we are in
	orbitTarget : undefined,
	orbitTargetDistance : undefined,
	
	displayMethod : 'None',
	
	//ugly ugly singletons
	starSystemsHasGotResults : false,

	slideDuration : 500,
};
cfg.integrateTimeStep = cfg.defaultIntegrateTimeStep;
cfg.planetInfluences = [];

function floatToGLSL(x) {
	x = ''+x;
	if (x.indexOf('.') == -1 && x.indexOf('e') == -1) x += '.';
	return x;
}

function mat3ToGLSL(m) {
	let s = 'mat3(';
	for (let i = 0; i < 9; ++i) {
		if (i > 0) s += ',';
		s += floatToGLSL(m[i]);
	}
	s += ')';
	return s;
}

cfg.geodeticPositionCode = `

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
`;

cfg.quatRotateCode = `
vec3 quatRotate(vec4 q, vec3 v){
	return v + 2. * cross(cross(v, q.xyz) - q.w * v, q.xyz);
}
`;

//used by milkyway and skycube
//which really both do the same thing: display the milky way in the foreground or background
cfg.coordinateSystemCode = 
`
const mat3 eclipticalToGalactic = `+mat3ToGLSL(eclipticalToGalacticTransform)+`;
const mat3 equatorialToEcliptical = `+mat3ToGLSL(equatorialToEclipticalTransform)+`;
`;

export {
	ids,
	urlparams,
	cfg,
	floatToGLSL,
}
