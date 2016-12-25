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

//TODO can I use glmatrix.js or does it still have math bugs? 
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

function solarSystemBarycentricToPlanetCartesian(destX, srcX, planet) {
	destX[0] = srcX[0] - planet.pos[0];
	destX[1] = srcX[1] - planet.pos[1];
	destX[2] = srcX[2] - planet.pos[2];
	var invAngle = [];
	quat.conjugate(invAngle, planet.angle);
	vec3TransformQuat(destX, destX, invAngle);
}

function planetGeodeticToSolarSystemBarycentric(destX, planet, lat, lon, height) {
	planet.geodeticPosition(destX, lat, lon, height);		// position relative to the planet center
	planetCartesianToSolarSystemBarycentric(destX, destX, planet);
}

function solarSystemBarycentricToPlanetGeodetic(planet, pos) {
	var tmp = [];
	solarSystemBarycentricToPlanetCartesian(tmp, pos, planet);
	return planet.fromGeodeticPosition(tmp);
}

//"Reconsidering the galactic coordinate system", Jia-Cheng Liu, Zi Zhu, and Hong Zhang, Oct 20, 2010 eqn 9
var eclipticalToGalacticTransform = mat3.create();	//column-major
eclipticalToGalacticTransform[0] = -0.054875539390;
eclipticalToGalacticTransform[1] = 0.494109453633;
eclipticalToGalacticTransform[2] = -0.867666135681;
eclipticalToGalacticTransform[3] = -0.873437104725;
eclipticalToGalacticTransform[4] = -0.444829594298;
eclipticalToGalacticTransform[5] = -0.198076389622;
eclipticalToGalacticTransform[6] = -0.483834991775;
eclipticalToGalacticTransform[7] = 0.746982248696;
eclipticalToGalacticTransform[8] = 0.455983794523;
var galacticToEclipticalTransform = mat3.create();
mat3.transpose(galacticToEclipticalTransform, eclipticalToGalacticTransform);

var equatorialToEclipticalTransform = mat3.create();
equatorialToEclipticalTransform[0] = 1;
equatorialToEclipticalTransform[1] = 0;
equatorialToEclipticalTransform[2] = 0;
equatorialToEclipticalTransform[3] = 0;
equatorialToEclipticalTransform[4] = 0.9177546256839811400496387250314000993967056274414;
equatorialToEclipticalTransform[5] = 0.3971478906347805648557880431326339021325111389160;
equatorialToEclipticalTransform[6] = 0;
equatorialToEclipticalTransform[7] = -0.3971478906347805648557880431326339021325111389160;
equatorialToEclipticalTransform[8] = 0.9177546256839811400496387250314000993967056274414;
var eclipticalToEquatorialTransform = mat3.create();
mat3.transpose(eclipticalToEquatorialTransform, equatorialToEclipticalTransform);

//units of Mpc
var milkyWayDistanceToCenterInKpc = 8.7;
var galaxyCenterInGalacticCoordsInMpc = [milkyWayDistanceToCenterInKpc / 1000, 0, 0];
var galaxyCenterInEclipticalCoordsInMpc = [];
vec3.transformMat3(galaxyCenterInEclipticalCoordsInMpc, galaxyCenterInGalacticCoordsInMpc, galacticToEclipticalTransform);
var galaxyCenterInEquatorialCoordsInMpc = [];
vec3.transformMat3(galaxyCenterInEquatorialCoordsInMpc, galaxyCenterInEclipticalCoordsInMpc, eclipticalToEquatorialTransform);
