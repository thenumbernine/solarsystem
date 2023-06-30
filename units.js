let metersPerUnits = {
	Mpc :	648000/Math.PI*149597870700*1000000,	//megaparsec
	Kpc :	648000/Math.PI*149597870700*1000,	//kiloparsec
	pc :	648000/Math.PI*149597870700,	//parsec
	lyr :	9460730472580800,	//light year
	AU :	149597870700,	//astronomical unit
	ls :	299792458,	//light second
	km :	1000,	//kilometer
	m :		1,		//meter
};

let speedOfLight = metersPerUnits.ls;	// m/s
let gravitationalConstant = 6.6738480e-11;	// m^3 / (kg * s^2)

//1 = c m/s  <-> c m = 1 s
let metersPerSecond = speedOfLight;

//1 = G m^3 / (kg s^2) <=> G m^3 / (c m)^2 = kg <=> G/c^2 = kg/m
const kilogramsPerMeter = gravitationalConstant / (metersPerSecond * metersPerSecond);

export {
	metersPerUnits,
	kilogramsPerMeter,
};
