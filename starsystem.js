/*
a planetary system orbitting a star / collection of stars
TODO rename to PlanetarySystem
*/
let StarSystem = makeClass({
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
		for (let i = 0; i < this.planets.length; ++i) {
			this.planets[i].index = i;
			let planet = this.planets[i];
			this.indexes[planet.name] = i;
			//while we're here...
			planet.starSystem = this;
		}
	},

	//map parent field from name to index (or should it be to object?)
	mapParents : function() {
		//convert parent from name to class (or undefined if no such name exists)
		for (let i = 0; i < this.planets.length; ++i) {
			let planet = this.planets[i];
			if (planet.parent !== undefined) {
				assert(typeof(planet.parent) === 'string');
				let index = assertExists(this.indexes, planet.parent);
				planet.parent = assertExists(this.planets, index);
			}
		}
	},

	//this is the old Planets behavior
	// which I might recreate
	// combination of Array and integration functions
	clonePlanets : function() {
		let planets = [];
		for (let i = 0; i < this.planets.length; ++i) {
			planets[i] = this.planets[i].clone();
		}
		return planets;
	},

	//static
	copyPlanets : function(dest, src) {
		assert(dest.length == src.length);
		for (let i = 0; i < src.length; ++i) {
			dest[i].copy(src[i]);
		}
	},

	updatePlanetsPos : function() {
		for (let i = 0; i < this.planets.length; ++i) {	//or do it for all systems?
			this.planets[i].updatePosVel();
		}
	}
});
