//our solar system
var SolarSystem = makeClass({
	super : StarSystem,
	name : 'Solar System',
	init : function() {
		SolarSystem.super.apply(this, arguments);

		//add our initial planets ...
		this.planets.push(mergeInto(new Planet(), {
			id : 10,
			name : 'Sun',
			mass : 1.9891e+30,
			radius : 6.960e+8,
			type : 'star'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 199,
			name : 'Mercury',
			parent : 'Sun',
			mass : 3.302e+23,
			radius : 2.440e+6,
			equatorialRadius : 2440e+3,
			rotationPeriod : 58.6462,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 299,
			name : 'Venus',
			parent : 'Sun',
			mass : 4.8685e+24,
			radius : 6.0518e+6,
			equatorialRadius : 6051.893e+3,
			rotationPeriod : -243.0185,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 3,
			name : 'Earth Barycenter',
			parent : 'Sun',
			type : 'barycenter',
			orbitalPeriod : 365.256363004
		}));
		/*
		hmm, technically the earth barycenter orbits the sun every 365 days
		and the earth and moon both orbit the barycenter ever 28 days ...
		*/
		this.planets.push(mergeInto(new Planet(), {
			id : 399,
			name : 'Earth',
			parent : 'Earth Barycenter',
			mass : 5.9736e+24,
			radius : 6.37101e+6,
			equatorialRadius : 6378.136e+3,
			inverseFlattening : 298.257223563,
			rotationPeriod : (23 + (56 + 4.09053083288 / 60) / 60) / 24,	//sidereal day
			orbitalPeriod : 27.321662,	//(tidally locked) moon's rotation period = orbit around barycenter
			rotationOffset : -2/24 * 2*Math.PI,	//looks about 2 hours off
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 301,
			name : 'Moon',
			parent : 'Earth Barycenter',
			mass : 7.349e+22,
			radius : 1.73753e+6,
			rotationPeriod : 27.321662,
			orbitalPeriod : 27.321662,
			rotationOffset : 180 * Math.PI / 360,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 4,
			name : 'Mars Barycenter',
			parent : 'Sun',
			type : 'barycenter'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 499,
			name : 'Mars',
			parent : 'Mars Barycenter',
			mass : 6.4185e+23,
			radius : 3.3899e+6,
			equatorialRadius : 3397e+3,
			inverseFlattening : 154.409,
			rotationPeriod : 24.622962/24,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 5,
			name : 'Jupiter Barycenter',
			parent : 'Sun',
			type : 'barycenter'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 599,
			name : 'Jupiter',
			parent : 'Jupiter Barycenter',
			mass : 1.89813e+27,
			radius : 6.9911e+7,
			equatorialRadius : 71492e+3,
			inverseFlattening : 1/0.06487,
			ringRadiusRange : [102200000,227000000],
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 6,
			name : 'Saturn Barycenter',
			parent : 'Sun',
			type : 'barycenter'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 699,
			name : 'Saturn',
			parent : 'Saturn Barycenter',
			mass : 5.68319e+26,
			radius : 5.8232e+7,
			equatorialRadius : 60268e+3,
			inverseFlattening : 1/0.09796,
			ringRadiusRange : [74510000,140390000],
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 7,
			name : 'Uranus Barycenter',
			parent : 'Sun',
			type : 'barycenter'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 799,
			name : 'Uranus',
			parent : 'Uranus Barycenter',
			mass : 8.68103e+25,
			radius : 2.5362e+7,
			equatorialRadius : 25559e+3,
			inverseFlattening : 1/0.02293,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 8,
			name : 'Neptune Barycenter',
			parent : 'Sun',
			type : 'barycenter'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 899,
			name : 'Neptune',
			parent : 'Neptune Barycenter',
			mass : 1.0241e+26,
			radius : 2.4624e+7,
			equatorialRadius : 24766e+3,
			inverseFlattening : 1/0.0171,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 9,
			name : 'Pluto Barycenter',
			parent : 'Sun',
			type : 'barycenter'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 999,
			name : 'Pluto',
			parent : 'Pluto Barycenter',
			mass : 1.314e+22,
			radius : 1.151e+6,
			rotationPeriod : -6.387230,
			type : 'planet'
		}));

		//sun is our only star
		this.stars.push(this.planets[0]);

		//start out with the current planets
		var currentIDs = {};
		for (var i = 0; i < this.planets.length; ++i) {
			currentIDs[this.planets[i].id] = this.planets[i];
		}

		//I need to fix up my export script ...
		if (horizonsDynamicData.coords.length != horizonsStaticData.length) throw 'static to dynamic data lengths differ: dynamic has '+horizonsDynamicData.coords.length+' while static has '+horizonsStaticData.length;
		for (var i = 0; i < horizonsDynamicData.coords.length; ++i) {
			var dynamicData = horizonsDynamicData.coords[i];
			var staticData = horizonsStaticData[i];
			if (dynamicData.id != staticData.id) {
				console.log('staticData',staticData);
				console.log('dynamicData',dynamicData);
				throw "got a mismatch in ids between static data and dynamic data";
			}
			if (dynamicData.name.match(/^L[1-5]/)) {
				//391-395 are Lagrangian points - skip
			} else if (dynamicData.name.match(/Barycenter$/)) {
				//1-9 are Barycenter points
			} else {
				if (currentIDs[dynamicData.id]) {
					//if the id already exists then skip it
					currentIDs[dynamicData.id].magnitude = staticData.magnitude;
				} else {
					//... finally, add the planet
					/*
					Things we are going to want:
					mass
					radius
					equatorialRadius (all but sun, pluto, moon)
					inverseFlattening (all but sun, mercury, venus, moon, pluto)
					magntidue
					*/
					var parentName = staticData.parent;
					var parentIndex = this.planets.findWithComparator(null, function(planet) {
						return planet.name == parentName;
					});
					var parent = this.planets[parentIndex];
					if (parent !== undefined &&
						parent.parent !== undefined)
					{
						if (parent.parent.match(/Barycenter$/g)) {
							parentName = parent.parent;
						}
					} else {
						console.log("couldn't read parent of ",staticData);
					}
					this.planets.push(mergeInto(new Planet(), {
						id:dynamicData.id,
						name:staticData.name,
						mass:staticData.mass,
						radius:staticData.radius,
						equatorialRadius:staticData.equatorialRadius,
						inverseFlattening:staticData.inverseFlattening,
						magnitude:staticData.magnitude,
						parent:parentName
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

		//calculate total mass of barycenter systems
		for (var i = this.planets.length-1; i >= 0; --i) {
			var planet = this.planets[i];
			planet.isBarycenter = planet.name.match(/Barycenter$/g);
			if (!planet.isBarycenter) continue;
			assert(planet.mass === undefined);
			planet.mass = 0;
			$.each(this.planets, function(j, childPlanet) {
				if (childPlanet.parent != planet.name) return;
				if (childPlanet.mass !== undefined) {
					planet.mass += childPlanet.mass;
				}
			});
		}



		// get texture name 
		for (var i = 0; i < this.planets.length; ++i) { 
			var planet = this.planets[i];
			if (planet.name in {
				Sun:1,
				Mercury:1,
				Venus:1,
				Earth:1,
				Moon:1,
				Mars:1,
					Phobos:1,
					Deimos:1,
				Jupiter:1,
					Io:1,
					Europa:1,
					Ganymede:1,
					Callisto:1,
				Saturn:1,
					Mimas:1,
					Enceladus:1,
					Tethys:1,
					Dione:1,
					Rhea:1,
					Titan:1,
					Iapetus:1,
					Phoebe:1,
				Uranus:1,
				Neptune:1,
				Pluto:1,
					Charon:1
			}) {
				planet.imgURL = 'textures/'+planet.name.toLowerCase()+'.png';
			}
	
		}

		//once we're done making planets, make a copy of the init
		this.initPlanets = this.clonePlanets();
	},

	// interface with smallbodies point cloud picking system
	// for creating/removing Planet objects into this StarSystem:

	//specific to solarsystem / smallbodies:

	removeSmallBody : function(row) {
		var name = row.name;

		//only add if it's already there
		if (solarSystem.indexes[name] === undefined) return;
			
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
	},

	addSmallBody : function(row) {
		var name = row.name;

		//only add if it's not there
		if (solarSystem.indexes[name] !== undefined) {
			return solarSystem.planets[solarSystem.indexes[name]];
		}

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

		planet.initColorSchRadiusAngle();
		planet.initSceneLatLonLineObjs();
		
		// if we pass the row info
		planet.getKOEFromSourceData();
		// if we just pass pos and vel
		//planet.calcKOEFromPosVel();

		planet.updatePosVel();

		return planet;
	}

});
