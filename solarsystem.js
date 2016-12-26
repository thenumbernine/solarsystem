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
		this.planets.push(mergeInto(new Planet(), {
			id : 399,
			name : 'Earth',
			parent : 'Sun',
			mass : 5.9736e+24,
			radius : 6.37101e+6,
			equatorialRadius : 6378.136e+3,
			inverseFlattening : 298.257223563,
			rotationPeriod : (23 + (56 + 4.09053083288 / 60) / 60) / 24,	//sidereal day
			orbitalPeriod : 365.256363004,
			rotationOffset : -4/24 * 2*Math.PI,	//looks about 4 hours off
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {
			id : 301,
			name : 'Moon',
			parent : 'Earth',
			mass : 7.349e+22,
			radius : 1.73753e+6,
			rotationPeriod : 27.321662,
			orbitalPeriod : 27.321662,
			rotationOffset : 180 * Math.PI / 360,
			type : 'planet'
		}));
		this.planets.push(mergeInto(new Planet(), {id:499, name:'Mars', parent:'Sun', mass:6.4185e+23, radius:3.3899e+6, equatorialRadius:3397e+3, inverseFlattening:154.409, rotationPeriod:24.622962/24, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:599, name:'Jupiter', parent:'Sun', mass:1.89813e+27, radius:6.9911e+7, equatorialRadius:71492e+3, inverseFlattening:1/0.06487, ringRadiusRange:[102200000,227000000], type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:699, name:'Saturn', parent:'Sun', mass:5.68319e+26, radius:5.8232e+7, equatorialRadius:60268e+3, inverseFlattening:1/0.09796, ringRadiusRange:[74510000,140390000], type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:799, name:'Uranus', parent:'Sun', mass:8.68103e+25, radius:2.5362e+7, equatorialRadius:25559e+3, inverseFlattening:1/0.02293, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {id:899, name:'Neptune', parent:'Sun', mass:1.0241e+26, radius:2.4624e+7, equatorialRadius:24766e+3, inverseFlattening:1/0.0171, type:'planet'}));
		this.planets.push(mergeInto(new Planet(), {
			id : 999,
			name : 'Pluto',
			parent : 'Sun',
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
				//1-9 are Barycenter points - skip
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
					this.planets.push(mergeInto(new Planet(), {
						id:dynamicData.id,
						name:staticData.name,
						mass:staticData.mass,
						radius:staticData.radius,
						equatorialRadius:staticData.equatorialRadius,
						inverseFlattening:staticData.inverseFlattening,
						magnitude:staticData.magnitude,
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
