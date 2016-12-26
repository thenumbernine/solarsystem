var Planet = makeClass({

	init : function(args) {
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
	},

	clone : function() {
		var p = new this.init();
		vec3.copy(p.pos, this.pos);
		vec3.copy(p.vel, this.vel);
		quat.copy(p.angle, this.angle);
		quat.copy(p.tiltAngle, this.tiltAngle);
		return p;
	},

	copy : function(src) {
		this.pos[0] = src.pos[0];
		this.pos[1] = src.pos[1];
		this.pos[2] = src.pos[2];
		this.vel[0] = src.vel[0];
		this.vel[1] = src.vel[1];
		this.vel[2] = src.vel[2];
	},

	//TODO the math ... don't be lazy
	fromGeodeticPosition : function(pos) {
		var x = pos[0];
		var y = pos[1];
		var z = pos[2];
		//if (this.inverseFlattening === undefined) {
		var r2 = Math.sqrt(x*x + y*y);
		var r = Math.sqrt(r2*r2 + z*z);
		var phi = Math.atan2(z, r2);
		var lambda = Math.atan2(y, x);

		var equatorialRadius = this.equatorialRadius !== undefined ? this.equatorialRadius : this.radius;
		if (equatorialRadius === undefined) {
			return {
				lat : 0,
				lon : 0,
				height : 0
			};
		}
		return {
			lat : Math.deg(phi),
			lon : Math.deg(lambda),
			height : r - equatorialRadius
		}
		//} else it is a nonlinear problem and I will think about it later 
	},
	
	//no longer used for mesh construction -- that all goes on in the shader
	//this is only used for getting positions for updating the tidal array calculations
	// and for (geosynchronously) orbiting geodetic locations (which is currently disabled)
	geodeticPosition : function(destX, lat, lon, height) {
		var phi = Math.rad(lat);
		var lambda = Math.rad(lon);
		var cosPhi = Math.cos(phi);
		var sinPhi = Math.sin(phi);

		var equatorialRadius = this.equatorialRadius !== undefined ? this.equatorialRadius : this.radius;
		if (equatorialRadius === undefined) {
			destX[0] = 0;
			destX[1] = 0;
			destX[2] = 0;
		} else if (this.inverseFlattening !== undefined) {
			var eccentricitySquared = (2 * this.inverseFlattening - 1) / (this.inverseFlattening * this.inverseFlattening);
			var sinPhiSquared = sinPhi * sinPhi;
			var N = equatorialRadius / Math.sqrt(1 - eccentricitySquared * sinPhiSquared);
			var NPlusH = N + height;
			destX[0] = NPlusH * cosPhi * Math.cos(lambda);
			destX[1] = NPlusH * cosPhi * Math.sin(lambda);
			destX[2] = (N * (1 - eccentricitySquared) + height) * sinPhi;
		} else {
			var NPlusH = equatorialRadius + height;
			destX[0] = NPlusH * cosPhi * Math.cos(lambda);
			destX[1] = NPlusH * cosPhi * Math.sin(lambda);
			destX[2] = NPlusH * Math.sin(phi);
		}
	},

	geodeticNormal : function(lat, lon) {
		var phi = Math.rad(lat)
		var lambda = Math.rad(lon)
		var cosPhi = Math.cos(phi)
		var sinPhi = Math.sin(phi)

		var equatorialRadius = this.equatorialRadius !== undefined
			? this.equatorialRadius
			: this.radius;
		if (equatorialRadius === undefined) {
			throw "don't know how to calculate this planet's surface "+this;
		}
		if (this.inverseFlattening !== undefined) {
			//mind you this is the planet shape eccentricity, not to be confused with the orbit path eccentricity
			var eccentricitySquared = (2 * this.inverseFlattening - 1) / (this.inverseFlattening * this.inverseFlattening);
			var sinPhiSquared = sinPhi * sinPhi;
			var N = equatorialRadius / Math.sqrt(1 - eccentricitySquared * sinPhiSquared);
			var oneMinusEccSq = 1 - eccentricitySquared;
			var NPlusH = N + height;
			var NPrime = sinPhi * eccentricitySquared / (equatorialRadius * equatorialRadius) * N * N * N;
			var sinPhi = Math.sin(phi);
			var cosPhi = Math.cos(phi);
			var sinLambda = Math.sin(lambda);
			var cosLambda = Math.cos(lambda);
			var nx = oneMinusEccSq * NPlusH * cosPhi * cosLambda * (NPrime * sinPhi + NPlusH * cosPhi);
			var nx = oneMinusEccSq * NPlusH * cosPhi * sinLambda * (NPrime * sinPhi + NPlusH * cosPhi);
			var nz = -NPlusH * cosPhi * (NPrime * cosPhi - NPlusH * sinPhi);
			var l = Math.sqrt(nx*nx + ny*ny + nz*nz);
			return [nx/l, ny/l, nz/l];
		} else {
			var x = cosPhi * Math.cos(lambda);
			var y = cosPhi * Math.sin(lambda);
			var z = Math.sin(phi);
			return [x, y, z]
		}
	}
});
