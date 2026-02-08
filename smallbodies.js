/*
A few options on how to handle small bodies:

1) downloads everything.  Names, KOE, pos/vel.  I think about 250 MB for double / 125MB for single precision. unreasonable.
	pro: you can do orbit evolution with all objects simultaneously.
	 	do it on GPU = realtime for all 1 million bodies?
	con: size

2) downloads pos/vel only, calculates KOE, and remote queries names.  22MB for float32 data.

3) downloads pos only, queries name and KOE.  11MB.

4) download visible nodes in an octree.
	?MB.  Hopefully less than above to just view a few small bodies.  But to view all it will be more than the flatfile of 250MB.
	Seems most octree nodes were visible, so this was just downloading more data than above.
	Also we have to issue multiple draw commands per octree node, so it will draw slower as well.

If we were going to simulate everything on the GPU, what would we need?
	- orbitType: elliptic vs hyperbolic vs parabolic
	elliptic:
		- orbitalPeriod			\
		- meanAnomalyAtEpoch 	 - these can be combined into meanAnomaly
		- epoch					/
	hyperbolic:
		- meanAnomaly
		- timeOfPeriapsisCrossing
	parabolic:
		haven't done this yet
	- eccentricity
	- A, B
	- sqrt(semiMajorAxis^3 / gravitationalParameter)
	- meanAnomaly is written
	- eccentricAnomaly is written
	- fractionOffset is written ... for the orbit path shader

xyz <= (const) A, B, (dyn) coeffA, coeffB
coeffA, coeffB <= (const) eccentricity, (dyn) pathEccentricAnomaly, dE_dt
dE_dt <= (const) sqrt(semiMajorAxis^3/gravitationalParameter), eccentricity, (dyn) pathEccentricAnomaly
pathEccentricAnomaly <= eccentricAnomaly + meanMotion * (extern) timeAdvanced


TODO don't even download this until requested
and another TODO - organize the external datasets and have them all downloaded on request
 kinda like the universe visualization does
*/
import {vec3, mat4} from '/js/gl-matrix-3.4.1/index.js';
import {ui} from './ui.js';
import {cfg} from './globals.js';

//used with isa in the cfg.orbitTarget detection
//instanciated when the user selects a node in the tree
class SmallBody {
	//callback from setOrbitTarget
	onSelect() {
		//add/removeSmallBody works with the UI controls to add/remove rows
		let planet = solarSystem.addSmallBody(this.row);

		//for the orbit-target popup menu on the right
		planet.type = this.row.bodyType;

		//why doesn't this match up with the point location?  float error?
		cfg.setOrbitTarget(planet);

		//hmm, this updates where the picking region goes
		//but even when mousing over that region, the pick highlight still shows up at the old position
		//let data = this.node.sceneObj.attrs.vertex.buffer.data;
		//data[0+3*this.nodeLocalIndex] = planet.pos[0];
		//data[1+3*this.nodeLocalIndex] = planet.pos[1];
		//data[2+3*this.nodeLocalIndex] = planet.pos[2];
		//instead of matching the pointfield pos to the small body
		// how about hiding the pointfield pos?
		//but then you'd have to restore it once the body was hidden
		//this.node.sceneObj.attrs.vertex.buffer.updateData();
	}
}

class SmallBodies {
	constructor() {
		const ModifiedDepthShaderProgram = ui.ModifiedDepthShaderProgram;
		this.shader = new ModifiedDepthShaderProgram({
			vertexCode : `
in vec3 vertex;
uniform mat4 viewMatInv;
uniform mat4 localMat;
uniform mat4 projMat;
uniform float pointSize;
void main() {
	vec4 vtx4 = vec4(vertex, 1.);
	vtx4 = localMat * vtx4;
	if (flatEarthCoeff != 0.) {
		vtx4 = flatEarthXForm(vtx4);
	}
	gl_Position = projMat * (viewMatInv * vtx4);
	gl_PointSize = pointSize / gl_Position.w;
	gl_PointSize = clamp(gl_PointSize, .25, 5.);
	gl_Position.z = depthfunction(gl_Position);
}
`,
			fragmentCode : `
uniform float alpha;
out vec4 fragColor;
void main() {
	fragColor = vec4(.75, .75, .75, alpha);
}
`,
		});


		let thiz = this;


		let xhr = new XMLHttpRequest();
		//xhr.open('GET', 'jpl-ssd-smallbody/posvel.f64', true);
		xhr.open('GET', 'jpl-ssd-smallbody/pos.f32', true);
		xhr.responseType = 'arraybuffer';
		// if we want a progress bar ...
		let lastProgress = 0;
		xhr.addEventListener('progress', e => {
			if (e.total) {
				let thisProgress = Math.floor(e.loaded / e.total * 10);
				if (thisProgress != lastProgress) {
					lastProgress = thisProgress;
					console.log("smallbody loaded "+(thisProgress*10));
				}
			}
		});

		this.mins = [Infinity, Infinity, Infinity];
		this.maxs = [-Infinity, -Infinity, -Infinity];
		this.added = 0;

		xhr.addEventListener('load', e => {
			let data = new DataView(xhr.response);
			//let len = data.byteLength / Float64Array.BYTES_PER_ELEMENT;
			let len = data.byteLength / Float32Array.BYTES_PER_ELEMENT;

			//units are in parsecs
			//don't forget velocity is not being rescaled (i'm not using it at the moment)
			let floatBuffer = new Float32Array(len);
			for (let j = 0; j < len; ++j) {
				//floatBuffer[j] = data.getFloat64(j * Float64Array.BYTES_PER_ELEMENT, true);
				floatBuffer[j] = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);
/*
				if (j % 6 == 5) {
let n = floatBuffer.length;
let x = floatBuffer[n-5];
let y = floatBuffer[n-4];
let z = floatBuffer[n-3];
if (isFinite(x) && isFinite(y) && isFinite(z)
	&& (x != 0 || y != 0 || z != 0)
) {
	thiz.mins[0] = Math.min(thiz.mins[0], x);
	thiz.maxs[0] = Math.max(thiz.maxs[0], x);
	thiz.mins[1] = Math.min(thiz.mins[1], y);
	thiz.maxs[1] = Math.max(thiz.maxs[1], y);
	thiz.mins[2] = Math.min(thiz.mins[2], z);
	thiz.maxs[2] = Math.max(thiz.maxs[2], z);
	++thiz.added;
//} else {
//	floatBuffer.splice(n-5, 6);
}
				}
*/
			}
			thiz.doneLoading(floatBuffer);
		});
		xhr.send();
	}

	doneLoading(float32data) {
		const gl = ui.gl;
		this.vertexArray = float32data;
		this.vertexBuffer  = new glutil.ArrayBuffer({
			dim : this.numElem,
			data : this.vertexArray
		});

		this.sceneObj = new glutil.SceneObject({
			mode : gl.POINTS,
			attrs : {
				vertex : new glutil.Attribute({
					buffer : this.vertexBuffer,
					size : 3,
					stride : this.numElem * Float32Array.BYTES_PER_ELEMENT,
					offset : 0
				}),
			},
			uniforms : {
				pointSize : this.pointSize,
				alpha : this.pointAlpha,
			},
			shader : this.shader,
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			pos : [0,0,0],
			parent : null,
			static : false
		});
	}

	draw(
		tanFovY,
		picking,
		viewPosInv,
		invRotMat,
		distFromSolarSystemInM
	) {
		let thiz = this;

		if (!cfg.showSmallBodies) return;
		if (picking && !cfg.allowSelectSmallBodies) return;
		if (this.sceneObj === undefined) return;

		this.doDraw(distFromSolarSystemInM, tanFovY, picking, viewPosInv, invRotMat);
	}

	doDraw(distInM, tanFovY, picking, viewPosInv, invRotMat) {
		//no geometry and no buffer if the points buffer is size zero
		if (!this.sceneObj) return;

		vec3.scale(viewPosInv, glutil.view.pos, -1);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);

		const canvas = ui.canvas;
		let pointSize =
			this.pointSize
			* canvas.width
			* Math.sqrt(distInM)
			/ tanFovY;
		this.sceneObj.uniforms.pointSize = pointSize;
		this.sceneObj.uniforms.alpha = this.pointAlpha;
		this.sceneObj.uniforms.julianDate = cfg.julianDate;

		if (cfg.orbitTarget !== undefined && cfg.orbitTarget.pos !== undefined) {
			this.sceneObj.pos[0] = -cfg.orbitTarget.pos[0];
			this.sceneObj.pos[1] = -cfg.orbitTarget.pos[1];
			this.sceneObj.pos[2] = -cfg.orbitTarget.pos[2];
		}

		this.sceneObj.uniforms.flatEarthCoeff = cfg.flatEarthCoeff;
		this.sceneObj.uniforms.earthPos = cfg.flatEarthRelativeEarthPos;
		this.sceneObj.uniforms.earthNorthDir = cfg.flatEarthRelativeEarthNorthDir;

		if (picking) {
/*
			//extra tough too-small threshold for picking
			if (this.visRatio < this.visRatioPickThreshold) return;
			let thiz = this;
			cfg.pickObject.drawPoints({
				sceneObj : this.sceneObj,
				targetCallback : function(i) {
					return thiz.onPick(thiz, i);
				},
				pointSize : pointSize,
				pointSizeScaleWithDist : true,
				//defined in shader
				pointSizeMin : .25,
				pointSizeMax : 5
			});
*/
		} else {
			this.sceneObj.draw({
				uniforms : {
					pointSize : this.pointSize,
					alpha : this.pointAlpha,
				}
			});
		}
	}

	onPick(node, index) {
		let data = this.vertexArray;
		let x = data[0+this.numElem*index];
		let y = data[1+this.numElem*index];
		let z = data[2+this.numElem*index];
		//let vx = data[3+this.numElem*index];
		//let vy = data[4+this.numElem*index];
		//let vz = data[5+this.numElem*index];

		//TODO toggle on/off orbit data if we're selecting on/off a small body
		//TODO even more - don't query this, but instead use the local keplar orbital elements
		//TODO even more - GPU update using keplar orbital elements to small body position
		//TODO even more - correct the eccentric anomaly but and update it every change in julian date
		//TODO even more - unify KOE systems of planets and small bodies
		//TODO even more - unify point cloud octree system of small bodies and starfields

		if (this.cache[index]) {
			return this.cache[index];
		}

		//called from PointOctree.onPick when it is creating a new object
		// either for mouse-over or for click selection
		// and stored in the cache

		//TODO a remote request to get the name

		let smallBody = mergeInto(new SmallBody(), {
			name : 'Small Body '+index,	//node.nameArray[index],
			pos : [x,y,z],
			//vel : [vx,vy,vz],
			radius : 1,
			smallBodyID : index,
		});
		this.cache[index] = smallBody;

		/*
		Now in setOrbitTarget there is special-case code looking for returned instances
		of SmallBody.  From there that code does an ajax query of jpl-ssd-smallbody/search.lua.
		But why bother when we have the info here already?
		*/

/*
		let row = {}
		for (let i = 0; i < this.rows.length; ++i) {
			let field = this.rows[i].name;
			row[field] = node[field+'Array'][nodeLocalIndex];
		}
		row.bodyType = ['comet', 'numbered asteroid', 'unnumbered asteroid'][row.bodyType];
*/
		smallBody.row = {
			name : smallBody.name
		};

		return smallBody;
	}
}

SmallBodies.prototype.pointSize = 1e+11;	//in m ... so maybe convert everything to AU
SmallBodies.prototype.pointAlpha = .75;
SmallBodies.prototype.showAllAtOnce = false;
SmallBodies.prototype.showWithDensity = false;	//don't change this after init
SmallBodies.prototype.visRatioThreshold = .03;
SmallBodies.prototype.visRatioPickThreshold = .1;
SmallBodies.prototype.numElem = 3;
SmallBodies.prototype.show = true;

export {SmallBody, SmallBodies}
