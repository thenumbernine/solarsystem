import {vec3} from '/js/gl-matrix-3.4.1/index.js';
import {cfg} from './globals.js';
import {ui} from './ui.js';
import {milkyWay} from './milkyway.js';
import {galaxyCenterInEquatorialCoordsInMpc} from './vec.js';
import {metersPerUnits} from './units.js';

//if we're orbiting at 1AU then we can only click things at 1000 AU
const ratioOfOrbitDistanceToAllowSelection = 10000;

class Galaxy {
	constructor(args) {
		for (k in args) {
			this[k] = args[k];
		}
	}
}

class Galaxies {
	constructor() {
		this.closestDistInM = undefined;	//for occlusion of all selection stuff
		this.cache = {};
	
		let thiz = this;
		let xhr = new XMLHttpRequest();
		xhr.open('GET', 'simbad/galaxies.f32', true);
		xhr.responseType = 'arraybuffer';
		xhr.addEventListener('load',  e => {
			thiz.onload(xhr.response);
		});
		xhr.send();
	}

	onload(arrayBuffer) {
		const gl = ui.gl;
		const ModifiedDepthShaderProgram = ui.ModifiedDepthShaderProgram;
		let data = new DataView(arrayBuffer);

		let floatBuffer = new Float32Array(data.byteLength / Float32Array.BYTES_PER_ELEMENT);
		let pos = [];
		for (let j = 0; j < floatBuffer.length; ++j) {
			//units in Mpc
			let x = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);
			if (Math.abs(x) > 1e+26) {
				console.log('galaxy '+Math.floor(j/3)+' has coordinate that position exceeds fp resolution');
			}
			//units still in Mpc
			x += galaxyCenterInEquatorialCoordsInMpc[j%3];
			//units converted to intergalactic render units
			x *= metersPerUnits.Mpc;
			pos[j%3] = x;
			if (j%3 == 2) {
				let len = vec3.length(pos);
				if (this.closestDistInM === undefined) {
					this.closestDistInM = len;
				} else {
					this.closestDistInM = Math.min(len, this.closestDistInM); 
				}
			}
			floatBuffer[j] = x / cfg.interGalacticRenderScale;
		}

		//now that we have the float buffer ...
		this.buffer = new glutil.ArrayBuffer({data : floatBuffer});
		this.sceneObj = new glutil.SceneObject({
			mode : gl.POINTS,
			shader : new ModifiedDepthShaderProgram({
				vertexCode : `
in vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;	// = constant sprite width / screen width, though I have a tapering function that changes size with scale
void main() {
	gl_Position = projMat * (mvMat * flatEarthXForm(vec4(vertex, 1.)));
	gl_PointSize = pointSize / gl_Position.w;
	gl_PointSize = max(1., gl_PointSize);
	gl_Position.z = depthfunction(gl_Position);
}
`,
				fragmentCode : `
uniform sampler2D tex;
out vec4 fragColor;
void main() {
	fragColor = texture(tex, gl_PointCoord);
}
`
			}),
			attrs : {
				vertex : new glutil.Attribute({
					buffer : this.buffer,
					size : 3,
					stride : 3 * Float32Array.BYTES_PER_ELEMENT
				}),
			},
			uniforms : {
				tex : 0
			},
			texs : [],
			useDepth : false,
			blend : [gl.SRC_ALPHA, gl.ONE],
			pos : [0,0,0],
			parent : null,
			static : false
		});
		
		console.log("loaded galaxies");
		
		this.loadNames();
	}

	//this is a 1mb json file ... maybe I should just remotely query it?
	//orrr just do another point octree, and load names with leaf nodes
	loadNames() {
		let thiz = this;
		let url = 'simbad/galaxyNames.json';
		fetch(url)
		.then(response => {
			if (!response.ok) return Promise.reject('not ok');
			response.json()
			.then(names => {
				thiz.names = names;
				console.log("loaded galaxy names");
			});
		}).catch(e => {
			console.log(e);
			console.log('failed to load url '+url+', trying again...');
			thiz.loadNames();
		});
	}

	draw(
		tanFovY,
		picking,
		distFromSolarSystemInMpc
	) {
		//wait for the milky way obj to load and grab its texture
		//TODO work out loading and what a mess it has become
		if (!milkyWay.sceneObj) return;
		if (!cfg.showGalaxies) return;
		if (!this.sceneObj) return;
			
		this.sceneObj.texs[0] = milkyWay.sceneObj.texs[0];	
	
		if (cfg.orbitTarget !== undefined && cfg.orbitTarget.pos !== undefined) {
			this.sceneObj.pos[0] = -cfg.orbitTarget.pos[0] / cfg.interGalacticRenderScale;
			this.sceneObj.pos[1] = -cfg.orbitTarget.pos[1] / cfg.interGalacticRenderScale;
			this.sceneObj.pos[2] = -cfg.orbitTarget.pos[2] / cfg.interGalacticRenderScale;
		}

		const canvas = ui.canvas;
		let pointSize = .02 * Math.sqrt(distFromSolarSystemInMpc) * canvas.width * metersPerUnits.Mpc / cfg.interGalacticRenderScale / tanFovY;
		if (picking) {
			if (cfg.allowSelectGalaxies &&
				this.closestDistInM < ratioOfOrbitDistanceToAllowSelection * cfg.orbitTargetDistance)
			//	&& dist < ratioOfOrbitDistanceToAllowSelection * cfg.orbitTargetDistance
			{
				let thiz = this;
				cfg.pickObject.drawPoints({
					sceneObj : this.sceneObj,
					targetCallback : function(index) {
						return thiz.getCached(index);
					},
					pointSize : pointSize,
					pointSizeScaleWithDist : true,
					//defined in shader:
					pointSizeMin : 1,
					pointSizeMax : Infinity
				});
			}
		} else {
			this.sceneObj.uniforms.pointSize = pointSize;
			this.sceneObj.draw();
		}
	}

	getCached(index) {
		if (this.cache[index]) return this.cache[index];
		let x = this.buffer.data[3*index+0] * cfg.interGalacticRenderScale;
		let y = this.buffer.data[3*index+1] * cfg.interGalacticRenderScale;
		let z = this.buffer.data[3*index+2] * cfg.interGalacticRenderScale;
		let galaxy = new Galaxy({
			name : this.names === undefined ? ('Galaxy #'+index) : this.names[index].id,
			index : index,
			pos : [x,y,z],
			radius : 10 * metersPerUnits.Kpc
		});
		this.cache[index] = galaxy;
		return galaxy;
	}
}

export {Galaxy, Galaxies};
