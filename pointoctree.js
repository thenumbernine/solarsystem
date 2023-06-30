import {vec3, mat4} from '/js/gl-matrix-3.4.1/index.js';
import {cfg} from './globals.js';

class PointOctreeNode {
	pointSize : 500,	//in m ... so maybe convert this to AU
	pointAlpha : .75,
	constructor(tree) {
		this.tree = tree;
		this.mins = [];
		this.maxs = [];
		this.center = [];
	}
	loadData() {
		if (this.sceneObj) return;		
		if (this.loadingData) return;
		this.unloaded = false;
		this.loadingData = true;
		let thiz = this;
		let url = this.tree.urlBase+'/nodes/'+this.nodeID+'.json';
		fetch(url)
		.then(response => {
			if (!response.ok) return Promise.reject('not ok');
			response.json()
			.then(data => {
				//if (thiz.unloaded) return;
				thiz.processData(data);
			})
		}).catch(e => {
			console.log('failed to get node '+thiz.nodeID+' from '+url);
		});
	}
	processData(data) {
		assert(!this.sceneObj);	//TODO use this.sceneObj instead of this.loaded

		//hold all parameters for each body
		let len = data.length;
		if (!len) return;

		let rows = this.tree.rows;
		for (let i = 0; i < rows.length; ++i) {
			let buffer;
			switch (rows[i].type) {
			case 'vec3': buffer = new Float32Array(3*len); break;
			case 'byte': buffer = new Uint8Array(len); break;
			case 'int': buffer = new Int32Array(len); break;
			case 'float': buffer = new Float32Array(len); break;
			case 'string': buffer = []; break;
			default:
				throw "don't know how to create array for type "+type;
			}
			this[rows[i].name+'Array'] = buffer;
		}

		for (let i = 0; i < data.length; ++i) {
			let e = 0;
			for (let j = 0; j < rows.length; ++j) {
				let buffer = this[rows[j].name+'Array'];
				switch (rows[j].type) {
				case 'vec3':
					for (let k = 0; k < 3; ++k) {
						buffer[k+3*i] = data[i][e++];
					}
					break;
				default:
					buffer[i] = data[i][e++];
					break;
				}
			}
		}
		
		// TODO put all the other attr arrays in *another* geometry object, or texture, 
		// and GPU render-to-buffer to update the positions
		this.sceneObj = new glutil.SceneObject({
			mode : gl.POINTS,
			attrs : {
				vertex : new glutil.Attribute(new glutil.ArrayBuffer({dim : 3, data : this.posArray})),
			},
			shader : this.tree.shader,
			uniforms : {
				pointSize : this.pointSize,
				alpha : this.pointAlpha
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			parent : null,
			static : true
		});
	}
	unloadData() {
//		if (!this.sceneObj) return;
		this.unloaded = true;	//tell any loading operations to stop
/*		this.loadedData = undefined;
		this.loadingData = undefined;

		if (this.sceneObj.attrs !== undefined) {
			this.sceneObj.attrs.forEach((attr,k) => {
				gl.deleteBuffer(attr.buffer.obj);
			});
		}
		delete this.sceneObj;
		this.sceneObj = undefined;
*/
	}
	prepDraw(drawList, tanFovY) {
		let radius = Math.max(
			this.maxs[0] - this.mins[0],
			this.maxs[1] - this.mins[1],
			this.maxs[2] - this.mins[2]);
		//from center 
		let dx = this.center[0] - glutil.view.pos[0] - cfg.orbitTarget.pos[0];
		let dy = this.center[1] - glutil.view.pos[1] - cfg.orbitTarget.pos[1];
		let dz = this.center[2] - glutil.view.pos[2] - cfg.orbitTarget.pos[2];
		let dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
		
		//TODO test occlusion
		
		this.visRatio = radius / (dist * tanFovY);
		if (this.visRatio < this.tree.visRatioThreshold 
|| this.unloaded		
		) {//too-small threshold
			//TODO remove any cached geometry
			this.unloadData();
			return;	
		}

		//TODO for all-at-once downloading, just look things up in the big binary structure
		//load data if needed
		this.loadData();

		//insert sorted by visRatio, lowest to highest, starting at the back (optimistic)
		for (let i = drawList.length-1; i >= 0; --i) {
			if (drawList[i].visRatio < this.visRatio) {
				drawList.splice(i+1, 0, this);
				return;
			}
		}
		drawList.splice(0, 0, this);
	}
	drawAndAdd(drawList, tanFovY, distInM, picking) {
		this.draw(distInM, tanFovY, picking);
		this.tree.drawnThisFrame.push(this);

		if (this.children !== undefined) {
			for (let i = 0; i < 8; ++i) {
				let ch = this.children[i];
				if (ch !== undefined) {
					ch.prepDraw(drawList, tanFovY);
				}
			}
		}
	}
	draw(distInM, tanFovY, picking) {
		//no geometry and no buffer if the points buffer is size zero
		if (!this.sceneObj) return;
		
		let pointSize = 
			this.pointSize 
			* canvas.width 
			* Math.sqrt(distInM) 
			/ tanFovY;
		this.sceneObj.uniforms.pointSize = pointSize;
		this.sceneObj.uniforms.alpha = this.pointAlpha;
		this.sceneObj.uniforms.julianDate = cfg.julianDate;
		
		if (picking) {
			//extra tough too-small threshold for picking
			if (this.visRatio < this.tree.visRatioPickThreshold) return;	
			let thiz = this;
			cfg.pickObject.drawPoints({
				sceneObj : this.sceneObj,
				targetCallback : function(i) {
					return thiz.tree.onPick(thiz, i);
				},
				pointSize : pointSize,
				pointSizeScaleWithDist : true,
				//defined in shader
				pointSizeMin : .25,
				pointSizeMax : 5
			});
		} else {
			this.sceneObj.draw();
		}
	}
	find(x,y,z) {
		if (this.children !== undefined) {
			for (let i = 0; i < 8; ++i) {
				let ch = this.children[i];
				if (ch === undefined) continue;
				if (x >= ch.mins[0] && x <= ch.maxs[0] &&
					y >= ch.mins[1] && y <= ch.maxs[1] &&
					z >= ch.mins[2] && z <= ch.maxs[2])
				{
					return ch.find(x,y,z);
				}
			}
		}
		return this;
	}
}

class PointOctree {
	maxDrawnNodes : 400,	//there are 2185 for the small bodies
	showAllAtOnce : false,
	showWithDensity : false,	//don't change this after init
	visRatioThreshold : .03,
	visRatioPickThreshold : .1,
	init : function() {
		const ModifiedDepthShaderProgram = ui.ModifiedDepthShaderProgram;
		assertExists(this, 'urlBase');
		assertExists(this, 'rows');
		//I could just have PointOctree provide these two fields:
		assertEquals(this.rows[0].name, 'pos');
		assertEquals(this.rows[0].type, 'vec3');
		assertEquals(this.rows[1].name, 'index');
		assertEquals(this.rows[1].type, 'int');
		
		
		this.allNodes = [];
		this.drawProcessing = [];
		this.drawnThisFrame = [];
		this.cache = {};
		
		this.shader = new ModifiedDepthShaderProgram({
			vertexCode : `
in vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;
void main() {
	gl_Position = projMat * (mvMat * vec4(vertex, 1.));
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


		if (this.showWithDensity) {
			this.overlayLogBase = 1;
			this.fboTexWidth = 2048;
			this.fboTexHeight = 2048;
			this.fboTex = new glutil.Texture2D({
				internalFormat : gl.RGBA,
				format : gl.RGBA,
				type : gl.FLOAT,
				width : this.fboTexWidth,
				height : this.fboTexHeight,
				magFilter : gl.LINEAR,
				minFilter : gl.LINEAR,
				wrap : {
					s : gl.CLAMP_TO_EDGE,
					t : gl.CLAMP_TO_EDGE
				}
			});
			this.fbo = new glutil.Framebuffer();
			gl.bindFramebuffer(gl.FRAMEBUFFER, this.fbo.obj);
			gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.fboTex.obj, 0);
			gl.bindFramebuffer(gl.FRAMEBUFFER, null);
			this.overlayShader = new glutil.ShaderProgram({
				vertexCode : `
in vec2 vertex;
out vec2 tc;
void main() {
	tc = vertex;
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
`,
				fragmentCode : `
uniform sampler2D tex;
uniform sampler2D hsvTex;
uniform float logBase;
uniform vec2 texSize;
in vec2 tc;
out vec4 fragColor;
void main() {
	float alpha = texture(tex, tc * texSize).x;
	alpha = log(alpha + 1.) / logBase;
	fragColor = texture(hsvTex, vec2(alpha, .5));
	fragColor.a = min(alpha * 100., 1.);
}
`,
				uniforms : {
					tex : 0,
					hsvTex : 1,
					logBase : 1,
					texSize : [1,1]
				}
			});
		} //this.showWithDensity
		
		this.load();
	},
	load : function() {
		let url = this.urlBase + '/octree.json';
		let thiz = this;
		fetch(url)
		.then(response => {
			if (!response.ok) return Promise.reject('not ok');
			response.json()
			.then(data => {
				thiz.processData(data);
			});
		}).catch(e => {
			console.log('failed to get octree info from '+url+' , trying again...');
			setTimeout(function() {
				thiz.load(url);
			}, 5000);
		});
	}
	
	/*
	data holds
		mins
		maxs
		indexes (all valid child nodes)
	*/
	processData(data) {
		this.data = data;
		this.nodeIDSet = {};
		for (let i = 0; i < data.nodes.length; ++i) {
			this.nodeIDSet[data.nodes[i]] = true;
		}

		this.root = new PointOctreeNode(this);
		this.root.nodeID = 0;
		this.root.depth = 0;
		this.root.levelID = 0;

		for (let j = 0; j < 3; ++j) {
			this.root.mins[j] = data.mins[j];
			this.root.maxs[j] = data.maxs[j];
			this.root.center[j] = .5 * (data.mins[j] + data.maxs[j]);
		}
		
		this.allNodes.push(this.root);

		this.processNode(this.root);
	}
	
	processNode(node) {
		//node.nodeID includes offsets into each level

		for (let ix = 0; ix < 2; ++ix) {
			for (let iy = 0; iy < 2; ++iy) {
				for (let iz = 0; iz < 2; ++iz) {
					let is = [ix,iy,iz];
					let childIndex = ix | ((iy | (iz<<1)) << 1);

					let childDepth = node.depth + 1;
					let childLevelStart = ((1 << (3*childDepth)) - 1) / 7;
					let childLevelID = node.levelID | (childIndex << 3*node.depth);
					let childNodeID = childLevelStart + childLevelID;
					//if we find it in the master list 
					// then create the node
					if (!this.nodeIDSet[childNodeID]) continue;
					if (!node.children) node.children = [];
					let child = new PointOctreeNode(this);

					child.nodeID = childNodeID;
					child.depth = childDepth;
					child.levelID = childLevelID;
					
					node.children[childIndex] = child;
					child.parent = node;
					this.allNodes.push(child);
					for (let j = 0; j < 3; ++j) {
						child.mins[j] = is[j] ? node.center[j] : node.mins[j];
						child.maxs[j] = is[j] ? node.maxs[j] : node.center[j];
						child.center[j] = .5 * (child.mins[j] + child.maxs[j]);
					}
					
					this.processNode(child);
				}
			}
		}
	}

	draw(
		tanFovY,
		picking,
		viewPosInv,
		invRotMat,
		distFromSolarSystemInM
	) {
		let thiz = this;

		vec3.scale(viewPosInv, glutil.view.pos, -1);
		vec3.sub(viewPosInv, viewPosInv, cfg.orbitTarget.pos);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);
	
		if (!this.showWithDensity) {
			//TODO adjust based on LOD node depth
			if (this.root) {
				if (this.showAllAtOnce) {
					for (let i = 0; i < this.allNodes.length; ++i) {
						let node = this.allNodes[i];
						node.draw(distFromSolarSystemInM, tanFovY, picking);
					}
				} else {	//good for selective rendering but bad for all rendering
					let drawList = this.drawProcessing;
					drawList.length = 0;
					this.drawnThisFrame.length = 0;
					this.root.prepDraw(drawList, tanFovY);
					for (let i = 0; i < this.maxDrawnNodes && drawList.length > 0; ++i) {
						let node = drawList.splice(drawList.length-1, 1)[0];
						node.drawAndAdd(drawList, tanFovY, distFromSolarSystemInM, picking);
					}
				}
			}
		} else { //this.showWithDensity
			if (!picking) {
				let viewport = gl.getParameter(gl.VIEWPORT);
				gl.viewport(0, 0, Math.min(canvas.width, this.fboTexWidth), Math.min(canvas.height, this.fboTexHeight));
				this.fbo.draw({
					callback : function() {
						gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
						for (let i = 0; i < thiz.maxDrawnNodes && i < thiz.allNodes.length; ++i) {
							let node = thiz.allNodes[i];
							node.draw(distFromSolarSystemInM, tanFovY, picking);
						}
					}
				});
				gl.viewport.apply(gl, viewport);
				gl.disable(gl.DEPTH_TEST);
				hsvTex
					.bind()
					.setWrap({s : gl.REPEAT, t : gl.REPEAT})
					.unbind();
				glutil.unitQuad.draw({
					shader : this.overlayShader,
					uniforms : {
						logBase : this.overlayLogBase,
						texSize : [
							Math.min(canvas.width / this.fboTexWidth, 1),
							Math.min(canvas.height / this.fboTexHeight, 1)]
					},
					texs : [this.fboTex, hsvTex],
					blend : [gl.SRC_ALPHA, gl.ONE]
				});
				gl.enable(gl.DEPTH_TEST);
			}
		} //this.showWithDensity
		
		vec3.scale(viewPosInv, glutil.view.pos, -1);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);
	}

	onPick(node, nodeLocalIndex) {
		let data = node.sceneObj.attrs.vertex.buffer.data;
		let x = data[0+3*nodeLocalIndex];
		let y = data[1+3*nodeLocalIndex];
		let z = data[2+3*nodeLocalIndex];
		let index = node.indexArray[nodeLocalIndex];

		//TODO toggle on/off orbit data if we're selecting on/off a small body
		//TODO even more - don't query this, but instead use the local keplar orbital elements
		//TODO even more - GPU update using keplar orbital elements to small body position
		//TODO even more - correct the eccentric anomaly but and update it every change in julian date
		//TODO even more - unify KOE systems of planets and small bodies
		//TODO even more - unify point cloud octree system of small bodies and starfields

		if (this.cache[index]) {
			return this.cache[index];
		}
	
		return this.createOrbitTarget(node,nodeLocalIndex,x,y,z,index);
	}
}
