var PointOctreeNode = makeClass({
	pointSize : 500,	//in m ... so maybe convert this to AU
	pointAlpha : .75,
	init : function(tree) {
		this.tree = tree;
		this.mins = [];
		this.maxs = [];
		this.center = [];
	},
	loadData : function() {
		if (this.sceneObj) return;		
		if (this.loadingData) return;
		this.unloaded = false;
		this.loadingData = true;
		var thiz = this;
		var url = this.tree.urlBase+'/nodes/'+this.nodeID+'.json';
		$.ajax({
			url : url,
			dataType : 'json',
		}).error(function() {
			console.log('failed to get node '+thiz.nodeID+' from '+url);
		}).done(function(data) {
//			if (thiz.unloaded) return;
			thiz.processData(data);
		});
	},
	processData : function(data) {
		assert(!this.sceneObj);	//TODO use this.sceneObj instead of this.loaded

		//hold all parameters for each body
		var len = data.length;
		if (!len) return;

		var rows = this.tree.rows;
		for (var i = 0; i < rows.length; ++i) {
			var buffer;
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

		for (var i = 0; i < data.length; ++i) {
			var e = 0;
			for (var j = 0; j < rows.length; ++j) {
				var buffer = this[rows[j].name+'Array'];
				switch (rows[j].type) {
				case 'vec3':
					for (var k = 0; k < 3; ++k) {
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
	},
	unloadData : function() {
//		if (!this.sceneObj) return;
		this.unloaded = true;	//tell any loading operations to stop
/*		this.loadedData = undefined;
		this.loadingData = undefined;

		if (this.sceneObj.attrs !== undefined) {
			$.each(this.sceneObj.attrs, function(k,attr) {
				gl.deleteBuffer(attr.buffer.obj);
			});
		}
		delete this.sceneObj;
		this.sceneObj = undefined;
*/
	},
	prepDraw : function(drawList, tanFovY) {
		var radius = Math.max(
			this.maxs[0] - this.mins[0],
			this.maxs[1] - this.mins[1],
			this.maxs[2] - this.mins[2]);
		//from center 
		var dx = this.center[0] - glutil.view.pos[0] - orbitTarget.pos[0];
		var dy = this.center[1] - glutil.view.pos[1] - orbitTarget.pos[1];
		var dz = this.center[2] - glutil.view.pos[2] - orbitTarget.pos[2];
		var dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
		
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
		for (var i = drawList.length-1; i >= 0; --i) {
			if (drawList[i].visRatio < this.visRatio) {
				drawList.splice(i+1, 0, this);
				return;
			}
		}
		drawList.splice(0, 0, this);
	},
	drawAndAdd : function(drawList, tanFovY, distInM, picking) {
		this.draw(distInM, tanFovY, picking);
		this.tree.drawnThisFrame.push(this);

		if (this.children !== undefined) {
			for (var i = 0; i < 8; ++i) {
				var ch = this.children[i];
				if (ch !== undefined) {
					ch.prepDraw(drawList, tanFovY);
				}
			}
		}
	},
	draw : function(distInM, tanFovY, picking) {
		//no geometry and no buffer if the points buffer is size zero
		if (!this.sceneObj) return;
		
		var pointSize = 
			this.pointSize 
			* canvas.width 
			* Math.sqrt(distInM) 
			/ tanFovY;
		this.sceneObj.uniforms.pointSize = pointSize;
		this.sceneObj.uniforms.alpha = this.pointAlpha;
		this.sceneObj.uniforms.julianDate = julianDate;
		
		if (picking) {
			//extra tough too-small threshold for picking
			if (this.visRatio < this.tree.visRatioPickThreshold) return;	
			var thiz = this;
			pickObject.drawPoints({
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
	},
	find : function(x,y,z) {
		if (this.children !== undefined) {
			for (var i = 0; i < 8; ++i) {
				var ch = this.children[i];
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
});

var PointOctree = makeClass({
	maxDrawnNodes : 400,	//there are 2185 for the small bodies
	showAllAtOnce : false,
	showWithDensity : false,	//don't change this after init
	visRatioThreshold : .03,
	visRatioPickThreshold : .1,
	init : function() {
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
			vertexCode : mlstr(function(){/*
attribute vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;
void main() {
	gl_Position = projMat * (mvMat * vec4(vertex, 1.));
	gl_PointSize = pointSize / gl_Position.w;
	gl_PointSize = clamp(gl_PointSize, .25, 5.);
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode : mlstr(function(){/*
uniform float alpha;
void main() {
	gl_FragColor = vec4(.75, .75, .75, alpha);
}
*/}),
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
				vertexPrecision : 'best',
				vertexCode : mlstr(function(){/*
attribute vec2 vertex;
varying vec2 tc;
void main() {
	tc = vertex;
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/}),
				fragmentPrecision : 'best',
				fragmentCode : mlstr(function(){/*
uniform sampler2D tex;
uniform sampler2D hsvTex;
uniform float logBase;
uniform vec2 texSize;
varying vec2 tc;
void main() {
	float alpha = texture2D(tex, tc * texSize).x;
	alpha = log(alpha + 1.) / logBase;
	gl_FragColor = texture2D(hsvTex, vec2(alpha, .5));
	gl_FragColor.a = min(alpha * 100., 1.);
}
*/}),
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
		var url = this.urlBase + '/octree.json';
		var thiz = this;
		$.ajax({
			url : url,
			dataType : 'json'
		}).error(function() {
			console.log('failed to get octree info from '+url+' , trying again...');
			setTimeout(function() {
				thiz.load(url);
			}, 5000);
		}).done(function(data) {
			thiz.processData(data);
		});
	},
	
	/*
	data holds
		mins
		maxs
		indexes (all valid child nodes)
	*/
	processData : function(data) {
		this.data = data;
		this.nodeIDSet = {};
		for (var i = 0; i < data.nodes.length; ++i) {
			this.nodeIDSet[data.nodes[i]] = true;
		}

		this.root = new PointOctreeNode(this);
		this.root.nodeID = 0;
		this.root.depth = 0;
		this.root.levelID = 0;

		for (var j = 0; j < 3; ++j) {
			this.root.mins[j] = data.mins[j];
			this.root.maxs[j] = data.maxs[j];
			this.root.center[j] = .5 * (data.mins[j] + data.maxs[j]);
		}
		
		this.allNodes.push(this.root);

		this.processNode(this.root);
	},
	
	processNode : function(node) {
		//node.nodeID includes offsets into each level

		for (var ix = 0; ix < 2; ++ix) {
			for (var iy = 0; iy < 2; ++iy) {
				for (var iz = 0; iz < 2; ++iz) {
					var is = [ix,iy,iz];
					var childIndex = ix | ((iy | (iz<<1)) << 1);

					var childDepth = node.depth + 1;
					var childLevelStart = ((1 << (3*childDepth)) - 1) / 7;
					var childLevelID = node.levelID | (childIndex << 3*node.depth);
					var childNodeID = childLevelStart + childLevelID;
					//if we find it in the master list 
					// then create the node
					if (!this.nodeIDSet[childNodeID]) continue;
					if (!node.children) node.children = [];
					var child = new PointOctreeNode(this);

					child.nodeID = childNodeID;
					child.depth = childDepth;
					child.levelID = childLevelID;
					
					node.children[childIndex] = child;
					child.parent = node;
					this.allNodes.push(child);
					for (var j = 0; j < 3; ++j) {
						child.mins[j] = is[j] ? node.center[j] : node.mins[j];
						child.maxs[j] = is[j] ? node.maxs[j] : node.center[j];
						child.center[j] = .5 * (child.mins[j] + child.maxs[j]);
					}
					
					this.processNode(child);
				}
			}
		}
	},

	draw : function(
		tanFovY,
		picking,
		viewPosInv,
		invRotMat,
		distFromSolarSystemInM
	) {
		var thiz = this;

		vec3.scale(viewPosInv, glutil.view.pos, -1);
		vec3.sub(viewPosInv, viewPosInv, orbitTarget.pos);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);
	
		if (!this.showWithDensity) {
			//TODO adjust based on LOD node depth
			if (this.root) {
				if (this.showAllAtOnce) {
					for (var i = 0; i < this.allNodes.length; ++i) {
						var node = this.allNodes[i];
						node.draw(distFromSolarSystemInM, tanFovY, picking);
					}
				} else {	//good for selective rendering but bad for all rendering
					var drawList = this.drawProcessing;
					drawList.length = 0;
					this.drawnThisFrame.length = 0;
					this.root.prepDraw(drawList, tanFovY);
					for (var i = 0; i < this.maxDrawnNodes && drawList.length > 0; ++i) {
						var node = drawList.splice(drawList.length-1, 1)[0];
						node.drawAndAdd(drawList, tanFovY, distFromSolarSystemInM, picking);
					}
				}
			}
		} else { //this.showWithDensity
			if (!picking) {
				var viewport = gl.getParameter(gl.VIEWPORT);
				gl.viewport(0, 0, Math.min(canvas.width, this.fboTexWidth), Math.min(canvas.height, this.fboTexHeight));
				this.fbo.draw({
					callback : function() {
						gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
						for (var i = 0; i < thiz.maxDrawnNodes && i < thiz.allNodes.length; ++i) {
							var node = thiz.allNodes[i];
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
	},

	onPick : function(node, nodeLocalIndex) {
		var data = node.sceneObj.attrs.vertex.buffer.data;
		var x = data[0+3*nodeLocalIndex];
		var y = data[1+3*nodeLocalIndex];
		var z = data[2+3*nodeLocalIndex];
		var index = node.indexArray[nodeLocalIndex];

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
});
