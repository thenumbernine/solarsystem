var showSmallBodies = true;
var allowSelectSmallBodies = true;

var SmallBody = makeClass({});

if (SHOW_ALL_SMALL_BODIES_AT_ONCE) { 

var smallBodyRootNode;
var allSmallBodyNodes = [];
var smallBodyMaxDrawnNodes = 400;	//there are 2185 for the small bodies
var showAllSmallBodiesAtOnce = false;
var smallBodyPointSize = 500;	//in m ... so maybe convert this to AU
var smallBodyPointAlpha = .75;
var smallBodyNodesDrawnThisFrame = []; 
var smallBodyShader;
var smallBodyOctreeData;

if (SHOW_ALL_SMALL_BODIES_WITH_DENSITY) {
var smallBodyFBOTexWidth = 2048;
var smallBodyFBOTexHeight = 2048;
var smallBodyOverlayLogBase = 1;
var smallBodyFBO;
var smallBodyFBOTex;
var smallBodyOverlayShader;
} // SHOW_ALL_SMALL_BODIES_WITH_DENSITY

var PointOctreeNode = makeClass({
	init : function() {
		this.mins = [];
		this.maxs = [];
		this.center = [];
	},
	loadData : function() {
		if (this.loadingData) return;
		this.loadingData = true;
		var thiz = this;
		var url = 'jpl-ssd-smallbody/nodes/'+this.nodeID+'.json';
		$.ajax({
			url : url,
			dataType : 'json',
			cache : false
		}).error(function() {
			console.log('failed to get small body '+thiz.nodeID+' from '+url);
		}).done(function(data) {
			thiz.processData(data);
		});
	},
	processData : function(data) {
		//hold all parameters for each body
		var len = data.length;
		if (!len) return;
	
		var vertexBuffer = new Float32Array(3*len);
		this.bodyTypeArray = new Uint8Array(len);
		this.idNumberArray = [];
		this.nameArray = [];
		this.semiMajorAxisArray = new Float32Array(len);
		this.longitudeOfAscendingNodeArray = new Float32Array(len);
		this.argumentOfPeriapsisArray = new Float32Array(len);
		this.inclinationArray = new Float32Array(len);
		this.eccentricityArray = new Float32Array(len);
		this.timeOfPerihelionPassageArray = new Float32Array(len);
		var orbitalPeriodArray = new Float32Array(len);
		this.meanAnomalyAtEpochArray = new Float32Array(len);
		this.epochArray = new Float32Array(len);
		this.perihelionDistanceArray = new Float32Array(len);
		this.absoluteMagnitudeArray = new Float32Array(len);
		this.magnitudeSlopeParameterArray = new Float32Array(len);
		this.orbitTypeArray = new Uint8Array(len);
		this.globalIndexArray = new Int32Array(len);
		this.orbitSolutionReferenceArray = [];

		for (var i = 0; i < data.length; ++i) {
			var e = 0;
			vertexBuffer[0+3*i] = data[i][e++];
			vertexBuffer[1+3*i] = data[i][e++];
			vertexBuffer[2+3*i] = data[i][e++];
			this.globalIndexArray[i] = data[i][e++]; 
			this.semiMajorAxisArray[i] = data[i][e++];
			this.longitudeOfAscendingNodeArray[i] = data[i][e++];
			this.argumentOfPeriapsisArray[i] = data[i][e++];
			this.inclinationArray[i] = data[i][e++];
			this.eccentricityArray[i] = data[i][e++];
			this.timeOfPerihelionPassageArray[i] = data[i][e++];
			orbitalPeriodArray[i] = data[i][e++];
			this.meanAnomalyAtEpochArray[i] = data[i][e++];
			this.epochArray[i] = data[i][e++];
			this.perihelionDistanceArray[i] = data[i][e++];
			this.absoluteMagnitudeArray [i] = data[i][e++];
			this.magnitudeSlopeParameterArray[i] = data[i][e++];
			this.bodyTypeArray[i] = data[i][e++];
			this.orbitTypeArray[i] = data[i][e++];
			this.idNumberArray[i] = data[i][e++];
			this.nameArray[i] = data[i][e++];
			this.orbitSolutionReferenceArray[i] = data[i][e++]; 
		}
		
		this.sceneObj = new glutil.SceneObject({
			mode : gl.POINTS,
			attrs : {
				vertex : new glutil.Attribute(new glutil.ArrayBuffer({dim : 3, data : vertexBuffer})),
				/* TODO put these all in *another* geometry object, or texture, and GPU render-to-buffer to update the positions
				semiMajorAxis : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.semiMajorAxisArray})),
				longitudeOfAscendingNode : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.longitudeOfAscendingNodeArray})),
				argumentOfPeriapsis : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.argumentOfPeriapsisArray})),
				inclination : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.inclinationArray})),
				eccentricity : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.eccentricityArray})),
				timeOfPerihelionPassage : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.timeOfPerihelionPassageArray})),
				orbitalPeriod : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : orbitalPeriodArray})),
				meanAnomalyAtEpoch : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.meanAnomalyAtEpochArray})),
				epoch : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.epochArray})),
				orbitType : new glutil.Attribute(new glutil.ArrayBuffer({dim : 1, data : this.orbitTypeArray}))
				*/
			},
			shader : smallBodyShader,
			uniforms : {
				pointSize : smallBodyPointSize,
				alpha : smallBodyPointAlpha
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			parent : null,
			static : true
		});
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
		if (this.visRatio < .03) {//too-small threshold
			//TODO remove any cached geometry
			return;	
		}

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
		smallBodyNodesDrawnThisFrame.push(this);

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
		
		var pointSize = smallBodyPointSize * canvas.width * Math.sqrt(distInM) / tanFovY;
		this.sceneObj.uniforms.pointSize = pointSize;
		this.sceneObj.uniforms.alpha = smallBodyPointAlpha;
		this.sceneObj.uniforms.julianDate = julianDate;
		
		if (picking) {
			if (this.visRatio < .1) return;	//extra tough too-small threshold for picking
			var thiz = this;
			pickObject.drawPoints({
				sceneObj : this.sceneObj,
				targetCallback : function(i) {
					return getCachedSmallBody(thiz, i);
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

function initSmallBodies() {
	
	smallBodyShader = new ModifiedDepthShaderProgram({
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


if (SHOW_ALL_SMALL_BODIES_WITH_DENSITY) {
	smallBodyFBOTex = new glutil.Texture2D({
		internalFormat : gl.RGBA,
		format : gl.RGBA,
		type : gl.FLOAT,
		width : smallBodyFBOTexWidth,
		height : smallBodyFBOTexHeight,
		magFilter : gl.LINEAR,
		minFilter : gl.LINEAR,
		wrap : {
			s : gl.CLAMP_TO_EDGE,
			t : gl.CLAMP_TO_EDGE
		}
	});
	smallBodyFBO = new glutil.Framebuffer();
	gl.bindFramebuffer(gl.FRAMEBUFFER, smallBodyFBO.obj);
	gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, smallBodyFBOTex.obj, 0);
	gl.bindFramebuffer(gl.FRAMEBUFFER, null);
	smallBodyOverlayShader = new glutil.ShaderProgram({
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
} //SHOW_ALL_SMALL_BODIES_WITH_DENSITY


	//holds the data used tree-wide
	smallBodyOctreeData = {};

	//holds the nodes, indexed by geometric position
	smallBodyNodeForID = {};

	var processSmallBodyOctree = function(data) {
		console.log('processing small body points...');
		var startTime = Date.now();
		/*
		holds
			mins
			maxs
			indexes (all valid child nodes)
		*/
		smallBodyOctreeData = data;
		
		var nodeExists = {};
		$.each(smallBodyOctreeData.nodes, function(i,nodeID) {
			nodeExists[nodeID] = true;
		});
	
		//some vertexes are infinite, so I'm just going to fix the octree bounds at Pluto's orbit distance
		var mins = smallBodyOctreeData.mins;
		var maxs = smallBodyOctreeData.maxs;

		smallBodyRootNode = new PointOctreeNode();
		smallBodyRootNode.nodeID = 0;
		smallBodyRootNode.depth = 0;
		smallBodyRootNode.levelID = 0;

		for (var j = 0; j < 3; ++j) {
			smallBodyRootNode.mins[j] = mins[j];
			smallBodyRootNode.maxs[j] = maxs[j];
			smallBodyRootNode.center[j] = .5 * (mins[j] + maxs[j]);
		}
		
		allSmallBodyNodes.push(smallBodyRootNode);
		smallBodyNodeForID[0] = smallBodyRootNode;

		var process;
		process = function(node) {
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
						if (!nodeExists[childNodeID]) continue;
						if (!node.children) node.children = [];
						var child = new PointOctreeNode();

						child.nodeID = childNodeID;
						child.depth = childDepth;
						child.levelID = childLevelID;
						
						node.children[childIndex] = child;
						child.parent = node;
						allSmallBodyNodes.push(child);
						for (var j = 0; j < 3; ++j) {
							child.mins[j] = is[j] ? node.center[j] : node.mins[j];
							child.maxs[j] = is[j] ? node.maxs[j] : node.center[j];
							child.center[j] = .5 * (child.mins[j] + child.maxs[j]);
						}
						
						process(child);
					}
				}
			}
		};

		process(smallBodyRootNode);

		var endTime = Date.now();
		console.log('done processing small bodies ', endTime - startTime, ' ms');
	};
		
	var smallBodyOctreeURL = 'jpl-ssd-smallbody/octree.json';
	var loadSmallBodyOctree;
	loadSmallBodyOctree = function() {
		$.ajax({
			url : smallBodyOctreeURL,
			dataType : 'json',
			cache : false,
			timeout : 30000
		}).error(function() {
			console.log('failed to get small body nodes from '+smallBodyNodeURL+' , trying again...');
			setTimeout(function() {
				loadSmallBodyOctree();
			}, 5000);
		}).done(function(data) {
			processSmallBodyOctree(data);
		});
	};
	loadSmallBodyOctree();

}	//SHOW_ALL_SMALL_BODIES_AT_ONCE

}


var cachedSmallBodies = {};
function getCachedSmallBody(node, localIndex) {
	var data = node.sceneObj.attrs.vertex.buffer.data;
	var x = data[3*localIndex+0];
	var y = data[3*localIndex+1];
	var z = data[3*localIndex+2];
	var globalIndex = node.globalIndexArray[localIndex];

	//TODO toggle on/off orbit data if we're selecting on/off a small body
	//TODO even more - don't query this, but instead use the local keplar orbital elements
	//TODO even more - GPU update using keplar orbital elements to small body position
	//TODO even more - correct the eccentric anomaly but and update it every change in julian date
	//TODO even more - unify KOE systems of planets and small bodies
	//TODO even more - unify point cloud octree system of small bodies and starfields

	if (cachedSmallBodies[globalIndex]) {
		return cachedSmallBodies[globalIndex];
	}
	
	var smallBody = mergeInto(new SmallBody(), {
		name : 'Small Body #'+globalIndex,
		pos : [x,y,z],
		radius : 1,
		smallBodyID : globalIndex
	});
	cachedSmallBodies[globalIndex] = smallBody;

	/*
	Now in setOrbitTarget there is special-case code looking for returned instances
	of SmallBody.  From there that code does an ajax query of jpl-ssd-smallbody/search.lua.
	But why bother when we have the info here already?
	*/
	//Here I'm convering the node body data to match the sql row returned from search.lua 
	var fields = [
		'bodyType',
		'idNumber',
		'name',
		'epoch',
		'perihelionDistance',
		'semiMajorAxis',
		'eccentricity',
		'inclination',
		'argumentOfPeriapsis',
		'longitudeOfAscendingNode',
		'meanAnomalyAtEpoch',
		'absoluteMagnitude',
		'magnitudeSlopeParameter',
		'timeOfPerihelionPassage',
		'orbitSolutionReference',
	];
	
	var row = {}
	for (var i = 0; i < fields.length; ++i) {
		var field = fields[i];
assert(node[field+'Array'], "failed to find "+field+"Array");		
		row[field] = node[field+'Array'][localIndex];
	}
	
	row.bodyType = ['comet', 'numbered asteroid', 'unnumbered asteroid'][row.bodyType];

	smallBody.row = row;

	return smallBody;
}

var smallBodies = new function() {
	this.draw = function(
		tanFovY,
		picking,
		viewPosInv,
		invRotMat,
		distFromSolarSystemInM
	) {
if (SHOW_ALL_SMALL_BODIES_AT_ONCE) {
		vec3.scale(viewPosInv, glutil.view.pos, -1);
		vec3.sub(viewPosInv, viewPosInv, orbitTarget.pos);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);
	
if (!SHOW_ALL_SMALL_BODIES_WITH_DENSITY) {
		//TODO adjust based on LOD node depth
		if (!picking || allowSelectSmallBodies) {
			if (smallBodyRootNode && showSmallBodies) {
				if (showAllSmallBodiesAtOnce) {
					for (var i = 0; i < allSmallBodyNodes.length; ++i) {
						var node = allSmallBodyNodes[i];
						node.draw(distFromSolarSystemInM, tanFovY, picking);
					}
				} else {	//good for selective rendering but bad for all rendering
					var drawList = [];
					smallBodyNodesDrawnThisFrame.length = 0;
					smallBodyRootNode.prepDraw(drawList, tanFovY);
					for (var i = 0; i < smallBodyMaxDrawnNodes && drawList.length > 0; ++i) {
						var node = drawList.splice(drawList.length-1, 1)[0];
						node.drawAndAdd(drawList, tanFovY, distFromSolarSystemInM, picking);
					}
				}
			}
		}
} else { //SHOW_ALL_SMALL_BODIES_WITH_DENSITY
		if (!picking) {
			var viewport = gl.getParameter(gl.VIEWPORT);
			gl.viewport(0, 0, Math.min(canvas.width, smallBodyFBOTexWidth), Math.min(canvas.height, smallBodyFBOTexHeight));
			smallBodyFBO.draw({
				callback : function() {
					gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
					for (var i = 0; i < smallBodyMaxDrawnNodes && i < allSmallBodyNodes.length; ++i) {
						var node = allSmallBodyNodes[i];
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
				shader : smallBodyOverlayShader,
				uniforms : {
					logBase : smallBodyOverlayLogBase,
					texSize : [
						Math.min(canvas.width / smallBodyFBOTexWidth, 1),
						Math.min(canvas.height / smallBodyFBOTexHeight, 1)]
				},
				texs : [smallBodyFBOTex, hsvTex],
				blend : [gl.SRC_ALPHA, gl.ONE]
			});
			gl.enable(gl.DEPTH_TEST);
		}
} //SHOW_ALL_SMALL_BODIES_WITH_DENSITY
		vec3.scale(viewPosInv, glutil.view.pos, -1);
		mat4.translate(glutil.scene.mvMat, invRotMat, viewPosInv);
} //SHOW_ALL_SMALL_BODIES_AT_ONCE
	};
};
