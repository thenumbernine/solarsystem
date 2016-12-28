var showGalaxies = true;
var allowSelectGalaxies = true;

//too small (1e+3) and the depth buffer precision gets destroyed - which is important for color-picking
//too big (1e+20) and the scene gets zfar clipped away
//used by milkyway.js and galaxies.js
var interGalacticRenderScale = 1e+15;

var galaxies = new function() {

	this.closestDistInM = undefined;	//for occlusion of all selection stuff
	
	this.init = function() {
		var thiz = this;
		var xhr = new XMLHttpRequest();
		xhr.open('GET', 'simbad/galaxies.f32', true);
		xhr.responseType = 'arraybuffer';
		xhr.onload = function(e) {
			thiz.onload(this.response);
		};
		xhr.send();
	};

	this.onload = function(arrayBuffer) {
		var data = new DataView(arrayBuffer);

		var floatBuffer = new Float32Array(data.byteLength / Float32Array.BYTES_PER_ELEMENT);
		var pos = [];
		for (var j = 0; j < floatBuffer.length; ++j) {
			//units in Mpc
			var x = data.getFloat32(j * Float32Array.BYTES_PER_ELEMENT, true);
			if (Math.abs(x) > 1e+26) {
				console.log('galaxy '+Math.floor(j/3)+' has coordinate that position exceeds fp resolution');
			}
			//units still in Mpc
			x += galaxyCenterInEquatorialCoordsInMpc[j%3];
			//units converted to intergalactic render units
			x *= metersPerUnits.Mpc;
			pos[j%3] = x;
			if (j%3 == 2) {
				var len = vec3.length(pos);
				if (this.closestDistInM === undefined) {
					this.closestDistInM = len;
				} else {
					this.closestDistInM = Math.min(len, this.closestDistInM); 
				}
			}
			floatBuffer[j] = x / interGalacticRenderScale;
		}

		//now that we have the float buffer ...
		this.buffer = new glutil.ArrayBuffer({data : floatBuffer});
		this.sceneObj = new glutil.SceneObject({
			mode : gl.POINTS,
			shader : new ModifiedDepthShaderProgram({
				vertexCode : mlstr(function(){/*
attribute vec3 vertex;
uniform mat4 mvMat;
uniform mat4 projMat;
uniform float pointSize;	// = constant sprite width / screen width, though I have a tapering function that changes size with scale
void main() {
	gl_Position = projMat * (mvMat * vec4(vertex, 1.));
	gl_PointSize = pointSize / gl_Position.w;
	gl_PointSize = max(1., gl_PointSize);
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
				fragmentCode : mlstr(function(){/*
uniform sampler2D tex;
void main() {
	gl_FragColor = texture2D(tex, gl_PointCoord);
}
*/})
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
	};

	//this is a 1mb json file ... maybe I should just remotely query it?
	//orrr just do another point octree, and load names with leaf nodes
	this.loadNames = function() {
		var thiz = this;
		var url = 'simbad/galaxyNames.json';
		$.ajax({
			url : url,
			dataType : 'json',
			timeout : 10000
		}).error(function() {
			console.log('failed to load url '+url+', trying again...');
			thiz.loadNames();
		}).done(function(names) {
			thiz.names = names;
			console.log("loaded galaxy names");
		});
	};

	this.draw = function(
		tanFovY,
		picking,
		distFromSolarSystemInMpc
	) {
		//wait for the milky way obj to load and grab its texture
		//TODO work out loading and what a mess it has become
		if (!milkyWay.sceneObj) return;
		if (!showGalaxies) return;
		if (!this.sceneObj) return;
			
		this.sceneObj.texs[0] = milkyWay.sceneObj.texs[0];	
	
		if (orbitTarget !== undefined && orbitTarget.pos !== undefined) {
			this.sceneObj.pos[0] = -orbitTarget.pos[0] / interGalacticRenderScale;
			this.sceneObj.pos[1] = -orbitTarget.pos[1] / interGalacticRenderScale;
			this.sceneObj.pos[2] = -orbitTarget.pos[2] / interGalacticRenderScale;
		}

		var pointSize = .02 * Math.sqrt(distFromSolarSystemInMpc) * canvas.width * metersPerUnits.Mpc / interGalacticRenderScale / tanFovY;
		if (picking) {
			if (allowSelectGalaxies &&
				this.closestDistInM < ratioOfOrbitDistanceToAllowSelection * orbitTargetDistance)
			//	&& dist < ratioOfOrbitDistanceToAllowSelection * orbitTargetDistance
			{
				var thiz = this;
				pickObject.drawPoints({
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
	};

	this.cache = {};

	this.getCached = function(index) {
		if (this.cache[index]) return this.cache[index];
		var x = this.buffer.data[3*index+0] * interGalacticRenderScale;
		var y = this.buffer.data[3*index+1] * interGalacticRenderScale;
		var z = this.buffer.data[3*index+2] * interGalacticRenderScale;
		var galaxy = new Galaxy({
			name : this.names === undefined ? ('Galaxy #'+index) : this.names[index].id,
			index : index,
			pos : [x,y,z],
			radius : 10 * metersPerUnits.Kpc
		});
		this.cache[index] = galaxy;
		return galaxy;
	}
};
