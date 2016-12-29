/*
this hold static data used by all the star systems
kind of like starfield
but that's just a point cloud
*/

var starSystems = [];
var starSystemForNames = {};
var solarSystem;	//the one and only.  don't construct until after WebGL init so we can populate our float tex for the planets

var planetSceneObj;
var planetLatLonObj;

//collection of all star systems
//kind of like starfield, but that's a point cloud
var StarSystems = makeClass({

	initSolarSystem : function() {
		solarSystem = new SolarSystem();
		starSystemForNames[solarSystem.name] = solarSystem;
		solarSystem.index = starSystems.length;
		starSystems.push(solarSystem);
		solarSystem.doneBuildingPlanets();
		solarSystem.magnitude = solarSystem.planets[solarSystem.indexes.Sun].magnitude;
	},

	initSolarSystemImgUrls : function() {
		for (var planetIndex_ = 0; planetIndex_ < solarSystem.planets.length; ++planetIndex_) { 
			var planetIndex = planetIndex_;
			var planet = solarSystem.planets[planetIndex];
			planet.initColorSchRadiusAngle();
			planet.initSceneLatLonLineObjs();

			// load texture
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
				Neptune:1,
				Uranus:1,
				Pluto:1,
					Charon:1
			}) {
				planet.imgURL = 'textures/'+planet.name.toLowerCase()+'.png';
			}
		}
	},

	initPlanetSceneObjs : function() {
		console.log('init planet shaders...');
		var quad = [[0,0],[0,1],[1,1],[1,1],[1,0],[0,0]];
		var triIndexArray = [];
		var latLonIndexArray = [];
		var vertexArray = [];

		for (var loni=0; loni <= longitudeDivisions; ++loni) {
			var lon = longitudeMin + loni * longitudeStep;
			for (var lati=0; lati <= latitudeDivisions; ++lati) {
				var lat = latitudeMin + lati * latitudeStep;

				vertexArray.push(lat);
				vertexArray.push(lon);

				if (loni < longitudeDivisions && lati < latitudeDivisions) {
					for (var j = 0; j < quad.length; ++j) {
						var ofs = quad[j];
						var index = (lati + ofs[0]) + (latitudeDivisions + 1) * (loni + ofs[1]);
						triIndexArray.push(index);
					}
					//if we're using 5 div step then every 6 will be 30 degrees
					if ((loni * longitudeStep) % 30 == 0) {
						latLonIndexArray.push(lati + (latitudeDivisions+1) * loni);
						latLonIndexArray.push(lati+1 + (latitudeDivisions+1) * loni);
					}
					if ((lati * latitudeStep) % 30 == 0) {
						latLonIndexArray.push(lati + (latitudeDivisions+1) * loni);
						latLonIndexArray.push(lati + (latitudeDivisions+1) * (loni+1));
					}
				}
			}
		}

		//these could be merged if you don't mind evaluating the cos() and sin() in-shader
		var vertexBuffer = new glutil.ArrayBuffer({dim : 2, data : vertexArray});

		planetSceneObj = new glutil.SceneObject({
			mode : gl.TRIANGLES,
			indexes : new glutil.ElementArrayBuffer({data : triIndexArray}),
			attrs : {
				vertex : vertexBuffer
			},
			uniforms : {
				color : [1,1,1,1],
				pos : [0,0,0],
				angle : [0,0,0,1],
				sunDir : [0,0,0, 0,0,0, 0,0,0, 0,0,0],
				ambient : .3
			},
			texs : [],
			parent : null,
			static : true
		});
		planetLatLonObj = new glutil.SceneObject({
			mode : gl.LINES,
			indexes : new glutil.ElementArrayBuffer({
				data : latLonIndexArray
			}),
			shader : latLonShader,
			attrs : {
				vertex : vertexBuffer
			},
			uniforms : {
				color : [1,1,1,.2]
			},
			blend : [gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA],
			pos : [0,0,0],
			angle : [0,0,0,1],
			parent : null
		});

	},

	planetShadersForNumStars : {},

	//request this per solar system.  rebuild if we need, return from cache if we don't.
	getPlanetShadersForNumberOfStars : function(numberOfStars) {
		if (numberOfStars <= 0) numberOfStars = 1;	//huh, I guess I have a star system with no stars ... "CFBDSIR2149 / CFBDSIR J214947.2-040308.9 / CFBDS J214947-040308"
		var shaders = this.planetShadersForNumStars[numberOfStars];
		if (shaders !== undefined) return shaders;

		shaders = {};
		shaders.colorShader = new ModifiedDepthShaderProgram({
			vertexCode :
'#define NUM_STARS '+numberOfStars+'\n'+
			mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs:
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
varying vec3 lightDir[NUM_STARS];		//light position
varying vec3 normal;		//surface normal

*/}) 
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
			+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
uniform vec4 color;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
uniform float ambient;
void main() {
	float litLum = 0.;
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
			+ mlstr(function(){/*
	float luminance = min(1., litLum);
	gl_FragColor = color * max(ambient, sqrt(luminance));
}
*/}),
			uniforms : {
				color : [1,1,1,1],
				pointSize : 4
			}
		});

		shaders.texShader = new ModifiedDepthShaderProgram({
			vertexCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;

*/}) 
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	vec3 modelVertex = geodeticPosition(vertex) * scaleExaggeration;
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 vtx3 = quatRotate(angle, modelVertex) + pos;
	normal = quatRotate(angle, normalize(modelVertex));
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - vtx3);')
			+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
uniform sampler2D tex;
uniform float ambient;
void main() {
	float litLum = 0.;
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
			+ mlstr(function(){/*
	float luminance = min(1., litLum);
	gl_FragColor = texture2D(tex, texCoordv) * max(ambient, sqrt(luminance));
}
*/}),
			uniforms : {
				tex : 0
			}
		});

		shaders.ringShadowShader = new ModifiedDepthShaderProgram({
			vertexCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ mlstr(function(){/*
attribute vec2 vertex;		//lat/lon pairs
uniform mat4 mvMat;			//modelview matrix
uniform mat4 projMat;		//projection matrix
uniform vec3 pos;			//offset to planet position
uniform vec4 angle;			//planet angle
uniform vec3 sunDir[NUM_STARS];		//sun pos, for lighting calculations
uniform float equatorialRadius;		//or use planet radius
uniform float inverseFlattening;	//default 1 if it does not exist
uniform float scaleExaggeration;	//exhaggerate planet sizes
//to fragment shader:
varying vec3 modelVertexv;
varying vec3 normal;
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];

*/}) 	
			+ geodeticPositionCode 
			+ quatRotateCode 
			+ mlstr(function(){/*

void main() {
	//vertex is really the lat/lon in degrees
	modelVertexv = geodeticPosition(vertex) * scaleExaggeration;
	texCoordv = vertex.yx / vec2(360., 180.) + vec2(.5, .5);
	vec3 worldVertex = quatRotate(angle, modelVertexv) + pos;
	normal = quatRotate(angle, normalize(modelVertexv));
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'lightDir[i] = normalize(sunDir[i] - worldVertex);')
			+ mlstr(function(){/*
	vec4 vtx4 = mvMat * vec4(worldVertex, 1.);
	gl_Position = projMat * vtx4;
	gl_Position.z = depthfunction(gl_Position);
}
*/}),
			fragmentCode :
'#define NUM_STARS '+numberOfStars+'\n'
			+ quatRotateCode
			+ mlstr(function(){/*
varying vec2 texCoordv;
varying vec3 lightDir[NUM_STARS];
varying vec3 normal;
varying vec3 modelVertexv;
uniform sampler2D tex;
uniform sampler2D ringTransparencyTex;
uniform float ambient;
uniform float ringMinRadius;
uniform float ringMaxRadius;
uniform vec4 angle;

float ringIntersect(vec3 startPos, vec3 dir) {
	if (dot(startPos, dir) < 0.) return -1.;	//occluded by planet
	float t = -startPos.z / dir.z;
	if (t < 0.) return -1.;	//trace intersects backwards
	vec2 intersect = startPos.xy + t * dir.xy;
	float r = length(intersect);
	return (r - ringMinRadius) / (ringMaxRadius - ringMinRadius);
}

void main() {
	float luminance = 1.;

	//inverse rotate lightDir[0]
	//TODO how to combine this per-light-source ... depenends on the intensity of each ...
	vec3 lightDirInModelSpace = quatRotate(vec4(-angle.xyz, angle.w), lightDir[0]);
	float intersectPos = ringIntersect(modelVertexv, lightDirInModelSpace);
	if (intersectPos >= 0. && intersectPos <= 1.) {
		luminance *= texture2D(ringTransparencyTex, vec2(intersectPos, .5)).r;
	}

	float litLum = 0.;
*/}) 
			+ unravelForLoop('i', 0, numberOfStars-1, 'litLum += max(0., dot(lightDir[i], normal));')
			+ mlstr(function(){/*
	luminance *= min(1., litLum);

	gl_FragColor.rgb = texture2D(tex, texCoordv).rgb * max(ambient, sqrt(luminance));
	gl_FragColor.a = 1.;
}
*/}),
			uniforms : {
				tex : 0,
				ringTransparencyTex : 1
			}
		});

		this.planetShadersForNumStars[numberOfStars] = shaders;
		return shaders;
	}
});

var starSystemsExtra = new StarSystems(); 
