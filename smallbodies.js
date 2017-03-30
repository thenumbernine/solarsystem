var showSmallBodies = true;
var allowSelectSmallBodies = true;

//used with isa in the orbitTarget detection
//instanciated when the user selects a node in the tree 
var SmallBody = makeClass({
	//callback from setOrbitTarget
	onSelect : function() {
		//add/removeSmallBody works with the UI controls to add/remove rows
		var planet = solarSystem.addSmallBody(this.row);
		
		//for the orbit-target popup menu on the right
		planet.type = this.row.bodyType;

		//why doesn't this match up with the point location?  float error?
		setOrbitTarget(planet);

		//hmm, this updates where the picking region goes
		//but even when mousing over that region, the pick highlight still shows up at the old position
		//var data = this.node.sceneObj.attrs.vertex.buffer.data;
		//data[0+3*this.nodeLocalIndex] = planet.pos[0];
		//data[1+3*this.nodeLocalIndex] = planet.pos[1];
		//data[2+3*this.nodeLocalIndex] = planet.pos[2];
		//instead of matching the pointfield pos to the small body
		// how about hiding the pointfield pos?
		//but then you'd have to restore it once the body was hidden
		//this.node.sceneObj.attrs.vertex.buffer.updateData();
	}
});

var SmallBodies = makeClass({
	super : PointOctree,
	urlBase : 'jpl-ssd-smallbody',

	//TODO get this from jpl-ssd-smallbody/parse.lua
	rows : smallBodyRows,
	
	show : true,
	init : function() {
		SmallBodies.superProto.init.apply(this, arguments);
	},
	draw : function(
		tanFovY,
		picking,
		viewPosInv,
		invRotMat,
		distFromSolarSystemInM
	) {
		if (!showSmallBodies) return;
		if (picking && !allowSelectSmallBodies) return;
		SmallBodies.superProto.draw.apply(this, arguments);
	},

	//called from PointOctree.onPick when it is creating a new object
	// either for mouse-over or for click selection
	// and stored in the cache
	createOrbitTarget : function(node, nodeLocalIndex, x,y,z,index) {
		var smallBody = mergeInto(new SmallBody(), {
			name : node.nameArray[nodeLocalIndex],
			pos : [x,y,z],
			radius : 1,
			smallBodyID : index,

			node : node,
			nodeLocalIndex : nodeLocalIndex
		});
		this.cache[index] = smallBody;

		/*
		Now in setOrbitTarget there is special-case code looking for returned instances
		of SmallBody.  From there that code does an ajax query of jpl-ssd-smallbody/search.lua.
		But why bother when we have the info here already?
		*/
		
		var row = {}
		for (var i = 0; i < this.rows.length; ++i) {
			var field = this.rows[i].name;
			row[field] = node[field+'Array'][nodeLocalIndex];
		}
		row.bodyType = ['comet', 'numbered asteroid', 'unnumbered asteroid'][row.bodyType];

		smallBody.row = row;
	
		return smallBody;
	}
});

var smallBodies;
