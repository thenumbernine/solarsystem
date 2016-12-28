var showSmallBodies = true;
var allowSelectSmallBodies = true;

//used with isa in the orbitTarget detection
//instanciated when the user selects a node in the tree 
var SmallBody = makeClass({
	//callback from setOrbitTarget
	onSelect : function() {
		var row = newTarget.row;
	
		//add/removeSmallBody works with the UI controls to add/remove rows
		var planet = addSmallBody(row);
		
		//TODO why doesn't this match up with the point location?  float error?
		setOrbitTarget(planet);
	}
});

var SmallBodies = makeClass({
	super : PointOctree,
	urlBase : 'jpl-ssd-smallbody',

	//TODO get this from jpl-ssd-smallbody/parse.lua
	rows : [
	//these need to be here for the PointOctree to work:
		{name:'vertex', type:'vec3'},
		{name:'globalIndex', type:'int'},	//index into the dense list of this point cloud. 
	//these are specific to SmallBodies
		{name:'semiMajorAxis', type:'float'},
		{name:'longitudeOfAscendingNode', type:'float'},
		{name:'argumentOfPeriapsis', type:'float'},
		{name:'inclination', type:'float'},
		{name:'eccentricity', type:'float'},
		{name:'timeOfPerihelionPassage', type:'float'},
		{name:'orbitalPeriod', type:'float'},
		{name:'meanAnomalyAtEpoch', type:'float'},
		{name:'epoch', type:'float'},
		{name:'perihelionDistance', type:'float'},
		{name:'absoluteMagnitudeArray', type:'float'},
		{name:'magnitudeSlopeParameter', type:'float'},
		{name:'bodyType', type:'byte'},
		{name:'orbitType', type:'byte'},
		{name:'idNumber', type:'string'},	//up to 6 bytes.  number in the horizons system.  sometimes a number, sometimes a letter, sometimes both
		{name:'name', type:'string'},	//up to 38 bytes
		{name:'orbitSolutionReference', type:'string'}	//up to 10 bytes
	],
	
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
	createOrbitTarget : function(node, nodeLocalIndex, x,y,z,globalIndex) {
		var smallBody = mergeInto(new SmallBody(), {
			name : node.nameArray[nodeLocalIndex],
			pos : [x,y,z],
			radius : 1,
			smallBodyID : globalIndex
		});
		this.cache[globalIndex] = smallBody;

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
