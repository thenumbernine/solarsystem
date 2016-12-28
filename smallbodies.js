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
		//Here I'm convering the node body data to match the sql row returned from search.lua  (found in coldesc.lua)
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
			row[field] = node[field+'Array'][nodeLocalIndex];
		}
		row.bodyType = ['comet', 'numbered asteroid', 'unnumbered asteroid'][row.bodyType];

		smallBody.row = row;
	
		return smallBody;
	}

});

var smallBodies;
