var slideDuration = 500;
var slideWidth = 300;
var currentOpenSidePanelID = undefined;
var showBodyInfo = false;
var allSidePanelIDs = [];
var displayConstellations = [];
var constellationIndexForName = {};

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;

	//fix side panel heights
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('height', window.innerHeight);
	});

	//fix info panel height
	if (showBodyInfo) {
		var infoDivDestTop = $('#timeControlDiv').offset().top + $('#timeControlDiv').height();
		$('#infoPanel').css('height', window.innerHeight - infoDivDestTop);
	}

	glutil.resize();

	//TODO fix mouseline to work with this
	//glutil.view.fovY = Math.min(1, canvas.height / canvas.width) * 90;

	/*
	var aspectRatio = canvas.width / canvas.height;
	var nearHeight = Math.tan(glutil.view.fovY * Math.PI / 360);
	var nearWidth = aspectRatio * nearHeight;
	*/

	/** /
	//setup infinite projection matrix
	//http://www.terathon.com/gdc07_lengyel.pdf
	//http://www.gamasutra.com/view/feature/131351/the_mechanics_of_robust_stencil_.php?page=2
	var epsilon = 0;
	mat4.identity(glutil.scene.projMat);
	glutil.scene.projMat[0] = glutil.view.zNear / nearWidth;
	glutil.scene.projMat[5] = glutil.view.zNear / nearHeight;
	glutil.scene.projMat[10] = epsilon-1;
	glutil.scene.projMat[11] = (epsilon - 2) * glutil.view.zNear;
	glutil.scene.projMat[14] = -1;
	glutil.scene.projMat[15] = 0;
	/**/

	/* or we could just linearly ramp the depth over the range we want, since the rest is going to full anyways * /
	//(-z - zNear) / (zFar - zNear) * 2 - 1
	//z * (-2 / (zFar - zNear)) + ( - 2 * zNear / (zFar - zNear) - 1)
	var resZNear = 1e+4;
	var resZFar = 1e+8;
	mat4.identity(glutil.scene.projMat);
	glutil.scene.projMat[0] = glutil.view.zNear / nearWidth;
	glutil.scene.projMat[5] = glutil.view.zNear / nearHeight;
	glutil.scene.projMat[10] = -2 * (resZFar - resZNear);
	glutil.scene.projMat[11] = -(1 + 2 * resZNear / (resZFar - resZNear));
	glutil.scene.projMat[14] = -1;
	glutil.scene.projMat[15] = 0;
	/**/
}

function showSidePanel(sidePanelID) {
	if (currentOpenSidePanelID !== undefined) {
		hideSidePanel(currentOpenSidePanelID, true);
	}
	$.each(allSidePanelIDs, function(i,sidePanelID) {
		$('#'+sidePanelID).css('z-index', 1);
	});
	var sidePanel = $('#'+sidePanelID);
	sidePanel.css('z-index', 2);
	sidePanel.show();
	$('#menu').animate(
		{
			left : slideWidth
		}, {
			duration : slideDuration,
			step : function(now, fx) {
				var degrees = now / slideWidth * 180;
				setCSSRotation($(this), degrees);
			}
		}
	);
	sidePanel.animate({left:0}, {duration:slideDuration});
	currentOpenSidePanelID = sidePanelID;
}

function hideSidePanel(sidePanelID, dontMoveOpenButton) {
	if (sidePanelID === currentOpenSidePanelID) currentOpenSidePanelID = undefined;
	var sidePanel = $('#'+sidePanelID);
	if (!dontMoveOpenButton) {
		$('#menu').animate(
			{
				left : 0
			},
			{
				duration : slideDuration,
				step : function(now, fx) {
					var degrees = (now / slideWidth) * 180;
					setCSSRotation($(this), degrees);
				}
			}
		);
	}
	sidePanel.animate({left:-slideWidth}, {
		duration:slideDuration,
		complete:function() {
			sidePanel.hide();
		}
	});
}

var ui = new function() {
	this.init = function() {
		allSidePanelIDs = [
			'mainSidePanel',
			'displayOptionsSidePanel',
			'overlaySidePanel',
			'solarSystemSidePanel',
			'smallBodiesSidePanel',
			'constellationsSidePanel',
			'starSystemsSidePanel',
			'controlsSidePanel'
		];
		mainSidePanel = $('#mainSidePanel');

		//keep track of what menu is open
		//if none are open, open the main menu
		//if any are ... if it's not main
		$('#menu').click(function() {
			if (currentOpenSidePanelID === undefined) {
				showSidePanel('mainSidePanel');
			} else if (currentOpenSidePanelID == 'mainSidePanel') {
				hideSidePanel(currentOpenSidePanelID);
			} else {
				showSidePanel('mainSidePanel');
			}
		});

		$.each([
			{buttonID:'mainButtonDisplayOptions', divID:'displayOptionsSidePanel'},
			{buttonID:'mainButtonOverlay', divID:'overlaySidePanel'},
			{buttonID:'mainButtonSolarSystem', divID:'solarSystemSidePanel'},
			{buttonID:'mainButtonSmallBodies', divID:'smallBodiesSidePanel'},
			{buttonID:'mainButtonConstellations', divID:'constellationsSidePanel'},
			{buttonID:'mainButtonStarSystems', divID:'starSystemsSidePanel'},
			{buttonID:'mainButtonControls', divID:'controlsSidePanel'},
		], function(i, info) {
			$('#'+info.buttonID).click(function() {
				showSidePanel(info.divID);
			});
		});

		$('#reset').click(function() {
			integrationPaused = true;
			integrateTimeStep = defaultIntegrateTimeStep;
			julianDate = initJulianDate;
			refreshCurrentTimeText();
		});
		$('#play').click(function() {
			integrationPaused = false;
			integrateTimeStep = Math.abs(integrateTimeStep);
		});
		$('#reverse').click(function() {
			integrationPaused = false;
			integrateTimeStep = -Math.abs(integrateTimeStep);
		});
		$('#pause').click(function() {
			integrationPaused = true;
		});
		//fast/slow vs ffwd/rewind?
		$('#ffwd').click(function() {
			integrateTimeStep *= 2;
		});
		$('#rewind').click(function() {
			integrateTimeStep /= 2;
		});

		canvas = $('<canvas>', {
			css : {
				left : 0,
				top : 0,
				position : 'absolute',
				background : 'red'
			}
		}).prependTo(document.body).get(0);

		$(canvas).disableSelection();

		try {
			glutil = new GLUtil({canvas:canvas});
			gl = glutil.context;
		} catch (e) {
			$('#menu').remove();
			$(canvas).remove();
			$('#webglfail').show();
			throw e;
		}

		$(document).keydown(function(e) {
			switch (e.keyCode) {
			case 32:	//space
				integrationPaused = !integrationPaused;
				break;
			case 38:	//up
				integrateTimeStep *= 2;
				break;
			case 40:	//down
				integrateTimeStep /= 2;
				break;
			case 39:	//right
				integrationPaused = false;
				integrateTimeStep = Math.abs(integrateTimeStep);
				break;
			case 37:	//left
				integrationPaused = false;
				integrateTimeStep = -Math.abs(integrateTimeStep);
				break;
			}
		});
	};

	this.initSidePanel = function() {
		var overlaySidePanelMetric = $('#overlaySidePanelMetric');
		$.each(metricInfos, function(metricIndex, metricInfo) {
			var radio = $('<input>', {
				type : 'radio',
				name : 'calculationMetric',
				value : metricIndex,
				click : function() {
					metric = new metricInfo.classObj();
					calcTides.invalidateForces();
				}
			})
				.attr('name', 'metric')
				.appendTo(overlaySidePanelMetric);
			if (metric.isa(metricInfo.classObj)) radio.attr('checked', 'checked');
			$('<span>', {text:metricInfo.name}).appendTo(overlaySidePanelMetric);
			$('<br>').appendTo(overlaySidePanelMetric);
		});

		var overlaySidePanelContents = $('#overlaySidePanelContents');
		$('<span>', {text:'Overlay:'}).appendTo(overlaySidePanelContents);
		$('<br>').appendTo(overlaySidePanelContents);
		$.each(displayMethods, function(displayMethodIndex,thisDisplayMethod) {
			var radio = $('<input>', {
				type : 'radio',
				name : 'displayMethods',
				value : displayMethodIndex,
				click : function() {
					displayMethod = thisDisplayMethod;
					calcTides.invalidateForces();
				}
			})
				.attr('name', 'display')
				.appendTo(overlaySidePanelContents);
			if (thisDisplayMethod == displayMethod) radio.attr('checked', 'checked');
			$('<span>', {text:thisDisplayMethod}).appendTo(overlaySidePanelContents);
			$('<br>').appendTo(overlaySidePanelContents);
		});
		$('<br>').appendTo(overlaySidePanelContents);
		$('<span>', {text:'Influencing Planets:'}).appendTo(overlaySidePanelContents);
		$('<br>').appendTo(overlaySidePanelContents);

		//add radio buttons hierarchically ...
		var overlayControlsForPlanets = {};

		var HierarchicalCheckboxControl = makeClass({
			/*
			args:
				title		<- used to identify this checkbox
				change		<- callback upon changing the checkbox value
				clickTitle	<- callback upon clicking the title
				isChecked
				... and anything else that can be referenced through this.args
			*/
			init : function(args) {
				this.args = args;
				this.childControls = [];
				this.div = $('<div>', {
					css : {paddingLeft:'5px'}
				});

				var thiz = this;
				this.checkbox = $('<input>', {
					type : 'checkbox',
					change : function() {
						//refresh all parent controls' moon checkboxes -- to whiteout or grey them
						for (var c = thiz.parentControls; c; c = c.parentControls) {
							c.recomputeMoonCheckbox();
						}

						args.change.call(thiz);
					}
				})
					.prop('checked', args.isChecked)
					.appendTo(this.div);

				$('<span>', {
					text : args.title,
					css : {
						textDecoration : 'underline',
						cursor : 'pointer',
					},
					click : function(e) {
						args.clickTitle.call(thiz);
					}
				}).appendTo(this.div);

				this.toggleChildDiv = $('<span>', {
					css : {
						paddingLeft : '10px',
						cursor : 'pointer'
					}
				}).appendTo(this.div);

				this.moonCheckbox = $('<input>', {
					type : 'checkbox',
					change : function() {
						//select / deselect all children
						thiz.setAllChildren($(this).prop('checked'));
					}
				})
					.prop('checked', 1)
					.appendTo(this.toggleChildDiv);

				$('<span>', {
					css : {
						cursor : 'pointer'
					},
					text : '...',
					click : function() {
						if (thiz.childDiv.css('display') == 'none') {
							thiz.childDiv.show();
						} else {
							thiz.childDiv.hide();
						}
					}
				}).appendTo(this.toggleChildDiv);

				$('<br>').appendTo(this.div);

				this.childDiv = $('<div>').appendTo(this.div);

			},
			addChild : function(childControl) {
				childControl.div.appendTo(this.childDiv);
				childControl.parentControls = this;
				this.childControls.push(childControl);
			},
			setAllChildren : function(checked) {
				for (var i = 0; i < this.childControls.length; ++i) {
					var ch = this.childControls[i].checkbox;
					if (checked) {
						if (!ch.prop('checked')) { ch.prop('checked', 1); ch.trigger('change'); }
					} else {
						if (ch.prop('checked')) { ch.prop('checked', 0); ch.trigger('change'); }
					}
					this.childControls[i].setAllChildren(checked);
				}
			},
			recomputeMoonCheckbox : function() {
				var numChecked = 0;
				var total = 0;
				for (var i = 0; i < this.childControls.length; ++i) {
					++total;
					//check the child
					if (this.childControls[i].checkbox.prop('checked')) {
						++numChecked;
					}
					//check the child's children if they exist
					if (this.childControls[i].childControls.length > 0) {
						if (this.childControls[i].moonCheckbox.prop('checked')) {
							if (this.childControls[i].moonCheckbox.prop('indeterminate')) {
								numChecked += .5;
							} else {
								++numChecked;
							}
						}
						++total;
					}
				}
				if (numChecked == 0) {
					this.moonCheckbox.prop('checked', 0)
						.prop('indeterminate', 0);
				} else if (numChecked == total) {
					this.moonCheckbox.prop('checked', 1)
						.prop('indeterminate', 0);
				} else {
					this.moonCheckbox.prop('checked', 1)
						.prop('indeterminate', 1);
				}
			}
		});

		//TODO add star systems
		$.each(solarSystem.planets, function(planetIndex, planet) {
			
			//if any other planet doesn't have recorded mass then skip it
			if (planet.mass === undefined) return;
			
			//and if it's a barycenter then skip it ... or at least process its children?
			if (planet.isBarycenter) return;

			var parentPlanet = planet.parent;
			if (parentPlanet !== undefined) {
				if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
				//skip past any parent planets that don't have controls (i.e. barycenters)
				while (parentPlanet && !overlayControlsForPlanets[parentPlanet]) {
					parentPlanet = parentPlanet.parent;
				}
			}

			var controls = new HierarchicalCheckboxControl({
				title : planet.name,
				isChecked : true,
				change : function() {
					planetInfluences[this.args.planetIndex] = this.checkbox.is(':checked');
					calcTides.invalidateForces();
				},
				clickTitle : function() {
					setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
				},
				planetIndex : planetIndex
			});

			//add to parent or to the side panel
			if (parentPlanet === undefined) {
				controls.div.appendTo(overlaySidePanelContents);
			} else {
				overlayControlsForPlanets[parentPlanet.index].addChild(controls);
			}

			planetInfluences[planetIndex] = true;
			overlayControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
		});

		for (var planetIndex in overlayControlsForPlanets) {
			var planetIndex = +planetIndex;
			var controls = overlayControlsForPlanets[planetIndex];

			controls.recomputeMoonCheckbox();

			//if a planet didn't get any children, hide its 'toggle child div'
			if (controls.childDiv.children().length == 0) {
				controls.toggleChildDiv.hide();
				continue;
			}

			//if a planet got children and it's not the sun or earth then hide by default
			if (planetIndex !== solarSystem.indexes.Sun &&
				planetIndex !== solarSystem.indexes.Earth)
			{
				controls.childDiv.hide();
			}
		}

		$('<br>').appendTo(overlaySidePanelContents);


		// display options side panel

		var radioGroups = [
			['gravityWellScaleNormalized', 'gravityWellScaleFixed'],
			['overlayShowOrbitTarget', 'overlayShowCurrentPosition']
		];

		$.each([
			'showLinesToOtherPlanets',
			'showVelocityVectors',
			'showRotationAxis',
			'showOrbitAxis',
			'showLatAndLonLines',
			'showGravityWell',
			'showNames',				//overlays.js
			'showOrbits',
			'showSmallBodies',			//smallbodies.js
			'allowSelectSmallBodies',	//smallbodies.js
			'showStars',				//in starfield.js
			'allowSelectStars',			//in starfield.js
			'showGalaxies',				//in galaxies.js
			'allowSelectGalaxies',
			'showPlanetsAsDistantPoints',
			//radio
			'gravityWellScaleNormalized',
			'gravityWellScaleFixed',
			//radio
			'overlayShowOrbitTarget',
			'overlayShowCurrentPosition'
		], function(_, toggle) {
			var checkbox = $('#'+toggle);
			if (window[toggle]) checkbox.attr('checked', 'checked');
			checkbox.change(function() {
				window[toggle] = checkbox.is(':checked');
				
				var found = false;
				for (var i = 0; i < radioGroups.length; ++i) {
					var group = radioGroups[i];
					for (var j = 0; j < group.length; ++j) {
						if (group[j] == toggle) {
							for (var k = 0; k < group.length; ++k) {
								if (k == j) continue;
								window[group[k]] = false;
							}
							found = true;
							break;
						}
					}
					if (found) break;
				}
			});
		});

		$.each([
			'gravityWellScaleFixedValue',
			'starsVisibleMagnitudeBias',	//in starfield.js
			'planetScaleExaggeration'
		], function(_, toggle) {
			(function(){
				var textfield = $('#'+toggle);
				textfield.val(window[toggle]);
				textfield.change(function() {
					window[toggle] = textfield.val();
				});
			})();
		});


		// celestial bodies side panel


		var solarSystemSidePanel = $('#solarSystemSidePanel');

		//add radio buttons hierarchically ...
		var celestialBodiesControlsForPlanets = {};

		var cometParent;

		$.each(solarSystem.planets, function(planetIndex,planet) {

			var parentPlanet = planet.parent;
			if (parentPlanet !== undefined) {
				if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
			}

			var controls = new HierarchicalCheckboxControl({
				title : planet.name,
				isChecked : !planet.hide,
				change : function() {
					solarSystem.planets[this.args.planetIndex].hide = !this.checkbox.is(':checked');
				},
				clickTitle : function() {
					setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
				},
				planetIndex : planetIndex
			});

			if (parentPlanet === undefined) {
				controls.div.appendTo($('#celestialBodiesVisiblePlanets'));
			} else {
				celestialBodiesControlsForPlanets[parentPlanet.index].addChild(controls);
			}

			celestialBodiesControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
		});

		if (cometParent) cometParent.recomputeMoonCheckbox();
		for (var planetIndex in celestialBodiesControlsForPlanets) {
			var planetIndex = +planetIndex;
			var controls = celestialBodiesControlsForPlanets[planetIndex];

			controls.recomputeMoonCheckbox();

			//if a planet didn't get any children, hide its 'toggle child div'
			if (controls.childDiv.children().length == 0) {
				controls.toggleChildDiv.hide();
				continue;
			}

			//if a planet got children and it's not the sun or earth then hide by default
			if (planetIndex !== solarSystem.indexes.Sun &&
				planetIndex !== solarSystem.indexes.Earth)
			{
				controls.childDiv.hide();
			}
		}

		$('<br>').appendTo($('#celestialBodiesVisiblePlanets'));

		//these are added to the end of the result
		//they should get greyed upon new query (search, prev, next click)
		//they should be regenerated upon new content
		//they should be re-enabled upon error
		var nextButton = undefined;
		var prevButton = undefined;

		var searchResults = [];

		var searchLastID = 0;	//uid to inc so if we search twice, the first can be invalidated

		//dataSource = 'remote' for remote queries, 'local' for the currently-selected planets
		var processSearch = function(pageIndex, dataSource) {
			var button = $('#celestialBodiesSearch');
			var searchText = $('#celestialBodiesSearchText');
			var searchStr = searchText.val();
			button.prop('disabled', 1);
			searchText.val('searching...');
			searchText.prop('disabled', 1);

			if (prevButton) prevButton.prop('disabled', 1);
			if (nextButton) nextButton.prop('disabled', 1);

			var searchID = ++searchLastID;

			var processResults = function(results) {
				if (searchID < searchLastID-1) return;	//someone else has searched

				searchText.val(searchStr);
				searchText.prop('disabled', 0);
				button.prop('disabled', 0);

				$('#celestialBodiesSearchToggleAll').prop('checked', 0);	//disable too?

				var pageSize = 20;	//fixed atm
				var pageMax = Math.floor((results.count-1) / pageSize);

				searchResults = [];

				var resultsDiv = $('#celestialBodiesSearchResults');
				resultsDiv.empty();
				$.each(results.rows, function(i,row) {
					var rowDiv = $('<div>');
					rowDiv.appendTo(resultsDiv);

					var name = row.name;

					var titleSpan;
					var checkbox = $('<input>', {
						type : 'checkbox',
						change : function() {
							if (!$(this).is(':checked')) {	//uncheck checkbox => remove planet
								solarSystem.removeSmallBody(row);
								titleSpan.css({textDecoration:'', cursor:''});
							} else {	//check checkbox => add planet
								solarSystem.addSmallBody(row);
								titleSpan.css({textDecoration:'underline', cursor:'pointer'});
							}
						}
					})
						.prop('checked', solarSystem.indexes[name] !== undefined)
						.appendTo(rowDiv);

					titleSpan = $('<span>', {
						text : name,
						click : function() {
							var targetPlanet = solarSystem.planets[solarSystem.indexes[name]];
							if (targetPlanet !== undefined) setOrbitTarget(targetPlanet);
						}
					}).appendTo(rowDiv);
					//TODO put an 'add' button next to each
					//on clicking it, add the body to the planet list
					//and repopulate the 'extra' div

					searchResults.push({
						checkbox : checkbox,
						div : rowDiv,
						data : row,
						name : name
					});
				});

				//and add prev/next/pages if there is
				if (pageMax > 0) {
					var changePage = function(dir) {
						//TODO remove or grey out results as well?
						if (nextButton) nextButton.prop('disabled', 1);
						if (prevButton) prevButton.prop('disabled', 1);
						processSearch(pageIndex+dir, dataSource);
					};
					if (pageIndex > 0) {
						prevButton = $('<button>', {
							text : 'Prev',
							click : function() {
								changePage(-1);
							}
						}).appendTo($('#celestialBodiesSearchResults'));
					}
					if (pageIndex < pageMax) {
						nextButton = $('<button>', {
							text : 'Next',
							click : function() {
								changePage(1);
							}
						}).appendTo($('#celestialBodiesSearchResults'));
					}
					$('<span>', {text:(pageIndex+1)+' of '+pageMax}).appendTo($('#celestialBodiesSearchResults'));
				}
			};

			if (dataSource == 'remote') {
				$.ajax({
					url : '/solarsystem/jpl-ssd-smallbody/search.lua',
					dataType : 'json',
					data : {
						comet : $('#celestialBodiesSearchComets').prop('checked')?1:0,
						numbered : $('#celestialBodiesSearchNumbered').prop('checked')?1:0,
						unnumbered : $('#celestialBodiesSearchUnnumbered').prop('checked')?1:0,
						text : searchStr,
						page : pageIndex+1	//1-based
					},
					cache : false
				}).error(function() {
					searchText.val(searchStr);
					searchText.prop('disabled', 0);
					button.prop('disabled', 0);
					//TODO animate background color of search text
					if (prevButton) prevButton.prop('disabled', 0);
					if (nextButton) nextButton.prop('disabled', 0);

					var warning = $('<div>', {text:'Connection Failed!', css:{color:'red'}});
					$('#celestialBodiesSearchWarning').after(warning);
					setTimeout(function() {
						//after five seconds, fade away
						warning.animate({
							height : 0,
							opacity : 0
						}, {
							duration : 500,
							complete : function() {
								warning.remove();
							}
						});
					}, 3000);
				}).done(processResults);
			} else if (dataSource == 'local') {
				var rows = [];
				var searchingComets = $('#celestialBodiesSearchComets').prop('checked');
				var searchingNumbered = $('#celestialBodiesSearchNumbered').prop('checked');
				var searchingUnnumbered = $('#celestialBodiesSearchUnnumbered').prop('checked');
				for (var i = 0; i < solarSystem.planets.length; ++i) {
					var planet = solarSystem.planets[i];
					var row = planet.sourceData;
					if (row) {
						if ((row.bodyType == 'comet' && searchingComets) ||
							(row.bodyType == 'numbered asteroid' && searchingNumbered) ||
							(row.bodyType == 'unnumbered asteroid' && searchingUnnumbered))
						{
							rows.push(row);
						}
					}
				}
				
/*hack for debugging* /
rows = [];
if (false) {
   //test case for hyperbolic
	rows.push({
		"perihelionDistance" : 209232954020.56,
		"inclination" : 2.2519519080902,
		"timeOfPerihelionPassage" : 2456955.82823,
		"argumentOfPeriapsis" : 0.042602963442406,
		"bodyType" : "comet",
		"orbitSolutionReference" : "JPL 101",
		"epoch" : 56691,
		"eccentricity" : 1.00074241,
		"idNumber" : "C",
		"name" : "2013 A1 (Siding Spring)",
		"longitudeOfAscendingNode" : 5.2530270564044
	});
//target C1/2013 Siding Spring data: on 2014-11-15 10:30:00
//position m: 	64031628815.774, -196629902235.24, 57392197050.865
//velocity m/day: -1916173862.297, -455287182.34414, 2315832701.4279
	
	//test case for elliptical
	rows.push({
		"perihelionDistance" : 87661077532.81,
		"inclination" : 2.8320181936429,
		"timeOfPerihelionPassage" : 2446467.39532,
		"argumentOfPeriapsis" : 1.9431185149437,
		"bodyType" : "comet",
		"orbitSolutionReference" : "JPL J863/77",
		"epoch" : 49400,
		"eccentricity" : 0.96714291,
		"idNumber" : "1P",
		"name" : "Halley",
		"longitudeOfAscendingNode" : 1.0196227452785
	});
}
if (true) {	//recent asteroid passing by
	rows.push({
		"meanAnomaly":6.1790425090414,
		"longitudeOfAscendingNode":2.2116873367796,
		"inclination":0.41440347267775,
		"name":"2004 BL86",
		"eccentricity":0.40307315,
		"idNumber":"357439",
		"bodyType":"numbered asteroid",
		"orbitSolutionReference":"JPL 34",
		"absoluteMagnitude":19.1,
		"argumentOfPeriapsis":5.4324242142291,
		"magnitudeSlopeParameter":0.15,
		"semiMajorAxis":224726235521.07,
		"epoch":57000
	});
}
/**/
				
				processResults({rows:rows, count:rows.length});
			} else {
				throw "got an unknown data source request " + dataSource;
			}
		};

		$('#celestialBodiesSearchText').keydown(function(e){
			if (e.keyCode == 13) {
				$('#celestialBodiesSearch').trigger('click');
			}
		});

		//change a check box, immediately update search results
		$.each([
			$('#celestialBodiesSearchComets'),
			$('#celestialBodiesSearchNumbered'),
			$('#celestialBodiesSearchUnnumbered'),
			$('#celestialBodiesSearchVisible')
		], function(_,checkbox) {
			checkbox.change(function() {
				$('#celestialBodiesSearch').trigger('click');
			});
		});

		$('#celestialBodiesSearchToggleAll').click(function() {
			var checked = $(this).is(':checked');
			$.each(searchResults, function(i,result) {
				if (result.checkbox.prop('checked') != checked) {
					result.checkbox.trigger('click');
				}
			});
		});

		$('#celestialBodiesSearch').click(function() {
			processSearch(0,
				$('#celestialBodiesSearchVisible').prop('checked') ? 'local' : 'remote'
			);
		});
		$('#celestialBodiesSearch').trigger('click');	//fire one off

		
		// constellations



		var constellationsResults = $('#constellationsResults');
		var sortedConstellationNames = [];
		$.each(constellations, function(i,con) {
			sortedConstellationNames.push(con.name);
		});
		sortedConstellationNames.sort();
		
		$.each(constellations, function(i,con) {
			constellationIndexForName[con.name] = i;
		});
		
		$.each(sortedConstellationNames, function(i,name) {
			$('<input>', {
				type : 'checkbox',
				change : function() {
					var index = constellationIndexForName[name];
					displayConstellations[index] = !displayConstellations[index];
				}
			})
				.appendTo(constellationsResults);
		
			$('<span>', {
				text : name,
			})
				.appendTo(constellationsResults);
			$('<br>').appendTo(constellationsResults);
		});


		// rest of the init


		$(window).resize(resize);
		resize();
	};
};
