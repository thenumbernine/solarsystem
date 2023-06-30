import {DOM, show, hide, toggleHidden, removeFromParent, animate} from '/js/util.js';
import {GLUtil} from '/js/gl-util.js';
import {ids, cfg} from './globals.js';
import {metricInfos} from './metric.js';
import {starSystemsExtra} from './starsystems.js';
import {makeGradient} from '/js/gl-util-Gradient.js';

let slideWidth = 300;
let currentOpenSidePanelID = undefined;
let allSidePanelIDs = [];
let displayConstellations = [];
let constellationIndexForName = {};
let canvas, gl, glutil;

const displayMethods = [
	'None',
	'Tangent Tidal',
	'Normal Tidal',
	'Total Tidal',
	'Tangent Gravitational',
	'Normal Gravitational',
	'Total Gravitational',
	'Tangent Total',
	'Normal Total',
	'Total'
];

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;

	//fix side panel heights
	allSidePanelIDs.forEach(sidePanelID => {
		ids[sidePanelID].style.height = window.innerHeight+'px';
	});

	//fix info panel height
	if (ui.showBodyInfo) {
		let infoDivDestTop = ids.timeControlDiv.offsetTop + ids.timeControlDiv.offsetHeight;
		ids.infoPanel.style.height = (window.innerHeight - infoDivDestTop)+'px';
	}

	glutil.resize();

	//TODO fix mouseline to work with this
	//glutil.view.fovY = Math.min(1, canvas.height / canvas.width) * 90;

	/*
	let aspectRatio = canvas.width / canvas.height;
	let nearHeight = Math.tan(glutil.view.fovY * Math.PI / 360);
	let nearWidth = aspectRatio * nearHeight;
	*/

	/** /
	//setup infinite projection matrix
	//http://www.terathon.com/gdc07_lengyel.pdf
	//http://www.gamasutra.com/view/feature/131351/the_mechanics_of_robust_stencil_.php?page=2
	let epsilon = 0;
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
	let resZNear = 1e+4;
	let resZFar = 1e+8;
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
	allSidePanelIDs.forEach(sidePanelID => {
		ids[sidePanelID].style.zIndex = 1;
	});
	let sidePanel = ids[sidePanelID];
	sidePanel.style.zIndex = 2;
	show(sidePanel);
	const startMenuLeft = ids.menu.offsetLeft;
	const startSidePanelLeft = sidePanel.offetLeft;
	animate({
		duration : cfg.slideDuration,
		callback : frac => {
			const degrees = frac * 180;
			ids.menu.style.transform = 'rotate('+degrees+'deg)';
			ids.menu.style.left = (startMenuLeft*(1-frac) + slideWidth*frac)+'px';
			sidePanel.style.left = (startSidePanelLeft*(1-frac) + 0*frac)+'px';
		},
	});
	currentOpenSidePanelID = sidePanelID;
}

function hideSidePanel(sidePanelID, dontMoveOpenButton) {
	if (sidePanelID === currentOpenSidePanelID) currentOpenSidePanelID = undefined;
	let sidePanel = ids[sidePanelID];
	if (!dontMoveOpenButton) {
		const startMenuLeft = ids.menu.offsetLeft;
		animate({
			duration : cfg.slideDuration,
			callback : frac => {
				const degrees = (1 - frac) * 180;
				ids.menu.style.transform = 'rotate('+degrees+'deg)';
				ids.menu.style.left = (startMenuLeft*(1-frac) + 0*frac)+'px';
			},
		});
	}
	const startSidePanelLeft = sidePanel.offetLeft;
	animate({
		duration : cfg.slideDuration,
		callback : frac => {
			sidePanel.style.left = (startSidePanelLeft*(1-frac) + -slideWidth*frac)+'px';
		},
		done : () => {
			hide(sidePanel);
		},
	});
}

let ui = new function() {
	this.init = function() {
		this.showBodyInfo = false;
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

		//keep track of what menu is open
		//if none are open, open the main menu
		//if any are ... if it's not main
		ids.menu.addEventListener('click', e => {
			if (currentOpenSidePanelID === undefined) {
				showSidePanel('mainSidePanel');
			} else if (currentOpenSidePanelID == 'mainSidePanel') {
				hideSidePanel(currentOpenSidePanelID);
			} else {
				showSidePanel('mainSidePanel');
			}
		});

		[
			{buttonID:'mainButtonDisplayOptions', divID:'displayOptionsSidePanel'},
			{buttonID:'mainButtonOverlay', divID:'overlaySidePanel'},
			{buttonID:'mainButtonSolarSystem', divID:'solarSystemSidePanel'},
			{buttonID:'mainButtonSmallBodies', divID:'smallBodiesSidePanel'},
			{buttonID:'mainButtonConstellations', divID:'constellationsSidePanel'},
			{buttonID:'mainButtonStarSystems', divID:'starSystemsSidePanel'},
			{buttonID:'mainButtonControls', divID:'controlsSidePanel'},
		].forEach(info => {
			ids[info.buttonID].addEventListener('click', e => {
				showSidePanel(info.divID);
			});
		});

		ids.reset.addEventListener('click', e => {
			cfg.integrationPaused = true;
			cfg.integrateTimeStep = cfg.defaultIntegrateTimeStep;
			cfg.julianDate = cfg.initJulianDate;
			refreshCurrentTimeText();
		});
		ids.play.addEventListener('click', e => {
			cfg.integrationPaused = false;
			cfg.integrateTimeStep = Math.abs(cfg.integrateTimeStep);
		});
		ids.reverse.addEventListener('click', e => {
			cfg.integrationPaused = false;
			cfg.integrateTimeStep = -Math.abs(cfg.integrateTimeStep);
		});
		ids.pause.addEventListener('click', e => {
			cfg.integrationPaused = true;
		});
		//fast/slow vs ffwd/rewind?
		ids.ffwd.addEventListener('click', e => {
			cfg.integrateTimeStep *= 2;
		});
		ids.rewind.addEventListener('click', e => {
			cfg.integrateTimeStep /= 2;
		});

		canvas = DOM('canvas', {
			css : {
				left : 0,
				top : 0,
				position : 'absolute',
				background : 'red',
				userSelect : 'none',
			},
			prependTo : document.body,
		});

		try {
			glutil = new GLUtil({canvas:canvas});
			gl = glutil.context;
		} catch (e) {
			removeFromParent(ids.menu);
			removeFromParent(canvas);
			show(ids.webglfail);
			throw e;
		}
		glutil.import('Gradient', makeGradient);
		ui.canvas = canvas;
		ui.gl = gl;
		ui.glutil = glutil;

		document.addEventListener('keydown', e => {
			switch (e.keyCode) {
			case 32:	//space
				cfg.integrationPaused = !cfg.integrationPaused;
				break;
			case 38:	//up
				cfg.integrateTimeStep *= 2;
				break;
			case 40:	//down
				cfg.integrateTimeStep /= 2;
				break;
			case 39:	//right
				cfg.integrationPaused = false;
				cfg.integrateTimeStep = Math.abs(cfg.integrateTimeStep);
				break;
			case 37:	//left
				cfg.integrationPaused = false;
				cfg.integrateTimeStep = -Math.abs(cfg.integrateTimeStep);
				break;
			}
		});
	};

	this.initSidePanel = function() {
		metricInfos.forEach((metricInfo, metricIndex) => {
			let radio = DOM('input', {
				type : 'radio',
				name : 'calculationMetric',
				value : metricIndex,
				click : e => {
					cfg.metric = new metricInfo.classObj();
					calcTides.invalidateForces();
				},
				attrs : {
					name : 'metric',
				},
				appendTo : ids.overlaySidePanelMetric,
			})
			if (cfg.metric instanceof metricInfo.classObj) radio.checked = true;
			DOM('span', {text:metricInfo.name, appendTo:ids.overlaySidePanelMetric});
			DOM('br', {appendTo:ids.overlaySidePanelMetric});
		});

		DOM('span', {text:'Overlay:', appendTo:ids.overlaySidePanelContents});
		DOM('br', {appendTo:ids.overlaySidePanelContents});
		displayMethods.forEach((thisDisplayMethod,displayMethodIndex) => {
			const radio = DOM('input', {
				type : 'radio',
				name : 'displayMethods',
				value : displayMethodIndex,
				click : e => {
					cfg.displayMethod = thisDisplayMethod;
					calcTides.invalidateForces();
				},
				attrs : {
					name : 'display',
				},
				appendTo : ids.overlaySidePanelContents,
			})
			if (thisDisplayMethod == cfg.displayMethod) radio.checked = true;
			DOM('span', {text:thisDisplayMethod, appendTo:ids.overlaySidePanelContents});
			DOM('br', {appendTo:ids.overlaySidePanelContents});
		});
		DOM('br', {appendTo:ids.overlaySidePanelContents});
		DOM('span', {text:'Influencing Planets:', appendTo:ids.overlaySidePanelContents});
		DOM('br', {appendTo:ids.overlaySidePanelContents});

		//add radio buttons hierarchically ...
		let overlayControlsForPlanets = {};

		class HierarchicalCheckboxControl {
			/*
			args:
				title		<- used to identify this checkbox
				change		<- callback upon changing the checkbox value
				clickTitle	<- callback upon clicking the title
				isChecked
				... and anything else that can be referenced through this.args
			*/
			constructor(args) {
				this.args = args;
				this.childControls = [];
				this.div = DOM('div', {
					css : {paddingLeft:'5px'},
				});

				let thiz = this;
				this.checkbox = DOM('input', {
					type : 'checkbox',
					change : e => {
						//refresh all parent controls' moon checkboxes -- to whiteout or grey them
						for (let c = thiz.parentControls; c; c = c.parentControls) {
							c.recomputeMoonCheckbox();
						}

						args.change.call(thiz);
					},
					checked : args.isChecked,
					appendTo : this.div,
				});

				DOM('span', {
					text : args.title,
					css : {
						textDecoration : 'underline',
						cursor : 'pointer',
					},
					click : e => {
						args.clickTitle.call(thiz);
					},
					appendTo : this.div,
				});

				this.toggleChildDiv = DOM('span', {
					css : {
						paddingLeft : '10px',
						cursor : 'pointer',
					},
					appendTo:this.div,
				});

				this.moonCheckbox = DOM('input', {
					type : 'checkbox',
					change : e => {
						//select / deselect all children
						thiz.setAllChildren(thiz.moonCheckbox.checked);
					},
					checked : true,
					appendTo : this.toggleChildDiv,
				});

				DOM('span', {
					css : {
						cursor : 'pointer'
					},
					text : '...',
					click : e => { toggleHidden(thiz.childDiv); },
					appendTo : this.toggleChildDiv,
				});

				DOM('br', {appendTo:this.div});

				this.childDiv = DOM('div', {appendTo:this.div});
			}
			
			addChild(childControl) {
				this.childDiv.appendChild(childControl.div);
				childControl.parentControls = this;
				this.childControls.push(childControl);
			}
			
			setAllChildren(checked) {
				for (let i = 0; i < this.childControls.length; ++i) {
					let ch = this.childControls[i].checkbox;
					if (checked) {
						if (!ch.checked) { ch.checked = true; ch.dispatchEvent(new Event('change')); }
					} else {
						if (ch.checked) { ch.checked = fale; ch.dispatchEvent(new Event('change')); }
					}
					this.childControls[i].setAllChildren(checked);
				}
			}
			
			recomputeMoonCheckbox() {
				let numChecked = 0;
				let total = 0;
				for (let i = 0; i < this.childControls.length; ++i) {
					++total;
					//check the child
					if (this.childControls[i].checkbox.checked) {
						++numChecked;
					}
					//check the child's children if they exist
					if (this.childControls[i].childControls.length > 0) {
						if (this.childControls[i].moonCheckbox.checked) {
							if (this.childControls[i].moonCheckbox.indeterminate) {
								numChecked += .5;
							} else {
								++numChecked;
							}
						}
						++total;
					}
				}
				if (numChecked == 0) {
					this.moonCheckbox.checked = false;
					this.moonCheckbox.indeterminate = false;
				} else if (numChecked == total) {
					this.moonCheckbox.checked = true;
					this.moonCheckbox.indeterminate = false;
				} else {
					this.moonCheckbox.checked = true;
					this.moonCheckbox.indeterminate = true;
				}
			}
		}

		//TODO add star systems
		const solarSystem = starSystemsExtra.solarSystem;
		solarSystem.planets.forEach((planet, planetIndex) => {
			
			//if any other planet doesn't have recorded mass then skip it
			if (planet.mass === undefined) return;
			
			//and if it's a barycenter then skip it ... or at least process its children?
			if (planet.isBarycenter) return;

			let parentPlanet = planet.parent;
			if (parentPlanet !== undefined) {
				if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
				//skip past any parent planets that don't have controls (i.e. barycenters)
				while (parentPlanet && !overlayControlsForPlanets[parentPlanet]) {
					parentPlanet = parentPlanet.parent;
				}
			}

			let controls = new HierarchicalCheckboxControl({
				title : planet.name,
				isChecked : true,
				change : function() {
					cfg.planetInfluences[this.args.planetIndex] = this.checkbox.is(':checked');
					calcTides.invalidateForces();
				},
				clickTitle : function() {
					setOrbitTarget(solarSystem.planets[solarSystem.indexes[planet.name]]);
				},
				planetIndex : planetIndex
			});

			//add to parent or to the side panel
			if (parentPlanet === undefined) {
				ids.overlaySidePanelContents.appendChild(controls.div);
			} else {
				overlayControlsForPlanets[parentPlanet.index].addChild(controls);
			}

			cfg.planetInfluences[planetIndex] = true;
			overlayControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
		});

		for (let planetIndex in overlayControlsForPlanets) {
			planetIndex = +planetIndex;
			let controls = overlayControlsForPlanets[planetIndex];

			controls.recomputeMoonCheckbox();

			//if a planet didn't get any children, hide its 'toggle child div'
			if (controls.childDiv.children.length == 0) {
				hide(controls.toggleChildDiv);
				continue;
			}

			//if a planet got children and it's not the sun or earth then hide by default
			if (planetIndex !== solarSystem.indexes.Sun &&
				planetIndex !== solarSystem.indexes.Earth)
			{
				hide(controls.childDiv);
			}
		}

		DOM('br', {appendTo:ids.overlaySidePanelContents});


		// display options side panel

		const radioGroups = [
			['gravityWellScaleNormalized', 'gravityWellScaleFixed'],
			['overlayShowOrbitTarget', 'overlayShowCurrentPosition']
		];

		[
			'showLinesToOtherPlanets',
			'showVelocityVectors',
			'showRotationAxis',
			'showOrbitAxis',
			'showEllipseAxis',
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
		].forEach(toggle => {
			let checkbox = ids[toggle];
			if (cfg[toggle]) checkbox.checked = true;
			checkbox.addEventListener('change', e => {
				cfg[toggle] = checkbox.checked;
				
				let found = false;
				for (let i = 0; i < radioGroups.length; ++i) {
					const group = radioGroups[i];
					for (let j = 0; j < group.length; ++j) {
						if (group[j] == toggle) {
							for (let k = 0; k < group.length; ++k) {
								if (k == j) continue;
								cfg[group[k]] = false;
							}
							found = true;
							break;
						}
					}
					if (found) break;
				}
			});
		});

		//checkbox but the reflected variable is a number not a bool
		ids.flatEarthMode.addEventListener('change', e => {
			if (ids.flatEarthMode.checked) {
				cfg.targetFlatEarthCoeff = 1;
			} else {
				cfg.targetFlatEarthCoeff = 0;
			}
		});

		[
			'gravityWellScaleFixedValue',
			'starPointSizeBias',	//in starfield.js
			'starPointSizeScale',	//in starfield.js
			'planetScaleExaggeration'
		].forEach(toggle => {
			const textfield = ids[toggle];
			textfield.value = cfg[toggle];
			textfield.addEventListener('change', e => {
				cfg[toggle] = textfield.value;
			});
		});


		// celestial bodies side panel


		let solarSystemSidePanel = ids.solarSystemSidePanel;

		//add radio buttons hierarchically ...
		let celestialBodiesControlsForPlanets = {};

		let cometParent;

		solarSystem.planets.forEach((planet, planetIndex) => {
			let parentPlanet = planet.parent;
			if (parentPlanet !== undefined) {
				if (parentPlanet.index >= planetIndex) throw "parent index should be < planet index or undefined";
			}

			let controls = new HierarchicalCheckboxControl({
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
				ids.celestialBodiesVisiblePlanets.appendChild(controls.div);
			} else {
				celestialBodiesControlsForPlanets[parentPlanet.index].addChild(controls);
			}

			celestialBodiesControlsForPlanets[planetIndex] = controls;	//JS only handles string keys, so get ready to typecast back to int
		});

		if (cometParent) cometParent.recomputeMoonCheckbox();
		for (let planetIndex in celestialBodiesControlsForPlanets) {
			planetIndex = +planetIndex;
			let controls = celestialBodiesControlsForPlanets[planetIndex];

			controls.recomputeMoonCheckbox();

			//if a planet didn't get any children, hide its 'toggle child div'
			if (controls.childDiv.children.length == 0) {
				hide(controls.toggleChildDiv);
				continue;
			}

			//if a planet got children and it's not the sun or earth then hide by default
			if (planetIndex !== solarSystem.indexes.Sun &&
				planetIndex !== solarSystem.indexes.Earth)
			{
				hide(controls.childDiv);
			}
		}

		DOM('br', {appendTo:ids.celestialBodiesVisiblePlanets});

		//these are added to the end of the result
		//they should get greyed upon new query (search, prev, next click)
		//they should be regenerated upon new content
		//they should be re-enabled upon error
		let nextButton = undefined;
		let prevButton = undefined;

		let searchResults = [];

		let searchLastID = 0;	//uid to inc so if we search twice, the first can be invalidated

		//dataSource = 'remote' for remote queries, 'local' for the currently-selected planets
		let processSearch = function(pageIndex, dataSource) {
console.log("small-bodies got a search request...");			
			let button = ids.celestialBodiesSearch;
			let searchText = ids.celestialBodiesSearchText;
			let searchStr = searchText.val();
			button.disabled = true;
			searchText.val('searching...');
			searchText.disabled = true;

			if (prevButton) prevButton.disabled = true;
			if (nextButton) nextButton.disabled = true;

			let searchID = ++searchLastID;

			let processResults = function(results) {
				if (searchID < searchLastID-1) return;	//someone else has searched

				searchText.val(searchStr);
				searchText.disabled = false;
				button.disabled = false;

				ids.celestialBodiesSearchToggleAll.checked = false;	//disable too?

				let pageSize = 20;	//fixed atm
				let pageMax = Math.floor((results.count-1) / pageSize);

				searchResults = [];

				let resultsDiv = ids.celestialBodiesSearchResults;
				resultsDiv.innerHTML = '';
				results.rows.forEach((row,i) => {
					let rowDiv = DOM('div', 
						{appendTo:resultsDiv}
					);

					let name = row.name;

					let titleSpan;
					let checkbox = DOM('input', {
						type : 'checkbox',
						change : e => {
							if (!checkbox.checked) {	//uncheck checkbox => remove planet
								solarSystem.removeSmallBody(row);
								titleSpan.css({textDecoration:'', cursor:''});
							} else {	//check checkbox => add planet
								solarSystem.addSmallBody(row);
								titleSpan.css({textDecoration:'underline', cursor:'pointer'});
							}
						},
						checked : solarSystem.indexes[name] !== undefined,
						appendTo : rowDiv,
					})

					titleSpan = DOM('span', {
						text : name,
						click : e => {
							let targetPlanet = solarSystem.planets[solarSystem.indexes[name]];
							if (targetPlanet !== undefined) setOrbitTarget(targetPlanet);
						},
						appendTo : rowDiv,
					});
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
					let changePage = function(dir) {
						//TODO remove or grey out results as well?
						if (nextButton) nextButton.disabled = true;
						if (prevButton) prevButton.disabled = true;
						processSearch(pageIndex+dir, dataSource);
					};
					if (pageIndex > 0) {
						prevButton = DOM('button', {
							text : 'Prev',
							click : e => { changePage(-1); },
							appendTo : ids.celestialBodiesSearchResults,
						});
					}
					if (pageIndex < pageMax-1) {
						nextButton = DOM('button', {
							text : 'Next',
							click : e => { changePage(1); },
							appendTo : ids.celestialBodiesSearchResults,
						});
					}
					DOM('span', {
						text : (pageIndex+1)+' of '+pageMax,
						appendTo : ids.celestialBodiesSearchResults,
					});
				}
			};

			if (dataSource == 'remote') {
				const fetchArgs = new URLSearchParams();
				fetchArgs.set('comet', ids.celestialBodiesSearchComets.checked?1:0);
				fetchArgs.set('numbered' , ids.celestialBodiesSearchNumbered.checked?1:0);
				fetchArgs.set('unnumbered', ids.celestialBodiesSearchUnnumbered.checked?1:0);
				fetchArgs.set('text', searchStr);
				fetchArgs.set('page', pageIndex+1);	//1-based
				//cache : false,
				//timeout : 30000
				fetch('/solarsystem/jpl-ssd-smallbody/search.lua?'+fetchArgs.toString())
				.then(response => {
					if (!response.ok) return Promise.reject('not ok');
					response.json()
					.then(obj => {
						processResults(obj);
					});
				}).catch(e => {
console.log("search error", arguments);
					searchText.val(searchStr);
					searchText.disabled = false;
					button.disabled = false;
					//TODO animate background color of search text
					if (prevButton) prevButton.setAttribute('disabled', 0);
					if (nextButton) nextButton.setAttribute('disabled', 0);

					let warning = DOM('div', {text:'Connection Failed!', css:{color:'red'}});
					ids.celestialBodiesSearchWarning.after(warning);
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
				});
			} else if (dataSource == 'local') {
				let rows = [];
				let searchingComets = ids.celestialBodiesSearchComets.checked;
				let searchingNumbered = ids.celestialBodiesSearchNumbered.checked;
				let searchingUnnumbered = ids.celestialBodiesSearchUnnumbered.checked;
				for (let i = 0; i < solarSystem.planets.length; ++i) {
					let planet = solarSystem.planets[i];
					let row = planet.sourceData;
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
		"horizonID" : "C",
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
		"horizonID" : "1P",
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
		"horizonID":"357439",
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

		ids.celestialBodiesSearchText.addEventListener('keydown', e => {
			if (e.keyCode == 13) {
				ids.celestialBodiesSearch.dispatchEvent(new Event('click'));
			}
		});

		//change a check box, immediately update search results
		[
			ids.celestialBodiesSearchComets,
			ids.celestialBodiesSearchNumbered,
			ids.celestialBodiesSearchUnnumbered,
			ids.celestialBodiesSearchVisible
		].forEach(checkbox => {
			checkbox.addEventListener('change', e => {
				ids.celestialBodiesSearch.dispatchEvent(new Event('click'));
			});
		});

		ids.celestialBodiesSearchToggleAll.addEventListener('click', e => {
			let checked = ids.celestialBodiesSearchToggleAll.checked;
			searchResults.forEach((result, i) => {
				if (result.checkbox.checked != checked) {
					result.checkbox.dispatchEvent(new Event('click'));
				}
			});
		});

		ids.celestialBodiesSearch.addEventListener('click', e => {
			processSearch(0,
				ids.celestialBodiesSearchVisible.checked ? 'local' : 'remote'
			);
		});
		
//		ids.celestialBodiesSearch.dispatchEvent(new Event('click'));	//fire one off

		
		// constellations



		let constellationsResults = ids.constellationsResults;
		
		
		//TODO - slider / min/max for filtering stars by ... (a) app. mag, or (b) abs mag
		// (or (c) app. mag from custom location?)
		DOM('span', {
			text : 'selected min:',
			appendTo : constellationsResults,
		});
		ids.constellationsSelMinMag = DOM('span', {
			id : 'constellationsSelMinMag',
			appendTo : constellationsResults,
		});
		DOM('br', {appendTo : constellationsResults});
		
		DOM('span', {
			text : 'selected max:',
			appendTo : constellationsResults,
		});
		ids.constellationsSelMaxMag = DOM('span', {
			id : 'constellationsSelMaxMag',
			appendTo : constellationsResults,
		});
		DOM('br', {appendTo : constellationsResults});
	
/* TODO finish me
		DOM('span', {
			text : 'filter min:',
			appendTo : constellationsResults,
		});
		ids.constellationsMinMag = DOM('span', {
			id : 'constellationsMinMag'
			appendTo : constellationsResults,
		});
		DOM('br', {appendTo : constellationsResults});
		
		DOM('span', {
			text : 'filter max:'
			appendTo : constellationsResults,
		});
		ids.constellationsMaxMag = DOM('span', {
			id : 'constellationsMaxMag',
			appendTo : constellationsResults,
		});
		DOM('br', {appendTo : constellationsResults});
*/


		let sortedConstellationNames = [];
		constellations.forEach((con,i) => {
			sortedConstellationNames.push(con.name);
		});
		sortedConstellationNames.sort();
		
		constellations.forEach((con,i) => {
			constellationIndexForName[con.name] = i;
		});

		let updateMagMinMax = function() {
			let magmin = Infinity;
			let magmax = -Infinity;
			for (let k = 0; k < constellations.length; ++k) {
				if (displayConstellations[k]) {
					magmin = Math.min(magmin, constellations[k].mag.min);
					magmax = Math.max(magmax, constellations[k].mag.max);
				}
			}
			ids.constellationsSelMinMag.innerText = magmin;
			ids.constellationsSelMaxMag.innerText = magmax;
		};

		sortedConstellationNames.forEach((name,i) => {
			DOM('input', {
				type : 'checkbox',
				change : function() {
					let index = constellationIndexForName[name];
					displayConstellations[index] = !displayConstellations[index];
					updateMagMinMax();
				},
				appendTo : constellationsResults,
			});
		
			DOM('span', {
				text : name,
				appendTo : constellationsResults,
			});
			DOM('br', {appendTo : constellationsResults});
		});

		updateMagMinMax();

		// rest of the init


		window.addEventListener('resize', resize);
		resize();
	};
};

export {ui}
