[![Donate via Stripe](https://img.shields.io/badge/Donate-Stripe-green.svg)](https://buy.stripe.com/00gbJZ0OdcNs9zi288)<br>
[![Donate via Bitcoin](https://img.shields.io/badge/Donate-Bitcoin-green.svg)](bitcoin:37fsp7qQKU8XoHZGRQvVzQVP8FrEJ73cSJ)<br>

![solar system](https://cdn.rawgit.com/thenumbernine/solarsystem/master/images/screenshot.png)

## Solar System in WebGL

Allows you see gravitational / tidal forces of various planets on various other planets.

Calculations are done using Newtonian coefficients rewritten in terms of Riemann metric tensor.  See Ch. 12 of "Gravitation" by Misner, Thorne, and Wheeler.
This is equivalent to linearized gravity when the velocity is assumed to be purely timelike.  I have a cartesian Schwarzschild metric in black-hole-skymap,
remind me to use that for this.  Until then, linearized-timelike/Newtonian will just have to do.

Somewhere in this repo is the original Lua version, designed to pull from astro-phys.com or from NASA's Ephemeris 406 data, and to subsequentially integrate forward or backwards in time with a RK45 integrator.
I've since been corrected that I should be using a Gauss Cubic integrator.  Expect to see the integrator, and a cross-domain pull from astro-phys.com.  Don't expect me to keep a copy of the Ephemeris data.  After parsing it was someting like 400mb.  Play with the Lua scripts if you want to mess with that yourself.

The textures are (c) NASA.  

- Sun is from http://commons.m.wikimedia.org/wiki/File:Map_of_the_full_sun.jpg
- Mercury from http://laps.noaa.gov/albers/sos/sos.html
- Venus without clouds from http://laps.noaa.gov/albers/sos/sos.html
- Venus with clouds from http://www.mmedia.is/~bjj/data/venus/
- Earth from http://m.earthobservatory.nasa.gov/Features/BlueMarble
- Mars, Phobos, Deimos from http://maps.jpl.nasa.gov
- Jupiter, Enceladus, Mimas, Tethys, Dione, Rhea, Titan, Iapetus, Phoebe from Cassini, images by http://www.ciclops.org/maps.php
	- Io, Europa from http://laps.noaa.gov/albers/sos/sos.html
	- Ganymede, Callisto from https://ksp.sarbian.com/nathankell/RealSolarSystem/Textures/4096/
- Jupiter rings from http://www.celestiamotherlode.net/catalog/jupiter.php
- Saturn from http://www.mmedia.is/~bjj/data/saturn/
- Saturn rings from http://www.mmedia.is/~bjj/data/s_rings/index.html
- Neptune from http://www.mmedia.is/~bjj/data/neptune/
- Pluto from http://en.spaceengine.org/forum/22-758-6
	- Charon from http://en.spaceengine.org/forum/22-758-6

Position and velocity data is from NASA Horizons telnet server horizons.jpl.nasa.gov
Planet characteristics are from NASA Horizons server and from http://solarsystem.nasa.gov/planets/index.cfm
Comets and small bodies are from NASA JPL SSD data files at http://ssd.jpl.nasa.gov/?sb_elem

See horizons/README for how to gather data from NASA Horizons telnet site
See jpl-ssd-smallbody/README for how to gather data from NASA JPL/SSD
See simbad/README for how to get intergalactic data from Simbad
See hyg/README for how to get the HYG star database

TODO

- planets:
	- put ephemeris data online, and re-enable the ffwd/rewind system
	- get keplerian parameters merged with integrator / with ephemeris data
		- then only integrate the main 11 bodies and update all other planets accordingly
		- option to use rk45 integrator?
- all objects: hierarchy of rendering for each planet and orbital children
	- give it a smaller bounding sphere for the planet mesh
	- give it a larger bounding sphere for children
	- occlusion test, traverse if passes 
	- give each node separate scale factor and coordinate system units, so we don't lose precision as we go from solar system to galaxy to supercluster 
	- when mouseover from a distance, if a planet and its moons are compressed to a single point in view, have the planet name show up first (related to hierarchical representation of planets)
- add a atmosphere shader to the Earth (and other planets) so it doesn't look so dull
- render screen to a float buffer and use HDR / post-processing / point spread function
- gravity / tidal force calculations:
	- make overlays computed by GPU so they can recompute realtime during integration.  will accuracy be good enough?
	- redo as Schwarzschild or Kerr calculations rather than newtonian limit calculations
- cron task for updating positions from NASA? fix bug in dynamic data program stalling after pluto?
- merge with larger star catalog
- merge with universe data
- small bodies:
	- add support for hyperbolas (Vinti-6)
	- make a million-point static buffer for rendering comets?
		- calculate via orbital elements 
		- download the current locations from NASA?  180 objects is slow, I'd imagine a million would take a while.
		- maybe make a script for calculating position and velocity from keplerian elements / GPU program / think about this one ...
- point source rendering
	- close-up render of correct sizes for stars in the HYG database (I don't have their radius data on file)
	- distance render of correct sizes / intensities / colors for both (a) extrasolar stars and (b) planets, when within the solar system
- extra feature: show stars proportional to our own sun
- merge planet toggle with asteroid/comet search/toggle with star search/select/toggle?
- hierarchy of exoplanets
- fix bugs in phone touch movement.  it gets a jump at the beginning or end of a rotation or pinch zoom.

Hierarchical structure:

starfield = new Starfield(); has Star() objects
	<- derived from the following:
		stardata.f32 = float buffer used for rendering stars as a point cloud
			contains star positions, absolute magnitude, color index, and maybe velocity later
		namedStars.json = data file for selecting individual stars and star systems

Star has the following:
	- planets ( should or shouldn't include moons?  no for hierarchy, yes for integrating n-body sim)
	- smallbodies 	
		<- local to our solar system.
		currently remote, because it contains a million objects
		might do point cloud with these?
		might do gpu-based kepler generation for time-varying coordinates ... if there is accuracy enough on the gpu ...
		doesn't seem intuitive to allow selection of all million at once, or visually sensible to show all orbits at once ...

planets has the following:
	- moons.  currently only for our solar system.  grabbed from NASA Horizons and another NASA site.
