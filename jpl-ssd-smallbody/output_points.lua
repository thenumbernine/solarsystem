local class = require 'ext.class'
local OutputMethod = require 'output'

local sunMassInKg = 1.9891e+30	-- kg
local gravitationalConstant = 6.6738480e-11		-- m^3 / (kg * s^2)
local julian = assert(loadfile('../horizons/julian.lua'))()

local vec3 = require 'vec.vec3'

local julianDate = julian.fromCalendar(os.date'!*t')

-- and for just outputting a cloud of points ...
local OutputToPoints = class(OutputMethod)

function OutputToPoints:init(args)
	self.currentBodyType = args.bodyType
end

function OutputToPoints:staticInit()
	self.bodies = table()
end

function OutputToPoints:processBody(body)
	local isComet = self.currentBodyType == 0
	local isAsteroid = self.currentBodyType > 0

	local semiMajorAxis
	local eccentricity
	local eccentricAnomaly
	local orbitType
	local meanAnomaly
	local A, B

	local eccentricity = body.eccentricity

	--[[
	local parabolicEccentricityEpsilon = 1e-7
	if math.abs(eccentricity - 1) < parabolicEccentricityEpsilon then
		orbitType = 'parabolic'
	else
	--]]
	if eccentricity > 1 then
		orbitType = 'hyperbolic'
	else
		orbitType = 'elliptic'
	end

	local pericenterDistance
	if isComet then
		pericenterDistance = assert(body.perihelionDistance)
		if orbitType ~= 'parabolic' then
			semiMajorAxis = pericenterDistance / (1 - eccentricity)
		else		--otherwise, if it is parabolic, we don't get the semi-major axis ...
			semiMajorAxis = 0/0
		end
	elseif isAsteroid then
		semiMajorAxis = assert(body.semiMajorAxis)
		pericenterDistance = semiMajorAxis * (1 - eccentricity)
	end

	local gravitationalParameter = gravitationalConstant * sunMassInKg	--assuming the comet mass is negligible, since the comet mass is not provided
	local semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis
	
	--orbital period is only defined for circular and elliptical orbits (not parabolic or hyperbolic)
	local orbitalPeriod
	if orbitType == 'elliptic' then
		orbitalPeriod = 2 * math.pi * math.sqrt(semiMajorAxisCubed / gravitationalParameter) / (60*60*24)	--julian day
	end

	local longitudeOfAscendingNode = assert(body.longitudeOfAscendingNode)
	local cosAscending = math.cos(longitudeOfAscendingNode)
	local sinAscending = math.sin(longitudeOfAscendingNode)

	local argumentOfPeriapsis = assert(body.argumentOfPeriapsis)
	local cosPericenter = math.cos(argumentOfPeriapsis)
	local sinPericenter = math.sin(argumentOfPeriapsis)

	local inclination = assert(body.inclination)
	local cosInclination = math.cos(inclination)
	local sinInclination = math.sin(inclination)

	local oneMinusEccentricitySquared = 1 - eccentricity * eccentricity
	--magnitude of A is a 
	A = vec3(semiMajorAxis * (cosAscending * cosPericenter - sinAscending * sinPericenter * cosInclination),
		 semiMajorAxis * (sinAscending * cosPericenter + cosAscending * sinPericenter * cosInclination),
		 semiMajorAxis * sinPericenter * sinInclination)
	--magnitude of B is a * sqrt(|1 - e^2|)
	B = vec3(semiMajorAxis * math.sqrt(math.abs(oneMinusEccentricitySquared)) * -(cosAscending * sinPericenter + sinAscending * cosPericenter * cosInclination),
		 semiMajorAxis * math.sqrt(math.abs(oneMinusEccentricitySquared)) * (-sinAscending * sinPericenter + cosAscending * cosPericenter * cosInclination),
		 semiMajorAxis * math.sqrt(math.abs(oneMinusEccentricitySquared)) * cosPericenter * sinInclination)
	--inner product: A dot B = 0

	local timeOfPeriapsisCrossing
	if isComet then
		timeOfPeriapsisCrossing = body.timeOfPerihelionPassage	--julian day
	end

	if orbitType == 'parabolic' then
		eccentricAnomaly = math.tan(argumentOfPeriapsis / 2)
		meanAnomaly = eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3 
	elseif orbitType == 'hyperbolic' then
		assert(timeOfPeriapsisCrossing)	--only comets are hyperbolic, and all comets have timeOfPeriapsisCrossing defined
		meanAnomaly = math.sqrt(-gravitationalParameter / semiMajorAxisCubed) * timeOfPeriapsisCrossing * 60*60*24	--in seconds
	elseif orbitType == 'elliptic' then
		--in theory I can say 
		--eccentricAnomaly = math.acos((eccentricity + math.cos(argumentOfPeriapsis)) / (1 + eccentricity * math.cos(argumentOfPeriapsis)))
		-- ... but math.acos has a limited range ...

		if isComet then
			timeOfPeriapsisCrossing = body.timeOfPerihelionPassage	--julian day
			local timeSinceLastPeriapsisCrossing = julianDate - timeOfPeriapsisCrossing
			meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * math.pi / orbitalPeriod
		elseif isAsteroid then
			meanAnomaly = body.meanAnomaly
		else
		--	error 'here'
		end
	else
		error 'here'
	end
	
	--solve Newton Rhapson
	--for elliptical orbits:
	--	f(E) = M - E + e sin E = 0
	--	f'(E) = -1 + e cos E
	--for parabolic oribts:
	--	f(E) = M - E - E^3 / 3
	--	f'(E) = -1 - E^2
	--for hyperbolic orbits:
	--	f(E) = M - e sinh(E) - E
	--	f'(E) = -1 - e cosh(E)
	eccentricAnomaly = meanAnomaly
	for i=1,10 do
		local func, deriv
		if orbitType == 'parabolic' then	--parabolic
			func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3
			deriv = -1 - eccentricAnomaly * eccentricAnomaly
		elseif orbitType == 'elliptic' then 	--elliptical
			func = meanAnomaly - eccentricAnomaly + eccentricity * math.sin(eccentricAnomaly)
			deriv = -1 + eccentricity * math.cos(eccentricAnomaly)	--has zeroes ...
		elseif orbitType == 'hyperbolic' then	--hyperbolic
			func = meanAnomaly + eccentricAnomaly - eccentricity  * math.sinh(eccentricAnomaly)
			deriv = 1 - eccentricity * math.cosh(eccentricAnomaly)
		else
			error 'here'
		end

		local delta = func / deriv
		if math.abs(delta) < 1e-15 then break end
		eccentricAnomaly = eccentricAnomaly - delta
	end

	local sinEccentricAnomaly = math.sin(eccentricAnomaly)
	local cosEccentricAnomaly = math.cos(eccentricAnomaly)
	
	--parabolas and hyperbolas don't define orbitalPeriod
	-- so no need to recalculate it
	if orbitalPeriod and meanAnomaly then
		timeOfPeriapsisCrossing = meanAnomaly * orbitalPeriod / (2 * math.pi) --if it is a comet then we're just reversing the calculation above ...
	end

	local dt_dE
	if orbitType == 'parabolic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 + eccentricAnomaly * eccentricAnomaly)
	elseif orbitType == 'elliptic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (1 - eccentricity * math.cos(eccentricAnomaly))
	elseif orbitType == 'hyperbolic' then
		dt_dE = math.sqrt(semiMajorAxisCubed / gravitationalParameter) * (eccentricity * math.cosh(eccentricAnomaly) - 1)
	end
	local dE_dt = 1/dt_dE
	--finally using http://en.wikipedia.org/wiki/Kepler_orbit like I should've in the first place ...
	local coeffA, coeffB
	local coeffDerivA, coeffDerivB
	if orbitType == 'parabolic' then
		--...?
	elseif orbitType == 'elliptic' then
		coeffA = math.cos(eccentricAnomaly) - eccentricity
		coeffB = math.sin(eccentricAnomaly)
		coeffDerivA = -math.sin(eccentricAnomaly) * dE_dt
		coeffDerivB = math.cos(eccentricAnomaly) * dE_dt
	elseif orbitType == 'hyperbolic' then
		coeffA = eccentricity - math.cosh(eccentricAnomaly)
		coeffB = math.sinh(eccentricAnomaly)
		coeffDerivA = -math.sinh(eccentricAnomaly) * dE_dt
		coeffDerivB = math.cosh(eccentricAnomaly) * dE_dt
	elseif orbitType == 'parabolic' then
		coeffA = 0/0
		coeffB = 0/0
	end

	local pos = A * coeffA + B * coeffB

	if not math.isfinite(pos[1]) or not math.isfinite(pos[2]) or not math.isfinite(pos[3]) then
		-- returning early and avoiding bad data is causing our max radius to report a higher value ?!
		return
	end
	
	local r = pos:length()
	if not math.isfinite(r) then return end
	
	OutputToPoints.maxR = math.max(OutputToPoints.maxR or r, r)

	body.pos = pos

	self.bodies:insert(body)
	body.index = #self.bodies
end

function OutputToPoints:staticDone()
	print('max radius:',self.maxR)

	local leafPointCount = 1000

	-- [[ writing octree
	-- now that we're done, construct the octree here
	-- store the individual node data and load it client-side
	-- then remotely grab names as leafs are mouse-over'd
	local mins = vec3(-6e+12, -6e+12, -6e+12)
	local maxs = vec3(6e+12, 6e+12, 6e+12)
	local size = maxs - mins

	local PointOctreeNode = class()
	function PointOctreeNode:init()
		self.bodies = table()
		self.mins = vec3()
		self.maxs = vec3()
		self.center = vec3()
	end
	
	local allNodes = table()

	local root = PointOctreeNode()
	allNodes:insert(root)
	root.index = #allNodes
	root.mins = vec3(mins:unpack())
	root.maxs = vec3(maxs:unpack())
	root.center = (mins + maxs) * .5

	for i=1,#self.bodies do
		local body = self.bodies[i]
		local x,y,z = body.pos:unpack()
		local node = root
		
		-- first add to leafmost until it passes a threshold, then split
		while node.children do
			local ix = x > node.center[1] and 1 or 0
			local iy = y > node.center[2] and 1 or 0
			local iz = z > node.center[3] and 1 or 0
			local childIndex = 1 + ix + 2 * (iy + 2 * iz)
			node = node.children[childIndex]
		end
	
		node.bodies:insert(body)
		if #node.bodies > leafPointCount then
			node.children = table()
			for ix=0,1 do
				for iy=0,1 do
					for iz=0,1 do
						local is = vec3(ix,iy,iz)
						local childIndex = 1 + ix + 2 * (iy + 2 * iz)
						local child = PointOctreeNode()
						allNodes:insert(child)
						child.index = #allNodes
						node.children[childIndex] = child
						child.parent = node
						for j=1,3 do
							child.mins[j] = is[j] == 1 and node.center[j] or node.mins[j]
							child.maxs[j] = is[j] == 1 and node.maxs[j] or node.center[j]
							child.center[j] = .5 * (child.mins[j] + child.maxs[j])
						end
					end
				end
			end
			for j=1,#node.bodies do
				local body = node.bodies[j]
				local x,y,z = body.pos:unpack()
				local ix = x > node.center[1] and 1 or 0
				local iy = y > node.center[2] and 1 or 0
				local iz = z > node.center[3] and 1 or 0
				local childIndex = 1 + ix + 2 * (iy + 2 * iz)
				local child = node.children[childIndex]
				child.bodies:insert(body)
			end
			node.bodies = table()
		end
	end
	
	-- now push upward based on random sampling of all children ...
	local function process(node)
		if not node.children then return end
		for i=1,leafPointCount do
			local child = node.children[math.random(8)]
			if child then
				if #child.bodies == 0 then
					process(child)
				end
				if #child.bodies > 0 then
					node.bodies:insert(child.bodies:remove(math.random(#child.bodies)))
				end
			end
		end
	end
	process(root)
	
	print('num nodes',#allNodes)

	local json = require 'dkjson'
	local pts = table()
	for _,body in ipairs(self.bodies) do
		for j=1,3 do
			pts:insert(body.pos[j])
		end
	end	
	for _,node in ipairs(allNodes) do
		setmetatable(node, nil)
		if node.parent then node.parent = node.parent.index end
		if node.children then
			for i=1,8 do
				if node.children[i] then
					node.children[i] = node.children[i].index
				else
					node.children[i] = -1	-- empty value -- for the sake of making the array dense
				end
			end
		end
		
		--[[
		-- now repackage points continuously with all nodes
		-- and store in each node the range information
		node.bodyStartIndex = #pts/3
		for _,body in ipairs(node.bodies) do
			for j=1,3 do
				pts:insert(body.pos[j])
			end
		end
		node.bodyEndIndex = #pts/3
		--]]

		node.bodies = node.bodies:map(function(body)
			return {
				name = body.name,
				-- either need one or the other : (a) repackage and keep range, or (b) keep individual node mapping
				index = body.index,	-- index in the total list 
			}
		end)
		file['nodes/'..node.index..'.json'] = json.encode(node.bodies):gsub('},{','},\n{')
		node.bodies = nil
	end
	for _,node in ipairs(allNodes) do
		node.index = nil	-- don't need this anymore
	end
	file['octree.json'] = json.encode(setmetatable(allNodes, nil)):gsub('},{','},\n{')

	local ffi = require 'ffi'
	local buffer = ffi.new('float[?]', #pts, pts)
	file['output.f32'] = ffi.string(ffi.cast('char*', buffer), ffi.sizeof(buffer))

	--]]

	--[[ writing all at once:
	local ffi = require 'ffi'
	require 'ffi.c.stdio'
	local fh = ffi.C.fopen('output.f32', 'wb')
	local v = ffi.new('float[3]')
	for i=1,#self.bodies do
		local p = self.bodies[i].pos
		v[0] = p[1]
		v[1] = p[2]
		v[2] = p[3]
		ffi.C.fwrite(v, ffi.sizeof(v), 1, fh)
	end
	ffi.C.fclose(fh)
	--]]

end

return OutputToPoints
