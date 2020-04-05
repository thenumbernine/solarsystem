require 'ext'

local sunMassInKg = 1.9891e+30	-- kg
local gravitationalConstant = 6.6738480e-11		-- m^3 / (kg * s^2)
local julian = assert(loadfile('../horizons/julian.lua'))()

local vec3 = require 'vec.vec3'
local json = require 'dkjson'
--<? if _VERSION ~= 'Lua 5.3' ?>
local bit = bit32 or require 'bit'
--<? else ?>
local bit = {
	bor = function(a,b) return a | b end,
	lshift = function(a,b) return a << b end,
}
--<? end ?>
local ffi = require 'ffi'

local julianDate = julian.fromCalendar(os.date'!*t')

local writeIndividualNodes = false
local writeOneBigFile = true

-- and for just outputting a cloud of points ...
local OutputToPoints = class()

function OutputToPoints:init(args)
	self.currentBodyType = args.bodyType
end

function OutputToPoints:staticInit()
	self.bodies = table()
end

-- http://www.bogan.ca/orbits/kepler/orbteqtn.html
-- https://www.mathworks.com/matlabcentral/fileexchange/13439-orbital-mechanics-library/content/randv.m
function OutputToPoints:processBody(body)
	local isComet = self.currentBodyType == 0
	local isAsteroid = self.currentBodyType > 0
	body.bodyType = self.currentBodyType

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
		
		-- if ffi is used then this will be true ... so I'll use nan as well
		pericenterDistance = assert(math.isfinite(body.perihelionDistance) and body.perihelionDistance)
		
		if orbitType ~= 'parabolic' then
			semiMajorAxis = pericenterDistance / (1 - eccentricity)
		else		--otherwise, if it is parabolic, we don't get the semi-major axis ...
			semiMajorAxis = math.nan
		end
	elseif isAsteroid then
		semiMajorAxis = assert(math.isfinite(body.semiMajorAxis) and body.semiMajorAxis)
		pericenterDistance = semiMajorAxis * (1 - eccentricity)
	end

	local gravitationalParameter = gravitationalConstant * sunMassInKg	--assuming the comet mass is negligible, since the comet mass is not provided
	local semiMajorAxisCubed = semiMajorAxis * semiMajorAxis * semiMajorAxis
	
	--orbital period is only defined for circular and elliptical orbits (not parabolic or hyperbolic)
	local orbitalPeriod = math.nan
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
		-- this must mean no hyperbolic orbits are numbered / unnumbered small bodies? 
		assert(timeOfPeriapsisCrossing)	--only comets are hyperbolic, and all comets have timeOfPeriapsisCrossing defined
		meanAnomaly = math.sqrt(-gravitationalParameter / semiMajorAxisCubed) * timeOfPeriapsisCrossing * 60*60*24	--in seconds
	elseif orbitType == 'elliptic' then
		--in theory I can say 
		--eccentricAnomaly = math.acos((eccentricity + math.cos(argumentOfPeriapsis)) / (1 + eccentricity * math.cos(argumentOfPeriapsis)))
		-- ... but math.acos has a limited range ...

		if isComet then
			-- comets should use perihelion distance and time of perihelion passage 
			-- instead of semimajor axis and mean anomaly
			-- so how to convert ...
			assert(timeOfPeriapsisCrossing)	--julian day
			local timeSinceLastPeriapsisCrossing = julianDate - timeOfPeriapsisCrossing
			meanAnomaly = timeSinceLastPeriapsisCrossing * 2 * math.pi / orbitalPeriod
		elseif isAsteroid then
			-- https://en.wikipedia.org/wiki/Mean_anomaly
			-- M = M0 + n (t - t0)
			-- M0 = mean anomaly at epoch
			-- t0 = epoch
			-- n = 2 pi / orbitalPeriod
			meanAnomaly = body.meanAnomalyAtEpoch + 2 * math.pi / orbitalPeriod * (julianDate - body.epoch)
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
			func = meanAnomaly - (eccentricAnomaly - eccentricity * math.sin(eccentricAnomaly))
			deriv = -(1 - eccentricity * math.cos(eccentricAnomaly))	--has zeroes ...
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

	-- E = eccentric anomaly
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
		coeffA = math.nan
		coeffB = math.nan
	end

	local pos = A * coeffA + B * coeffB

	local r = pos:length()
	if math.isfinite(r) then 
		OutputToPoints.maxR = math.max(OutputToPoints.maxR or r, r)
	end

	body.pos = pos

	self.bodies:insert(body)
	body.index = #self.bodies

	body.semiMajorAxis = semiMajorAxis
	body.eccentricity = eccentricity
	body.eccentricAnomaly = eccentricAnomaly
	body.longitudeOfAscendingNode = longitudeOfAscendingNode
	body.argumentOfPeriapsis = argumentOfPeriapsis
	body.inclination = inclination
	body.timeOfPeriapsisCrossing = timeOfPeriapsisCrossing
	body.meanAnomaly = meanAnomaly
	body.orbitType = ({
		elliptic = 0,
		hyperbolic = 1,
		parabolic = 2,
	})[orbitType]
	body.orbitalPeriod = orbitalPeriod	--only exists for elliptical orbits
	body.A = A
	body.B = B
end

function OutputToPoints:staticDone()
	print('max radius:',self.maxR)
	
	print('building octree leafs')
	local leafPointCount = 1000

	-- [[ writing octree
	-- now that we're done, construct the octree here
	-- store the individual node data and load it client-side
	-- then remotely grab names as leafs are mouse-over'd
	
	--[[
	nodeID's can be omitted if we use a hash to find them instead
	let node -1 be the root
	then nodes 0-7 represent bits 000..111 that correspond with the zL/R, yL/R, xL/R sides of the spatial median for where the child goes 
	then nodes 000xxx..111xxx represent the next layer's children
	etc

	so for N layers you need 3N bits
	so for 32 bit indexes, we get 10 layers in our tree 

	and the 0-2 bits are always the level0 position
	3-5 bits are always the level1 position
	etc

	that'll overlap everything with parent 000
	so we need to add previous geometric sums ...

	L0 = 0
	L1 = [1,9) = 1 + 0-7
	L2 = [1+8,1+8+64)
	L3 = [1+8+64,1+8+64+128)

	geom sum = (8^(n+1) - 1) / (8 - 1)
	n=0: sum = 1
	n=1: sum = 1+8 = 9
	n=2: sum = 1+8*(1+8) = 1+8+64 = 73

	so, for the n+1'th depth, add the sum to level n
	for root, add 0
	for depth=1, add the sum for n=0 
	
	--]]

	--some vertexes are infinite, so I'm just going to fix the octree bounds at Pluto's orbit distance
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
	root.nodeID = 0
	root.mins = vec3(mins:unpack())
	root.maxs = vec3(maxs:unpack())
	root.center = (mins + maxs) * .5

local function long(x)
	return ffi.cast('long', x)
end

	for i=1,#self.bodies do
		local body = self.bodies[i]
		local x,y,z = body:getPos()
		local node = root
		local levelID = long(0)

		if math.isfinite(x)
		and math.isfinite(y)
		and math.isfinite(z)
		then
			-- first add to leafmost until it passes a threshold, then split
			-- 0 level is the root node, involves no bit writes	into nodeID
			local depth = 0
			
			while node.children do
				local ix = x > node.center[1] and 1 or 0
				local iy = y > node.center[2] and 1 or 0
				local iz = z > node.center[3] and 1 or 0

--<? if _VERSION ~= 'Lua 5.3' then ?>
				local childIndex = bit.bor(
					ix,
					bit.lshift(iy, 1),
					bit.lshift(iz, 2))
				
				node = node.children[childIndex+1]
			
				levelID = bit.bor(
					levelID,
					bit.lshift(
						long(childIndex), 
						long(3*depth)
					))
				depth = depth + 1
				
				-- depth==0 has no bits
				-- depth==1 writes to bits 0..2
			end

			-- childDepth == 21 writes to bits 60..62
			-- so if depth >= 21 then we are writing past the end
			-- then again, doubles hold 0..52 bits of precision
			-- then again, javascript bit shifts only work on 32bit integers ...
			local maxDepth = 9	-- inclusive
			assert(depth < maxDepth, "looks like you need a bigger data size, or a cap on the max depth")

			node.bodies:insert(body)
			if #node.bodies > leafPointCount
			and depth < maxDepth-1
			then
				node.children = table()
				for ix=0,1 do
					for iy=0,1 do
						for iz=0,1 do
							local is = vec3(ix,iy,iz)
							local childIndex = bit.bor(
								ix,
								bit.lshift(iy, 1),
								bit.lshift(iz, 2))
							local child = PointOctreeNode()

							local childDepth = depth + 1
--[[
childDepth 	start	end	size
1			1		9	8
2			9		73	64
3			73		
--]]

						
							-- this is the size of the level
							-- child depth 0 means our level size is 1 = 8^0
							-- child depth 1 means our level size is 8 = 8^1
							local levelSize = bit.lshift(long(1), long(3*childDepth))
				
							-- this is where the level starts
							-- for depth = 1, levelStart = 1, levelEnd = 9
							-- level start (inclusive) = (8^depth-1)/(8-1)
							-- level end (exclusive) = (8^(depth+1)-1)/(8-1)
							local levelStart = (bit.lshift(long(1), long(3*childDepth)) - 1) / 7
						
							-- this is the levelID in the current level
							local childLevelID = bit.bor(
								levelID,
								bit.lshift(
									long(childIndex),
									long(3*depth)
								))

							-- offset by sum from n=0 to depth of 8^n
							-- depth is the parent's depth, so depth == 0 when the child's parent is root
							-- and the offset is (8^1 - 1) / 7 = 1, so we offset past the root only
							
							child.nodeID = childLevelID + levelStart
							child.childIndex = childIndex
							child.depth = childDepth

--print()
--print('parent depth', depth)
--print('parent nodeID', nodeID)
--print('childIndex', childIndex)
--print('childLevelID', childLevelID)
--print('childNodeID', child.nodeID)
--print('level size',levelSize)
--print('level start',levelStart)

							assert(childLevelID >= 0 and childLevelID < levelSize, 
								"got an oob child nodeID "..tostring(childLevelID)
								.. " should be between 0 and "..tostring(levelSize))
							
							
							-- TODO count bits instead of using math.log
							local check_depth = math.floor(math.log(7*tonumber(child.nodeID)+1,8))
							assert(check_depth == childDepth, 
								"for child.nodeID "..tostring(child.nodeID)
								.." with levelStart="..tostring(levelStart)
								.." expected child.depth to be "..tostring(childDepth)
								.." got "..tostring(check_depth))
							local check_levelStart = (bit.lshift(long(1), long(3*check_depth)) - 1) / 7
							local levelEnd = (bit.lshift(long(1), long(3*(check_depth+2))) - 1) / 7
							assert(
								check_levelStart <= child.nodeID
								and child.nodeID < levelEnd,
								"expected "..tostring(child.nodeID)
								.." on depth "..tostring(check_depth)
								.." to be between "..tostring(check_levelStart)
								.." and "..tostring(levelEnd))
							local check_levelID = child.nodeID - check_levelStart 
							assert(check_levelID == childLevelID, 
								"expected "..tostring(check_levelID)
								.." == "..tostring(levelNodeID)
								.." for depth "..tostring(check_depth)
								.." level start "..tostring(check_levelStart)
								.." level end "..tostring(levelEnd))

							node.children[childIndex+1] = child
							child.parent = node
							for j=1,3 do
								child.mins[j] = is[j] == 1 and node.center[j] or node.mins[j]
								child.maxs[j] = is[j] == 1 and node.maxs[j] or node.center[j]
								child.center[j] = .5 * (child.mins[j] + child.maxs[j])
							end
			
							for _,node in ipairs(allNodes) do
								assert(node.nodeID ~= child.nodeID, 
									"found matching nodeIDs: "..tostring(node.nodeID)
									.." childIndex "..tostring(node.childIndex)
									.." depth "..tostring(node.depth)
									.." and "..tostring(child.nodeID)
									.." childIndex "..tostring(child.childIndex)
									.." depth "..tostring(child.depth))
							end
							allNodes:insert(child)
						end
					end
				end
				for j=1,#node.bodies do
					local body = node.bodies[j]
					local x,y,z = body:getPos()
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
	end

	print('pushing points rootward')
	-- now push upward based on random sampling of all children ...
	local function process(node)
		if not node.children then return end
		for i=1,8 do
			local child = node.children[i]
			if child then
				if #child.bodies == 0 then 
					process(child) 
				end
				if #child.bodies > 0 then
					child.bodies:sort(function(a,b)
						return a.absoluteMagnitude < b.absoluteMagnitude
					end)
				end
			end
		end
		for i=1,leafPointCount do
			local bestAbsoluteMagnitude = math.huge 
			local bestChildren = {}
			for j=1,8 do
				local child = node.children[j]
				if child 
				and #child.bodies > 0 
				and child.bodies[1].absoluteMagnitude <= bestAbsoluteMagnitude
				then
					if bestAbsoluteMagnitude == child.bodies[1].absoluteMagnitude then
						table.insert(bestChildren,child)
					else
						bestAbsoluteMagnitude = child.bodies[1].absoluteMagnitude
						bestChildren = {child}
					end
				end
			end
			if not bestChildren then break end
			-- TODO pick the brightest instead of at random
			node.bodies:insert(bestChildren[math.random(1,#bestChildren)].bodies:remove(1))
		end
	end
	process(root)
	
	print('num nodes',#allNodes)

	local function tostringnum(x)
		return (tostring(x):match('(%d+)'))
	end	

	do
		print('building node dict...')
		local nodeDict = table()
		for _,node in ipairs(allNodes) do
			for localIndex,body in ipairs(node.bodies) do
				nodeDict:insert{
					name = ffi.string(body.name),
					nodeID = node.nodeID,
					index = localIndex,
				}
			end
		end
		print('sorting...')
		nodeDict:sort(function(a,b)
			return assert(a.name) < assert(b.name)
		end)
		print('serializing and writing...')
		local f = assert(io.open('node-dict.csv', 'w'))
		f:write'#name,nodeID,index\n'
		for i=1,#nodeDict do
			local entry = nodeDict[i]
			f:write('"'..entry.name..'",'..tostringnum(entry.nodeID)..','..entry.index..'\n')
		end
		f:close()
		
		print'done with node dictionary!'
	end

	local rowDesc = json.decode(file['row-desc.json']:match('.-=(.*)'))

	-- TODO: automatically make dirs for file[] access?
	if writeIndividualNodes then
		os.execute('mkdir nodes')
	end
	print('writing out nodes')

	local allNodeIDs = setmetatable(allNodes:map(function(node)
		return node.nodeID
	end):sort():map(function(nodeID)
		assert(tostring(tonumber(nodeID)) == tostringnum(nodeID))
		return tonumber(nodeID)
	end), nil)
	
	local longp = ffi.new('long[1]')
	local allNodesFile 
	if writeOneBigFile then	
		-- I can't json.encode all at once without luajit crapping out, so I will do it piecewise ...
		--local allNodesFile = io.open('nodes.json', 'w')
		--allNodesFile:write'{\n' 
		require 'ffi.c.stdio'
		-- make sure body_t is defined
		ffi.sizeof'body_t'
		allNodesFile = assert(ffi.C.fopen('nodes.raw', 'wb'))

		ffi.cdef[[
typedef struct header_s {
	real mins[3];
	real maxs[3];
	long numNodes;
} header_t;
]]
		local header = ffi.new('header_t[1]')
		for i=0,2 do
			header[0].mins[i] = mins[i+1]
			header[0].maxs[i] = maxs[i+1]
		end
	
		header[0].numNodes = #allNodeIDs
		ffi.C.fwrite(header, ffi.sizeof'header_t', 1, allNodesFile)
		for _,nodeID in ipairs(allNodeIDs) do
			longp[0] = nodeID
			ffi.C.fwrite(longp, ffi.sizeof'long', 1, allNodesFile)
		end
	end

	for _,node in ipairs(allNodes) do
		setmetatable(node, nil)
		node.parent = nil
		node.children = nil
		if writeIndividualNodes then
			local jsonBodies = setmetatable(node.bodies:map(function(body)
				-- TODO store max length of idnumber, name in octree.json
				-- then write out binary data
				local row = table()
				for _,desc in ipairs(rowDesc) do
					local name = desc.name
					if desc.type == 'vec3' then
						row:insert(body[name][0])
						row:insert(body[name][1])
						row:insert(body[name][2])
					elseif desc.type == 'int' 
					or desc.type == 'byte'
					or desc.type == 'float'
					then
						row:insert(tonumber(body[name]))
					elseif desc.type == 'string' then
						row:insert(ffi.string(body[name]))
					else
						error("got unknown type "..desc.type)
					end
				end
				return setmetatable(row, nil)
			end), nil)
			local nodeJSONData = json.encode(jsonBodies):gsub('%],%[','%],\n%[')
			file['nodes/'..tostringnum(node.nodeID)..'.json'] = nodeJSONData
		end
		if writeOneBigFile then
			--allNodesFile:write('\t"',tostring(node.nodeID),'" = ',nodeJSONData,',\n')
			longp[0] = node.nodeID
			ffi.C.fwrite(longp, ffi.sizeof'long', 1, allNodesFile)
			longp[0] = #node.bodies
			ffi.C.fwrite(longp, ffi.sizeof'long', 1, allNodesFile)
			for i,body in ipairs(node.bodies) do
				ffi.C.fwrite(body._ptr, ffi.sizeof'body_t', 1, allNodesFile)
			end
		end
		node.bodies = nil
	end
	--allNodesFile:write'}\n'
	--allNodesFile:close()
	ffi.C.fclose(allNodesFile)

	if writeIndividualNodes then
		print('writing octree info')
		file['octree.json'] = json.encode({
			mins = {mins:unpack()},
			maxs = {maxs:unpack()},
			nodes = allNodeIDs,
		}, {indent=true})
	end

	print('done')
end

return OutputToPoints
