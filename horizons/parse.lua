#!/usr/bin/env lua
--[[
I was smart enough to capture the output of the telnet session
but not smart enough to save the data in individual files per-object
oh well, here I will attempt to split the data out of the one big file
--]]

require 'ext'
local json = require 'dkjson'

local data = io.readfile('horizons.txt'):trim()
data = data:gsub('\r\n', '\n'):gsub('\r', '\n')
local lines = data:split('\n')

local lineIndex
local readline = coroutine.wrap(function()
	for i,line in ipairs(lines) do
		lineIndex = i
		coroutine.yield(line)
	end
end)

-- single element lookahead buffer
local nextline
function getline()
	nextline = readline()
end

function canbe(pattern)
	local res = {nextline:match(pattern)}
	if #res > 0 then getline() end
	return table.unpack(res)
end

function mustbe(pattern)
	local res = {canbe(pattern)}
	if #res == 0 then error('failed to match pattern '..pattern) end
	return table.unpack(res)
end

local unknownDatas = {}
local db = {}	-- our database
local allkeys = table()

planetInfos = table()
xpcall(function()
	getline()
	while nextline do 
		local id = canbe('^Horizons>%[%[%[my input:(%d+)%]%]%]$')
		if id then
			mustbe('^ '..id..'$')
			mustbe('^%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*$')
			local name = canbe('^ ?Revised ?: ... %d%d, %d%d%d%d%s+(.*)%s+     ')
						or canbe('^JPL/HORIZONS%s+(.*)%s+     ')
			name = name:trim()
			--[[
			name could be in one of the following forms:
				name
				name Barycenter
				name / (parent)
			--]]
			local parent
			local a, b = name:match('^(.*)%s*/%s*%((.*)%)$')
			if a and b then
				name = a:trim()
				parent = b
			end

			local vars = table()
			
			planetInfos:insert{
				name = name,
				parent = parent,
				id = id,
				vars = vars,
			}

			if name:match('^L%d ') then
				-- L-n lagrangian points don't have any data
			else
				canbe('^ Description: ... %d%d, %d%d%d%d')	-- Hegemone / (Jupiter) and Kore (Jupiter) have an extra description line 
				while canbe('^%s+http://') do end	-- Mimas / (Saturn) and Helene / (Saturn)
				while canbe('^%s*$') do end	-- maybe some spaces 
				
				if canbe('^ Solution fit to Voyager and Cassini data%.$') then
					while canbe('^%s*$') do end
				end

				if canbe('^ GEOPHYSICAL DATA')
				or canbe('^ PHYSICAL DATA')
				or canbe('^ ?PHYSICAL PROPERTIES[^\n]*:$')
				or canbe('^ SATELLITE PHYSICAL PROPERTIES:')
				then
					local varlines = table()
					while true do
					
						if canbe('^%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*?$') then break end
						
						canbe('^%[1%]')	-- skip footnotes
						
						while canbe('^%s*$') do end	-- stop at empty lines
				
						-- multiple headers ... i could organize this parser better, I know
						if canbe('^ GEOPHYSICAL DATA')
						or canbe('^ PHYSICAL DATA')
						or canbe('^ ?PHYSICAL PROPERTIES:$')
						or canbe('^ SATELLITE PHYSICAL PROPERTIES:')
						or canbe('^ ?DYNAMICAL CHARACTERISTICS:$')
						or canbe('^ SATELLITE ORBITAL DATA')
						then
	
						-- special variable formatting
						elseif nextline:match('^%s*Motn%.') then
							vars:insert{nextline}
							getline()
							vars:insert{nextline}
						else
							varlines:insert(nextline)
						end
						getline()
					end

					local mid = 1
					for _,line in ipairs(varlines) do
						for _,w in ipairs{'Mass', 'Density'} do
							local pos = line:find(w)
							if pos then
								mid = math.max(mid, pos)
							end
						end
					end

					for _,line in ipairs(varlines) do
						
						local cols = table()
						cols:insert(line:sub(1,mid-1))
						cols:insert(line:sub(mid))
						for _,col in ipairs(cols) do
							local equalsLoc = col:find('=',nil,true)
							if equalsLoc then	-- otherwise got a non-parameter line
								local key = col:sub(1, equalsLoc-1):trim()
								local value = col:sub(equalsLoc+1):trim()
								vars[key] = value
							else
								vars:insert(col)
							end
						end
					end
				else
					-- new discovery / no data
					unknownDatas[nextline] = true
				end
			end

			print(id, name, parent)
			for k,v in pairs(vars) do
				print('',k,v)
				if type(k) == 'string' then
					allkeys[k] = true
				elseif type(v) == 'string' then
					allkeys[v] = true
				end
			end
		else
			getline()
		end
	end
end, function(err)
	io.stderr:write('line '..lineIndex..'\n')
	io.stderr:write(err..'\n'..debug.traceback()..'\n')
end)

print('got unknown datas')
for line,_ in pairs(unknownDatas) do
	print(line)
end

print('keys')
for k,v in pairs(allkeys) do print(k,v) end
print('sorted')
for _,k in ipairs(allkeys:keys():sort()) do
	print(k)
end


--[[
extract the variables we want
Mass
Radius
Equatorial Radius
Flattening
--]]
for _,planetInfo in ipairs(planetInfos) do
	planetInfo.id = tonumber(planetInfo.id)
	local vars = planetInfo.vars
	--print()
	for k,v in pairs(planetInfo.vars) do
		--print('(vars)','',k,v)
		if type(v) == 'string' then
			local plusminus = v:find('%(%d*%+%-',nil)
			if plusminus then
				v = v:sub(1,plusminus-1):trim()
			end
			local plusminus = v:find('+-',nil,true)
			if plusminus then
				v = v:sub(1,plusminus-1):trim()
			end
		end
		if type(k) == 'string' and v ~= '' then
			if k:match('^Mass %(10%^%d+ k?g ?%)$') then
				local scale, kilo = k:match('^Mass %(10%^(%d+) (k?)g ?%)$')
				scale = tonumber(scale)
				-- our scale is kg, so if it's in g's then bring it down
				if kilo ~= 'k' then scale = scale - 3 end
				
				local v2, scale2 = v:match('^([%.%d]+) %(10%^(%-?%d+)%)$')
				if v2 and scale2 then
					scale2 = tonumber(scale2)
					planetInfo.mass = tonumber(v2) * 10^(scale + scale2)
				else
					planetInfo.mass = tonumber(v) * 10^scale
				end
			elseif k:match('^Mass, 10%^%d+ kg$') then
				local scale = k:match('^Mass, 10%^(%d+) kg$')
				scale = tonumber(scale)
				planetInfo.mass = tonumber(v) * 10^scale
			end
			if k == 'Radius (photosphere)' then
				local v2, scale = v:match('^([.%d]+)%(10%^(%d+)%) km')
				v2 = tonumber(v2)
				scale = tonumber(scale)
				planetInfo.radius = v2 * 10^(scale+3)
			end
			if k:lower() == 'mean radius (km)' 
			or k == 'Mean radius, km'
			or k == 'Radius (km)'
			or k == 'Radius, km'
			or k == 'Volumetric mean radius'
			then
				local x,y,z = v:match('^([.%d]+) ?x ?([.%d]+) ?x ?([.%d]+)$')
				x = tonumber(x)
				y = tonumber(y)
				z = tonumber(z)
				if x and y and z then
					x = x * 1000
					y = y * 1000
					z = z * 1000
					planetInfo.radius = (x + y + z) / 3	-- can't think of a better way to do this 
				else
					local x,y = v:match('^([.%d]+) ?x ?([.%d]+)')
					x = tonumber(x)
					y = tonumber(y)
					if x and y then
						x = x * 1000
						y = y * 1000
						planetInfo.radius = (x + y) / 2
					else
						planetInfo.radius = tonumber(v) * 1000
					end
				end
			end
			if k:match('^Equatorial Radius') then
				local v2 = v:match('^(.*) km$')
				v2 = tonumber(v2)
				v2 = v2 * 1000	--km
				planetInfo.equatorialRadius = tonumber(v2)
			end
			
			--[[ I'm not sure if this is orbital or surface eccentricity
			-- eccentricity^2 = (2 * inverse flattening - 1) / inverse flattening^2
			-- eccentricity^2 inverse flattening^2 = 2 * inverse flattening - 1
			-- eccentricity^2 inverse flattening^2 - 2 * inverse flattening + 1 = 0
			-- inverse flattening = (2 +- sqrt(4 - 4 eccentricity^2)) / (2 eccentricity^2)
			-- inverse flattening = 1/eccentricity^2  +-  sqrt(1 - eccentricity^2) / eccentricity^2
			-- inverse flattening = eccentricity^-2  +-  sqrt(1 - eccentricity^2) / eccentricity^2
			if k:match('^Eccentricity') then
				if v:sub(1,1) == '~' 
				or v:sub(1,1) == '<'
				then 
					v = v:sub(2):trim() 
				end
				if v:sub(#v) == 'o' then v = v:sub(1,#v-1) end
				if v:sub(#v) == 'R' then v = v:sub(1,#v-1) end
				local e = tonumber(v)
				planetInfo.inverseFlattening = (1 - math.sqrt(1 - e^2)) / e^2
			end
			--]]
			if k:match('^Flattening') then
				local inv = v:match('^1/(.*)$')
				if inv then
					planetInfo.inverseFlattening = tonumber(inv)
				else
					planetInfo.inverseFlattening = 1 / tonumber(v)
				end
			end
		end
	end
	--[[
	for k,v in pairs(planetInfo) do
		if k ~= 'vars' then
			print(k,v)
		end
	end
	--]]
end

--remove that straggler that doesn't show up in the dynamic data
if planetInfos[#planetInfos].id == 1000041 then
	planetInfos:remove(#planetInfos)
end

local planetInfoForName = {}
for _,planetInfo in ipairs(planetInfos) do
	planetInfoForName[planetInfo.name] = planetInfo
end
-- append some facts I found on another site
-- http://solarsystem.nasa.gov/planets/profile.cfm?Object=Sun
planetInfoForName.Sun.radius= 695508000
planetInfoForName.Sun.mass= 1.9891e+30
planetInfoForName.Mercury.radius = 2.4397e+6
planetInfoForName.Mercury.mass = 3.30104e+23
planetInfoForName.Venus.radius = 6.0581e+6
planetInfoForName.Venus.mass = 4.86732e+24
planetInfoForName.Earth.radius = 6.37101e6
planetInfoForName.Earth.mass = 5.97291e+24
planetInfoForName.Moon.radius = 1.73753e+6
planetInfoForName.Moon.mass = 7.34767309245735e+22
planetInfoForName.Mars.radius = 3.3895e+6
planetInfoForName.Mars.mass = 6.41693e+23
--planetInfoForName.Phobos.radius =
planetInfoForName.Amalthea.mass = 2068162437674130000 
planetInfoForName.Himalia.mass = 6744007948937370000 
planetInfoForName.Elara.mass = 869227691196372000 
planetInfoForName.Pasiphae.mass = 299733686619439000
planetInfoForName.Sinope.mass = 74933421654859700
planetInfoForName.Lysithea.mass = 62944074190082100

print('parsed '..#planetInfos..' planets')
io.writefile('static-vars.json', 'horizonsStaticData = '..json.encode(planetInfos, {indent=true})..';')

