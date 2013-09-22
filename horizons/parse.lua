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
	return unpack(res)
end

function mustbe(pattern)
	local res = {canbe(pattern)}
	if #res == 0 then error('failed to match pattern '..pattern) end
	return unpack(res)
end

local unknownDatas = {}
local db = {}	-- our database

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
				name = a
				parent = b
			end

			local vars = table()

			if name:match('^L%d ') then
				-- L-n lagrangian points don't have any data
			else
				canbe('^ Physical: ... %d%d, %d%d%d%d$')	-- sun has an extra date line 
				canbe('^ Description: ... %d%d, %d%d%d%d')	-- Hegemone / (Jupiter) and Kore (Jupiter) have an extra description line 
				while canbe('^%s+http://') do end	-- Mimas / (Saturn) and Helene / (Saturn)
				mustbe('^%s*$') 	-- then always comes a newline
				while canbe('^%s*$') do end	-- sometimes two
				
				if canbe('^ Solution fit to Voyager and Cassini data%.$') then
					while canbe('^%s*$') do end
				end

				if canbe('^ GEOPHYSICAL DATA')
				or canbe('^ PHYSICAL DATA')
				or canbe('^ ?PHYSICAL PROPERTIES:$')
				or canbe('^ SATELLITE PHYSICAL PROPERTIES:')
				then
					-- variables at 3 and 42
					while true do
						if canbe('^%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*%*$') then break end
						canbe('^%[1%]')	-- skip footnotes
						while canbe('^%s*$') do end	-- stop at empty lines
						local cols = table()
						local numEquals = select(2, nextline:gsub('=', '=')) == 2
						if numEquals  then	-- two columns
							-- TODO Deimos needs to be shifted one way, Venus needs to be shifted the other ...
							cols:insert(nextline:sub(1,41))
							cols:insert(nextline:sub(42))
						elseif numEquals == 1 then
							cols:insert(nextline)
						end
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
						getline()
					end
				else
					-- new discovery / no data
					unknownDatas[nextline] = true
				end
			end
			table.insert(db, {
				id=id,
				name=name,
				parent=parent,
				vars=vars,
			})
			--[[
			print(id, name, parent)
			for k,v in pairs(vars) do
				print('',k,v)
			end
			--]]
		else
			getline()
		end
	end
end, function(err)
	io.stderr:write('line '..lineIndex..'\n')
	io.stderr:write(err..'\n'..debug.traceback()..'\n')
end)

--[[
print('got unknown datas')
for line,_ in pairs(unknownDatas) do
	print(line)
end
--]]

print(json.encode(db))
