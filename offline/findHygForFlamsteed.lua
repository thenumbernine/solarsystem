require 'ext'

local CSV = require 'csv'
print'reading hyg'
local hyg = CSV.file'../hyg/hygdata_v3.csv'
print'done reading hyg'
hyg:setColumnNames(hyg.rows:remove(1))

local num, con = ...
if not con then
	con = num
	num = nil
else
	assert(tonumber(num), "couldn't parse number "..tostring(num))
end
assert(con, "exected [flamsteed#] [constellation-con] , or just the con")


local rows = table()
for _,row in ipairs(hyg.rows) do
	if row.con == con then
		if not num
		or row.flam == num
		then
			rows:insert(row)
		end
	end
end

rows:sort(function(a,b)
	return (tonumber(a.mag) or math.huge) 
		> (tonumber(b.mag) or math.huge) 
end)
for _,row in ipairs(rows) do
	print(tolua(table.mapi(hyg.columns, function(col)
		local value = row[col]
		if value ~= '' then
			return row[col], col
		end
	end)))
end
print('found '..#rows..' rows')

--[[ maybe this is what 'findIndexForExoplanet' is now
local exoplanets = fromlua(file'../exoplanet/openExoplanetCatalog.lua':read())
if namedStars then
	local exoplanetSystemsMatched = 0
	local starIndexForName = table.map(namedStars, function(name,index)
		return index, name
	end):setmetatable(nil)
	for _,system in ipairs(exoplanets) do
		local name = system.name
		-- "14 Her" ... is where in the HYG database?
		-- "24 Sex" ... HYG has 23 and 25, but not 24
		-- "4 UMa"
		-- "1RXS1609" ... ?
		-- "WASP-*" ... ?
		-- "2MASS *" ... ?
		-- "2M *" ... ?
		-- "TOI-*" ... ?
		-- "YSES *" ... ?
		for abbrev,name in pairs(constellationNamesForAbbrevs) do
			name = name:gsub(name, abbrev)
		end
		-- off-spellings of the names:
		name = name
			:gsub('Aquarii', 'Aqr')				-- Aquarius
			:gsub('Cancri', 'Cnc')				-- Cancer
			:gsub('Coronae Borealis', 'CrB')	-- Corona Borealis
			:gsub('Cygni', 'Cyg')				-- Cygnus
			:gsub('Eridani', 'Eri')				-- Eridanus
			:gsub('Indi', 'Ind')				-- Indus
			:gsub('Leonis', 'Leo')				-- Leo I think? and not LeoMinor?
			:gsub('Serpentis', 'Ser')			-- Serpens
			:gsub('Ursae Minoris', 'UMi')
		-- Greek letters
		name = name
			:gsub('Epsilon ', 'Eps ')
			:gsub('eps ', 'Eps ')
			:gsub('gamma ', 'Gam ')
			:gsub('kappa ', 'Kap ')
			:gsub('Omega ', 'Ome ')
			:gsub('alf ', 'Alp ')
			:gsub('eta ', 'Eta ')
			:gsub('mu ', 'Mu ')
			:gsub('nu ', 'Nu ')
			:gsub('tau ', 'Tau ')
			:gsub('xi ', 'Xi ')
			:gsub('omi ', 'Omi ')
		
		local index = starIndexForName[system.name]
		if index then
			exoplanetSystemsMatched = exoplanetSystemsMatched + 1
		end
		print(name, index)
	end
	print('matched '..exoplanetSystemsMatched..' of '..#exoplanets..' openExoplanetDatabase starsystems with HYG stars')
end
--]]
