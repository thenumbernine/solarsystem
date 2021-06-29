#!/usr/bin/env luajit
require 'ext'

--[[
build constellation lines from descriptions of stars and from hyg database
--]]


local constellationNameForAbbrev = fromlua(file['constellationNamesForAbbrevs.lua'])


-- [[ debugging - what info we have for
local CSV = require 'csv'
print'reading hyg'
local hyg = CSV.file'../../solarsystem/hyg/hygdata_v3.csv'
print'done reading hyg'
hyg:setColumnNames(hyg.rows:remove(1))


local hygForCons = {}
local rowForConsAndBayer = {}
local rowForProper = {}

-- namedStars is full of either proper, or bf, whichever exists
-- flam == some number, idk what
-- bayer == greek letter coinciding with order of brightness
-- bf = the three combined
-- proper = name, separate of constellation
for _,row in ipairs(hyg.rows) do
	if row.cons ~= '' then
		if row.bayer ~= '' then
			local bayer = row.bayer:lower()
			rowForConsAndBayer[row.con] = rowForConsAndBayer[row.con] or {}
			rowForConsAndBayer[row.con][bayer] = row
		end
		hygForCons[row.con] = hygForCons[row.con] or table()
		hygForCons[row.con]:insert(row)
	end
	if row.proper ~= '' then
		rowForProper[row.proper] = row
	end
end
for _,k in ipairs(table.keys(hygForCons)) do
	hygForCons[k]:sort(function(a,b)
		return (tonumber(a.mag) or math.huge) > (tonumber(b.mag) or math.huge)
	end)
end

for _,row in ipairs(hygForCons.And) do
	print(table{row.mag, row.bayer, row.bf, row.proper, row.flam, row.con}:mapi(function(x)
		return ('%q'):format(x)
	end):concat'\t')
end


local constellationLines = {
	And = {
		"alpha delta beta gamma-1",
		"beta mu nu phi",	-- , some unknown star},
		"beta pi",
		"nu zeta epsilon delta pi iota kappa lambda",
		"iota omicron",
	},
	Ant = {
		"iota alpha epsilon",
	},
	Aps = {
		"alpha delta beta gamma",
	},
	Aqr = {
		"epsilon mu beta alpha pi zeta gamma alpha theta rho lambda phi psi-2",	-- something
		-- psi-2 something2
		"iota beta",
		"zeta eta",
		"lambda tau delta psi-2",
	},
	Aql = {
		"beta alpha gamma delta lambda iota theta eta delta zeta lambda",
		"zeta epsilon",
		-- lambda something 
	},
	Ara = {
		"beta alpha kappa epsilon-1 zeta eta delta gamma",
	},
	Ari = {
		"gamma beta alpha",	-- something
	},
	Aur = {
		"alpha beta pi delta alpha epsilon zeta",
		"alpha eta iota",	-- something maybe in another constellation
		"beta theta", 	-- and then same star as above
	},
	Boo = {
		"zeta alpha eta tau",
		"lambda theta kappa lambda gamma rho alpha epsilon delta beta gamma",
	},
	Cae = {
		"delta alpha beta gamma",
	},
	Cam = {
		"beta alpha gamma",
		-- beta something
		-- alpha something2 something3 
		-- gamma something4 something5
	},
	Cnc = {
		"alpha delta gamma iota",
		"beta delta",
	},
	CVn = {
		"alpha beta",
	},
	CMa = {
		"alpha beta nu-2 omicron-1 epsilon sigma delta omicron-2 pi alpha iota gamma theta iota",
		"delta eta",
	},
	CMi = {
		"alpha beta",
	},
	Cap = {
		"alpha-1 beta rho omicron psi omega zeta epsilon delta gamma iota theta tau nu alpha-1",
	},
	Car = {
		"alpha beta omega theta",
		"theta eta", 	-- theta connects to unknown then to eta
		"theta iota",	-- same unknown connects to iota
		-- iota connects to some others out of constellation
		"iota epsilon chi",	-- then so does chi
		-- and that connects back to alpha
	},
}
-- TODO just steal from stellarium, i guess that's what the other cool kids do
for constellationAbbrev,lines in pairs(constellationLines) do
	local rowForBayer = assert(rowForConsAndBayer[constellationAbbrev])
	for i=1,#lines do
		line = lines[i]:split'%s+'
		for j=1,#line do
			local name = line[j]
			local suffix = name:match'(%-%d)$'
			name = name:sub(1,3):lower() .. (suffix or '')
			local row = assert(rowForBayer[name], "failed to find "..name)
			line[j] = row.id
		end
		lines[i] = line
	end
	print(constellationAbbrev)
	print(tolua(lines))
end
