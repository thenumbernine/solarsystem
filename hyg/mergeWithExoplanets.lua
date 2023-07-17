local fromlua = require 'ext.fromlua'
local hygStars = fromlua(path'simbad-data.lua':read())
local hygForOID = {}
for _,hygStar in pairs(hygStars) do
	if hygStar.simbadOIDRef then
		hygForOID[hygStar.simbadOIDRef] = hygStar
	end
end

local json = require 'dkjson'
local exoplanets = json.decode(path'../exoplanet/openExoplanetCatalog.json':read())
local exoplanetForOID = {}
for _,exoplanet in ipairs(exoplanets) do
	if exoplanet.simbadOIDRef then
		exoplanetForOID[exoplanet.simbadOIDRef] = exoplanet

		if hygForOID[exoplanet.simbadOIDRef] then
			print('overlap',exoplanet.simbadOIDRef)
		end
	end
end
