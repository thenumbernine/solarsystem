local file = require 'ext.file'
local json = require 'dkjson'

local hygStars = assert(load('return '..file['simbad-data.lua']))()
local hygForOID = {}
for _,hygStar in pairs(hygStars) do
	if hygStar.simbadOIDRef then
		hygForOID[hygStar.simbadOIDRef] = hygStar
	end
end

local exoplanets = json.decode(file['../exoplanet/openExoplanetCatalog.json'])
local exoplanetForOID = {}
for _,exoplanet in ipairs(exoplanets) do
	if exoplanet.simbadOIDRef then
		exoplanetForOID[exoplanet.simbadOIDRef] = exoplanet

		if hygForOID[exoplanet.simbadOIDRef] then
			print('overlap',exoplanet.simbadOIDRef)
		end
	end
end
