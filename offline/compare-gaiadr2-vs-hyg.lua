#!/usr/bin/env lua
--[[
not sure if i should put this in universe/offline/ or in solarsystem/hyg/

from https://gea.esac.esa.int/archive/ ADQL query
select * from gaiadr2.hipparcos2_best_neighbour order by original_ext_source_id asc;
then compare 'original_ext_source_id' field to the hyg_v37.csv 'hip' field

summary of results:

119614 stars in HYG
117956 stars in HYG and Hipparcos
83034 stars in Gaia-DR2 and Hipparcos
83033 stars in Gaia-DR2, Hipparcos, and HYG
(Looks like only Hipparcos #54976 is missing from HYG ...
 I hope I didn't delete that line accidentally)
145 named stars in HYG and Hipparcos
46 named stars in Gai-DR2, Hipparcos, and HYG

now what else is interesting is that the log10lum range of hyg is much dimmer than the log10lum range of gaia-dr2
gaia-dr2 seems to cut off just below -1.5

I wonder if the luminosity data between hyg and gaia compares?
--]]
local table = require 'ext.table'
local CSV = require 'csv'

-- TODO HERE
-- load the gaia-dr2 fits file
-- and extract the rows in it associated with (named/all) hipparcos stars
-- and compare stats

print'reading hyg'
local hyg = CSV.file'../hyg/hyg_v37.csv'
hyg:setColumnNames(hyg.rows:remove(1))

print'building hyg lookup for hip'
local hygForHip = table()
local hygForHipCount = 0
for _,row in ipairs(hyg.rows) do
	if row.hip ~= ''
	and row.proper ~= ''
	then
		hygForHip[tonumber(row.hip)] = row
		hygForHipCount = hygForHipCount + 1
	end
end
print('found '..hygForHipCount..' named stars shared between the hipparcos and hyg databases')

print'reading gaia-cross-hipparcos'
local gaia = CSV.file'gaiadr2.hipparcos2_best_neighbour.csv'
gaia:setColumnNames(gaia.rows:remove(1))

-- for each gaia-cross row find the hyg row
-- and then find
print'searching'
gaia_hyg_count = 0
for _,gaiarow in ipairs(gaia.rows) do
	local hip = tonumber(gaiarow.original_ext_source_id)
	hip = assert(tonumber(hip), "failed to convert hipparcos id "..hip.." to number")
	-- looks like not all hipparcos stars (such as #54976) are in the hyg database
	local hygrow = hygForHip[hip]
	if hygrow
	then
		gaia_hyg_count = gaia_hyg_count + 1
		print(hygrow.proper, hygrow.hip, gaiarow.source_id)
	end
end
print('found '..gaia_hyg_count..' named stars shared between the gaia, hipparcos, and hyg databases')
