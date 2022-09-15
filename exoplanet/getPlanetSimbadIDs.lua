require 'ext'
local querySimbad = require 'query-simbad'
local json = require 'dkjson'
local catalogue = json.decode(file'openExoplanetCatalog.json':read())
for i,system in ipairs(catalogue.systems) do
	if not system.simbadID then
		local names = system.name:split' / '
		for _,body in ipairs(system.bodies) do names:append(body.name:split' / ') end
		-- we should be able to just use simbad's entry for the system (right)
		for _,name in ipairs(names) do
			-- catalogue has Tau Boötis, simbad has Tau Bootis	
			name = name:gsub('ö', 'o')
			
			-- the catalogue uses "Corona" in place of "Coronae".  I should fix this in the parsing.
			-- 'Epsilon Coronae Borealis' and 'Rho Coronae Borealis' are catalgued correctly (and match simbad and wiki)
			-- 'Omicron Corona Borealis' and 'Kappa Corona Borealis' don't match either (until you add the 'e' at the end)
			name = name:gsub('Corona([^e])', 'Coronae%1')
			
			-- simbad recognizes "UZ Fornacis" but not "UZ Fornacis (ab)"
			if name == 'UZ Fornacis (ab)' then name = 'UZ Fornacis' end

			local results = querySimbad("select oidref,id from ident where id='"..name:gsub("'", "''").."'")
			if results and #results.data > 0 then
				system.simbadOIDRef, system.simbadID = table.unpack(results.data[1])
			end
			if system.simbadOIDRef and system.simbadID then break end
		end
		print(i..'/'..#catalogue.systems..': '..tolua{
			system=system.name,
			oidref=system.simbadOIDRef,
			simbadID=system.simbadID,
		})
	end
end
file'openExoplanetCatalog.json':write(json.encode(catalogue, {indent=true}))
