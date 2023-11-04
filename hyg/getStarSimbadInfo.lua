#!/usr/bin/env lua
require 'ext'
local querySimbad = require 'query-simbad'

local csv = require'csv'
print('reading csv file...')
local stars = csv.file'hyg_v37.csv'
print('...done reading csv file')
stars:setColumnNames(stars.rows:remove(1))

print('parsing simbad info...')
local simbadInfos = fromlua((path'simbad-data.lua':read())) or {}
local simbadInfosForHipID = {}
local simbadInfosForHDID = {}
print('...done parsing simbad info')

local hipIDsToRequest = table()
local hdIDsToRequest = table()
for _,star in ipairs(stars.rows) do
	local id = assert(tonumber(star.id))
	local simbadInfo = simbadInfos[id]
	if not simbadInfo then
		simbadInfo = {hygID=id}
		simbadInfos[id] = simbadInfo
	end
	
	if #star.hip > 0 then
		local hipID = assert(tonumber(star.hip), "failed to get number value of HIP id "..tostring(star.hip))
		simbadInfo.hipID = hipID
		simbadInfosForHipID[hipID] = simbadInfo 
		if not simbadInfo.simbadID then
			hipIDsToRequest:insert(hipID)
		end
	end

	if #star.hd > 0 then
		local hdID = assert(tonumber(star.hd), "failed to get number for HD id "..tostring(star.hd))
		simbadInfo.hdID = hdID
		simbadInfosForHDID[hdID] = simbadInfo
		if not simbadInfo.simbadID then
			hdIDsToRequest:insert(hdID)
		end
	end
end

print('hipIDs to query:', table.unpack(hipIDsToRequest))
local batchSize = 250
for i=1,#hipIDsToRequest,batchSize do
	print('querying simbad...')
	local query = "select oidref,id from ident where id in ("
		..hipIDsToRequest:sub(i,i+batchSize-1):map(function(id)
			return "'HIP "..id.."'"
		end):concat', '..")"
	local results = assert(querySimbad(query))
	print('...done querying simbad')
	for _,row in ipairs(results.data) do
		local simbadOIDRef, simbadID = table.unpack(row) 
		if simbadOIDRef and simbadID then
			local spaces, hipID = simbadID:match('^HIP(%s+)(%d+)$')
			assert(spaces and hipID)
			assert(#spaces + #hipID == 7)
			hipID = assert(tonumber(hipID))
			local simbadInfo = assert(simbadInfosForHipID[hipID])
			simbadInfo.simbadOIDRef = simbadOIDRef
			simbadInfo.simbadID = simbadID
		else
			--print("query",query,"failed to find info for row???")	-- can this even happen?  it's not like we queried nothing.
		end
	end
end

print('hdIDs to query:', table.unpack(hdIDsToRequest))
local batchSize = 250
for i=1,#hdIDsToRequest,batchSize do
	print('querying simbad...')
	local results = assert(querySimbad("select oidref,id from ident where id in ("
		..hdIDsToRequest:sub(i,i+batchSize-1):map(function(id)
			return "'HD "..id.."'"
		end):concat', '..")"))
	print('...done querying simbad')
	for _,row in ipairs(results.data) do
		local simbadOIDRef, simbadID = table.unpack(row) 
		if simbadOIDRef and simbadID then
			local hdID = assert(simbadID:match('^HD%s+(%d+)$'))
			hdID = assert(tonumber(hdID))
			local simbadInfo = assert(simbadInfosForHDID[hdID])
			simbadInfo.simbadOIDRef = simbadOIDRef
			simbadInfo.simbadID = simbadID
		else
			--print("query",query,"failed to find info for row???")	-- can this even happen?  it's not like we queried nothing.
		end
	end
end

path'simbad-data.lua':write('{' .. table.map(simbadInfos, function(v,k,t) return '['..k..']='..tolua(v), #t+1 end):concat',\n' .. '}')
