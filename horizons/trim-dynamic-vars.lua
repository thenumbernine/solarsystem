#!/usr/bin/lua
require 'ext'
local json = require 'myjson'

local staticIntro = 'horizonsStaticData = '
local static = json.decode(file['static-vars.json']:match(staticIntro..'(.*)'))

local dynamicIntro = 'horizonsDynamicData = '
local dynamic = json.decode(file['dynamic-vars.json']:match(dynamicIntro..'(.*)'))


print('#static', #static)
print('#dynamic', #dynamic.coords)

-- TODO duplicate names?
local dynamicKeys = {}
for _,d in ipairs(dynamic.coords) do
	dynamicKeys[d.id] = true
end
local common = {}
for _,s in ipairs(static) do
	if dynamicKeys[s.id] then
		common[s.id] = true
	end
end
local commonList = table.keys(common):sort()
print('#common:', #commonList)
print('common: '..commonList:concat', ')


-- remove SEMB-L (lagrangian points?)
-- these are the id's of the SEMB-L* points
common[31] = nil
common[32] = nil
common[34] = nil
common[35] = nil

for i=#dynamic.coords,1,-1 do
	if not common[dynamic.coords[i].id] then
		print('removing dynamic '..dynamic.coords[i].id..' '..dynamic.coords[i].name)
		table.remove(dynamic.coords, i)
	end
end
for i=#static,1,-1 do
	if not common[static[i].id] then
		print('removing static '..static[i].id..' '..static[i].name)
		table.remove(static, i)
	end
end

file['static-vars.json' ] = staticIntro..json.encode(static, {indent=true})..';'
file['dynamic-vars.json'] = dynamicIntro..json.encode(dynamic, {indent=true})..';'
