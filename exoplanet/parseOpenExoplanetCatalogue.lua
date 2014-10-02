#!/usr/bin/env lua -lluarocks.require
require 'ext'
local json = require 'dkjson'
local htmlparser = require 'htmlparser'
require 'htmlparser.xpath'
htmlparser.htmlnonclosing = {}
local tree = htmlparser.parse(io.readfile('systems.xml'))
local function buildResults(src)
	-- text node
	if src.child
	and #src.child == 1
	and type(src.child[1]) == 'string'
	then
		if not src.attr then
			return src.tag, assert(src.child[1])
		end
		if #src.attr == 2 then
			return src.tag, {
				[src.attrs[1].name] = src.attrs[1].value,
				[src.attrs[2].name] = src.attrs[2].value,
				value = src.child[1],
			}
		end
		error("unknown attributes on child node")
	end

	-- object/array node
	local results = {}
	if src.attr then
		for _,attr in ipairs(src.attr) do
			assert(not attr.name)
			results[attr.name] = attr.value
		end
	end
	local multiples = {}
	if src.child then
		for _,child in ipairs(src.child) do
			local key, value = buildResults(child)
			if results[key] then	-- multiple results
				multiples[key] = {}
				table.insert(multiples[key], value)
			else
				results[key] = value
			end
		end
	end
	for key,value in pairs(multiples) do
		table.insert(value, 1, results[key])	-- insert first up front
		results[key] = value		-- replace 1st with array
	end
	return src.tag, assert(results)
end
local key, value = buildResults(tree[1])
local results = {}
results[key] = value

io.writefile('openExoplanetCatalog.json', json.encode(results, {indent=true}))
