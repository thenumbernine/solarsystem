local f = assert(loadfile('output_points.template.lua'))
local r = assert(f())
return r
--return assert(load(require 'template'(require 'ext.file'['output_points.lua'])))()
