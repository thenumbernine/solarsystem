local file = require 'ext.file'

return assert(loadfile('output_points.template.lua'))()
--return assert(load(require 'template'(file['output_points.lua'])))()
