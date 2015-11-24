local class = require 'ext.class'
local OutputMethod = require 'output'

local OutputToJSON = class(OutputMethod)

--[[
args:
	filename
	variableName
--]]
function OutputToJSON:init(args)
	self.filename = assert(args.filename)
	self.variableName = assert(args.variableName)

	self.outputLines = table()
	self.outputLines:insert(self.variableName..' = [')
end

function OutputToJSON:processBody(body)
	self.outputLines:insert('\t'..json.encode(body)..',')
end

function OutputToJSON:done()
	self.outputLines:insert('];')
	file[self.filename] = self.outputLines:concat('\n')
end

return OutputToJSON
