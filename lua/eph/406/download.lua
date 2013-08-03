require 'ext'
require 'socket'
require 'socket.ftp'
require 'ltn12'
local url = 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de406/'
for l in io.lines('remotels') do
	local ws = l:split('%s+')
	local fn = ws[1]
	print(fn)
	socket.ftp.get{
		url = url..fn,
		sink = ltn12.sink.file(io.open(fn, 'wb')),
	}
end
