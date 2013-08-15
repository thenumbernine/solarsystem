require 'socket'
if not json then json = pcall(require,'json') end
if not json then json = pcall(require,'dkjson') end
if not json then error 'failed to find json module' end

-- try to match the astro-phys date on file

print('connecting...')
local conn = assert(socket.connect('horizons.jpl.nasa.gov', 6775))
conn:settimeout(0, 'b')
print('connected!')
socket.sleep(2)	-- stop two seconds for terminal negotiation ... or implement it

local f = io.open('horizons-output.txt', 'w')
function run()
	local function readall()
		socket.sleep(1)
		while true do
			local result, reason = conn:receive('*l')
			if reason == 'closed' then error 'closed' end
			if reason == 'timeout' then return reason end
			print(result, reason)
			if result then
				print(result)
				f:write(result..'\n')
				f:flush()
			end
		end
	end

	local function sendAndRead(cmd)
		conn:send(cmd..'\n')
		readall()
	end

	readall()
	sendAndRead('page')	-- get rid of interactive output

	sendAndRead('sun')

	--[[
	sendAndRead('mb')	-- list major bodies
	sendAndRead('')		-- end final prompt
	--]]

	sendAndRead('q')	-- quit
end
xpcall(run, function(err)
	if err == 'closed' then return end
	io.stderr:write(err..'\n'..debug.traceback()..'\n')
end)

f:close()
conn:close()
