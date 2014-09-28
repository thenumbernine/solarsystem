local julian = {}

-- testing according to http://aa.usno.navy.mil/cgi-bin/aa_jdconv.pl

function julian.toCalendar(julian)
-- [[ http://www.astro-phys.com/js/astro/api.js
	local jd = 0.5 + julian
	local I = math.floor(jd)
	local F = jd - I
	local B
	if I > 2229160 then
		local A = math.floor((I - 1867216.25) / 36524.25)
		B = I + 1 + A - math.floor(A / 4.0)
	else
		B = I
	end
	local C = B + 1524
	local D = math.floor((C - 122.1) / 365.25)
	local E = math.floor(365.25 * D)
	local G = math.floor((C - E) / 30.6001)
	local d = C - E + F - math.floor(30.6001 * G)
	local month
	if G < 13.5 then
		month = G - 1
	else
		month = G - 13
	end
	local year
	if month > 2.5 then
		year = D - 4716
	else
		year = D - 4715
	end
	local day = math.floor(d)
	local h = (d - day) * 24
	local hour = math.floor(h)
	local mind = math.abs(h - hour) * 60
	local min = math.floor(mind)
	local sec = ((mind - min) * 60)	--math.floor
	return {year=year, month=month, day=day, hour=hour, min=min, sec=sec}
--]]
--[[ http://www.tondering.dk/claus/cal/julperiod.php
	local jd = julian
	local a = jd + 32044
	local b = math.floor((4 * a + 3) / 146097)
	local c = a - math.floor((146097 * b) / 4)
	local d = math.floor((4 * c + 3) / 1461)
	local e = c - math.floor((1461 * d) / 4)
	local m = math.floor((5 * e + 2) / 153)
	local day = e - math.floor((153 * m + 2) / 5) + 1
	local month = m + 3 - 12 * math.floor(m / 10)
	local year = 100 * b + d - 4800 + math.floor(m / 10)
	return {year=year, month=month, day=day, hour=0, min=0, sec=0}
--]]
end

--[[
date has year, month, day, and maybe hour, min, sec
--]]
function julian.fromCalendar(date)
	local Y = assert(date.year)
	local M = assert(date.month)
	local D = assert(date.day)
	local hour = date.hour or 0
	local min = date.min or 0
	local sec = date.sec or 0
	if Y < 0 then Y = Y + 1 end	-- account for 1BC->1AD missing zero year	
	if M == 1 or M == 2 then
		Y = Y - 1
		M = M + 12
	end
	local A = math.floor(Y/100)
	local B = math.floor(A/4)
	local C = 2-A+B
	local E = math.floor(365.25 * (Y + 4716))
	local F = math.floor(30.6001 * (M + 1))
	local JD = C + D + E + F - 1524.5
	JD = JD + (hour + (min + sec / 60) / 60) / 24 
	return JD
end

return julian

--http://www.madore.org/~david/misc/time.html
