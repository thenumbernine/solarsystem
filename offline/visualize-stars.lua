#!/usr/bin/env luajit
-- simple pointcloud visualizer
-- TODO move this to solarsystem project, so i can keep building off of it, and comparing it to the solarsystem renderer, and eventually add in things like the exoplanets
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local string = require 'ext.string'
local math = require 'ext.math'
local file = require 'ext.file'
local fromlua = require 'ext.fromlua'
local ffi = require 'ffi'
local template = require 'template'
local Image = require 'image'
local gl = require 'gl'
local glreport = require 'gl.report'
local GLElementArrayBuffer = require 'gl.elementarraybuffer'
local GLProgram = require 'gl.program'
local GLFBO = require 'gl.fbo'
local GLTex2D = require 'gl.tex2d'
local GLHSVTex = require 'gl.hsvtex'
local GLArrayBuffer = require 'gl.arraybuffer'
local CLEnv = require 'cl.obj.env'
local vec3f = require 'vec-ffi.vec3f'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'
local vector = require 'ffi.cpp.vector'
local matrix_ffi = require 'matrix.ffi'
matrix_ffi.real = 'float'	-- default matrix_ffi type

--[[
set = directory to use if you don't specify individual files:

pointfile = filename for our point data
namefile = lua file mapping indexes to names of stars
consfile = lua file containing constellation info for star points ... not used anymore ...

format = whether to use xyz or our 9-col format ... not used anymore ...

-- TODO do this per-dataset, and allow multiple datasets (like gaia + hyg + exo)
lummin = set this to filter out for min value in solar luminosity
lumhicount = show only this many of the highest-luminosity stars
appmaghicount = show only this many of the highest apparent magnitude stars
rlowcount = show only this many of the stars with lowest r
rmax = filter out points further than this value
print = set this to print out all the remaining points
nolumnans = remove nan luminosity ... not needed, this was a mixup of abs mag and lum, and idk why abs mag had nans (maybe it too was being calculated from lums that were abs mags that had negative values?)

addsun = add our sun to the dataset.  I don't see it in Gaia.
addexo = add exoplanets.  this requires testing by the hyg id, which is the hyg 9th param, so this will not work with gaia right now.
nounnamed = remove stars that don't have entries in the namefile
buildnbhds = build neighbhoods
buildvels = build velocity field lines

showStarNames = whether to show star names
--]]
local cmdline = require 'ext.cmdline'(...)

--local set = cmdline.set or 'gaia'
local set = cmdline.set or 'hyg'

-- assign defaults first, then override with cmdline
local pointfile = '../'..set..'/stardata.f32'
local namefile = '../'..set..'/namedStars.lua'

if cmdline.pointfile then
	pointfile = cmdline.pointfile
	-- AND invalidate any default values (like the star names etc)
	namefile = nil
end
if cmdline.namefile then
	namefile = cmdline.namefile
end

-- lua table mapping from index to string
-- TODO now upon mouseover, determine star and show name by it
local namedStars
if namefile then
	local namedata = file[namefile]
	if namedata then
		namedStars = fromlua(namedata)
	end
end

local constellationNamesForAbbrevs = fromlua(file['../constellations/constellationNamesForAbbrevs.lua'])
local constellationAbbrevsForNames = table.map(constellationNamesForAbbrevs, function(name, abbrev)
	return abbrev, name
end):setmetatable(nil)

local constellations = table.map(constellationNamesForAbbrevs, function(name, abbrev, t)
	return {
		name = name,
		abbrev = abbrev,
		enabled = false,
	}, #t+1
end)

-- TODO make sure nothing is indexing constellations (like the HYG-generated star data was)
constellations = constellations:sort(function(a,b)
	return a.name:lower() < b.name:lower()
end)

local constellationForAbbrev = table.mapi(constellations, function(cons)
	return cons, (cons.abbrev or 'unnamed')
end):setmetatable(nil)

-- constellation lines are in hip, so use this to convert them to hyg
local indexForHip = fromlua(file['../hyg/index-for-hip.lua'])

for name, lines in pairs(fromlua(file['../constellations/constellation-lines.lua'])) do
	local abbrev = table.find(constellationNamesForAbbrevs, nil, function(name2)
		return name2:gsub(' ', '') == name
	end)
	assert(abbrev, "couldn't find abbrev for lines name "..name)
	local constellation = assert(constellationForAbbrev[abbrev], "failed to find constellation for abbrev "..abbrev)
	
	-- remap lines from hipparcos index to stardata.f32 index using index-for-hip.lua
	constellation.lines = lines
	for _,line in ipairs(lines) do
		for j=#line,1,-1 do
			local index = indexForHip[line[j]]
			if not index then
				table.remove(line, j)
			else
				line[j] = index
			end
		end
	end
end


local App = class(require 'glapp.orbit'(require 'imguiapp'))
local ig = require 'imgui'		-- windows bug, gotta include ig after imguiapp (or after imgui?)
	
local _1_log_10 = 1 / math.log(10)

App.title = 'pointcloud visualization tool'
App.viewDist = 5e-4

--[[
my gaia data:
	float pos[3]
	float vel[3]
	float luminosity
	float temperature
	float radius
	TODO uint64_t source_id
my hyg data:
	float pos[3]
	float vel[3]
	float luminosity
	float temperature
	float constellationIndex
--]]
ffi.cdef[[
typedef struct {
	vec3f_t pos;
} pt_3col_t;

typedef struct {
	vec3f_t pos;		// in Pc
	vec3f_t vel;
	
	/*
	solar luminosity
	https://en.wikipedia.org/wiki/Solar_luminosity
	3.828e+26 watts
	total emitted over the whole surface
	L_sun = 4 pi k I_sun A^2
		for L_sum solar luminosity (watts)
		and I_sun solar irradiance (watts/meter^2)
	
	https://en.wikipedia.org/wiki/Luminosity
	
	flux = luminosity / area
	area of sphere = 4 pi r^2
	area normal with observer = 2 pi r
	flux = luminosity / (2 pi r)
	
	absolute magnitude:
	MBol = -2.5 * log10( LStar / L0 ) ~ -2.5 * log10 LStar + 71.1974
	MBolStar - MBolSun = -2.5 * log10( LStar / LSun )
	absolute bolometric magnitude of sun: MBolSun = 4.74
	apparent bolometric magnitude of sun: mBolSun = -26.832
	(negative = brighter)

	*/
	float lum;

	/*
	temperature, in Kelvin
	*/
	float temp;

	//TODO
	//in HYG I'm actually outputting the constellation index here
	// which I'm rebuilding later from the 'constellations' file ...
	/// sooo ... dont' rebuild it? just make this 'constellation'
	//and as far as radius ... HYG has r1 and r2 for spheroid sizes
	//Gaia has 'radius' as well
	//how many valid columns in both?
	// I should add on 'r1' and 'r2' .. or just 'radius' ... for both
	// and make this 11 cols
	//radius, probably in sun-radii
	//alright, I'm using the hyg id for the hyg database here.  no need for constellation.
	//and maybe i'm using it to cross-correlate with the exoplanet database
	float radius;

} pt_9col_t;
]]

assert(ffi.sizeof'pt_3col_t' == 3 * ffi.sizeof'float')
assert(ffi.sizeof'pt_9col_t' == 9 * ffi.sizeof'float')

--[[
local format = cmdline.format or pointfile:match'%-(.*)%.f32$' or '3col'
local pt_t = ({
	['3col'] = 'pt_3col_t',
	['9col'] = 'pt_9col_t',
})[format] or error("couldn't deduce point type from format")
--]]
local pt_t = 'pt_9col_t'

-- 2x
local gpuVelLineBuf
local cpuVelLineBuf


ffi.cdef[[
typedef struct {
	vec3f_t pos;	//pos of each endpoint of the line
	float dist;		//distance between.  will be doubled up in pairs of vertexes for each line
} nbhd_t;
]]

local cpuNbhdLineBuf
local gpuNbhdLineBuf

local env	-- I'm not using CL at the moment
local drawIDShader
local accumStarPointShader
local accumStarLineShader
local renderAccumShader
local drawNbhdLineShader

local fbo
local fbotex

-- black body color table from http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color_D58.html
local colorTempMin = 1000
local colorTempMax = 40000
local tempTex

local LSun = 3.828e+26 	-- Watts
local L0 = 3.0128e+28	-- Watts
local LSunOverL0 = LSun / L0

-- set to true to put some rough connections between neighboring stars
-- which TODO looks like soup until I dist-atten it somehow
local buildNeighborhood = cmdline.buildnbhds
local buildVelocityBuffers = cmdline.buildvels

-- _G so that ig.luatableSliderFloat can use them
starPointAlpha = 1
starPointSizeScale = 3
starPointSizeBias = -3
lineVelAlpha = .5
hsvRange = .2
hdrScale = .001
hdrGamma = 1
bloomLevels = 0
showDensity = false
lineVelScalar = 1
drawPoints = true
drawVelLines = false
normalizeVel = false
showPickScene = false
showInformation = true			-- show information on the current star being hovered over
drawGrid = true			-- draw a sphere with 30' segments around our orbiting star
tiltGridToPlanet = false
gridRadius = 100
showNeighbors = false	-- associated with the initialization flag buildNeighborhood
nbhdLineAlpha = 1
showConstellations = cmdline.showConstellations or false
sliceRMin = 0
sliceRMax = math.huge
sliceLumMin = 0
sliceLumMax = math.huge
showStarNames = true
if cmdline.showStarNames ~= nil then showStarNames = cmdline.showStarNames end

--[[
picking is based on point size drawn
point size is based on apparent magnitude
but it's still depth-sorted
so weaker closer stars will block brighter further ones
a better way would be drawing the pick ID + app mag to a buffer
and then post-process flooding it outward to the appropriate size based on its app mag
that way brighter stars would overbear closer weaker neighbors
--]]
pickSizeBias = 4


local starTex

local numPts	-- number of points
local gpuPointBuf, cpuPointBuf

function App:initGL(...)
	App.super.initGL(self, ...)

	gl.glDisable(gl.GL_DEPTH_TEST)

	self.view.znear = 1e-3
	self.view.zfar = 1e+6

--	self.view.angle = (quatd():fromAngleAxis(0, 0, 1, 90) * self.view.angle):normalize()

	local data = file[pointfile]
	numPts = #data / ffi.sizeof(pt_t)
print('loaded '..numPts..' stars...')
--numPts = math.min(numPts, 100000)
	
	local src = ffi.cast(pt_t..'*', ffi.cast('char*', data))
	
	--[[ I would just cast, but luajit doesn't seem to refcount it, so as soon as 'data' goes out of scope, 'cpuPointBuf' deallocates, and when I use it outside this function I get a crash ....
	cpuPointBuf = ffi.cast(pt_t..'*', s)
	--]]
	-- [[
	cpuPointBuf = ffi.new(pt_t..'[?]', numPts)
	for i=0,numPts-1 do
		cpuPointBuf[i] = ffi.new(pt_t, src[i])
	end
	--]]

	if cmdline.lummin
	or cmdline.lumhicount
	or cmdline.appmaghicount
	or cmdline.rlowcount
	or cmdline.rmax
	or cmdline.nolumnans
	or cmdline.addsun
	or (cmdline.nounnamed and namedStars)
	or cmdline.addexo
	then
print'filtering some of our data...'
print'converting from binary blob to lua table...'
		local pts = table()
		for i=0,numPts-1 do
			-- keep track of the original index, for remapping the name file
			pts:insert{obj=ffi.new(pt_t, cpuPointBuf[i]), index=i}
		end
print'...done converting'

		-- do this first, while namedStars <-> objs, before moving any objs
		if cmdline.nounnamed then
print'filtering out unnamed...'
			assert(namedStars, "you can't filter out unnamed stars if you don't have a name file")
			for i=numPts-1,0,-1 do
				if not namedStars[i] then
					pts:remove(i+1)
				end
			end
			print('nounnamed filtered down to '..#pts)
		end

		if cmdline.nolumnans then
print'filtering out luminosity nans...'
			pts = pts:filter(function(pt)
				return math.isfinite(pt.obj.lum)
			end)
			print('nolumnans filtered down to '..#pts)
		end

		if cmdline.addsun then
print'adding sun point...'
			local sun = ffi.new(pt_t)
			sun.pos:set(0,0,0)
			-- vel is in Pc/year?  double check plz.
			-- Gaia velocity is in solar reference frame, so the sun's vel will be zero
			sun.vel:set(0,0,0)
			sun.lum = 1
			sun.temp = 5772	-- well, the B-V is 0.63.  maybe I should just be storing that?
			sun.radius = 1	-- in solar radii
			
			pts:insert{obj=sun, index=numPts}
			numPts = numPts + 1	-- use a unique index.  don't matter about modifying numPts, we will recalculate numPts soon
		end

		if cmdline.lummin then
print('filtering out luminosity min '..cmdline.lummin)
			pts = pts:filter(function(pt)
				return pt.obj.lum >= cmdline.lummin
			end)
			print('lummin filtered down to '..#pts)
		end
		if cmdline.lumhicount then
print('filtering out only the highest '..cmdline.lumhicount..' luminous objects...')
			pts = pts:sort(function(a,b)
				return a.obj.lum > b.obj.lum
			end):sub(1, cmdline.lumhicount)
			print('lumhicount filtered down to '..#pts)
		end
		if cmdline.appmaghicount then
print('filtering out only the highest '..cmdline.appmaghicount..' apparent magnitude objects...')
			pts = pts:sort(function(a,b)
				return  a.obj.lum / a.obj.pos:lenSq()
						> b.obj.lum / b.obj.pos:lenSq()
			end):sub(1, cmdline.appmaghicount)
			print('appmaghicount filtered down to '..#pts)
		end
		if cmdline.rlowcount then
print('filtering out only the closest '..cmdline.rlowcount..' objects...')
			pts = pts:sort(function(a,b)
				return  a.obj.pos:lenSq() < b.obj.pos:lenSq()
			end):sub(1, cmdline.rlowcount)
			print('rlowcount filtered down to '..#pts)
		end
		if cmdline.rmax then
print('filtering out only objects of distance '..cmdline.rmax..' or less...')
			pts = pts:filter(function(pt)
				return pt.obj.pos:length() <= cmdline.rmax
			end)
			print('rmin filtered down to '..#pts)
		end
	
		if cmdline.addexo then
print('adding exoplanets not done yet')
			--[[
			now load the open exoplanet catalog parsed data
			map from exoplanet system name to hyg name
			if there's no hyg number associated then add the system
			notice, the exoplanet catalog provides
			visualMagnitude per-body, sometimes
			temperature per-body sometimes
			radius per-body sometimes, maybe more often than not
			J K H magnitudes.
			in absense of visual magnitude, should I derive it from the J K H magnitudes?
			or should I use the radius+temp equation:
			luminosity = sigma area temperature^4, for area = 4 pi radius^2, sigma = Stefan-Boltzmann constant 5.670374419e-8 W / (m^2 K^4)
			--]]
		end

		if namedStars then
print('remapping old named stars...')
			local newnames = namedStars and {} or nil
			for i,pt in ipairs(pts) do
				-- pts will be 0-based
				if namedStars then
					newnames[i-1] = namedStars[pt.index]
				end
			end
			if namedStars then
				namedStars = newnames
			end
print('...done remapping old named stars')
		end

		if cmdline.print then
print('printing radius and luminosity...')
			for _,pt in ipairs(pts) do
				local r = pt.obj:length()
				print('r='..r..' lum='..pt.obj.lum)
			end
print('...done printing radius and luminosity')
		end

		numPts = #pts
print('allocating new binary blob...')
		cpuPointBuf = ffi.new(pt_t..'[?]', numPts)
print('copying from lua table to new binary blob...')
		for i=0,numPts-1 do
			cpuPointBuf[i] = pts[i+1].obj
		end
print('...done copying from lua table to new binary blob')
	end

	env = CLEnv{precision='float', size=numPts, useGLSharing=false}
	local real3code = template([[
<?
for _,t in ipairs{'float', 'double'} do
?>
typedef union {
	<?=t?> s[3];
	struct { <?=t?> s0, s1, s2; };
	struct { <?=t?> x, y, z; };
} _<?=t?>3;
typedef union {
	<?=t?> s[4];
	struct { <?=t?> s0, s1, s2, s3; };
	struct { <?=t?> x, y, z, w; };
} _<?=t?>4;
<?
end
?>
typedef _<?=env.real?>3 real3;
typedef _<?=env.real?>4 real4;
]], {
		env = env,
	})

	ffi.cdef(real3code)
	env.code = env.code .. real3code
	
	-- get range
	local s = require 'stat.set'('x','y','z','r', 'v', 'lum', 'temp', 'log10lum')
--[[
HYG:
r = {min = 0, max = 990.09945066522, avg = 245.82425693605, sqavg = 96753.244535951, stddev = 190.58772058502, count = 109399},
v = {min = 0, max = 0.0017771109717193, avg = 3.619610793339e-05, sqavg = 2.6188731301195e-09, stddev = 3.6176164813229e-05, count = 109399},
lum = {min = 1.2257446542208e-06, max = 67483.875, avg = 65.549295330161, sqavg = 231284.10516467, stddev = 476.43194167308, count = 109399},
temp = {min = 1499.3383789063, max = 21707.421875, avg = 6182.1188766896, sqavg = 42211879.557415, stddev = 1998.3207329886, count = 109399},
log10lum = {min = -5.9115999921446, max = 4.8292000123107, avg = 1.0959642811746, sqavg = 2.1518952720112, stddev = 0.97506798039962, count = 109399},

Gaia:
r = {min = 3.2871054403464, max = 2121788012.1646, avg = 3722.3396439795, sqavg = 913550045428.65, stddev = 955790.87127688, count = 7112113},
v = {min = 1.2907013896885e-07, max = 86.409967259077, avg = 0.00012526101936957, sqavg = 0.0010910844398197, stddev = 0.033031329817262, count = 7112113},
lum = {min = 0.030576450750232, max = 96606.6328125, avg = 51.305569407229, sqavg = 29367.297608038, stddev = 163.50852013225, count = 6081418},
temp = {min = 3229, max = 9715.6669921875, avg = 4822.2554139669, sqavg = 23954513.485383, stddev = 836.87884896789, count = 7102644},
log10lum = {min = -1.5146129279835, max = 4.9850069452041, avg = 1.0599897393166, sqavg = 1.8131122068599, stddev = 0.83038181543398, count = 6081418},

Seems with Gaia I am getting only a small subset of the luminosity range that I am with the HYG database.  I wonder why?  Maybe luminosity is one of the things they calculated last?  And I should just be deriving it from the magnitude?
--]]

-- [=[ why is this crashing intermittantly for gaia data?
-- maybe something to with how, when it doesn't crash, luajit has a 'not enough memory' error shortly after
	--local log10lumbin = require 'stat.bin'(-10, 10, 200)
	-- TODO better way of setting these flags ... in ctor maybe?
	for _,s in ipairs(s) do
		s.onlyFinite = true
	end
	--[[
	in Mpc:
	r min = 3.2871054551763e-06
	r max = 2121.7879974726 Pc
	r avg = 0.0037223396472864
	r stddev = 0.95579087454634
	so 3 sigma is ... 3
	
	SagA* = 8178 Pc from us, so not in the HYG or Gaia data

	Proxima Centauri is 1.3 Pc from us

	lum stddev is 163, so 3 stddev is ~ 500
	--]]
print'calculating stats on data...'
	for i=0,numPts-1 do
		local pt = cpuPointBuf[i]
		local pos = pt.pos
		local x,y,z = pos:unpack()
		local r = pos:length()
		local v = pt.vel:length()
		local lum = pt.lum
		local temp = pt.temp
		local log10lum = math.log(lum) * _1_log_10		-- multiply this by -2.5 to get the abs mag ... but how is that intuitive?
		s:accum(x, y, z, r, v, lum, temp, log10lum)
	end
	print("data range (Pc):")
	print(s)
--]=]

	if cmdline.makeLog10LumBins then
		-- TODO why not just do this in the getstats or another offline tool?
		local log10lumbincount = 200
		local log10lumbins = require 'stat.bin'(
			s.log10lum.min + (s.log10lum.min - s.log10lum.max) * .5 / log10lumbincount,
			s.log10lum.max + (s.log10lum.max - s.log10lum.min) * .5 / log10lumbincount,
			log10lumbincount)
		
		for i=0,numPts-1 do
			log10lumbins:accum(math.log(cpuPointBuf[i].lum) * _1_log_10)
		end
		do	--if normalizebins then
			local total = 0
			local dx = log10lumbins:getDx()
			for i,x in ipairs(log10lumbins) do	-- I'd use table.sum but it sums non-natural keys
				total = total + x * dx
			end
			for i=1,#log10lumbins do	-- I'd use table.sum but it sums non-natural keys
				log10lumbins[i] = log10lumbins[i] / total
			end
		end
		
		-- plot it:
		local plotdatafn = 'log10lum-dist-'..set..'.txt'
		file[plotdatafn] = log10lumbins:getTextData()
		require 'gnuplot'{
			output = 'log10lum-dist-'..set..'.png',
			style = 'data linespoints',
			xlabel = 'log_{10} lum',
			ylabel = '%age',
			{datafile=plotdatafn, using='1:2', title=set},
		}
	end


	gpuPointBuf = GLArrayBuffer{
		size = numPts * ffi.sizeof(pt_t),
		data = cpuPointBuf,
	}

	local gpuPointBuf_attrs = {
		pos = {
			buffer = gpuPointBuf,
			size = 3,
			type = gl.GL_FLOAT,
			stride = ffi.sizeof(pt_t),
			offset = ffi.offsetof(pt_t, 'pos'),
		},
		vel = pt_t == 'pt_9col_t' and {
			buffer = gpuPointBuf,
			size = 3,
			type = gl.GL_FLOAT,
			stride = ffi.sizeof(pt_t),
			offset = ffi.offsetof(pt_t, 'vel'),
		} or nil,
		lum = pt_t == 'pt_9col_t' and {
			buffer = gpuPointBuf,
			size = 1,
			type = gl.GL_FLOAT,
			stride = ffi.sizeof(pt_t),
			offset = ffi.offsetof(pt_t, 'lum'),
		} or nil,
		temp = pt_t == 'pt_9col_t' and {
			buffer = gpuPointBuf,
			size = 1,
			type = gl.GL_FLOAT,
			stride = ffi.sizeof(pt_t),
			offset = ffi.offsetof(pt_t, 'temp'),
		} or nil,
		radius = pt_t == 'pt_9col_t' and {
			buffer = gpuPointBuf,
			size = 1,
			type = gl.GL_FLOAT,
			stride = ffi.sizeof(pt_t),
			offset = ffi.offsetof(pt_t, 'radius'),
		} or nil,
	}


	-- TODO show orbit diagrams for all stars around Sgt A*?
	local gpuVelLineBuf_attrs
	if buildVelocityBuffers then
		cpuVelLineBuf = ffi.new(pt_t..'[?]', numPts*2)
		for i=0,numPts-1 do
			cpuVelLineBuf[0+2*i] = ffi.new(pt_t, cpuPointBuf[i])
			cpuVelLineBuf[1+2*i] = ffi.new(pt_t, cpuPointBuf[i])
		end
		
		gpuVelLineBuf = GLArrayBuffer{
			size = ffi.sizeof(pt_t) * 2 * numPts,
			data = cpuVelLineBuf,
		}
	
		gpuVelLineBuf_attrs = {
			pos = {
				buffer = gpuVelLineBuf,
				size = 3,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof(pt_t),
				offset = ffi.offsetof(pt_t, 'pos'),
			},
			vel = pt_t == 'pt_9col_t' and {
				buffer = gpuVelLineBuf,
				size = 3,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof(pt_t),
				offset = ffi.offsetof(pt_t, 'vel'),
			} or nil,
			lum = pt_t == 'pt_9col_t' and {
				buffer = gpuVelLineBuf,
				size = 1,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof(pt_t),
				offset = ffi.offsetof(pt_t, 'lum'),
			} or nil,
			temp = pt_t == 'pt_9col_t' and {
				buffer = gpuVelLineBuf,
				size = 1,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof(pt_t),
				offset = ffi.offsetof(pt_t, 'temp'),
			} or nil,
			radius = pt_t == 'pt_9col_t' and {
				buffer = gpuVelLineBuf,
				size = 1,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof(pt_t),
				offset = ffi.offsetof(pt_t, 'radius'),
			} or nil,
		}



	end

	local nbhd_t_attrs
	if buildNeighborhood then
		-- maybe i can quick sort and throw in some random compare and that'll do some kind of rough estimate
print('looking for max')
		cpuNbhdLineBuf = vector'nbhd_t'
		-- heuristic with bins
		-- r stddev is 190, so 3 is 570
		local threeSigma = 570
		local min = -570
		local max = 570
		local nodeCount = 0
		local box3 = require 'vec.box3'
		local root = {
			index = nodeCount,
			box = box3{
				min = {min, min, min},
				max = {max, max, max},
			},
			pts = table(),
		}
		root.mid = (root.box.min + root.box.max) * .5
		nodeCount = nodeCount + 1
		local nodemax = 10
		local function addToTree(node, i)
			local pos = cpuPointBuf[i].pos
			
			-- have we divided?  pay it forward.
			if node.children then
				local childIndex = bit.bor(
					(pos.x > node.mid[1]) and 1 or 0,
					(pos.y > node.mid[2]) and 2 or 0,
					(pos.z > node.mid[3]) and 4 or 0)
				return addToTree(node.children[childIndex], i)
			end
		
			-- not divided yet?  push into leaf until it gets too big, then divide.
			node.pts:insert(i)
			if #node.pts >= nodemax then
				-- make children
				node.children = {}
				for childIndex=0,7 do
					local xL = bit.band(childIndex,1) == 0
					local yL = bit.band(childIndex,2) == 0
					local zL = bit.band(childIndex,4) == 0
					local child = {
						index = nodeCount,
						box = box3{
							min = {
								xL and node.box.min[1] or node.mid[1],
								yL and node.box.min[2] or node.mid[2],
								zL and node.box.min[3] or node.mid[3]
							},
							max = {
								xL and node.mid[1] or node.box.max[1],
								yL and node.mid[2] or node.box.max[2],
								zL and node.mid[3] or node.box.max[3]
							},
						},
						pts = table()
					}
					child.mid = (child.box.min + child.box.max) * .5
					child.parent = node
					node.children[childIndex] = child
					nodeCount = nodeCount + 1
				end
				-- split the nodes up into the children
				for _,i in ipairs(node.pts) do
					local pos = cpuPointBuf[i].pos
					local childIndex = bit.bor(
						(pos.x > node.mid[1]) and 1 or 0,
						(pos.y > node.mid[2]) and 2 or 0,
						(pos.z > node.mid[3]) and 4 or 0)
					addToTree(node.children[childIndex], i)
				end
				node.pts = nil
			end
		end
print'pushing into bins'
		for i=0,numPts-1 do
			addToTree(root, i)
		end
print('created '..nodeCount..' nodes')
print'searching bins'
	
		-- when connecting each star to all closest stars:
		-- nbhdThrehsold = 5 <=> 248422 lines
		-- nbhdThrehsold = 7 <=> 689832 lines
		-- nbhdThrehsold = 10 <=> too many
		local nbhdThreshold = 7	-- in Pc
		-- targetApparentMagnitude = -1 <=> 103454 lines (incl dups)
		-- targetApparentMagnitude = 0 <=> 412688 lines
		-- but this method looks dumb
		--local targetApparentMagnitude = 0

		local ai = 1
		local lastTime = os.time()
		local function searchTree(node)
			if node.children then
				assert(not node.pts)
				for childIndex=0,7 do
					searchTree(node.children[childIndex])
				end
			else
				assert(node.pts)

				-- TODO for each point in the leaf
				-- search all neighbors and all parents, and all parents' neighbors
				-- which means everything except siblings-of-siblings

				local pts = node.pts
				local n = #pts
				for i=1,n-1 do
					local pi = cpuPointBuf[pts[i]]

					--[[
					what if dist threshold of connection varies with its lum?
					like, dist threshold is how far away an apparent magnitude of X would be
					
					m = -5/2 log10(LStar/L0) - 5 + 5 log10(d)
					m/5 + 1/2 log10(LStar/L0) + 1 = log10(d)
					d = 10^(m/5 + 1) sqrt(LStar/L0)
					ex: for LSun, the dist at which it has
					app mag = -1 <=> d = 0.7 pc
					app mag = 0 <=> d = 1.1 pc
					app mag = 1 <=> d = 1.8 pc
					app mag = 2 <=> d = 2.8 pc
					app mag = 3 <=> d = 4.5 pc
					app mag = 4 <=> d = 7.1 pc
					app mag = 6 <=> d = 17.9 pc

					--]]
					--local nbhdThreshold = 10^(targetApparentMagnitude  / 5 + 1) * math.sqrt(pi.lum * LSunOverL0)

					-- TODO how about nbhd threshold based on abs mag? log10 of distance from neighboring stars?
					-- so you see more connectiosn for more visible stars?
					-- so the lines represent what can see the star?
					local touchbnd = box3{
						min = {
							pi.pos.x - nbhdThreshold,
							pi.pos.y - nbhdThreshold,
							pi.pos.z - nbhdThreshold,
						},
						max = {
							pi.pos.x + nbhdThreshold,
							pi.pos.y + nbhdThreshold,
							pi.pos.z + nbhdThreshold,
						},
					}

					local touchingNodes = table()
					local function search2(node2)
						-- if we don't contain the node, or if we don't touch the node, then bail
						if not node2.box:touches(touchbnd) then return end
						if node2.children then
							assert(not node2.pts)
							for ch=0,7 do
								search2(node2.children[ch])
							end
						else
							-- if touchbnd contains node2.box then just add all
							assert(node2.pts)
							touchingNodes:insert(node2)
						end
					end
					-- search through all touching nodes of this node
					search2(root)

					-- now compare up the tree
					-- but there are no points in non-leaf nodes
					-- so we have to find all leafs that are neighboring this node
					-- so ... i guess that means traverse the whole tree
					-- and look at whose bounding points are contained by our bounding points
					for _,node2 in ipairs(touchingNodes) do
						for j = (node == node2 and i+1 or 1),#node2.pts do
							-- don't process point pairs twice
							if node.pts[i] > node2.pts[j] then
								-- TODO Prevent duplicates, though now distance tests are asymmetric
								local pj = cpuPointBuf[node2.pts[j]]
								local distSq = (pi.pos - pj.pos):lenSq()
								if distSq < nbhdThreshold * nbhdThreshold then
									local dist = math.sqrt(distSq)
									local na = ffi.new'nbhd_t'
									na.pos:set(pi.pos:unpack())
									na.dist = dist
									cpuNbhdLineBuf:push_back(na)
									local nb = ffi.new'nbhd_t'
									nb.pos:set(pj.pos:unpack())
									nb.dist = dist
									cpuNbhdLineBuf:push_back(nb)
								end
							end
						end
					end

					--[[
					local thisTime = os.time()
					if lastTime ~= thisTime then
						lastTime = thisTime
						print((100 * (ai / nkbins)) .. '% done')
					end
					--]]
				end
			end
		end
		searchTree(root)
print('created '..#cpuNbhdLineBuf..' nbhd lines')

		gpuNbhdLineBuf = GLArrayBuffer{
			size = ffi.sizeof(cpuNbhdLineBuf.type) * #cpuNbhdLineBuf,
			data = cpuNbhdLineBuf.v,
		}
	
		nbhd_t_attrs = {
			pos = {
				buffer = gpuNbhdLineBuf,
				size = 3,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof'nbhd_t',
				offset = ffi.offsetof('nbhd_t', 'pos'),
			},
			dist = {
				buffer = gpuNbhdLineBuf,
				size = 1,
				type = gl.GL_FLOAT,
				stride = ffi.sizeof'nbhd_t',
				offset = ffi.offsetof('nbhd_t', 'dist'),
			},
		}
	end

	--refreshPoints()

	hsvtex = GLHSVTex(1024, nil, true)

	local starTexSize = 256
	starTex = GLTex2D{
		image = Image(starTexSize, starTexSize, 3, 'unsigned char', function(i,j)
			local u = ((i + .5) / starTexSize) * 2 - 1
			local v = ((j + .5) / starTexSize) * 2 - 1
			
			--[[
			local l = math.exp(-5 * (u*u + v*v))
			--]]
			-- [[
			local r = math.sqrt(u*u + v*v)
			local l = math.exp(-50 * (r - .75)^2)
			--]]
			
			l = math.floor(l * 255)
			return l,l,l
		end),
		wrap = {s = gl.GL_CLAMP_TO_EDGE, t = gl.GL_CLAMP_TO_EDGE},
		magFilter = gl.GL_LINEAR,
		minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
		generateMipmap = true,
	}

	-- black body color table from http://www.vendian.org/mncharity/dir3/blackbody/UnstableURLs/bbr_color_D58.html
	do
--[[ writing colorForTemp.png from the black body data:
		local rgbs = table()
		for l in io.lines'bbr_color_D58.txt' do
			if l ~= '' and l:sub(1,1) ~= '#' then
				local cmf = l:sub(11,15)
				if cmf == '10deg' then
					local temp = tonumber(string.trim(l:sub(2,6)))
					-- TODO instead of colorTempMin/colorTempMax determining the texture bounds, have the texture determine the colorTempMin/colorTempMax bounds
					-- and maybe remap it - logarithmically - or based on abs mag (log of lum)
					if temp >= colorTempMin and temp <= colorTempMax then
						local r = tonumber(l:sub(82,83), 16)
						local g = tonumber(l:sub(84,85), 16)
						local b = tonumber(l:sub(86,87), 16)
						rgbs:insert{r,g,b}
					end
				end
			end
		end
		local tempImg = Image(#rgbs, 1, 3, 'unsigned char', function(i,j)
			return table.unpack(rgbs[i+1])
		end)
		tempImg:save'../colorForTemp.png'
--]]
-- [[ or just reading it:
		local tempImg = Image'../colorForTemp.png'
--]]
		tempTex = GLTex2D{
			image = tempImg,
			wrap = {s = gl.GL_CLAMP_TO_EDGE, t = gl.GL_CLAMP_TO_EDGE},
			magFilter = gl.GL_LINEAR,
			minFilter = gl.GL_NEAREST,
		}
		glreport'here'
	end


local calcPointSize = template([[
<?
local clnumber = require 'cl.obj.number'
?>
#define M_1_LOG_10	<?=clnumber(1/math.log(10))?>
#define LSunOverL0	<?=clnumber(LSunOverL0)?>

	//how to calculate this in fragment space ...
	// coordinates are in Pc
	float distInPcSq = dot(vmv.xyz, vmv.xyz);
	
	//log(distInPc^2) = 2 log(distInPc)
	//so log10(distInPc) = .5 log10(distInPc^2)
	//so log10(distInPc) = .5 / log(10) * log(distInPc^2)
	float log10DistInPc = (.5 * M_1_LOG_10) * log(distInPcSq);

	//MBolStar - MBolSun = -2.5 * log10( LStar / LSun)
	float LStarOverLSun = lum;
	float LStarOverL0 = LSunOverL0 * LStarOverLSun;
	float absoluteMagnitude = (-2.5 * M_1_LOG_10) * log(LStarOverL0);	// abs magn

	/*
	apparent magnitude:
	M = absolute magnitude
	m = apparent magnitude
	d = distance in parsecs
	m = M - 5 + 5 * log10(d)
	*/
	float apparentMagnitude = absoluteMagnitude - 5. + 5. * log10DistInPc;
	
	/*
	ok now on to the point size ...
	and magnitude ...
	hmm ...
	HDR will have to factor into here somehow ...
	*/

	gl_PointSize = (6.5 - apparentMagnitude) * starPointSizeScale + starPointSizeBias;
]], {
		LSunOverL0 = LSunOverL0,
	})



	-- since the point renderer varies its gl_PointSize with the magnitude, I gotta do that here as well
	drawIDShader = GLProgram{
		vertexCode = template([[
#version 460

in vec3 pos;
in float lum;

out float discardv;
out vec3 color;

uniform float starPointSizeScale;
uniform float starPointSizeBias;
uniform float pickSizeBias;

uniform float sliceRMin, sliceRMax;

uniform float sliceLumMin, sliceLumMax;

uniform vec3 orbit;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

void main() {
	if (lum < sliceLumMin ||
		lum > sliceLumMax)
	{
		discardv = 1.;
		return;
	}

	vec4 vtx = vec4(pos.xyz, 1.);
	vec4 vmv = modelViewMatrix * vtx;
	gl_Position = projectionMatrix * vmv;

	<?=calcPointSize?>
	gl_PointSize = max(0., gl_PointSize);
	gl_PointSize += pickSizeBias;

	discardv = 0.;
	vec3 orbitToPos = pos - orbit;
	float distFromOrbitSq_in_Pc2 = dot(orbitToPos, orbitToPos);
	if (distFromOrbitSq_in_Pc2 < sliceRMin * sliceRMin ||
		distFromOrbitSq_in_Pc2 > sliceRMax * sliceRMax)
	{
		discardv = 1.;
		return;
	}

	float i = gl_VertexID;
	color.r = mod(i, 256.);
	i = (i - color.r) / 256.;
	color.r *= 1. / 255.;
	color.g = mod(i, 256.);
	i = (i - color.g) / 256.;
	color.g *= 1. / 255.;
	color.b = mod(i, 256.);
	i = (i - color.b) / 256.;
	color.b *= 1. / 255.;
}
]], {calcPointSize = calcPointSize}),
		fragmentCode = template[[
#version 460

in vec3 color;
in float discardv;
out vec4 fragColor;

void main() {
	if (discardv > 0.) {
		discard;
	}

	fragColor = vec4(color, 1.);
}
]],
		attrs = gpuPointBuf_attrs,
	}
	glreport'here'
	
	drawIDShader:useNone()
	glreport'here'


	-- i've neglected the postprocessing, and this has become the main shader
	-- which does the fake-postproc just with a varying gl_PointSize
	accumStarPointShader = GLProgram{
		vertexCode = template([[
<?
local clnumber = require 'cl.obj.number'
?>
#version 460

in vec3 pos;
in float lum;
in float temp;

out float lumv;
out vec3 tempcolor;
out float discardv;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float starPointSizeScale;
uniform float starPointSizeBias;
uniform sampler2D tempTex;

//how to occlude, based on radial distance, in Parsecs -- *FROM THE ORBIT*
uniform float sliceRMin, sliceRMax;
uniform float sliceLumMin, sliceLumMax;

//current orbiting position, used for occlusion
uniform vec3 orbit;

void main() {
	discardv = 0.;
	
	if (lum < sliceLumMin ||
		lum > sliceLumMax)
	{
		discardv = 1.;
		return;
	}

	vec4 vtx = vec4(pos.xyz, 1.);
	vec4 vmv = modelViewMatrix * vtx;
	gl_Position = projectionMatrix * vmv;

	<?=calcPointSize?>

	vec3 orbitToPos = pos - orbit;
	float distFromOrbitSq_in_Pc2 = dot(orbitToPos, orbitToPos);
	if (distFromOrbitSq_in_Pc2 < sliceRMin * sliceRMin ||
		distFromOrbitSq_in_Pc2 > sliceRMax * sliceRMax)
	{
		discardv = 1.;
		return;
	}
	
	lumv = 1.;

	// if the point size is < .5 then just make the star dimmer instead
	const float pointSizeMin = .5;
	float dimmer = gl_PointSize - pointSizeMin;
	if (dimmer < 0) {
		gl_PointSize = pointSizeMin;
		lumv *= pow(2., dimmer);
	}

	float tempfrac = (temp - <?=clnumber(colorTempMin)?>) * <?=clnumber(1/(colorTempMax - colorTempMin))?>;
	tempcolor = texture(tempTex, vec2(tempfrac, .5)).rgb;
}
]], 	{
			calcPointSize = calcPointSize,
			colorTempMin = colorTempMin,
			colorTempMax = colorTempMax,
		}),
		fragmentCode = template([[
#version 460

in float lumv;
in vec3 tempcolor;
in float discardv;

uniform float starPointAlpha;
uniform sampler2D starTex;

out vec4 fragColor;

void main() {
	if (discardv > 0.) {
		discard;
	}

	float lumf = lumv;

#if 0
	// make point smooth
	vec2 d = gl_PointCoord.xy * 2. - 1.;
	float rsq = dot(d,d);
	lumf *= 1. / (10. * rsq + .1);
#endif

	fragColor = vec4(tempcolor * lumf * starPointAlpha, 1.)
		* texture(starTex, gl_PointCoord)
	;
}
]]),
		uniforms = {
			tempTex = 0,
			starTex = 1,
		},
		attrs = gpuPointBuf_attrs,
	}
	

	accumStarPointShader:useNone()
	glreport'here'

	accumStarLineShader = GLProgram{
		vertexCode = [[
#version 460

in vec3 pos;
in vec3 vel;

uniform float lineVelScalar;
uniform bool normalizeVel;

uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

void main() {
	vec4 vtx = vec4(pos.xyz, 1.);

	vec3 velv = vel;
	if (normalizeVel) velv = normalize(velv);
	
	float end = float(gl_VertexID & 1);

	vtx.xyz += velv * (end * lineVelScalar);
	
	gl_Position = projectionMatrix * (modelViewMatrix * vtx);
}
]],
		fragmentCode = [[
#version 460

uniform float lineVelAlpha;

out vec4 fragColor;

void main() {
	fragColor = vec4(.1, 1., .1, lineVelAlpha);
}
]],
		attrs = gpuVelLineBuf_attrs,
	}
	accumStarLineShader:useNone()
	glreport'here'

	renderAccumShader = GLProgram{
		vertexCode = [[
#version 460

in vec3 pos;
out vec2 texcoord;

void main() {
	texcoord = pos.xy;
	gl_Position = vec4(pos.x * 2. - 1., pos.y * 2. - 1., 0., 1.);
}
]],
		fragmentCode = template[[
<?
local clnumber = require 'cl.obj.number'
?>
#version 460

in vec2 texcoord;

uniform sampler2D fbotex;
uniform sampler1D hsvtex;
uniform float hdrScale;
uniform float hdrGamma;
uniform float hsvRange;
uniform bool showDensity;
uniform float bloomLevels;

out vec4 fragColor;

void main() {
	fragColor = vec4(0., 0., 0., 0.);
<?
local maxLevels = 8
for level=0,maxLevels-1 do
?>	if (bloomLevels >= <?=clnumber(level)?>) fragColor += texture(fbotex, texcoord, <?=clnumber(level)?>);
<?
end
?>
	fragColor *= hdrScale * <?=clnumber(1/maxLevels)?>;

	if (showDensity) {
		fragColor = texture(hsvtex, log(dot(fragColor.rgb, vec3(.3, .6, .1)) + 1.) * hsvRange);
	} else {
		//tone mapping, from https://learnopengl.com/Advanced-Lighting/HDR
		//fragColor.rgb = fragColor.rgb / (fragColor.rgb + vec3(1.));
		fragColor.rgb = pow(fragColor.rgb, vec3(1. / hdrGamma));
		fragColor.rgb = log(fragColor.rgb + vec3(1.));
	}

	fragColor.a = 1.;
}
]],
		uniforms = {
			fbotex = 0,
			hsvtex = 1,
			showDensity = false,
		},
	}
	glreport'here'
	renderAccumShader:useNone()
	glreport'here'

	if buildNeighborhood then
		drawNbhdLineShader = GLProgram{
			vertexCode = [[
#version 460

in vec3 pos;
in float dist;

out float lumv;

uniform float nbhdLineAlpha;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

void main() {
	vec4 vtx = vec4(pos, 1.);
	vec4 vmv = modelViewMatrix * vtx;
	gl_Position = projectionMatrix * vmv;
	float viewDist = length(vmv.xyz);
	lumv = nbhdLineAlpha / (dist * dist
		* (viewDist + 1.)
);
}
]],
			fragmentCode = [[
#version 460

in float lumv;
out vec4 fragColor;

void main() {
	float lumf = lumv;
	vec3 color = vec3(.1, 1., .1) * lumf;
	fragColor = vec4(color, 1.);
}
]],
			attrs = nbhd_t_attrs,
		}
		glreport'here'
		drawNbhdLineShader:useNone()
		glreport'here'
	end
end

ffi.cdef[[
typedef union {
	uint8_t rgba[4];
	uint32_t i;
} pixel4_t;
]]

local pixel = ffi.new'pixel4_t'

local selectedIndex

local modelViewMatrix = matrix_ffi.zeros(4,4)
local projectionMatrix = matrix_ffi.zeros(4,4)

function App:drawPickScene()
	gl.glClearColor(1,1,1,1)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	gl.glEnable(gl.GL_DEPTH_TEST)

	-- TODO
	-- 1) smaller viewpoint of 3x3 pixels
	-- 2) pick matrix to zoom in on the specific mouse location

	gl.glEnable(gl.GL_PROGRAM_POINT_SIZE)
	
	drawIDShader:use()
	drawIDShader:setUniforms{
		starPointSizeScale = starPointSizeScale,
		starPointSizeBias = starPointSizeBias,
		pickSizeBias = pickSizeBias,
		sliceRMin = sliceRMin,
		sliceRMax = sliceRMax,
		sliceLumMin = sliceLumMin,
		sliceLumMax = sliceLumMax,
		orbit = {self.view.orbit:unpack()},	-- TODO
		modelViewMatrix = modelViewMatrix.ptr,
		projectionMatrix = projectionMatrix.ptr,
	}

	drawIDShader.vao:use()
	gl.glDrawArrays(gl.GL_POINTS, 0, numPts)
	drawIDShader.vao:useNone()

	drawIDShader:useNone()

	gl.glDisable(gl.GL_PROGRAM_POINT_SIZE)

	gl.glFlush()

	gl.glReadPixels(
		self.mouse.ipos.x,
		self.height - self.mouse.ipos.y - 1,
		1,
		1,
		gl.GL_RGB,
		gl.GL_UNSIGNED_BYTE,
		pixel)

	selectedIndex = tonumber(pixel.i)
	gl.glDisable(gl.GL_DEPTH_TEST)
end

function App:drawScene()
	gl.glClearColor(0,0,0,0)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	
	if drawPoints then
		gl.glEnable(gl.GL_PROGRAM_POINT_SIZE)
		gl.glEnable(gl.GL_POINT_SPRITE)			-- i thought this was perma-on after 3.2?  guess not
		
		accumStarPointShader:use()
		accumStarPointShader:setUniforms{
			starPointAlpha = starPointAlpha,
			starPointSizeScale = starPointSizeScale,
			starPointSizeBias = starPointSizeBias,
			sliceRMin = sliceRMin,
			sliceRMax = sliceRMax,
			sliceLumMin = sliceLumMin,
			sliceLumMax = sliceLumMax,
			modelViewMatrix = modelViewMatrix.ptr,
			projectionMatrix = projectionMatrix.ptr,
			--orbit = self.view.orbit.s,		-- TODO orbit is a vec3 is a glsl float type
			orbit = {self.view.orbit:unpack()},	-- but view.orbit is double[3], so we can't pass it as a ptr of the type associated with the glsl type(float)
		}

		tempTex:bind(0)
		starTex:bind(1)

		accumStarPointShader.vao:use()
		gl.glDrawArrays(gl.GL_POINTS, 0, numPts)
		accumStarPointShader.vao:useNone()
	
		starTex:unbind(1)
		tempTex:unbind(0)
		
		accumStarPointShader:useNone()
		
		gl.glDisable(gl.GL_POINT_SPRITE)
		gl.glDisable(gl.GL_PROGRAM_POINT_SIZE)
	end
	if drawVelLines and buildVelocityBuffers then
		accumStarLineShader:use()
		accumStarLineShader:setUniforms{
			lineVelAlpha = lineVelAlpha,
			lineVelScalar = lineVelScalar,
			normalizeVel = normalizeVel and 1 or 0,
			modelViewMatrix = modelViewMatrix.ptr,
			projectionMatrix = projectionMatrix.ptr,
		}

		accumStarLineShader.vao:use()
		gl.glDrawArrays(gl.GL_LINES, 0, 2 * numPts)
		accumStarLineShader.vao:useNone()
		accumStarLineShader:useNone()
	end

	gl.glDisable(gl.GL_BLEND)
end

local lastWidth, lastHeight
function App:drawWithAccum()
	if self.width ~= lastWidth or self.height ~= lastHeight then
		lastWidth = self.width
		lastHeight = self.height

		fbo = GLFBO{
			width=self.width,
			height=self.height,
		}
		
		fbotex = GLTex2D{
			width = fbo.width,
			height = fbo.height,
			format = gl.GL_RGBA,
			type = gl.GL_FLOAT,
			internalFormat = gl.GL_RGBA32F,
			minFilter = gl.GL_LINEAR_MIPMAP_LINEAR,
			magFilter = gl.GL_LINEAR,
		}
	end
	
	fbo:draw{
		viewport = {0,0,fbo.width,fbo.height},
		dest = fbotex,
		callback = function()
			self:drawScene()
		end,
	}

	fbotex:bind()
	gl.glGenerateMipmap(fbotex.target)
	fbotex:unbind()

	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glViewport(0, 0, self.width, self.height)
	
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glPushMatrix()
	gl.glLoadIdentity()
	gl.glOrtho(0,1,0,1,-1,1)
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glPushMatrix()
	gl.glLoadIdentity()

	renderAccumShader:use()
	renderAccumShader:setUniforms{
		hdrScale = hdrScale,
		hdrGamma = hdrGamma,
		hsvRange = hsvRange,
		bloomLevels = bloomLevels,
		showDensity = showDensity and 1 or 0,
	}
	fbotex:bind(0)
	hsvtex:bind(1)

	gl.glBegin(gl.GL_QUADS)
	gl.glVertex2f(0,0)
	gl.glVertex2f(0,1)
	gl.glVertex2f(1,1)
	gl.glVertex2f(1,0)
	gl.glEnd()

	hsvtex:unbind(1)
	fbotex:unbind(0)
	renderAccumShader:useNone()

	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glPopMatrix()
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glPopMatrix()

end

local function sphericalToCartesian(r,theta,phi)
	local ct = math.cos(theta)
	local st = math.sin(theta)
	local cp = math.cos(phi)
	local sp = math.sin(phi)
	return vec3d(
		r * cp * st,
		r * sp * st,
		r * ct
	)
end

function App:update()
	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, projectionMatrix.ptr)
	
	self:drawPickScene()

	if not showPickScene then
		-- [[
		self:drawScene()
		--]]
		--[[
		self:drawWithAccum()
		--]]
	end


	-- TODO inv square reduce this.... by inv square of one another, and by inv square from view
	if buildNeighborhood and showNeighbors then
		gl.glEnable(gl.GL_BLEND)
		drawNbhdLineShader:use()
		drawNbhdLineShader:setUniforms{
			nbhdLineAlpha = nbhdLineAlpha,
			modelViewMatrix = modelViewMatrix.ptr,
			projectionMatrix = projectionMatrix.ptr,
		}

		drawNbhdLineShader.vao:use()
		gl.glDrawArrays(gl.GL_LINES, 0, cpuNbhdLineBuf.size)
		drawNbhdLineShader.vao:useNone()
		drawNbhdLineShader:useNone()
		gl.glDisable(gl.GL_BLEND)
	end


	-- TODO draw around origin?  or draw around view orbit?
	if drawGrid then
		if tiltGridToPlanet then
			gl.glPushMatrix()
			gl.glRotatef(23.4365472133, -1, 0, 0)
		end
		gl.glEnable(gl.GL_BLEND)
		gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
		gl.glColor3f(.25, .25, .25)
		gl.glBegin(gl.GL_LINES)
		local idiv = 24
		local dphi = 2 * math.pi / idiv
		local jdiv = 12
		local dtheta = math.pi / jdiv
		for i=0,idiv-1 do
			local phi = 2 * math.pi * i / idiv
			for j=0,jdiv-1 do
				local theta = math.pi * j / jdiv
				gl.glVertex3f((sphericalToCartesian(gridRadius, theta, phi) + self.view.orbit):unpack())
				gl.glVertex3f((sphericalToCartesian(gridRadius, theta + dtheta, phi) + self.view.orbit):unpack())
				if j > 0 then
					gl.glVertex3f((sphericalToCartesian(gridRadius, theta, phi) + self.view.orbit):unpack())
					gl.glVertex3f((sphericalToCartesian(gridRadius, theta, phi + dphi) + self.view.orbit):unpack())
				end
			end
		end
		gl.glEnd()
		gl.glDisable(gl.GL_BLEND)
		if tiltGridToPlanet then
			gl.glPopMatrix()
		end
	end

	if showConstellations then
		gl.glColor3f(1,1,0)
		for _,constellation in ipairs(constellations) do
			if constellation.enabled then
				for _,line in ipairs(constellation.lines) do
					gl.glBegin(gl.GL_LINE_STRIP)
					for _,i in ipairs(line) do
						gl.glVertex3f(cpuPointBuf[i].pos:unpack())
					end
					gl.glEnd()
				end
			end
		end
	end

	glreport'here'
	App.super.update(self)

	-- with the gaia data, when filtering/rebuilding the data, i get 'out of memory' errors that hopefully this will help
	collectgarbage()
end


local function checkboxTooltipTable(title, ...)
	ig.igPushID_Str(title)
	local result = ig.luatableCheckbox('', ...)
	if ig.igIsItemHovered(ig.ImGuiHoveredFlags_None) then
		ig.igBeginTooltip()
		ig.igText(title)
		ig.igEndTooltip()
	end
	ig.igPopID()
	return result
end

local nameWithAppMagLastPos
local namesWithAppMag
local function guiShowStars(self)
	-- do this before any other tooltip, so it will be on bottom
	-- there's usually just 5000 or so of these
	-- and if we filter by apparent magnitude then there can't be many visible at once
	if showStarNames
	and namedStars
	-- and showAllNamedStarsAtOnce
	then
		ig.igPushID_Str('star names')
	
		-- global / persist: nameWithAppMagLastPos
		if not nameWithAppMagLastPos
		or (nameWithAppMagLastPos - self.view.pos):lenSq() > .01	-- greater than some epsilon of how far to move, squared.  make the dist less than the closest stars in the dataset
		then
			nameWithAppMagLastPos = nameWithAppMagLastPos or vec3d()
			nameWithAppMagLastPos:set(self.view.pos:unpack())
			-- now sort all named indexes basedon their apparent magnitude
			-- global / persist: namesWithAppMag
			namesWithAppMag = table.map(namedStars, function(name, index, t)
				local pt = cpuPointBuf[index]
				local distSq = (pt.pos - self.view.pos):lenSq()
				local log10DistInPc = .5 * _1_log_10 * math.log(distSq)
				local LStarOverL0 = pt.lum * LSunOverL0
				local absmag = -2.5 * _1_log_10 * math.log(LStarOverL0)
				local appmag = absmag - 5. + 5. * log10DistInPc

				return {
					name = name,
					index = index,
					appmag = appmag,
				}, #t+1
			end)
			namesWithAppMag:sort(function(a,b)
				return a.appmag < b.appmag
			end)
		end
		

		local windowCount = 0
		local function addWindowForStar(x, y, name)
			ig.igPushID_Str(name)
			ig.igSetNextWindowPos(
				ig.ImVec2(x,y), 	-- ImVec2 pos
				0,					-- ImGuiCond cond
				ig.ImVec2()			-- ImVec2 pivot
			)
			-- there's got to be a better way to have floating text in imgui
			-- imgui can't handle more than one tooltip
			-- and popups seem to eat up the input and i can't drag the screen any more.
			-- even if they have all the 'no input' and 'no nav' flags possible set
			--ig.igOpenPopup('star name', 0)
			-- igBeginPopup vs igBeginPopupEx, the only difference is normal takes str for title, Ex takes ImGuiID which is an int ... ?
--					if ig.igBeginPopup('star name',
--						bit.bor(
--							ig.ImGuiWindowFlags_NoFocusOnAppearing,
--							ig.ImGuiWindowFlags_NoBringToFrontOnFocus,
--						)
--					) then
				ig.igBegin(name, nil, bit.bor(
--						ig.ImGuiWindowFlags_NoSavedSettings,
--						ig.ImGuiWindowFlags_NoNav,
--						ig.ImGuiWindowFlags_NoInputs,	-- crashes?
				
				ig.ImGuiWindowFlags_NoDecoration,
--						ig.ImGuiWindowFlags_NoResize,
--						ig.ImGuiWindowFlags_NoScrollbar,
--						ig.ImGuiWindowFlags_NoCollapse

				ig.ImGuiWindowFlags_Tooltip
			))
			ig.igText(name)
			ig.igEnd()
			ig.igPopID()
			windowCount = windowCount + 1
			return windowCount > 10
		end

		local visibleNamedStars
--[[
		for index,name in pairs(namedStars) do
--]]
-- [[
		for _,info in ipairs(namesWithAppMag) do
			local index, name = info.index, info.name
--]]
			assert(index >= 0 and index < numPts)
			local pt = cpuPointBuf[index]
			local vpt = modelViewMatrix * require 'matrix'{pt.pos.x, pt.pos.y, pt.pos.z, 1}
--print('pt.pos', pt.pos)
--print('pt.pos.x', pt.pos.x)
--print('pt.pos.y', pt.pos.y)
--print('pt.pos.z', pt.pos.z)
--print("require 'matrix'{pt.pos.x, pt.pos.y, pt.pos.z, 1}", require 'matrix'{pt.pos.x, pt.pos.y, pt.pos.z, 1})
--print('modelViewMatrix ', modelViewMatrix)
--print('vpt', vpt)
--print('vpt', require 'ext.tolua'(vpt))
			local mz = vpt[3]
			if -mz > self.view.znear
			and -mz < self.view.zfar
			then
				-- find the point in screen coords
				local spt = projectionMatrix * vpt
				spt = spt / spt[4]

				-- only in bounds in the projection
				local x = spt[1]
				local y = spt[2]
				if x >= -1 and x <= 1
				and y >= -1 and y <= 1
				then
					x = (1 + x) * .5 * self.width
					y = (1 - y) * .5 * self.height

					local mx = vpt[1]
					local my = vpt[2]

					-- [[ handle immediately
					if addWindowForStar(x, y, name) then break end
					--]]
					--[=[ insert and sort now or later:
					local star = {
						x = x,
						y = y,
						name = name,
					}
					local distSq = mx*mx + my*my + mz*mz
					local log10DistInPc = .5 * _1_log_10 * math.log(distSq)
					local LStarOverL0 = pt.lum * LSunOverL0
					local absmag = -2.5 * _1_log_10 * math.log(LStarOverL0)
					local appmag = absmag - 5. + 5. * log10DistInPc
					star.appmag = appmag
					--[[ insert and sort later?
					visibleNamedStars = visibleNamedStars or table()
					visibleNamedStars:insert(star)
					--]]
					--[[ insert in-order.
					-- TODO insert in-order in a tree?
					visibleNamedStars = visibleNamedStars or table()
					local found
					for j=1,#visibleNamedStars do
						if visibleNamedStars[j].appmag > star.appmag then
							visibleNamedStars:insert(j, star)
							found = true
							break
						end
					end
					if not found then
						visibleNamedStars:insert(star)
					end
					--]]
					--]=]
					-- TODO how about re-sorting the stars based on apparent magnitude from the oribiting star
					-- and only recalculating them when the orbit moves?
				end
			end
		end
		
		--[[ insert and sort later
		visibleNamedStars:sort(function(a,b)
			return a.appmag > b.appmag
		end)
		--]]

		--[[ don't handle immediately
		if visibleNamedStars then
			for _,star in ipairs(visibleNamedStars) do
				if addWindowForStar(star.x, star.y, star.name) then break end
			end
		end
		--]]
		
		ig.igPopID()
	end
end

showAllConstellations = false
local search = {
	orbit = '',
	lookat = '',
}
function App:updateGUI()

	ig.luatableCheckbox('draw points', _G, 'drawPoints')
	ig.luatableSliderFloat('point size scale', _G, 'starPointSizeScale', -10, 10)
	ig.luatableSliderFloat('point size bias', _G, 'starPointSizeBias', -10, 10)
	ig.luatableSliderFloat('pick size', _G, 'pickSizeBias', 0, 20)
	ig.luatableInputFloat('point alpha value', _G, 'starPointAlpha')

--[[ not used atm
	ig.luatableCheckbox('show density', _G, 'showDensity')
	ig.luatableSliderFloat('hdr scale', _G, 'hdrScale', 0, 1000, '%.7f', 10)
	ig.luatableSliderFloat('hdr gamma', _G, 'hdrGamma', 0, 1000, '%.7f', 10)
	ig.luatableSliderFloat('hsv range', _G, 'hsvRange', 0, 1000, '%.7f', 10)
	ig.luatableSliderFloat('bloom levels', _G, 'bloomLevels', 0, 8)
--]]
	
	ig.luatableInputFloat('slice r min', _G, 'sliceRMin')
	ig.luatableInputFloat('slice r max', _G, 'sliceRMax')
	
	ig.luatableInputFloat('slice lum min', _G, 'sliceLumMin')
	ig.luatableInputFloat('slice lum max', _G, 'sliceLumMax')

	
	ig.luatableCheckbox('show pick scene', _G, 'showPickScene')
	ig.luatableCheckbox('show grid', _G, 'drawGrid')
	ig.luatableCheckbox('tilt to planet', _G, 'tiltGridToPlanet')
	ig.luatableInputFloat('grid radius', _G, 'gridRadius')

	-- draw nhbd lines stuff which looks dumb right now
	if buildNeighborhood then
		ig.luatableCheckbox('show nbhd', _G, 'showNeighbors')
		ig.luatableInputFloat('nbhd alpha', _G, 'nbhdLineAlpha')
	end

	-- view stuff
	if ig.igButton'reset view' then
		-- making some assumptions here
		self.view.orbit:set(0,0,0)
		self.viewDist = self.viewDist
		self.view.pos:set(0, 0, self.viewDist)
		self.view.angle:set(0, 0, 0, 1)
	end
	ig.igText('dist (Pc) '..(self.view.pos - self.view.orbit):length())
	ig.luatableInputFloat('znear', self.view, 'znear')
	ig.luatableInputFloat('zfar', self.view, 'zfar')
	ig.luatableSliderFloat('fov y', self.view, 'fovY', 0, 180)

	if namedStars then
		if ig.luatableInputText('orbit', search, 'orbit', ig.ImGuiInputTextFlags_EnterReturnsTrue) then
			for i,v in pairs(namedStars) do
				if v == search.orbit then
					assert(i >= 0 and i < numPts, "oob index in name table "..i)
					local pt = cpuPointBuf[i]
					self.view.orbit:set(pt.pos:unpack())
				end
			end
		end
		if ig.luatableInputText('look at', search, 'lookat', ig.ImGuiInputTextFlags_EnterReturnsTrue) then
			for i,v in pairs(namedStars) do
				if v == search.lookat then
					local orbitDist = (self.view.pos - self.view.orbit):length()
					local fwd = -self.view.angle:zAxis()
					
					assert(i >= 0 and i < numPts, "oob index in name table "..i)
					local pt = cpuPointBuf[i]
					local to = (pt.pos - self.view.pos):normalize()

					local angle = math.acos(math.clamp(fwd:dot(to), -1, 1))
					local axis = fwd:cross(to):normalize()

					local rot = quatd():fromAngleAxis(axis.x, axis.y, axis.z, math.deg(angle))
					quatd.mul(rot, self.view.angle, self.view.angle)
					
					self.view.pos = self.view.orbit + self.view.angle:zAxis() * orbitDist
				end
			end
		end
	end


	-- draw vel lines
	if buildVelocityBuffers then
		ig.luatableCheckbox('draw vel lines', _G, 'drawVelLines')
		ig.luatableInputFloat('line alpha value', _G, 'lineVelAlpha')
		ig.luatableInputFloat('line vel scalar', _G, 'lineVelScalar')
		ig.luatableCheckbox('normalize velocity', _G, 'normalizeVel')
	end

	-- gui stuff
	ig.luatableCheckbox('show mouseover info', _G, 'showInformation')
	ig.luatableCheckbox('show star names', _G, 'showStarNames')

	ig.luatableCheckbox('show constellations', _G, 'showConstellations')
	if showConstellations then
		if checkboxTooltipTable('All', _G, 'showAllConstellations') then	-- TODO infer
			for _,constellation in ipairs(constellations) do
				constellation.enabled = showAllConstellations
			end
		end
		ig.igSameLine()

		local count = 0
		for _,constellation in ipairs(constellations) do
			if constellation.name then
				if count > 0
				and count % 10 ~= 9
				then
					ig.igSameLine()
				end
				checkboxTooltipTable(constellation.name, constellation, 'enabled')
				count = count + 1
			end
		end
	end

	ig.igImage(
		ffi.cast('void*', ffi.cast('intptr_t', tempTex.id)),
		ig.ImVec2(128, 24),
		ig.ImVec2(0,0),
		ig.ImVec2(1,1))

	if showInformation
	and selectedIndex < 0xffffff
	and selectedIndex >= 0
	and selectedIndex < numPts
	then
		local s = table()
		s:insert('index: '..('%06x'):format(selectedIndex))

		local name = namedStars and namedStars[selectedIndex] or nil
		if name then
			s:insert('name: '..tostring(name))
		end
		local pt = cpuPointBuf[selectedIndex]
		local dist = (pt.pos - self.view.orbit):length()
			
		local LStarOverLSun = pt.lum
		local absmag = (-2.5 / math.log(10)) * math.log(LStarOverLSun * LSunOverL0)
		local appmag = absmag - 5 + (5 / math.log(10)) * math.log(dist)

		s:insert('dist (Pc): '..dist)
		s:insert('lum (LSun): '..LStarOverLSun)
		s:insert('temp (K): '..pt.temp)
		s:insert('abs mag: '..absmag)
		s:insert('app mag: '..appmag)

		ig.igBeginTooltip()
		ig.igText(s:concat'\n')
		ig.igEndTooltip()
	end

	-- hmm, whether i call this first or last, it always shows up over the tooltips ..
	-- maybe because of the type of popup i'm using?
	guiShowStars(self)
end

--[[
watershed of velocity functions ...

dphi/dx(p1) = v1_x
dphi/dy(p1) = v1_y
dphi/dz(p1) = v1_z

dphi/dx(p2) = v1_x
dphi/dy(p2) = v1_y
dphi/dz(p2) = v1_z

...

dphi/dx(pn) = v1_x
dphi/dy(pn) = v1_y
dphi/dz(pn) = v1_z

becomes

(phi(p[i] + h e_x) - phi(p[i] - h e_x)) / (2h) ~= v[i]_x
(phi(p[i] + h e_y) - phi(p[i] - h e_y)) / (2h) ~= v[i]_y
(phi(p[i] + h e_z) - phi(p[i] - h e_z)) / (2h) ~= v[i]_z

is not invertible for odd # of elements ...


--]]

App():run()
