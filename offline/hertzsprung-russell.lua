#!/usr/bin/env luajit
local ffi = require 'ffi'
local gl = require 'gl'
local glreport = require 'gl.report'
local vec2f = require 'vec-ffi.vec2f'
local vec3f = require 'vec-ffi.vec3f'
local vector = require 'stl.vector'
local math = require 'ext.math'
local table = require 'ext.table'

local _1_log10 = 1 / math.log(10)

local lpts = table()

local stats = require 'stat.set'('temp', 'log10lum')
print'loading hyg'
local hyg = require 'csv'.file'../hyg/hyg_v38.csv'
hyg:setColumnNames(hyg.rows:remove(1))
print'reading hyg'
for _,row in ipairs(hyg.rows) do
	if row.ci ~= ''
	and row.lum ~= ''
	then
		local colorIndex = tonumber(row.ci)
		local lum = tonumber(row.lum)
		if colorIndex and lum then
			local log10lum = math.log(lum) * _1_log10
			local temp = 4600 * (1 / (.92 * colorIndex + 1.7) + 1 / (.92 * colorIndex + .62))
			lpts:insert(vec3f(temp, log10lum, 0))
			stats:accum(temp, log10lum)
		end
	end
end
local n = #lpts
print('read '..n..' points')

local pts = ffi.new('vec3f[?]', n)
for i,pt in ipairs(lpts) do
	pts[i-1] = lpts[i]
end

local App = require 'imgui.appwithorbit'()
function App:initGL()
	App.super.initGL(self)
	self.view.ortho = true

	local vtxdata = ffi.new('vec2f[?]', n)
	for i=0,n-1 do
		local temp, log10lum = pts[i]:unpack()
		vtxdata[i]:set(
			(temp - stats.temp.min) / (stats.temp.max - stats.temp.min) * 2 - 1,
			(log10lum - stats.log10lum.min) / (stats.log10lum.max - stats.log10lum.min) * 2 - 1
		)
	end
	drawPoints = require 'gl.sceneobject'{
		program = {
			version = 'latest',
			precision = 'best',
			vertexCode = [[
in vec2 vertex;
uniform mat4 mvProjMat;
void main() {
	gl_Position = mvProjMat * vec4(vertex, 0., 1.);
}
]],
			fragmentCode = [[
out vec4 fragColor;
uniform vec4 color;
void main() {
	fragColor = color;
}
]],
		},
		uniforms = {
			color = {1,1,1,.1},
			mvProjMat = self.view.mvProjMat.ptr,
		},
		geometry = {
			mode = gl.GL_POINTS,
		},
		vertexes = {
			dim = 2,
			size = ffi.sizeof'vec2f' * n,
			count = n,
			data = vtxdata,
		},
	}
end
function App:update()
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)
	
	gl.glPointSize(2)
	
	drawPoints:draw()

	glreport'here'
	App.super.update(self)
end

return App():run()
