#!/usr/bin/env luajit
-- quick point visualization of the dynamic vars
local path = require 'ext.path'
local string = require 'ext.string'
local vec3f = require 'vec-ffi.vec3f'
local gl = require 'gl'
local App = require 'imgui.appwithorbit'()
App.viewUseGLMatrixMode = true
App.title = 'NASA Horizons Data Viewer'
function App:initGL(...)
	App.super.initGL(self, ...)
	self.bodies = require 'dkjson'.decode((path'dynamic-vars.json':read():match'.-=(.*)')).coords
	for i,body in ipairs(self.bodies) do
		for j=1,3 do
			body.pos[j] = assert(tonumber(string.trim(body.pos[j]:lower()))) * 1e-7
			body.vel[j] = assert(tonumber(string.trim(body.vel[j]:lower()))) * 1e-7
		end
		body.color = vec3f():map(function() return math.random() end):normalize()
	end
	self.view.znear = .1
	self.view.zfar = 10000
end
function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))
	gl.glBegin(gl.GL_LINES)
	gl.glColor3f(1,0,0) gl.glVertex3f(0,0,0) gl.glVertex3f(1,0,0)
	gl.glColor3f(0,1,0) gl.glVertex3f(0,0,0) gl.glVertex3f(0,1,0)
	gl.glColor3f(0,0,1) gl.glVertex3f(0,0,0) gl.glVertex3f(0,0,1)
	gl.glEnd()
	gl.glPointSize(3)
	gl.glBegin(gl.GL_POINTS)
	for _,body in ipairs(self.bodies) do
		gl.glColor3fv(body.color.s)
		gl.glVertex3f(table.unpack(body.pos))
	end
	gl.glEnd()
	gl.glPointSize(1)
	App.super.update(self, ...)
end
return App():run()
