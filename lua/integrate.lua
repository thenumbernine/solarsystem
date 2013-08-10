local integrate = {}
setmetatable(integrate, integrate)


function integrate.euler(t,x,dt,f)
	x = x + f(t, x) * dt
	return x
end

function integrate.midpoint(t,x,dt,f)
	local k = f(t, x) * dt
	x = x + f(t + dt * .5, x + k * .5) * dt
	return x
end

function integrate.heun(t,x,dt,f)
	local fAtX = f(t, x)
	local xTilde = x + fAtX * dt
	x = x + (fAtX + f(t + dt, xTilde)) * (dt * .5)
	return x
end

function integrate.rk2alpha(t,x,dt,f,args)
	local fAtX = f(t, x)
	local k = fAtX * dt
	local alpha = args and args.alpha or .5		-- alpha = .5 <=> midpoint, alpha = 1 <=> Heun
	local frac = 1 / (2 * alpha)
	x = x + (fAtX * (1 - frac) + f(t + alpha * dt, x + k * alpha) * frac) * dt
	return x
end

function integrate.rk4(t,x,dt,f)
	local k1 = f(t, x) * dt
	local k2 = f(t + dt * .5, x + k1 * .5) * dt
	local k3 = f(t + dt * .5, x + k2 * .5) * dt
	local k4 = f(t + dt, x + k3) * dt
	x = x + (k1 + k2 * 2 + k3 * 2 + k4) * (1 / 6)
	return x
end

function integrate.rkf45(t,x,dt,f,args)
	local accuracy = args and args.accuracy or 1e-5
	local tEndFinal = t + dt
	local tEnd = tEndFinal
	local norm = args and args.norm or math.abs
	local maxiterations = args and args.iterations or 100
	repeat
		local currentDt = tEnd - t
		for i=1,maxiterations do
			local k1 = f(t, x) * currentDt
			local k2 = f(t + currentDt * .25, x + k1 * .25) * currentDt
			local k3 = f(t + currentDt * (3/8), x + k1 * (3/32) + k2 * (9/32)) * currentDt
			local k4 = f(t + currentDt * (12/13), x + k1 * (1932/2197) - k2 * (7200/2197) + k3 * (7296/2197)) * currentDt
			local k5 = f(t + currentDt, x + k1 * (439/216) - k2 * 8 + k3 * (3680/513) - k4 * (845/4104)) * currentDt
			local k6 = f(t + currentDt * .5, x - k1 * (8/27) + k2 * 2 - k3 * (3544/2565) + k4 * (1859/4104) - k5 * (11/40)) * currentDt
			local xHi = x + k1 * (16/135) + k3 * (6656/12825) + k4 * (28561/56430) - k5 * (9/50) + k6 * (2/55)
			local xLo = x + k1 * (25/216) + k3 * (1408/2565) + k4 * (2197/4104) - k5 * (1/5)
			local xErr = xHi - xLo
			-- here's the test: error threshold
			-- depends on calculating a magnitude of x, which depends on normalizing its values (so all contribute equally)
			xErr = norm(xErr)
			--print('err',xErr,'vs',accuracy)
			if xErr < accuracy then
				x = xHi
				t = tEnd
				tEnd = tEndFinal
				break
			end
			--print('error threshold won!')
			currentDt = currentDt * .5
			tEnd = t + currentDt
		end
	until t == tEndFinal
	return x
end



integrate.methods = {
	euler = integrate.euler,
	midpoint = integrate.midpoint,
	heun = integrate.heun,
	rk2alpha = integrate.rk2Alpha,
	rk4 = integrate.rk4,
	rkf45 = integrate.rkf45,
}

--[[
arguments:
	t = parameter
	x = initial integral function value, F(t)
	dt = change in parameter to integrate
	f(t,x) = function to integrate
	methodName = name of method to use.  default: euler
		options: euler, midpoint, heun, rk2alpha, rk4, rkf45
	args = (optional) extra arguments used by some integrators

extra arguments:
	rk2alpha:
		alpha = alpha parameter
		
	rkf45:
		norm(x) = norm of x.  by default this is math.abs
		accuracy = threshold of subdivision.  default is 1e-5
		iterations = maximum iterations.  default is 100
	
operators used: (assuming f() returns an object of metatable x)
	+(x,x) vector-vector addition
	*(x,t) vector-scalar product
	
	-(x,x) vector-vector subtraction is used by rkf45
	
returns:
	F(t+dt)
--]]
function integrate:__call(t, x, dt, f, methodName, args)
	methodName = methodName or 'euler'
	local method = self.methods[methodName]
	return method(t,x,dt,f,args)
	
	-- Minkowski: post-integration, re-normalize velocity:
	--x.u[0] = math.sqrt(x.u[1] * x.u[1] + 1)	
end

return integrate
