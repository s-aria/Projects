import sympy as sp
from sympy import *
import numpy as np
import matplotlib.pyplot as plt



def least_squares_orth(f, psi, Omega):
	"""
	Same as least_squares, but for orthogonal
	basis such that one avoids calling up standard
	Gaussian elimination.
	"""
	N = len(psi) - 1
	A = [0]*(N+1) # plain list to hold symbolic expressions
	b = [0]*(N+1)
	x = sp.Symbol('x')
	#print '...evaluating matrix...',
	for i in range(N+1):
	#	print '(%d,%d)' % (i, i)
		# Assume orthogonal psi can be be integrated symbolically...
		A[i] = sp.integrate(psi[i]**2, (x, Omega[0], Omega[1]))
		# Fallback on numerical integration if f*psi is too difficult
		# to integrate
		integrand = psi[i]*f
		I = sp.integrate(integrand, (x, Omega[0], Omega[1]))
		if isinstance(I, sp.Integral):
	#		print 'numerical integration of', integrand
			integrand = sp.lambdify([x], integrand)
			I = sp.mpmath.quad(integrand, [Omega[0], Omega[1]])
		b[i] = I
	#print 'A:\n', A, '\nb:\n', b
	c = [b[i]/A[i] for i in range(len(b))]
	#print 'coeff:', c
	u = 0
	for i in range(len(psi)):
		u += c[i]*psi[i]
	# Alternative:
	# u = sum(c[i,0]*psi[i] for i in range(len(psi)))
	#print 'approximation:', u
	return u, c

######################################

from matplotlib import animation

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

fig = plt.figure()
ax = plt.axes(xlim=(0, 2*np.pi), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

#N=4
#x = sp.Symbol('x')
#psi = []
#for n in xrange(N):
#		psi.append(sp.sin((2*n+1)*x))

psi = []
def animate(i):
	Omega = [0,2*np.pi]
	x = sp.Symbol('x')
	f = sp.tanh(20*(x-np.pi))
	
	
	psi.append(sp.sin((2*i+1)*x))
	
	u, c = least_squares_orth(f=f, psi = psi, Omega=Omega)
	##################################
	

	f = sp.lambdify([x], f, modules="numpy")
	u = sp.lambdify([x], u, modules="numpy")
	
	xcoor = np.linspace(Omega[0], Omega[1], 1001)
	exact = f(xcoor)
	approx = u(xcoor)
	#plt.plot(xcoor, approx, '-')
	#plt.hold('on')

	

	line.set_data(xcoor,approx)
	return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=20, blit=True)
#plt.plot(xcoor, exact, '--')

#print len(exact), len(approx), len(xcoor)

plt.xlabel('x')
plt.show()
"""


def least_squares_orth(f, psi, Omega):
	""
	Same as least_squares, but for orthogonal
	basis such that one avoids calling up standard
	Gaussian elimination.
	"
	N = len(psi) - 1
	A = [0]*(N+1) # plain list to hold symbolic expressions
	b = [0]*(N+1)
	x = sp.Symbol('x')
	#print '...evaluating matrix...',
	for i in range(N+1):
	#	print '(%d,%d)' % (i, i)
		# Assume orthogonal psi can be be integrated symbolically...
		A[i] = sp.integrate(psi[i]**2, (x, Omega[0], Omega[1]))
		# Fallback on numerical integration if f*psi is too difficult
		# to integrate
		integrand = psi[i]*f
		I = sp.integrate(integrand, (x, Omega[0], Omega[1]))
		if isinstance(I, sp.Integral):
	#		print 'numerical integration of', integrand
			integrand = sp.lambdify([x], integrand)
			I = sp.mpmath.quad(integrand, [Omega[0], Omega[1]])
		b[i] = I
	#print 'A:\n', A, '\nb:\n', b
	c = [b[i]/A[i] for i in range(len(b))]
	#print 'coeff:', c
	u = 0
	for i in range(len(psi)):
		u += c[i]*psi[i]
	# Alternative:
	# u = sum(c[i,0]*psi[i] for i in range(len(psi)))
	#print 'approximation:', u
	return u, c

######################################

from matplotlib import animation

psi_list = []

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

fig = plt.figure()
ax = plt.axes(xlim=(0, 2*np.pi), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)

N=4
x = sp.Symbol('x')
psi = []
for n in xrange(N):
		psi.append([sp.sin((2*n+1)*x)])


def animate(i):
	Omega = [0,2*np.pi]
	x = sp.Symbol('x')
	f = sp.tanh(20*(x-np.pi))
	#,sp.sin((2*1+1)*x),sp.sin((2*2+1)*x)]

	psiva = psi

	u, c = least_squares_orth(f=f, psi = psiva[i], Omega=Omega)
	##################################3333


	f = sp.lambdify([x], f, modules="numpy")
	u = sp.lambdify([x], u, modules="numpy")

	xcoor = np.linspace(Omega[0], Omega[1], 1001)
	exact = f(xcoor)
	approx = u(xcoor)
         
        psi_list.append(approx)

        for delta in range(0, len(psi_list)):
            value += psi_list[delta] 

	line.set_data(xcoor, value)
	return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=20, blit=True)
#plt.plot(xcoor, exact, '--')

#print len(exact), len(approx), len(xcoor)

plt.xlabel('x')
plt.show()
"""