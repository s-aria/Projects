"""
Approximation of functions by linear combination of basis functions in
function spaces and the least squares method or the collocation method
for determining the coefficients.
"""
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

#import scitools.std as plt


def least_squares(f, psi, Omega):
	"""
	Given a function f(x) on an interval Omega (2-list)
	return the best approximation to f(x) in the space V
	spanned by the functions in the list psi.
	"""
	N = len(psi) - 1
	A = sp.zeros((N+1, N+1))
	b = sp.zeros((N+1, 1))
	x = sp.Symbol('x')
	print '...evaluating matrix...',
	for i in range(N+1):
		for j in range(i, N+1):
			print '(%d,%d)' % (i, j)
			integrand = psi[i]*psi[j]
			I = sp.integrate(integrand, (x, Omega[0], Omega[1]))
			if isinstance(I, sp.Integral):
				# Could not integrate symbolically, use numerical int.
				print 'numerical integration of', integrand
				integrand = sp.lambdify([x], integrand)
				I = sp.mpmath.quad(integrand, [Omega[0], Omega[1]])
			A[i,j] = A[j,i] = I
		integrand = psi[i]*f
		I = sp.integrate(integrand, (x, Omega[0], Omega[1]))
		if isinstance(I, sp.Integral):
			# Could not integrate symbolically, use numerical int.
			print 'numerical integration of', integrand
			integrand = sp.lambdify([x], integrand)
			I = sp.mpmath.quad(integrand, [Omega[0], Omega[1]])
		b[i,0] = I
	print
	print 'A:\n', A, '\nb:\n', b
	c = A.LUsolve(b) # symbolic solve
	print 'coeff:', c
	# c is a sympy Matrix object, numbers are in c[i,0]
	u = sum(c[i,0]*psi[i] for i in range(len(psi)))
	print 'approximation:', u
	return u, [c[i,0] for i in range(len(c))]



def comparison_plot(f, u, Omega, filename='tmp', plot_title='', ymin=None, ymax=None, u_legend='approximation'):
	"""Compare f(x) and u(x) for x in Omega in a plot."""
	x = sp.Symbol('x')
	print 'f:', f
	f = sp.lambdify([x], f, modules="numpy")
	u = sp.lambdify([x], u, modules="numpy")
	if len(Omega) != 2:
		raise ValueError('Omega=%s must be an interval (2-list)' % str(Omega))
	# When doing symbolics, Omega can easily contain symbolic expressions,
	# assume .evalf() will work in that case to obtain numerical
	# expressions, which then must be converted to float before calling
	# linspace below
	if not isinstance(Omega[0], (int,float)):
		Omega[0] = float(Omega[0].evalf())
	if not isinstance(Omega[1], (int,float)):
		Omega[1] = float(Omega[1].evalf())
	resolution = 601 # no of points in plot
	xcoor = np.linspace(Omega[0], Omega[1], resolution)
	# Vectorized functions expressions does not work with
	# lambdify'ed functions without the modules="numpy"
	exact = f(xcoor)
	approx = u(xcoor)
	plt.plot(xcoor, approx, '-')
	plt.hold('on')
	plt.plot(xcoor, exact, '--')
	plt.legend([u_legend, 'exact'])
	plt.title(plot_title)
	plt.xlabel('x')
	if ymin is not None and ymax is not None:
		plt.axis([xcoor[0], xcoor[-1], ymin, ymax])
	plt.savefig(filename + '.png')
	plt.show()


x = sp.Symbol('x')
f = sp.exp(-x)

N2 = [1,x,x**2] # N = 2
N4 = [1,x,x**2, x**3, x**4] #N = 4
N6 = [1,x,x**2, x**3, x**4, x**5, x**6] #N = 6


u, c = least_squares(f=f, psi = N2, Omega = [0, 8])
print u,c
comparison_plot(f,u,Omega=[0,8])