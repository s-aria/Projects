from main import *

def test_constant():
	u0 =  Constant(1.0)
	def alpha(u):
		return 1.0
	f = Constant(0.0)
	rho = 1.0
	
	for d in range(1,4):
		for p in xrange(1,3):
			u = solver(T=1.0, dt=0.3, mesh_points=8, u0=u0, alpha=alpha, f=f, rho=rho, u_exact=None, d=d, p=p)
			u = u.vector().array()
			#exact solution i.e our initial condition, but we will just write here once again
			#instead of converting the dolfin Constant(1.0)
			u_ex = 1.0 
			error = abs(u_ex - u[1])
			#uncomment the plot if you wish to see the plots, but do comment out the part .write_png(...)
			#if you do not wish to save the figures.
			#plot(u,interactive=True)#.write_png("dolfin_plot_%.i" %d)
			print  "Error %.15E for P%.i elemnt in %.dD" %(error, p, d)


if __name__ == '__main__':
	test_constant()