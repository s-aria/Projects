from main import *
"""
Here dim1 corresponds to 1D, dim2 to 2D and error corresponds to finding the error or the K
The boolean expressions can be set to either False or True depending on what we are after.
"""

def test_exact(dim1=False, error=False, dim2=True):
	u0 =  Expression('cos(pi*x[0])',t=0)
	def alpha(u):
		return 1.0
	f = Constant(0.0)
	rho = 1.0
	T = 0.8
	dt = 0.003
	u_exact = Expression('exp(-pi*pi*t)*cos(pi*x[0])', t=T)
	
	if dim2==True:
		u, u_e = solver(T=T, dt=dt, mesh_points=32, u0=u0, alpha=alpha, f=f, rho=rho, u_exact=u_exact, d=2, p=1)
		plot(u,title="FEM")
		plot(u_e,title="Exact")
		interactive()

	if error == True:
		Ns = range(2,32)
		#list containing the errors for the various meshpoints
		E_list = []
		for N in Ns:
			dt = (1.0/N)**2	
			u, u_e = solver(T=T, dt=dt, mesh_points=N, u0=u0, alpha=alpha, f=f, rho=rho, u_exact=u_exact, d=2, p=1)
			e = u_e.vector().array() - u.vector().array()
			E = std.sqrt(std.sum(e**2)/u.vector().array().size)
			E_list.append(E)
		K_list = []
		#making another loop to divide by dt and contain them all in a list
		for i in xrange(len(E_list)):
			K = E_list[i]/dt
			K_list.append(K)
		
		plt.plot(Ns,K_list)
		plt.title("Error plot")
		plt.xlabel('$N$')
		plt.ylabel('E/h')
		plt.show()

	#in case of wanting to plot 1D
	if dim1 == True:
		u, u_e = solver(T=T, dt=dt, mesh_points=32, u0=u0, alpha=alpha, f=f, rho=rho, u_exact=u_exact, d=1, p=1)
		#converting dolfin expression to vectors
		u, u_e = u.vector().array(), u_e.vector().array() 
		x_mesh = std.linspace(0,1,len(u_e)) #making x axis
		line1, =plt.plot(x_mesh, u[::-1],'r') #the numerical solution
	

		plt.hold('on')
	
		line2, =plt.plot(x_mesh,u_e[::-1],'b',label="Exact")
		
		plt.xlabel('$x$')
		plt.ylabel('u(x)')
		
		plt.title('dt=%.5f , T=%.5f' %(dt,T))
		# Create a legend.
		plt.legend([line1, line2], ["FEM", "Exact"])
		plt.show()





if __name__ == '__main__':
	test_exact()