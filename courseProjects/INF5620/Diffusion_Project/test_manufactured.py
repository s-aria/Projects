from main import *


def test_manufactured(graph=False):
	u0 =  Constant(0.0)
	
	def alpha(u):
		return 1 + u**2

	rho = 1.0
	
	T = [0.08,0.01,0.04,2.8]

	for i in T:
		f = Expression('-(rho*pow(x[0],3))/3 + (rho*pow(x[0],2))/2.0 + \
					   (8.0*pow(t,3)*pow(x[0],7))/9.0 - (28.0*pow(t,3)*pow(x[0],6))/9.0 + \
					   (7.0*pow(t,3)*pow(x[0],5))/2.0 - (5.0*pow(t,3)*pow(x[0],4))/4.0 +\
					    2.0*t*x[0] - t', t=i,rho=rho)
		dt = 0.0001
		u_exact = Expression('t*x[0]*x[0]*(0.5 - x[0]/3)', t=i)
		u, u_e = solver(T=i, dt=dt, mesh_points=32, u0=u0, alpha=alpha, f=f, rho=rho, u_exact=u_exact, d=1, p=1)
		e = abs(u_e.vector().array() - u.vector().array())
		print "The error: [%.10E] for T=%.4f and dt=%.4f" %(e[0],i,dt)

	if graph==True:
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
	test_manufactured()