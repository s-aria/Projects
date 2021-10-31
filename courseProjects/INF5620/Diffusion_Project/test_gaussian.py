from main import *

"""
Here we have only written the expressions, 
please go to main.py and uncomment the animation part to see the plot in action.
"""
def test_gaussian():
	sigma = 1.0 #the width of the function
	u0 =  Expression('exp(-(1/(2*sigma*sigma))*(pow(x[0],2) + pow(x[1],2)))',sigma=sigma)
	
	def alpha(u):
		beta = 2.0
		return 1 + beta*u**2
	f = Constant(0.0)
	rho = 1.0
	T = 0.4
	dt = 0.01
	u_exact = None
	u = solver(T=T, dt=dt, mesh_points=16, u0=u0, alpha=alpha, f=f, rho=rho, u_exact=u_exact, d=2, p=1)


if __name__ == '__main__':
	test_gaussian()