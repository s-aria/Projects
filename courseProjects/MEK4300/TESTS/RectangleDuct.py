from dolfin import *
from numpy import cosh, cos
set_log_active(False)

dpdx = Constant(-0.01)
mu = Constant(0.01)
a = 2
b = 1

factor = 16.*a*a/mu(0)/pi**3*(-dpdx(0))
class u_exact(Expression):
        
    def eval(self, value, x):
        u = 0
        for i in range(1, 500, 2):
            u += ((-1)**((i-1)/2)*(1-cosh(i*pi*x[1]/2./a)/
                   cosh(i*pi*b/2/a))*cos(i*pi*x[0]/2./a)/i**3)
        value[0] = u * factor        

# Much faster C++ version
ue_code = '''
class U : public Expression
{
  public:

    double a, b, mu, dpdx;

  void eval(Array<double>& values, const Array<double>& x) const
    {
      double u = 0.;
      double factor = 16.*a*a/mu/pow(DOLFIN_PI, 3)*(-dpdx);
      for (std::size_t i=1; i<600; i=i+2)
        u += pow(-1, (i-1)/2 % 2)*(1.-cosh(i*DOLFIN_PI*x[1]/2./a)/
           cosh(i*DOLFIN_PI*b/2./a))*cos(i*DOLFIN_PI*x[0]/2./a)/pow(i, 3);
      values[0] = u*factor;      
    }
};'''
u_c = Expression(ue_code)
u_c.a = float(a); u_c.b = float(b)
u_c.mu = float(mu(0)); u_c.dpdx = float(dpdx(0))

def main(N, degree=1):
    mesh = RectangleMesh(-a, -b, a, b, N, N)
    V = FunctionSpace(mesh, 'CG', degree)
    u = TrialFunction(V)
    v = TestFunction(V)
    F = inner(grad(u), grad(v))*dx + 1/mu*dpdx*v*dx
    bc = DirichletBC(V, Constant(0), DomainBoundary())
    u_ = Function(V)
    solve(lhs(F) == rhs(F), u_, bcs=bc)

    #u_e = interpolate(u_exact(), V)
    u_e = interpolate(u_c, V)
    bc.apply(u_e.vector())
    u_error = errornorm(u_e, u_, degree_rise=0)
    return u_error, mesh.hmin()

E = []; h = []; degree = 1
for n in [5, 10, 20, 40, 80]:
    ei, hi = main(n, degree=degree)
    E.append(ei)
    h.append(hi)
    
from math import log as ln
for i in range(1, len(E)):
    r = ln(E[i]/E[i-1])/ln(h[i]/h[i-1])
    print "h=%2.2E E=%2.2E r=%.2f" %(h[i], E[i], r)



