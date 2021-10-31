from dolfin import *
from numpy import cosh, cos
set_log_active(False)
from math import log as ln
from mshr import *

dpdx = Constant(-0.01)
mu = Constant(0.01)
a = 2
b = 1
c = 0
# Much faster C++ version
ue_code = '''
class U : public Expression
{
  public:

    double a, b, mu, dpdx;

  void eval(Array<double>& values, const Array<double>& x) const
    {
      double u = 0.;
      double factor = 1.0/(2*mu)*(-dpdx)*(a*a*b*b)/(a*a + b*b);
      u = 1.0-(x[0]*x[0])/(a*a) - (x[1]*x[1])/(b*b);
      values[0] = u*factor;      
    }
};'''

# assigning the parametres into the c++ code?
u_c = Expression(ue_code)
u_c.a = float(a); u_c.b = float(b)
u_c.mu = float(mu(0)); u_c.dpdx = float(dpdx(0))

mesh = generate_mesh(Ellipse(Point(c),a ,b, 32), 32)

def main(N, degree=1):
  mesh = generate_mesh(Ellipse(Point(c),a ,b, N), N)
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

  if N==5 or N==20 or N==80:
    plot(u_, title="Numerical")
    plot(u_e, title="Exact")
    interactive()
  return u_error, mesh.hmin()

E = []; h = []; degree = 3
for n in [5, 10, 20, 40, 80]:
  ei, hi = main(n, degree=degree)
  E.append(ei)
  h.append(hi)

for i in range(1, len(E)):   
    r = ln(E[i]/E[i-1])/ln(h[i]/h[i-1])
    print "h=%2.2E E=%2.2E r=%.2f" %(h[i], E[i], r)

