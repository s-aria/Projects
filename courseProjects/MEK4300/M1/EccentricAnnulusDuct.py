from dolfin import *
from numpy import log, sqrt, cosh, sinh, cos, sin
set_log_active(False)
from mshr import *

dpdx = Constant(-0.01)
mu = Constant(0.01)
a = 2
b = 0.5
c = 1

F = (a*a - b*b + c*c)/(2*c)
M = sqrt((F*F - a*a))
alpha = 0.5*log((F + M)/(F - M))
beta = 0.5 * log((F - c + M)/(F - c - M))

# Much faster C++ version
ue_code = '''
class U : public Expression
{
  public:

    double a, b, c, mu, dpdx, F, M, alpha, beta;
  void eval(Array<double>& values, const Array<double>& x) const
    {
      double Q = 0.;
      double factor = DOLFIN_PI/(8*mu)*(-dpdx);
      for (std::size_t n=1; n<600; n++)
        Q += (n * exp(-n*(beta + alpha)))/(sinh(n*beta - n*alpha));

      values[0] = factor * (pow(a,4) - pow(b,4) - (4 * c*c * M*M)/(beta - alpha) - 8 * c*c * M*M * Q);      
    }
};'''

u_c = Expression(ue_code)
u_c.a = float(a); u_c.b = float(b)
u_c.mu = float(mu(0)); u_c.dpdx = float(dpdx(0))
u_c.c = float(c); u_c.F = float(F); u_c.M = float(M)
u_c.alpha = float(alpha); u_c.beta = float(beta)

def main(N, degree=1):
    mesh = generate_mesh(Circle(Point(0.0,0.0), a, N) - Circle(Point(a-c), b, N), N)
    V = FunctionSpace(mesh, 'CG', degree)
    u = TrialFunction(V)
    v = TestFunction(V)
    F = inner(grad(u), grad(v))*dx + 1/mu*dpdx*v*dx
    bc = DirichletBC(V, Constant(0), DomainBoundary())
    u_ = Function(V)
    solve(lhs(F) == rhs(F), u_, bcs=bc)

    # u_e = interpolate(u_exact(), V)
    u_e = interpolate(u_c, V)
    bc.apply(u_e.vector())
    u_error = errornorm(u_e, u_, degree_rise=0)

    if N==5 or N==20 or N==80:
      plot(u_,title="Numerical")
      plot(u_e,title="Exact")
      interactive()
    return u_error, mesh.hmin()

E = []; h = []; degree = 2
for n in [5, 10, 20, 40, 80]:
    ei, hi = main(n, degree=degree)
    E.append(ei)
    h.append(hi)

from math import log as ln
for i in range(1, len(E)):
    r = ln(E[i]/E[i-1])/ln(h[i]/h[i-1])
    print "h=%2.2E E=%2.2E r=%.2f" %(h[i], E[i], r)



