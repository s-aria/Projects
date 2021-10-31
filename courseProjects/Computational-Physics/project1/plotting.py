from matplotlib.pyplot import *
from numpy import *


f = open("/home/Ace/Documents/Computational-Physics/data")

lines = f.readlines()
f.close()

x = []
v = [] #numerical values
u = [] #analytical values

#Error for exercise c:
eps = zeros(len(x))
#--------------

for line in lines:
	p = line.split()
	x.append(float(p[0]))
	v.append(float(p[1]))
	u.append(float(p[2]))


#--------------c---------------
x = array(x)
v = array(v)
u = array(u)

epsi=max(log(abs(v[:]-u[:])/u[:]))
print epsi
#------------------------------

plot(x,v, label="Numerical solution")
plot(x, u, label="Analytical solution")
legend(bbox_to_anchor=(0,0,1,1))
title('Comparison graph')
xlabel('x')
ylabel('Function values')
show()