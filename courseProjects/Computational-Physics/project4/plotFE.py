from matplotlib.pyplot import *
from numpy import *

f = open("/home/Ace/Documents/Computational-Physics/P4/data")

lines = f.readlines()
f.close()
columns = lines.T
x = []
print columns
for line in lines:
	p = line.split()
	x.append(float(p[0]))
	



x = array(x)
t = array(t)

plot(x,t)
#hold('ON')
#plot(rho,v2) #, label="Numerical solution")
#title('Comparison graph')
#xlabel('v1')
#ylabel('rho')
show()
