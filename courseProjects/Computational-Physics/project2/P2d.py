from matplotlib.pyplot import *
from numpy import *


f = open("/home/Ace/Documents/Computational-Physics/P5/data")

lines = f.readlines()
f.close()

v1 = []
#v2 = []
rho = [] 

for line in lines:
	p = line.split()
	#print p
	rho.append(float(p[0]))
	v1.append(float(p[1]))
	#v2.append(float(p[2]))



#--------------c---------------
#print v1
v1 = array(v1)
#v2 = array(v2)
rho = array(rho)
#u = array(u)
#print rho
#epsi=mav1(log(abs(rho[:]-u[:])/u[:]))
#print epsi
#------------------------------
#print rho[-1], v1[-1]

plot(rho,v1)
#hold('ON')
#plot(rho,v2) #, label="Numerical solution")
#title('Comparison graph')
#xlabel('v1')
#ylabel('rho')
show()