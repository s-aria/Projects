from matplotlib.pyplot import *
from numpy import *


f = open("/home/Ace/Documents/Computational-Physics/P3/data")

lines = f.readlines()
f.close()

x = []
y = []


x1 = []
y1 = []

x2 = []
y2 = []

x3 = []
y3 = []

x4 = []
y4 = []

x5 = []
y5 = []

x6 = []
y6 = []

x7 = []
y7 = []

x8 = []
y8 = []
#--------------

for line in lines:
	p = line.split()
	x.append(float(p[0]))
	y.append(float(p[1]))
	
	x1.append(float(p[2]))
	y1.append(float(p[3]))
	
	x2.append(float(p[4]))
	y2.append(float(p[5]))

	x3.append(float(p[6]))
	y3.append(float(p[7]))

	x4.append(float(p[8]))
	y4.append(float(p[9]))

	x5.append(float(p[10]))
	y5.append(float(p[11]))

	x6.append(float(p[12]))
	y6.append(float(p[13]))

	x7.append(float(p[14]))
	y7.append(float(p[15]))

	x8.append(float(p[16]))
	y8.append(float(p[17]))


#--------------c---------------
x = array(x)
y = array(y)

x1 = array(x1)
y1 = array(y1)

x2 = array(x2)
y2 = array(y2)

x3 = array(x3)
y3 = array(y3)

x4 = array(x4)
y4 = array(y4)

x5 = array(x5)
y5 = array(y5)

x6 = array(x6)
y6 = array(y6)

x7 = array(x7)
y7 = array(y7)

x8 = array(x8)
y8 = array(y8)

#------------------------------

plot(x,y, label="Earth's orbit")

plot(x1, y1, label="Jupiter's orbit")

plot(x2, y2, label="Mars's orbit")

plot(x3, y3, label="Venus's orbit")
plot(x4, y4, label="Saturn's orbit")
plot(x5, y5, label="Mercury's orbit")
plot(x6, y6, label="Uranus's orbit")
plot(x7, y7, label="Neptune's orbit")
plot(x8, y8, label="Pluto's orbit")

legend(bbox_to_anchor=(0,0,1,1))
title('All planets + Pluto')
xlabel('x [AU]')
ylabel('y [AU]')
show()