from scitools.std import *
from matplotlib.pyplot import *


#u = T/theta
N = 1000
b = linspace(0.01,3,N)
zi = zeros(N)



for i in range(N):
	for j in range(N):
		zi[i] += sum((2*j+1)*exp(-j*(j+1)/b[i]))

plot(b,zi)
show()
		