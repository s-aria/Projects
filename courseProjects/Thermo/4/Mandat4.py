from matplotlib.pyplot import *
from scitools.std import zeros,exp


n = 15

x = 3 #linspace(0,n,1000)

z = zeros(n)



for j in range(n):
	z[j] = (2*j+1)*exp((-j*(j+1))/(x))

plot(range(n),z)
show()