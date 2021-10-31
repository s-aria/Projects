from scitools.std import *
import random as rd
from scitools.numpyutils import compute_histogram

M = 10000
N = 50

S = zeros((M,N))

for i in range(M):
	for k in range(N):
		S[i,k] = rd.randint(0,1)
		if S[i,k] == 0:
			S[i,k] = -1
			
E = zeros(M)
for l in range(M):
	E[l] = sum(S[l,:])

x = linspace(0,1,M)

#plot(x, E)
sample = E
x,y=compute_histogram(sample, nbins=101, piecewise_constant=True)
plot(x,y)
"""sigma = max(y)#0.25
s=linspace(-25,25,len(x))
stirling=sigma*exp(-2*s**2/N)
hold("on")
plot(x,stirling)"""