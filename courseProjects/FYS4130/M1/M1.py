# from numpy import *
# from matplotlib.pyplot import *

# N = 1000
# p = [0.5,0.1,0.01]

# samples = 10**5

# X = zeros(samples)


# def createRW_1(N,p,samples,X):
# 	x = zeros(N)
# 	x[0] = 1

# 	for j in range(samples):
# 		ens = random.random(N)
# 		for i in range(len(x)-1):
# 			if ens[i] < p:
# 				x[i+1] = -x[i]
# 			else:
# 				x[i+1] = x[i]

# 		X[j] = sum(x)
	
# 	return X



# def createRW_2(N,p,samples,X):
# 	x = zeros(samples)
# 	x[0] = 1

# 	X2mean = zeros(N) # <X^2>
#  	Xmean2 = zeros(N) # <X>^2
#  	variance = zeros(N)
#  	X2mean[0] = 1
#  	Xmean2[0] = 1
# 	for j in range(1,N):
# 		ens = random.random(samples)
# 		for i in range(samples-1):
# 			if ens[i] < p:
# 				x[i+1] = -x[i]
# 			else:
# 				x[i+1] = x[i]

# 		X[j] = sum(x)
# 		X2mean[j] = mean(X**2)
# 		Xmean2[j] = mean(X)**2
# 		variance[j] = X2mean[j] - Xmean2[j]
# 	return X2mean, Xmean2, variance

# def calculations(N, p, samples, X,find_mean_square=False, D_value=False, plotting=False):
# 	for i in p:
# 		if find_mean_square==True:
# 			Xsum = createRW_1(N,i,samples,X)
# 			X2mean = mean(Xsum*Xsum)
# 			print "For p = %.3f we get <X^2> = %.2f" %(i, X2mean)
# 		#Finding the diffusion coefficient:
# 		if D_value==True and p==0.5:
# 			t = linspace(1,N,len(X2))
# 			D = X2/(2*t)
# 			print mean(D)
# 		if plotting==True:
# 			X2, X_2, DX = createRW_2(N,i,samples,X)
# 			# print len(X2), len(DX)
# 			t = linspace(1,N,len(X2))
# 			fig = figure()
# 			ax1 = fig.add_subplot(111)
# 			ax1.plot(log10(t),log10(DX))
# 			show()



# def updateX(X1,samples,p,x):
# 	ens = random.random(samples)
	
# 	for i in range(samples):
# 		if ens[i] < p:
# 			x[i] = -x[i] 
# 		else:
# 			x[i] = x[i] 
# 		X1[i] += x[i]
# 	return X1

# def buildRW(N,samples,p):
# 	x = zeros(samples)+1
# 	X1 = zeros(samples) + 1.0
# 	MeanX_sqr = zeros(N) #<X>**2
# 	Mean_Xsqr = zeros(N) #<X**2>
# 	MeanX_sqr[0] = 1
# 	Mean_Xsqr[0] = 1

# 	for j in range(N-1):
# 		X1 = updateX(X1,samples,p,x)
# 		MeanX_sqr[j+1] = mean(X1)**2
# 		Mean_Xsqr[j+1] = mean((X1)**2)
# 	var = Mean_Xsqr - MeanX_sqr
# 	return var

# def Dfactor(xvar,t,N):
# 	D=zeros(N)
# 	for i in range(0,len(D)):
# 		D[i]=xvar[i]/(2*t[i])
# 	return mean(D)

# X1 = createRW_1(N,0.5,samples,X)

# fig = figure(1)
# ax1 = fig.add_subplot(211)
# subplots_adjust(hspace=0.4)
# title("$p=0.5$")
# ax1.hist(X1, bins=30)

# X1 = createRW_1(N,0.9,samples,X)

# ax2 = fig.add_subplot(212)
# title("$p=0.9$")
# ax2.hist(X1, bins=30)
# show()

"""
#Uncomment for the plotting
t = linspace(1,N,N)
xvar = buildRW(N,samples,p)
#print Dfactor(xvar,t,N)
loglog(t,xvar,label="p=0.5")
# title("p = %g, N = %g, sample = %g" % (p,N, samples))
xlabel("log time") 
ylabel("log delta X^2")
legend()
hold('on')

p = 0.1
xvar = buildRW(N,samples,p)
loglog(t,xvar,label="p=0.1")
# title("p = %g, N = %g, sample = %g" % (p,N, samples))
xlabel("log time") 
ylabel("log delta X^2") 
legend()
hold('on')

p=0.01
xvar = buildRW(N,samples,p)
loglog(t,xvar,label="p=0.01")
xlabel("log time") 
ylabel("log delta X^2") 
legend()

show()
"""
from numpy import *
from matplotlib.pyplot import *
import matplotlib.mlab as mlab

N = 1000
p = 0.5 #[0.5,0.01,0.1]
samples = 10**5

def updateX(N,p,x):
	ens = random.random(N)
	
	for i in range(N-1):
		if ens[i] < p:
			x[i+1] = -x[i] 
		else:
			x[i+1] = x[i] 
	return sum(x)

def buildRW(N,samples,p):
	x0 = zeros(N)
	x0[0] = 1
	X = zeros(samples)

	for j in range(samples):
		X[j] = updateX(N,p,x0)

	return X

t = linspace(1,N,samples)
XN = buildRW(N,samples,p)
mu, sigma = 0, std(XN)
#n, x, patches = hist(XN, 35, normed=1, facecolor='blue', alpha=.60)


PX = lambda x: 1.0/(sigma*sqrt(2*pi))*exp(-x**2/(2*sigma**2))
"""
y = PX(x)

plot(x, y, 'r--', linewidth=2)
xlabel('$\mathrm{X_{N}}$')
ylabel('Probability P(X)')
title(r'$\mathrm{Histogram\ of\ X_{N}:}\ \mu=%g,\ \sigma=%g\ p = %g$'%(mu, sigma,p))
figure()"""

xb = sort(XN)
yb = PX(xb)

plot(sort(xb**2),log(yb))
xlabel('$\mathrm{X^2}$')
ylabel('log(P(X)')
print(sigma**2) #sigma/2 is where the maxima occurs)
show()