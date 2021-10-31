from scitools.std import *


T_hat = [1.15, 1.0, 0.85]

#V_hat = linspace(0.4, 20, 1000)
rho_hat = linspace(0,2,1000)

T1 = []
T2 = []
T3 = []

for i in T_hat:
	P_hat = 8*i*rho_hat/(3-rho_hat) - 3*rho_hat**2	
	#P_hat = 8*i/(3*V_hat - 1) - 3/V_hat**2
	if i == 1.15:
		T1.append(P_hat)
	if i == 1.0:
		T2.append(P_hat)
	if i == 0.85:
		T3.append(P_hat)

T1 = array(T1)
T2 = array(T2)
T3 = array(T3)

plot(rho_hat,T1,rho_hat,T2,rho_hat,T3,xlabel="rho",ylabel="T")