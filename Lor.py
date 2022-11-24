# import libraries
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Dt = 0.005              # timestep 

t_start = 0             # starttime
t_end = 4              # endtime 4 secondi
n_steps = int(round((t_end-t_start)/Dt))   # number of timesteps

X_arr = np.zeros(n_steps+1)    # create an array of zeros 
Y_arr = np.zeros(n_steps+1)    # create an array of zeros
Z_arr = np.zeros(n_steps+1)    # create an array of zeros 
t_arr = np.zeros(n_steps+1)    # create an array of zeros 
t_arr[0] = t_start              # add starttime to array
X_arr[0] = 9             	# add initial value 
Y_arr[0] = 10              	# add initial value 
Z_arr[0] = 20             	# add initial value

#parameters       
sigma=10
b=8/3
r=28

# Euler's method
for j in range (1, n_steps+1):
   
   X = X_arr[j-1] 
   Y = Y_arr[j-1]
   Z = Z_arr[j-1]
   t = t_arr[j-1]

   dXdt = sigma*(Y-X)               # calculate the derivative of X
   dYdt = (r*X) -(X*Z)-Y            # calculate the derivative of Y
   dZdt = (X*Y) - (b*Z)             # calculate the derivative of Z

   X_arr[j] = X + Dt*dXdt           # calc. X at next timestep,add to array
   Y_arr[j] = Y + Dt*dYdt           # calc. Y at next timestep,add to array
   Z_arr[j] = Z + Dt*dZdt           # calc. Z at next timestep,add to array
   t_arr[j] = t + Dt                # add new value of t to array


epsilon=np.zeros(100)
for k in range (1,100):				#k=100 valori di epsilon
	epsilon[k]=np.random.uniform(-0.75,0.75)


x=np.zeros((100,n_steps+1), dtype=np.float64)
y=np.zeros((100,n_steps+1), dtype=np.float64)
z=np.zeros((100,n_steps+1), dtype=np.float64)

R= np.zeros((100,n_steps+1))
y[k,0]= 10 
z[k,0]= 20
x[k,0]= 9+epsilon[k]
R[k,0]=epsilon[k]
for k in range (1,100):				#k=100 valori di epsilon
	y[k,0]= 10 
	z[k,0]= 20
	x[k,0]= 9+epsilon[k]
	R[k,0]=epsilon[k]
	 
	for j in range (1, n_steps+1):   		#800 steps-> stop at time=4

			x[k,j] = x[k,j-1]+ Dt*(sigma*(y[k,j-1]-x[k,j-1])) 
			y[k,j] = y[k,j-1]+ Dt*((r*x[k,j-1]) -(x[k,j-1]*z[k,j-1])-y[k,j-1])
			z[k,j] = z[k,j-1]+ Dt*((x[k,j-1]*y[k,j-1]) - (b*z[k,j-1]))
			t_arr[j] = t_arr[j-1] + Dt 
			x_sum_sq= np.power((X_arr[j]-x[k,j]),2)
			y_sum_sq= np.power((Y_arr[j]-y[k,j]),2)
			z_sum_sq= np.power((Z_arr[j]-z[k,j]),2)

			R[k,j] =np.sqrt(x_sum_sq + y_sum_sq + z_sum_sq)
 
S=np.zeros(n_steps+1)
XP=np.zeros(n_steps+1)
YP=np.zeros(n_steps+1)
ZP=np.zeros(n_steps+1)
RMSE=np.zeros(n_steps+1)
	
for j in range (n_steps+1):
	for k in range (100):
		S[j] += R[k,j]/100 
		XP[j] += x[k,j]/100 
		YP[j] += y[k,j]/100 
		ZP[j] += z[k,j]/100 

	x_sum_sq= np.power((X_arr[j]-XP[j]),2)
	y_sum_sq= np.power((Y_arr[j]-YP[j]),2)
	z_sum_sq= np.power((Z_arr[j]-ZP[j]),2)

	RMSE[j]= np.sqrt(x_sum_sq + y_sum_sq + z_sum_sq)
 
S_4=np.zeros(n_steps+1)
g=np.zeros((100, n_steps+1))
for j in range (n_steps+1):
	g=0
	for k in range (100):
		g += np.power((x[k,j]-XP[j]),2)/100
	S_4[j]= np.sqrt(g)

#RMSE CAMPO MEDIO
plt.figure() 		
plt.plot(t_arr,S, color="green", label= "media dei RMSE") 
plt.plot(t_arr,RMSE, color="blue", label= "RMSE del campo medio") 
plt.plot(t_arr,S_4, color="red", label= "RMSE terzo caso") 
plt.xlabel("t")
plt.grid()
plt.legend(loc= "upper left")

plt.tight_layout()
plt.show()

















