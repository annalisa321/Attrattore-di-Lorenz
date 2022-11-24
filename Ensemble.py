# import libraries
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

Dt = 0.005              # timestep 

t_start = 0             # starttime
t_end = 60              # endtime
n_steps = int(round((t_end-t_start)/Dt))   # number of timesteps

X_arr = np.zeros(n_steps +1)    # create an array of zeros 
Y_arr = np.zeros(n_steps + 1)   # create an array of zeros 
Z_arr = np.zeros(n_steps +1)    # create an array of zeros 
t_arr = np.zeros(n_steps + 1)   # create an array of zeros 
t_arr[0] = t_start              # add start time to array
X_arr[0] = 9             	# add initial value 
Y_arr[0] = 10              	# add initial value 
Z_arr[0] = 20             	# add initial value

X_arr2 = np.zeros(n_steps +1) 
X1 = np.zeros(n_steps +1)
X2 = np.zeros(n_steps +1)
X3 = np.zeros(n_steps +1)
X_arr2[0] = 9 + 1e-10 
X1[0] = 9+1e-5
X2[0] = 9+1e-3
X3[0]= 10

Y_arr2 = np.zeros(n_steps + 1)
Y1 = np.zeros(n_steps +1)
Y2 = np.zeros(n_steps +1)
Y3 = np.zeros(n_steps +1)
Y1[0] = Y2[0] = Y3[0]= Y_arr2[0] = 10

Z_arr2 = np.zeros(n_steps +1)
Z1 = np.zeros(n_steps +1)
Z2 = np.zeros(n_steps +1)
Z3 = np.zeros(n_steps +1)
Z1[0] = Z2[0] = Z3[0]= Z_arr2[0]= 20
 
#parameters       
sigma=10
b=8/3
r=28

#Metodo di Eulero
for i in range (1, n_steps + 1):
   
   X = X_arr[i-1] 
   Y = Y_arr[i-1]
   Z = Z_arr[i-1]
   t = t_arr[i-1]

   dXdt = sigma*(Y-X)               # calculate the derivative of X
   dYdt = (r*X) -(X*Z)-Y            # calculate the derivative of Y
   dZdt = (X*Y) - (b*Z)             # calculate the derivative of Z

   X_arr[i] = X + Dt*dXdt           # calc. X at next timestep,add to array
   Y_arr[i] = Y + Dt*dYdt           # calc. Y at next timestep,add to array
   Z_arr[i] = Z + Dt*dZdt           # calc. Z at next timestep,add to array
   t_arr[i] = t + Dt                # add new value of t to array

   #CALCOLO DEI CAMPI PERTURBATI
   X1[i] = X1[i-1]+ Dt*(sigma*(Y1[i-1]-X1[i-1])) 
   X2[i] = X2[i-1]+ Dt*(sigma*(Y2[i-1]-X2[i-1])) 
   X3[i] = X3[i-1]+ Dt*(sigma*(Y3[i-1]-X3[i-1]))

   Y1[i] = Y1[i-1]+ Dt*((r*X1[i-1]) -(X1[i-1]*Z1[i-1])-Y1[i-1])
   Y2[i] = Y2[i-1]+ Dt*((r*X2[i-1]) -(X2[i-1]*Z2[i-1])-Y2[i-1])
   Y3[i] = Y3[i-1]+ Dt*((r*X3[i-1]) -(X3[i-1]*Z3[i-1])-Y3[i-1])

   Z1[i] = Z1[i-1]+ Dt*((X1[i-1]*Y1[i-1]) - (b*Z1[i-1]))
   Z2[i] = Z2[i-1]+ Dt*((X2[i-1]*Y2[i-1]) - (b*Z2[i-1]))
   Z3[i] = Z3[i-1]+ Dt*((X3[i-1]*Y3[i-1]) - (b*Z3[i-1]))             

   X_arr2[i] = X_arr2[i-1] + Dt*(sigma*(Y_arr2[i-1]-X_arr2[i-1]))           
   Y_arr2[i] = Y_arr2[i-1] + Dt*((r*X_arr2[i-1]) -(X_arr2[i-1]*Z_arr2[i-1])-Y_arr2[i-1])         
   Z_arr2[i] = Z_arr2[i-1] + Dt*((X_arr2[i-1]*Y_arr2[i-1]) - (b*Z_arr2[i-1]) )            

            
#RMSE  
R_arr=np.zeros(n_steps + 1) 		
R_arr=((X_arr-X_arr2)**2) +((Y_arr-Y_arr2)**2)+((Z_arr-Z_arr2)**2)**(0.5)
R1 = np.zeros(n_steps + 1)
R1=((X_arr-X1)**2) +((Y_arr-Y1)**2)+((Z_arr-Z1)**2)**(0.5) 
R2 = np.zeros(n_steps + 1) 
R2=((X_arr-X2)**2) +((Y_arr-Y2)**2)+((Z_arr-Z2)**2)**(0.5) 
R3 = np.zeros(n_steps + 1) 
R3=((X_arr-X3)**2) +((Y_arr-Y3)**2)+((Z_arr-Z3)**2)**(0.5)    

a=0.5
ta=np.zeros(4)
for i in range(n_steps):
	if R_arr[i]>a:
	   ta[0]=t_arr[i]
	   break
for i in range(n_steps):
	if R1[i]>a:
	   ta[1]=t_arr[i]
	   break
for i in range(n_steps):
	if R2[i]>a:
	   ta[2]=t_arr[i]
	   break
for i in range(n_steps):
	if R3[i]>a:
	   ta[3]=t_arr[i]
	   break

d = [ [1e-10, ta[0]],
      [1e-5,  ta[1]],
      [1e-3,  ta[2]],
      [1,     ta[3]]]		

#Plot con r=28
fig,(ax)=plt.subplots() 		
ax.plot(X_arr,Z_arr, 'k') 
ax.set_title('Plot con r=28')
ax.set_ylabel('z')
ax.set_xlabel('x')
ax.grid()

#Figura 3D
fig = plt.figure() 	
bx = plt.axes(projection='3d')
bx.grid()
bx.plot3D(X_arr,Y_arr,Z_arr,t_arr, 'k')
bx.set_title('3D Plot')
bx.set_xlabel('x', labelpad=60)
bx.set_ylabel('y', labelpad=60)
bx.set_zlabel('t', labelpad=60)


#X(t) difference plotting
fig,(cx)=plt.subplots()
cx.plot(t_arr,(X_arr2-X_arr)) 		
cx.set_title('differenza tra i valori di x(t)') 
cx.set_xlabel('t')
cx.set_ylabel('differenza tra x(t)')
cx.grid()

#RMSE a confronto
fig,(ax, bx)=plt.subplots(2,1) 
ax.plot(t_arr,R_arr, color="green", label= "RMSE 1e-10") 
bx.plot(t_arr,R1, color="blue", label= "RMSE 1e-5")
ax.legend()
bx.legend()
ax.set_xlabel("t")
bx.set_xlabel("t")
fig,(cx, dx)=plt.subplots(2,1)
cx.plot(t_arr,R2, color="red", label= "RMSE 1e-3")
dx.plot(t_arr,R3, color="violet", label= "RMSE 1e+0")
cx.set_xlabel("t")
cx.legend()
dx.legend()
dx.set_xlabel("t")
plt.tight_layout()

#x(t) a confronto
plt.figure() 
plt.plot(t_arr,X_arr, color="green", label= "X(t)") 
plt.plot(t_arr,X_arr2, color="blue", label= "X(t) perturbata")
plt.xlabel("t")
plt.legend(bbox_to_anchor =(0.75, 1.15))

#RMSE con diverse perturbazioni
plt.figure()
plt.semilogy(t_arr,R_arr, "green", label= "RMSE 1e-10") 
plt.xlabel("t")
plt.legend(bbox_to_anchor =(0.75, 1.15))

plt.figure()	
plt.semilogy(t_arr,R1,"blue" ,label= "RMSE 1e-5")
plt.xlabel("t")
plt.legend(bbox_to_anchor =(0.75, 1.15))

plt.figure()
plt.semilogy(t_arr,R2, color="red", label= "RMSE 1e-3")
plt.xlabel("t")
plt.legend(bbox_to_anchor =(0.75, 1.15))

plt.figure()
plt.semilogy(t_arr,R3, color="violet", label= "RMSE 1e+0")
plt.xlabel("t")
plt.legend(bbox_to_anchor =(0.75, 1.15))
plt.show()

#tabella
df = pd.DataFrame(d, columns = ['perturbazione','tempo con R>0.5'])
print(df)

plt.tight_layout()
plt.show()
