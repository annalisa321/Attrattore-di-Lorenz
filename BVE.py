
###################################################
##############PARAMETERS###########################
################################################### 
DAYLEN  = 1   #  Forecast length in days.
DtHours = .5   #  Timestep in hours.
###################################################
##############IMPORT LIBRARIES#####################
################################################### 
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d
from matplotlib import ticker, cm   #COLORE CAMPO
from mpl_toolkits.axes_grid1 import make_axes_locatable  #COLORBAR
###################################################
##############DEFINE FUNCTIONS#####################
###################################################     

def make_Laplacian(Z):
    #  Compute the Laplacian 
    #of the geopotential height
    # within the boundary
    #(boundary excluded)
    M       = Z.shape[0]
    N       = Z.shape[1]
    Zxx  = np.zeros([M,N])     #  second x derivative of Z
    Zyy  = np.zeros([M,N])     #  second y derivative of Z
    # Compute within the domain (no boundaries) 
    # Second x-derivative of Z
    Zxx[1:M-1,:] = (Z[2:M,:]+Z[0:M-2,:]-2*Z[1:M-1,:])/(736e+3**2)
    # Second y-derivative of Z
    Zyy[:,1:N-1] = (Z[:,2:N]+Z[:,0:N-2]-2*Z[:,1:N-1])/(736e+3**2)    
    ##  Laplacian of height in the inner domain
    L0in = Zxx[1:M-1,1:N-1]+Zyy[1:M-1,1:N-1]
    return L0in

def make_Jacobian(Z,ABS_VOR):
    # Compute the Jacobian within 
    #the domain (boundary excluded)
    M       = Z.shape[0]
    N       = Z.shape[1]
    Zx    = np.zeros([M,N])     #  x derivative of Z
    Zy    = np.zeros([M,N])     #  y derivative of Z
    ABS_VORx  = np.zeros([M,N])     #  x derivative of ABS_VOR
    ABS_VORy  = np.zeros([M,N])     #  y derivative of ABS_VOR
    # x-derivative of Z
    Zx[1:M-1,:] = (Z[2:M,:]-Z[0:M-2,:])/(2*736e+3)
    # y-derivative of Z
    Zy[:,1:N-1] = (Z[:,2:N]-Z[:,0:N-2])/(2*736e+3)
    # x-derivative of the absolute vorticity 
    ABS_VORx[1:M-1,:] = (ABS_VOR[2:M,:]-ABS_VOR[0:M-2,:])/(2*736e+3)
    # y-derivative of the absolute vorticity 
    ABS_VORy[:,1:N-1] = (ABS_VOR[:,2:N]-ABS_VOR[:,0:N-2])/(2*736e+3)
    ##  Compute the Jacobian J(ABS_VOR,Z)
    Jacobi = ABS_VORx * Zy - ABS_VORy * Zx
    return Jacobi[1:M-1,1:N-1] 

def Poisson_solver(Jacobi):
    # Compute the height tendency from the L tendency
    # (that is from the Jacobian)
    # only within the domain (boundary excluded)    
    M       = Jacobi.shape[0]+2
    N       = Jacobi.shape[1]+2
    SM=np.zeros([M-2,M-2])
    SN=np.zeros([N-2,N-2])
    EIGEN=np.zeros([M-2,N-2])    
    ##  Coefficients for x-transformation
    for m1 in range(0,M-2):
     for m2 in range(0,M-2):
      SM[m1,m2] = np.sin(np.pi*(m1+1)*(m2+1)/(M-1))       
    ##  Coefficients for y-transformation
    for n1 in range(0,N-2):
     for n2 in range(0,N-2):
      SN[n1,n2] = np.sin(np.pi*(n1+1)*(n2+1)/(N-1))        
    ##  Eigenvalues of Laplacian operator
    for mm in range(0,M-2):
     for nn in range(0,N-2):
      eigen = (np.sin(np.pi*(mm+1)/(2*(M-1))))**2 +(np.sin(np.pi*(nn+1)/(2*(N-1))))**2
      EIGEN[mm,nn] = (-4/736e+3**2) * eigen
    #  Tendency values in interior.
    Ldot = Jacobi
    #  Compute the transform of the solution
    LDOT = np.dot(SM,np.dot(Ldot,SN))
    #  Convert transform of d(xi)/dt to transform of d(Z)/dt
    ZDOT = LDOT / EIGEN 
    #  Compute inverse transform to get the height tendency.
    Zdot = (4/((M-1)*(N-1))) *np.dot(SM,np.dot(ZDOT,SN))
    return Zdot
      
 
def make_f_and_h(N,M,Xp,Yp):
    FCOR=np.zeros([M,N])
    h=np.zeros([M,N])    
    a = (4*10**7)/(2*np.pi)      #  Radius of the Earth
    grav = 9.80665           #  Gravitational acceleration
    Omega = 2*np.pi/(24*60*60)  #  Angular velocity of Earth.
    ##  Compute Coriolis Parameter and Map Factor
    ##  and parameter h = g*m**2/f used in the BVE
    for ny in range(0,N):
     for nx in range(0,M):
      xx = (nx-Xp)*736e+3
      yy = (ny-Yp)*736e+3
      rr = np.sqrt(xx**2+yy**2)
      phi = 2*((np.pi/4)-np.arctan(rr/(2*a)))#Latitude
      mapPS = 2 / (1+np.sin(phi))
      f = 2*Omega*np.sin(phi)
      FCOR[nx,ny] = f
      h[nx,ny] = grav * mapPS**2 / f
    return FCOR,h

###################################################
##############FIXED PARAMETERS###########################
###################################################   
M  = 19 # Points in x direction
N  = 16 # Points in y direction
Xp = 8 # Coord. of North Pole
Yp = 12 # Coord. of North Pole

###################################################
###############COORDINATES AND TIME################
###################################################
daylen = DAYLEN                   #  Integration time (in days)
seclen = int(daylen*24*60*60)     #  Integration time (in secon736e+3)
Dt = DtHours*60*60                #  Timestep in secon736e+3
nt = int(seclen//Dt )             #  Total number of time-steps.
# Define the (X,Y) grid (for plotting)
X, Y  = np.meshgrid(np.linspace(1,M,M),np.linspace(1,N,N))
X = np.transpose(X)
Y = np.transpose(Y)
#Coriolis and map factor
FCOR,h=make_f_and_h(N,M,Xp,Yp)

###################################################
###############DEFINE WORKING ARRAYS###############
###################################################
Zout=np.zeros([nt+1,M,N])   
L0=np.zeros([nt+1,M,N])
###################################################
###############READ INPUT DATA##########################
###################################################
###   Read and plot the initial and verification height data
## Input files
File1 = 'Case1-1949010503.z00' #The initial value 
File2 = 'Case1-1949010603.z00' #The final value 
Z0  = np.genfromtxt(File1)
Z24 = np.genfromtxt(File2)

Zout[0,:,:]  = Z0      #  Copy initial height field

##  Compute the Laplacian of the height/psi
#################################################
#################################################
############      MAIN LOOP     #################
#################################################
#################################################
###                                           ###
###      Integrate the BVE in time            ###
###                                           ###
#################################################
#################################################
######## Start of the time-stepping loop  #######

L0[0,1:M-1,1:N-1]=make_Laplacian(Z0) #calcolo L0 da Z0 

#boundary conditions
L0[0,0,1:N-2]=2*L0[0,1,1:N-2]-L0[0,2,1:N-2] 
L0[0,M-1,1:N-2]=2*L0[0,M-2,1:N-2]-L0[0,M-3,1:N-2] 
L0[0,0:M-1,0]=2*L0[0,0:M-1,1]-L0[0,0:M-1,2] 
L0[0,0:M-1,N-1]=2*L0[0,0:M-1,N-2]-L0[0,0:M-1,N-3]

Ldot=np.zeros([M-2,N-2])
Zdot=np.zeros([M-2,N-2])

for n in range(0,nt): 
	  
 	if n==0: #eulero  
  		Ldot[:,:]=make_Jacobian(Z0,h*L0[0,:,:]+FCOR) #DERIVATA_T L
  		Zdot[:,:]=Poisson_solver(Ldot) #DERIVATA_T Z
  		Zout[1,:,:]=Zout[0,:,:] #CONTORNO INCREMENTATO AL TEMPO 1
  		Zout[1,1:M-1,1:N-1]=Zout[0,1:M-1,1:N-1]+Dt*Zdot[:,:] #EULERO 
  		L0[1,1:M-1,1:N-1]=make_Laplacian(Zout[1,:,:])  

	#leapfrog method 
 	else:
  		Ldot[:,:]=make_Jacobian(Zout[n,:,:],h*L0[n,:,:]+FCOR)
  		Zdot[:,:]=Poisson_solver(Ldot)
  		Zout[n+1,:,:]=Zout[n,:,:] #INCREMENTO I BORDI AL TSTEP SUCCESSIVO
  		Zout[n+1,1:M-1,1:N-1]=Zout[n-1,1:M-1,1:N-1]+2*Dt*Zdot[:,:] #LEAPFROG METHOD
  		L0[n+1,:,:]=L0[0,:,:] #INCREMENTO I BORDI
  		L0[n+1,1:M-1,1:N-1]=make_Laplacian(Zout[n+1,:,:])   
		
######## End of the time-stepping loop  #######

Znew=np.zeros([M,N])
Znew=Zout[nt,:,:] #Z at 24 hours

#Compute the RMSE for the computed forecast and for a persistence forecast defined as a stationary state:
differenza=np.zeros([M,N])
differenza_obs=np.zeros([M,N])

differenza=Znew-Z0  #DIFFERENZA 
differenza_obs=Z24-Z0  #DIFFERENZA

Rcomp=np.sqrt(np.sum((Znew[2:M-2,2:N-2]-Z24[2:M-2,2:N-2])**2)/((N-4)*(M-4)-1)) #RMSE
Rcomp_2=np.sqrt(np.sum((Z24[2:M-2,2:N-2]-Z0[2:M-2,2:N-2])**2)/((N-4)*(M-4)-1)) #RMSE

print ("RMSE Z(t)-Z24:",Rcomp) #STAMPO VALORI RMSE
print ("RMSE Z24-Z0:",Rcomp_2) #STAMPO VALORI RMSE


## Plot the analysis 
fig,ax=plt.subplots()
levels = np.linspace(480,600,9)
c=ax.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],Z24[2:M-2,2:N-2]/10,levels=levels,colors='green')
bb=ax.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],Znew[2:M-2,2:N-2]/10,levels=levels, colors='blue')
ax.clabel(c,fmt='%1.0f')
ax.scatter(X[Xp,Yp],Y[Xp,Yp])
ax.set_title("osservati vs calcolati")
ax.set_xlabel("x")
ax.set_ylabel("y")

fig,qqx=plt.subplots()
levels = np.linspace(480,600,9)
b=qqx.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],Z24[2:M-2,2:N-2]/10,levels=levels, colors="green")
bbb=qqx.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],Z0[2:M-2,2:N-2]/10,levels=levels, colors='violet')
qqx.clabel(b,fmt='%1.0f')
qqx.scatter(X[Xp,Yp],Y[Xp,Yp])
qqx.set_title("Z24 and Z0")
qqx.set_xlabel("x")
qqx.set_ylabel("y")


fig,dx=plt.subplots()
levels = np.linspace(480,600,9)
bb=dx.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],Znew[2:M-2,2:N-2]/10,levels=levels, colors='blue')
bbb=dx.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],Z0[2:M-2,2:N-2]/10,levels=levels, colors='violet')
dx.clabel(bb,fmt='%1.0f')
dx.scatter(X[Xp,Yp],Y[Xp,Yp])
dx.set_title("Znew and Z0")
dx.set_xlabel("x")
dx.set_ylabel("y")


#plot differenze tendency
fig = plt.figure(figsize=(8,8))
bx = fig.add_subplot(121)
levels = np.linspace(480,600,9)
e=bx.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],differenza_obs[2:M-2,2:N-2]/10,10,colors="k")
d=bx.contourf(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],differenza_obs[2:M-2,2:N-2]/10,10,cmap="winter")
bx.scatter(X[Xp,Yp],Y[Xp,Yp])
bx.set_title('Z24-Z0')
bx.set_xlabel("x")
bx.set_ylabel("y")
divider = make_axes_locatable(bx)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(d, cax=cax, orientation='vertical')
zx = fig.add_subplot(122)
levels = np.linspace(480,600,9)
ee=zx.contour(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],differenza[2:M-2,2:N-2]/10,10,colors="k")
kk=zx.contourf(X[2:M-2,2:N-2],Y[2:M-2,2:N-2],differenza[2:M-2,2:N-2]/10,10,cmap="winter")
zx.scatter(X[Xp,Yp],Y[Xp,Yp])
zx.set_title('Z(t)-Z0')
zx.set_xlabel("x")
zx.set_ylabel("y")
divider = make_axes_locatable(zx)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(kk, cax=cax, orientation='vertical')

plt.tight_layout()
plt.show()



