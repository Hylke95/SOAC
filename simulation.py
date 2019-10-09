import numpy as np
import matplotlib.pyplot as plt

def gauss(x,a,b,c):
	return a*np.e**(-(x-b)**2/(2*c**2))

def div(f):
    return np.ufunc.reduce(np.add,np.gradient(f,dL))

#simulation parameters
dL = 1000	#pixel size [m pix^-1]
dt = 365.25*24*3600/10	#time step [s]

#environmental parameters
L = 256	#grid size [pix]
coords = np.meshgrid(np.arange(L),np.arange(L))
b = np.zeros((L,L))	#bed elevation [m]
bstag = (b+np.roll(b,1,0)+np.roll(b,1,1)+np.roll(np.roll(b,1,0),1,))/4
h = np.zeros((L,L))	#ice thickness [m]

#surface mass balance [m yr^-1]
r = np.sqrt((coords[0]-L/2)**2+(coords[1]-L/2)**2)
adot = np.zeros((L,L))	#ice thickness [m]
dlt = r < 32
adot[dlt] = 0.3
adot /= 365.25*24*3600	#[m s^-1]

#constants
A = 1	#ice flow parameter [Pa^-3 yr^-1]
A /= 365.25*24*3600	#[Pa^-3 s^-1]
rho = 910	#ice density [kg m^-3]
g = 9.81	#gravitational acceleration [m s^-2]
n = 3	#Glen index
cnst = -2*(rho*g)**n*A/(n+2)	#combined constant [(kg m^-2 s^-2)^3]

plt.ion()

m = 0

plt.figure()
while True:
	#q = np.abs(div(elev))**3*cnst*h**5	#[kg^3 m^-1 s^-6] <- weird units
	#h += dt*(adot-np.abs(div(q)))
	hstag = (h+np.roll(h,1,0)+np.roll(h,1,1)+np.roll(np.roll(h,1,0),1,))/4
	elev = hstag+bstag
	maxes = np.array([np.max(np.abs(div(elev))),np.max(h)])
	print(maxes)
	if all(np.isfinite(maxes)):
		d = cnst*np.abs(div(elev))**(n-1)*hstag**(n+2)
		dstag = (d+np.roll(d,1,0)+np.roll(d,1,1)+np.roll(np.roll(d,1,0),1,))/4
		h += dt*(adot+dstag*div(elev))
		dlt = h < 0
		h[dlt] = 0
		
		m += 1
		
		plt.clf()
		plt.imshow(elev)
		plt.colorbar()
		plt.title('t = '+str(np.around(n*dt/(365.25*24*3600),decimals=1))+' yr')
		plt.pause(10**-6)
	else:
		break