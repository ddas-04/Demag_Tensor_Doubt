###################  Fong Read write failure paper simulation ####################

import numpy as np
import scipy.special as sc

def elliptical_incomplete_integral_1st_kind_F(k, v):
	# numerical result is verified with https://keisan.casio.com/exec/system/1244989948
	return (sc.ellipkinc(v, k*k))  
	
def elliptical_incomplete_integral_2nd_kind_E(k, v):
	# numerical result is verified with https://keisan.casio.com/exec/system/1244989948
	return (sc.ellipeinc(v, k*k))
	
def Demag_Tensor(a,b,c):
	a=a/2.0
	b=b/2.0
	c=c/2.0
	
	cos_v=c/a
	v=np.arccos(cos_v)
	sin_v=np.sin(v)
	
	cos_psi=b/a
	psi=np.arccos(cos_psi)
	sin_psi=np.sin(psi)
	
	k=(np.sin(psi))/(np.sin(v))
	alpha=np.arcsin(k)
	sin_alpha=np.sin(alpha)
	cos_alpha=np.cos(alpha)
	
	F = elliptical_incomplete_integral_1st_kind_F(k, v)
	E = elliptical_incomplete_integral_2nd_kind_E(k, v)
	
	prefac_Nx = (cos_psi*cos_v)/(sin_v**3 * sin_alpha**2)
	Nxx = prefac_Nx*(F-E)

	prefac_Ny = (cos_psi*cos_v)/(sin_v**3 * sin_alpha**2 * cos_alpha**2)
	t3=(sin_alpha**2 * sin_v* cos_v)/cos_psi
	Nyy = prefac_Ny*(E-cos_alpha**2*F - t3)

	prefac_Nz= (cos_psi*cos_v)/(sin_v**3 * cos_alpha**2)
	t1 = (sin_v*cos_psi)/cos_v
	Nzz = prefac_Nz*(t1-E)
	
	return Nxx, Nyy, Nzz

######## Universal constant #################
mu_0=4*np.pi*1e-7           # free space permeability
k_B=1.380649*1e-23    # Boltzmann constant
T=300                              # Temperature in Kelvin


Major_axis=float(input('Enter major axis dimension(in nm) = ' ))
Minor_axis = float(input('Enter minor axis dimension(in nm) = ' ))
thickness = float(input('Enter thickness(in nm) = ' ))

Major_axis=Major_axis*1e-9
Minor_axis=Minor_axis*1e-9
thickness=thickness*1e-9

V_FL=(np.pi/4.0)*Major_axis*Minor_axis*thickness

[Nx, Ny, Nz] = Demag_Tensor(Major_axis, Minor_axis, thickness)

print('Nxx = ' + str(Nx))
print('Nyy = ' + str(Ny))
print('Nzz = ' + str(Nz))

if Major_axis==Minor_axis:   # PMA
	Ms=1080e3
	KB=float(input('Enter KB (in KJ/m3) = ' ))
	KB=KB*1e3
	Kdemag=0.5*mu_0*(Nz-Nx)*Ms**2
	Ku2=(KB-Kdemag)
	Delta = (Ku2*V_FL)/(k_B*T)
	print('Delta = ' + str(Delta))
else: # IMA
	Ms=1257.3e3
	Ku2=0.5*mu_0*(Ny-Nx)*Ms**2
	Delta=(Ku2*V_FL)/(k_B*T)
	print('Delta = ' + str(Delta))



