
# coding: utf-8

# In[141]:

#Importing libraries
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


# In[11]:

#Importing the parameters of the box
Rbox=20.0
Mesh=0.1
Npoints=int(Rbox/Mesh)
r = np.arange(Mesh,Rbox+Mesh,Mesh)

#Plots of the potentials
PotentialN = np.loadtxt('neutron_potential.dat')
PotentialP = np.loadtxt('proton_potential.dat')
fig,ax = plt.subplots(3)
ax[0].plot(r,PotentialN[:,1])
ax[0].set(xlabel='r[fm]',ylabel='V$_{WS}$(r)')
ax[1].plot(r,PotentialN[:,2])
ax[1].set(xlabel='r[fm]',ylabel='V$_{so}$(r)')
ax[2].plot(r,PotentialN[:,1]+PotentialN[:,2])
ax[2].set(xlabel='r[fm]',ylabel='V$_{tot}$(r)')
plt.subplots_adjust()
plt.savefig('PotentialN.pdf')
fig,ax = plt.subplots(4)
ax[0].plot(r,PotentialP[:,1])
ax[0].set(xlabel='r[fm]',ylabel='V$_{WS}$(r)')
ax[1].plot(r,PotentialP[:,2])
ax[1].set(xlabel='r[fm]',ylabel='V$_{so}$(r)')
ax[2].plot(r,PotentialP[:,3])
ax[2].set(xlabel='r[fm]',ylabel='V$_{coul}$(r)')
ax[3].plot(r,PotentialP[:,1]+PotentialP[:,2]+PotentialP[:,3])
ax[3].set(xlabel='r[fm]',ylabel='V$_{tot}$(r)')
plt.subplots_adjust()
plt.savefig('PotentialP.pdf')

#Plots of the densities
DensityN = np.loadtxt('DensityN.dat')
DensityP = np.loadtxt('DensityP.dat')
fig,ax = plt.subplots(3)
ax[0].plot(r,DensityN[:,1])
ax[0].set(xlabel='r[fm]',ylabel='$\\rho [fm^{-3}]$')
ax[1].plot(r,DensityN[:,2])
ax[1].set(xlabel='r[fm]',ylabel='$ \\tau [fm^{-5}]$')
ax[2].plot(r,DensityN[:,3])
ax[2].set(xlabel='r[fm]',ylabel='J $[fm^{-4}]$')
plt.subplots_adjust()
plt.savefig('DensityN.pdf')
fig,ax = plt.subplots(3)
ax[0].plot(r,DensityP[:,1])
ax[0].set(xlabel='r[fm]',ylabel='$\\rho [fm^{-3}]$')
ax[1].plot(r,DensityP[:,2])
ax[1].set(xlabel='r[fm]',ylabel='$ \\tau [fm^{-5}]$')
ax[2].plot(r,DensityP[:,3])
ax[2].set(xlabel='r[fm]',ylabel='J $[fm^{-4}]$')
plt.subplots_adjust()
plt.savefig('DensityP.pdf')


#Plots of the eigenvalues
EigN = np.loadtxt('neutron_singleparticles.dat')
EigP = np.loadtxt('proton_singleparticles.dat')
fig,ax = plt.subplots(1,2)
ax[0].hlines(EigN[:,0],0,1)
plt.setp(ax[0].get_xticklabels(), visible=False)
ax[1].hlines(EigP[:,0],0,1)
plt.setp(ax[1].get_xticklabels(), visible=False)
plt.savefig('eigen.pdf')


#Plots of the wave-functions
Wfn = np.loadtxt('neutron_wfs.dat')
Wfp = np.loadtxt('protons_wfs.dat')
Wfn2=Wfn[:,1].reshape(int(len(Wfn[:,1])/Npoints),Npoints)
Wfp2=Wfp[:,1].reshape(int(len(Wfp[:,1])/Npoints),Npoints)
fig,ax = plt.subplots(6)
for i in range(6):
    ax[i].plot(r,Wfn2[i,:])
    ax[i].set(ylabel='n=$i')
plt.savefig('wfneutron.pdf')
fig,ax = plt.subplots(10)
for i in range(6):
    ax[i].plot(r,Wfp2[i,:])
plt.savefig('wfproton.pdf')
