# =====================================
# show Pi^z from 3 files in subplots
# =====================================

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
import os
import math
import sys
font = {'family' : 'Arial',
        'weight' : 'normal',
        'size'   : 20}
rc('font', **font)
matplotlib.use('tkagg')

mass = 1.115683

xFactor = 2.0   # !!!!!!!!!! extra factor, 1.0 to plot S, 2.0 to plot Pi

global_title = r'(standard + shear contributions), no T grad, $\epsilon_{fo}=0.5$ GeV/fm3'
figWidth = 5.33  # in inches
figHeight = 4.0 # in inches
file = sys.argv[1]
file_dim = file+'.dim'

dimP = 16
dimPhi = 40

rapShow = 0.0

dimRap, dimP, dimPhi = np.loadtxt(file_dim, dtype=int)
print('rapidity is', rapShow)
print(dimP)
#dndp = np.array([0.0] * dimP * dimPhi)
#Pizu = []
#for i in range(len(files)):
 #Pizu.append(np.array([0.0] * dim))
#print type(Pizu[0])

#Piz = []
#for i in range(len(files)):
 #Piz.append(np.array([0.0] * dim))

a = np.loadtxt(file, unpack=True)
inds, = np.where(a[0]==rapShow)
pT = a[1][inds]
phi = a[2][inds]
dndp = a[3][inds]
Pixu =  a[5][inds] #+ a[9]  # 4 or 8
Piyu =  a[6][inds] #+ a[10]  # 5 or 9
Pizu =  a[7][inds] #+ a[11]  # 6 or 10
pxE = pT * np.cos(phi)  # exact values (before the shift for better visual rep)
pyE = pT * np.sin(phi)  # exact values (before the shift for better visual rep)
dphi = phi[1] - phi[0]
phi = np.add(phi, -0.5*dphi)
# getting Cartesian components
px = pT * np.cos(phi)
py = pT * np.sin(phi)

# boost to Lambda rest frame
Pixurf = np.array([0.0] * len(Pixu))
Piyurf = np.array([0.0] * len(Pixu))
for i in range(len(Pixu)):
 e = math.sqrt(mass*mass + pT[i]*pT[i])
 Pixurf[i] = Pixu[i] - (Pixu[i]*pxE[i] + Piyu[i]*pyE[i])/(e*(e+mass))*pxE[i]  # 
 Piyurf[i] = Piyu[i] - (Pixu[i]*pxE[i] + Piyu[i]*pyE[i])/(e*(e+mass))*pyE[i]  # 

Pi = []
Pi.append(np.divide(Pixurf, dndp))
Pi.append(-np.divide(Piyurf, dndp))   # P^z --> P_J
Pi.append(np.divide(Pizu, dndp))

#z, = np.where(pT>2.9)
#print type(Pizu[z]), len(Pizu[z])
#print np.divide(Pizu[z], dndp[z])

pT = pT.reshape((dimP,dimPhi)).T
phi = phi.reshape((dimP,dimPhi)).T
px = px.reshape((dimP,dimPhi)).T
py = py.reshape((dimP,dimPhi)).T
#px = np.vstack((px,px[0]))
#py = np.vstack((py,py[0]))

Piyurf2d = Piyurf.reshape((dimP,dimPhi)).T
Pizurf2d = Pizu.reshape((dimP,dimPhi)).T
dndp2d = dndp.reshape((dimP,dimPhi)).T
Py_phi = np.array([0.0] * dimPhi)
Pz_phi = np.array([0.0] * dimPhi)
dn_phi = np.array([0.0] * dimPhi)
meanPi = 0.0
dndpInt = 0.0
for iphi in range(dimPhi):
 for ipt in range(dimP):
  if pT[iphi,ipt]>0.4 and pT[iphi,ipt]<2.0:          ### pT range!
   Py_phi[iphi] += Piyurf2d[iphi,ipt] * pT[iphi,ipt]
   Pz_phi[iphi] += Pizurf2d[iphi,ipt] * pT[iphi,ipt]
   dn_phi[iphi] += dndp2d[iphi,ipt] * pT[iphi,ipt]
   meanPi += Piyurf2d[iphi,ipt] * pT[iphi,ipt]
   dndpInt += dndp2d[iphi,ipt] * pT[iphi,ipt]
print('** average polarization = ', xFactor*meanPi/dndpInt)


#===========plotting:
# ================== P^z
fig = plt.figure(1, facecolor="white", figsize=(figWidth,figHeight))
plt.suptitle(global_title, color='blue')
#plt.title(r'$P^z(\phi)$')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$P^z$($\phi$)')
PzPhiPlot = xFactor*np.divide(Pz_phi,dn_phi)
PyPhiPlot = xFactor*np.divide(Py_phi,dn_phi)
plt.xlim(0., 1.6)
ticks=[0., 0.5*math.pi, math.pi]
tlabels=['0', r'$\pi/2$', r'$\pi$']
plt.xticks(ticks, tlabels)
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#plt.ylim(np.min(PzPhiPlot), 0.)
plt.plot(np.add(phi[:,0],0.5*dphi), PzPhiPlot)
plt.grid()

# ================== P_J
plt.figure(2, facecolor="white", figsize=(figWidth,figHeight))
#plt.title(r'$S_J(\phi)$')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$S_J$($\phi$)')
plt.xlim(0., 1.6)
plt.ylim(0., 1.1*np.max(-PyPhiPlot))
ticks=[0., 0.5*math.pi, math.pi]
tlabels=['0', r'$\pi/2$', r'$\pi$']
plt.xticks(ticks, tlabels)
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
#plt.ylim(np.min(PzPhiPlot), 0.)
plt.plot(np.add(phi[:,0],0.5*dphi), -PyPhiPlot)    # !!! plotting -S^y
plt.grid()

# =================== 2D plots
#title=[r'$P^x$',r'$P^J_\xi$',r'$P^z_\xi$']
#title=[r'$P^x$',r'$P^J_\varpi+P^J_\xi$',r'$P^z_\varpi+P^z_\xi$']
title=[r'$P^x$',r'$P^J_{\rm ILE}$',r'$P^z_{\rm ILE}$']
fig = plt.figure(3, facecolor="white", figsize=(1.7*figWidth, figHeight))
#plt.suptitle(global_title, color='blue')
for i in range(1,3):
 subplot_command = 121 + i - 1
 plt.subplot(subplot_command, aspect='equal')
 plt.title(title[i], pad=10)
 plt.xlabel(r'$p_x$ [GeV]')
 plt.ylabel(r'$p_y$ [GeV]')
 gridPi = (xFactor*Pi[i]).reshape((dimP,dimPhi)).T
 img=plt.pcolor(px, py, gridPi, cmap='seismic',
   vmax=np.max(gridPi), vmin=-np.max(gridPi))
 print(title[i], np.min(gridPi), np.max(gridPi))
 cb = plt.colorbar(img)
 cb.formatter.set_powerlimits((-2,2))
 cb.ax.yaxis.set_offset_position('left')
 cb.update_ticks()
 #plt.ticklabel_format(axis='z', style='sci', scilimits=(-2,2))
 plt.tight_layout()
 plt.subplots_adjust(top=0.86, bottom=0.2, left=0.09, right=0.98)
 #ax = plt.gca()
 #plt.text(0.75, 0.92, r'$\xi$', transform=ax.transAxes)
 #plt.text(0.75, 0.82 , r'no $\nabla T$', transform=ax.transAxes)

#plt.subplot(224)
#plt.title('dn/(dpx dpy), arb. units')
#plt.xlabel('px [GeV]')
#plt.ylabel('py [GeV]')
#gridDndp = dndp.reshape((size,size)).T
#img=plt.imshow(gridDndp, extent=(px.min(), px.max(), py.min(), py.max()),
           #interpolation='nearest', origin='lower', cmap='seismic')
#plt.colorbar(img)

#plt.suptitle(global_title)
plt.show()

