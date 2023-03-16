import  matplotlib.pyplot as plt 
import numpy as np
import scipy as sp
import sys
from matplotlib import gridspec
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import SymLogNorm
from matplotlib.colors import LogNorm

from openpmd_viewer import OpenPMDTimeSeries
from openpmd_viewer.addons import LpaDiagnostics

FontSize=18
FontSizeLabelAxis=18
rc('legend',fontsize=FontSize)
rc('xtick',labelsize=FontSize)
rc('ytick',labelsize=FontSize)
font = {'family' : 'normal',
        'size'   : FontSize}

rc('font', **font)
rc('text', usetex=True)

# scaling factors
mm=1e3
MV=1e-6 
kV=1e-3

cms = sp.constants.speed_of_light

arg = sys.argv 
argnum = len(sys.argv)
filedir = arg[1]
Nitera  = int(arg[2])

# generate a plot of the Ez wake in (z,x) with lineout of Ez
fig=plt.figure(9999)
gs1 = gridspec.GridSpec(11, 13)
ax1 = fig.add_subplot(gs1[0:7,1:10])
ax2 = fig.add_subplot(gs1[7:10,1:10], sharex =ax1)
ax3 = fig.add_subplot(gs1[0:7,11:12])

ts = LpaDiagnostics(filedir)

N_iterations = len(ts.iterations)

if Nitera==-1:
   it = ts.iterations[N_iterations-1]
else:
   it= ts.iterations[Nitera]
   
Ez, info_Ez = ts.get_field( iteration=it,  field='E', coord='z', 
                            slicing_dir='x', plot=False)
Ex, info_Ex = ts.get_field( iteration=it,  field='E', coord='x', 
                            slicing_dir='x', plot=False)
Ey, info_Ey = ts.get_field( iteration=it,  field='E', coord='y', 
                            slicing_dir='x', plot=False)
Bx, info_Bx = ts.get_field( iteration=it,  field='B', coord='x', 
                            slicing_dir='x', plot=False)
By, info_By = ts.get_field( iteration=it,  field='B', coord='y', 
                            slicing_dir='x', plot=False)

shapeEz= np.shape(Ez)

absmax=MV*max(np.max(np.max(np.abs(Ez))),np.abs(np.min(np.min(np.abs(Ez)))))

print ('------------- INFO_EZ')
print(info_Ez.__dict__)
print ('axes.....',info_Ez.axes)
print ('shape....', shapeEz) 

# record Ez in the (z,y) plane for x=0
Ezslice=Ez[:,:,int(shapeEz[2]/2)].transpose()*MV

# plot slice
imField=ax1.imshow(Ezslice,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.ymin*mm,info_Ez.ymax*mm], 
                   cmap='RdBu_r', vmin=-absmax, vmax=absmax, aspect='auto')
cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2) 
cbar1  = fig.colorbar(imField, cax=cbaxes, orientation='horizontal').set_label(label=r'$E_z$ (MV/m)',size=16)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)

# record particle and plot a histogram
zf, xf = ts.get_particle( ['z','x'],  species='myparticle', iteration=it)
particles=ax1.hexbin(zf*mm, xf*mm, gridsize=50, cmap='copper', alpha=0.5,  mincnt=5)
ax3.axis('off')
cbaxes2 = inset_axes(ax3, width="20%", height="75%", loc=2) 
cbar2 = fig.colorbar(particles, cax=cbaxes2).set_label(label=r'$n$ ',size=16)

zsc=np.linspace(info_Ez.zmin*1e3, info_Ez.zmax*1e3, shapeEz[0])
ax2.plot (zsc,Ezslice[int(shapeEz[1]/2),:])	
#ax2.plot (zsc,Ey[16,48,:]*1e-6)	
#ax2.plot (zsc,Ex[16,48,:]*1e-6)	
ax2.set_xlabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis) 
ax2.set_ylabel(r'$E_z$ (MV/m) ', fontsize=FontSizeLabelAxis, color="C0")
ax22 = ax2.twinx()
ax22.hist(zf*mm,512, alpha=0.3, color='C1')
ax22.set_ylabel(r'population', fontsize=FontSizeLabelAxis, color="C1")

plt.setp(ax1.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax2.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False) 
# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
#plt.subplots_adjust(hspace=.0)

figX=plt.figure(9998)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figX.add_subplot(gs1[1:10,1:10])
Exslice=(Ex[:,:,int(shapeEz[2]/2)].transpose() \
               -cms*By[:,:,int(shapeEz[2]/2)].transpose())*kV
imField=ax1.imshow(Exslice,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.ymin*mm,info_Ez.ymax*mm], 
                   cmap='RdBu_r', aspect='auto')
cbar1  = figX.colorbar(imField).set_label(label=r'$E_x-cB_y$ (kV/m)',size=16)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis) 

figX=plt.figure(9998)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figX.add_subplot(gs1[1:10,1:10])
Exslice=(Ex[:,:,int(shapeEz[2]/2)].transpose() \
               -cms*By[:,:,int(shapeEz[2]/2)].transpose())*kV
imField=ax1.imshow(Exslice,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.ymin*mm,info_Ez.ymax*mm], 
                   cmap='RdBu_r', aspect='auto')
cbar1  = figX.colorbar(imField).set_label(label=r'$E_x-cB_y$ (kV/m)',size=16)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis) 


figSliceZ=plt.figure(9997)  # transverse slice at a given z
gs1 = gridspec.GridSpec(11, 11)
ax1 = figSliceZ.add_subplot(gs1[1:10,1:10])

Ex, info_Ex = ts.get_field( iteration=it,  field='E', coord='x', 
                            slicing_dir='z', plot=False)
Ey, info_Ey = ts.get_field( iteration=it,  field='E', coord='y', 
                            slicing_dir='z', plot=False)
Bx, info_Bx = ts.get_field( iteration=it,  field='B', coord='x', 
                            slicing_dir='z', plot=False)
By, info_By = ts.get_field( iteration=it,  field='B', coord='y', 
                            slicing_dir='z', plot=False)
			    
print ('xySlice::',np.shape(Bx))
slideInd=78
Bx=Bx[slideInd,:,:]
Ex=Ex[slideInd,:,:]
Ey=Ey[slideInd,:,:]
By=By[slideInd,:,:]

# F/e computed for an electron
Fx=-(Ex-cms*By)
Fy=-(Ey+cms*Bx)
	       
imField=ax1.imshow(np.sqrt(Fx**2+Fy**2)*kV, cmap='gray', norm=LogNorm(vmin=0.01,vmax=1e4),
           extent=[info_Ex.xmin*mm,info_Ex.xmax*mm, info_Ex.ymin*mm,info_Ex.ymax*mm],
           aspect='auto') 
figSliceZ.colorbar(imField).set_label(label=r'$F_{\perp}$ (kV/m)',size=16)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$x$ (mm)', fontsize=FontSizeLabelAxis) 


figSliceY=plt.figure(9996)  # transverse slice at a given x in (y,z)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figSliceY.add_subplot(gs1[1:10,1:10])
Ex, info_Ex = ts.get_field( iteration=it,  field='E', coord='x', 
                            slicing_dir='y', plot=False)
Ey, info_Ey = ts.get_field( iteration=it,  field='E', coord='y', 
                            slicing_dir='y', plot=False)
Bx, info_Bx = ts.get_field( iteration=it,  field='B', coord='x', 
                            slicing_dir='y', plot=False)
By, info_By = ts.get_field( iteration=it,  field='B', coord='y', 
                            slicing_dir='y', plot=False)
			    
print ('xySlice::',np.shape(Bx))
slideInd=32
Bx=Bx[:,slideInd,:]
Ex=Ex[:,slideInd,:]
Ey=Ey[:,slideInd,:]
By=By[:,slideInd,:]

# F/e computed for an electron
Fx=-(Ex-cms*By)
Fy=-(Ey+cms*Bx)
	       
imField=ax1.imshow(np.sqrt(Fx**2+Fy**2)*kV, cmap='gray', norm=LogNorm(vmin=0.01,vmax=1e4),
           aspect='auto', origin='lower') 
figSliceY.colorbar(imField).set_label(label=r'$F_{\perp}$ (kV/m)',size=16)
ax1.set_ylabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$x$ (mm)', fontsize=FontSizeLabelAxis) 
ax1.set_title (r'$y=0$', fontsize=FontSizeLabelAxis) 

plt.show()

