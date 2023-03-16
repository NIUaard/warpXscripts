import  matplotlib.pyplot as plt 
import numpy as np
import sys
from matplotlib import gridspec
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

plt.show()

