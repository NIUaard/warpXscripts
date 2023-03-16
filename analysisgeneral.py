import  matplotlib.pyplot as plt 
import numpy as np
from scipy.fft import fft, fftfreq
import os
import shutil
import scipy as sp
import sys
from matplotlib import gridspec
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import SymLogNorm
from matplotlib.colors import LogNorm

from openpmd_viewer import OpenPMDTimeSeries
from openpmd_viewer.addons import LpaDiagnostics

FontSize=25
FontSizeLabelAxis=25
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
filedi = arg[1]
filedir = filedi + '/'
Nitera  = int(arg[2])


ts = LpaDiagnostics(filedir)

N_iterations = len(ts.iterations)

if Nitera==-1:
   it = ts.iterations[N_iterations-1]
else:
   it = ts.iterations[Nitera]

plotspath = filedi+'Plots_ts'+str(it)

# figure names to be used throughout
fig9999 = 'Fig9999_Fullimagexz_ts'+str(it)+'.png'
fig9998 = 'Fig9998_Fullimageyz_ts'+str(it)+'.png'
fig9997 = 'Fig9997_freq_ts'+str(it)+'.png'
fig9996 = 'Fig9996_Ex_cBy_ts'+str(it)+'.png'
fig9995 = 'Fig9995_Ftrans_atz_ts'+str(it)+'.png'
fig9994 = 'Fig9994_Ftrans_aty0_ts'+str(it)+'.png'
fig9993 = 'Fig9993_Ftrans_atx0_ts'+str(it)+'.png'
fig9992 = 'Fig9992_Permittivity_ts'+str(it)+'.png'

fig9999p = 'Fig9999_Fullimagexz_ts'+str(it)+'.pdf'
fig9998p = 'Fig9998_Fullimageyz_ts'+str(it)+'.pdf'
fig9997p = 'Fig9997_freq_ts'+str(it)+'.pdf'
fig9996p = 'Fig9996_Ex_cBy_ts'+str(it)+'.pdf'
fig9995p = 'Fig9995_Ftrans_atz_ts'+str(it)+'.pdf'
fig9994p = 'Fig9994_Ftrans_aty0_ts'+str(it)+'.pdf'
fig9993p = 'Fig9993_Ftrans_atx0_ts'+str(it)+'.pdf'
fig9992p = 'Fig9992_Permittivity_ts'+str(it)+'.pdf'



if len(arg)==4:
   save_dir = arg[3]
else:
   if os.path.isdir(plotspath):
      if os.path.isfile(plotspath+'/'+fig9999):
         os.remove(plotspath+'/'+fig9999)
      if os.path.isfile(plotspath+'/'+fig9998):
         os.remove(plotspath+'/'+fig9998)
      if os.path.isfile(plotspath+'/'+fig9997):
         os.remove(plotspath+'/'+fig9997)
      if os.path.isfile(plotspath+'/'+fig9996):
         os.remove(plotspath+'/'+fig9996)
      if os.path.isfile(plotspath+'/'+fig9995):
         os.remove(plotspath+'/'+fig9995)
      if os.path.isfile(plotspath+'/'+fig9994):
         os.remove(plotspath+'/'+fig9994)
      if os.path.isfile(plotspath+'/'+fig9993):
         os.remove(plotspath+'/'+fig9993)
      if os.path.isfile(plotspath+'/'+fig9992):
         os.remove(plotspath+'/'+fig9992)
      save_dir = plotspath+'/'
   else:
      os.mkdir(plotspath)
      save_dir = plotspath+'/'

# retrieve the data from the hdf5 file for the electric and magnetic fields
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

print('x, y, z grid size: ',np.flip(shapeEz))

absmax=MV*max(np.max(np.max(np.abs(Ez))),np.abs(np.min(np.min(np.abs(Ez)))))


# ----------Start to initialize figures and plotting -----------

# generate a plot of the Ez wake in (z,y) with lineout of Ez
fig=plt.figure(9999)
gs1 = gridspec.GridSpec(11, 13)
ax1 = fig.add_subplot(gs1[0:7,1:10])
ax2 = fig.add_subplot(gs1[7:10,1:10], sharex =ax1)
ax3 = fig.add_subplot(gs1[0:7,11:12])


# record Ez in the (z,x) plane for y=0
Ezslice=Ez[:,int(shapeEz[1]/2),:].transpose()*MV # gets the y=0 slice
shapeEzslice = np.shape(Ezslice)

# plot slice (upper plot with red and blue lines)
imField=ax1.imshow(Ezslice,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.xmin*mm,info_Ez.xmax*mm], 
                   cmap='RdBu_r', vmin=-absmax, vmax=absmax, aspect='auto')
cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2) 
cbar1  = fig.colorbar(imField, cax=cbaxes, orientation='horizontal').set_label(label=r'$E_z$ (MV/m)',size=20)
ax1.set_ylabel(r'$x$ (mm)', fontsize=FontSizeLabelAxis)

# record particle and plot a histogram on top of 3d of Ez - NEED TO ADJUST AXES FOR THIS
zf, xf = ts.get_particle( ['z','x'],  species='myparticle', iteration=it)
particles=ax1.hexbin(zf*mm, xf*mm, gridsize=200,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.xmin*mm,info_Ez.xmax*mm], cmap='copper', alpha=0.5,  mincnt=5)
ax3.axis('off')
cbaxes2 = inset_axes(ax3, width="30%", height="75%", loc=2) 
cbar2 = fig.colorbar(particles, cax=cbaxes2).set_label(label=r'$n$ ',size=20)

# Line plot of Ez plus histogram of the population on lower plot
zsc=np.linspace(info_Ez.zmin*1e3, info_Ez.zmax*1e3, shapeEz[0])
ax2.plot (zsc,Ezslice[int(shapeEzslice[0]/2),:]) # gets the y=0 slice
ax2.set_xlabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis) 
ax2.set_ylabel(r'$E_z$ (MV/m) ', fontsize=FontSizeLabelAxis, color="C0")
lenz = info_Ez.zmax*mm-info_Ez.zmin*mm
xticklist = np.arange(info_Ez.zmin*mm,info_Ez.zmax*mm+lenz/2,lenz/3)
#print(xticklist)
ax2.set_xticks(xticklist)
ax22 = ax2.twinx()
ax22.hist(zf*mm,512, alpha=0.3, color='C1')
ax22.set_ylabel(r'population', fontsize=FontSizeLabelAxis, color="C1")

# Set up the plots
plt.setp(ax1.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax2.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False) 
# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
#plt.subplots_adjust(hspace=.0)

# Save the plot to be viewed 
plt.savefig(save_dir+fig9999)
plt.savefig(save_dir+fig9999p,dpi=600)

# record Ez in the (z,y) plane for x=0
fig=plt.figure(9998)
gs1 = gridspec.GridSpec(11, 13)
ax1 = fig.add_subplot(gs1[0:7,1:10])
ax2 = fig.add_subplot(gs1[7:10,1:10], sharex =ax1)
ax3 = fig.add_subplot(gs1[0:7,11:12])

# record Ez in the (z,y) plane for x=0
Ezslice=Ez[:,:,int(shapeEz[2]/2)].transpose()*MV # gets the x=0 slice
shapeEzslice = np.shape(Ezslice)

# plot slice (upper plot with red and blue lines)
imField=ax1.imshow(Ezslice,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.ymin*mm,info_Ez.ymax*mm], 
                   cmap='RdBu_r', vmin=-absmax, vmax=absmax, aspect='auto')
cbaxes = inset_axes(ax1, width="30%", height="3%", loc=2) 
cbar1  = fig.colorbar(imField, cax=cbaxes, orientation='horizontal').set_label(label=r'$E_z$ (MV/m)',size=20)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)

# record particle and plot a histogram on top of 3d of Ez
zf, yf = ts.get_particle( ['z','y'],  species='myparticle', iteration=it)
particles=ax1.hexbin(zf*mm, yf*mm, gridsize=200,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.ymin*mm,info_Ez.ymax*mm], cmap='copper', alpha=0.5,  mincnt=5)
ax3.axis('off')
cbaxes2 = inset_axes(ax3, width="30%", height="75%", loc=2) 
cbar2 = fig.colorbar(particles, cax=cbaxes2).set_label(label=r'$n$ ',size=20)

# Line plot of Ez plus histogram of the population on lower plot
zsc=np.linspace(info_Ez.zmin*mm, info_Ez.zmax*mm, shapeEz[0])
ax2.plot (zsc,Ezslice[int(shapeEzslice[0]/2),:]) # Ez slice has shape (ylen,zlen)
#ax2.plot (zsc,Ey[16,48,:]*1e-6)	
#ax2.plot (zsc,Ex[16,48,:]*1e-6)	
ax2.set_xlabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis) 
ax2.set_ylabel(r'$E_z$ (MV/m) ', fontsize=FontSizeLabelAxis, color="C0")
lenz = info_Ez.zmax*mm-info_Ez.zmin*mm
xticklist = np.arange(info_Ez.zmin*mm,info_Ez.zmax*mm+lenz/2,lenz/3)
ax2.set_xticks(xticklist)
ax22 = ax2.twinx()
ax22.hist(zf*mm,512, alpha=0.3, color='C1')
ax22.set_ylabel(r'population', fontsize=FontSizeLabelAxis, color="C1")

# Set up the plots
plt.setp(ax1.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax2.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False) 
# remove vertical gap between subplots
plt.subplots_adjust(hspace=.0)
#plt.subplots_adjust(hspace=.0)

# Save the plot to be viewed 
plt.savefig(save_dir+fig9998)
plt.savefig(save_dir+fig9998p,dpi=600)

# Doing FFT to get the wavelength of the accelerated bunch for (y,z) slice
Ezslice =Ez[:,:,int(shapeEz[2]/2)].transpose()*MV # takes the z slice
zint = zsc/cms # convert z axis to a time axis
Ezslicelen = len(Ezslice[int(shapeEz[2]/2),:]) # get length of Ezslice for N

samplerate = Ezslicelen/((info_Ez.zmax*1e3-info_Ez.zmin*1e3)/cms) # N/(t_f - t_i) = N/((zmax(mm) - zmin(mm))/cms) so that end wavelength is mm

Fy = fft(Ezslice[int(shapeEz[1]/2)]) # do fft
Fx = fftfreq(Ezslicelen,1 / samplerate) # get the fft frequencies

# get the maximum of Fy
freqmax = np.amax(np.abs(Fy))
# which frequency does the max of Fy correspond to?
indexfreq = np.where(np.abs(Fy) == freqmax)
flatind = indexfreq[0]
print('max freq ',freqmax, ' located at ',Fx[flatind[0]])
# convert the frequency to wavelength using f = c/lambda
wavelength = cms/np.abs(Fx[flatind[0]])
print('wavelength = ', wavelength, ' mm')

# Plotting to see the frequencies to check with wavelength on plot
figfreq = plt.figure(9997)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figfreq.add_subplot(gs1[1:10,1:10])
ax1.text(Fx[flatind[0]]*2+Fx[flatind[1]]/2,freqmax*3/4,'lambda = '+str(round(wavelength,8))+'mm',bbox=dict(facecolor='red',alpha = 0.5))
imFreqs = ax1.plot(Fx,np.abs(Fy))
plt.xlim(0,(np.amax(Fx)/Fx[flatind[1]])*(3/5)*Fx[flatind[1]])
ax1.set_xlabel("Freq (Hz)")
plt.savefig(save_dir+fig9997)
plt.savefig(save_dir+fig9997p,dpi=600)


figX=plt.figure(9996)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figX.add_subplot(gs1[1:10,1:10])
Exslice=(Ex[:,:,int(shapeEz[2]/2)].transpose() \
               -cms*By[:,:,int(shapeEz[2]/2)].transpose())*kV
imField=ax1.imshow(Exslice,extent=[info_Ez.zmin*mm,info_Ez.zmax*mm, info_Ez.ymin*mm,info_Ez.ymax*mm], 
                   cmap='RdBu_r', aspect='auto')
cbar1  = figX.colorbar(imField).set_label(label=r'$E_x-cB_y$ (kV/m)',size=16)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis) 

plt.savefig(save_dir+fig9996)
plt.savefig(save_dir+fig9996p,dpi=600)


# transverse slice at a given z
figSliceZ=plt.figure(9995)  
gs1 = gridspec.GridSpec(11, 11)
ax1 = figSliceZ.add_subplot(gs1[1:10,1:10])

#print ('xySlice::',np.shape(Bx))
Bxshape = np.shape(Bx)
sliceInd=int(Bxshape[0]/2) # HOW IS THIS THOUGHT OF: NEEDS TO CHANGE TO GET THE MODES
Bxslice=Bx[sliceInd,:,:]
Exslice=Ex[sliceInd,:,:]
Eyslice=Ey[sliceInd,:,:]
Byslice=By[sliceInd,:,:]

# F/e computed for an electron
Fx=-(Exslice-cms*Byslice)
Fy=-(Eyslice+cms*Bxslice)

imField=ax1.imshow(np.sqrt(Fx**2+Fy**2)*kV, cmap='hot', norm=LogNorm(vmin=0.01,vmax=1e4),
           extent=[info_Ex.xmin*mm,info_Ex.xmax*mm, info_Ex.ymin*mm,info_Ex.ymax*mm],
           aspect='auto') 
figSliceZ.colorbar(imField).set_label(label=r'$F_{\perp}$ (kV/m)',size=16)
ax1.set_ylabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$x$ (mm)', fontsize=FontSizeLabelAxis) 
ax1.set_title(r'$z=%f$'%zsc[sliceInd], fontsize=FontSizeLabelAxis)
plt.savefig(save_dir+fig9995)


# transverse slice at a given y in (x,z)
figSliceX=plt.figure(9994)  # transverse slice at a given y in (x,z)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figSliceX.add_subplot(gs1[1:10,1:10])

Bxshape = np.shape(Bx)
sliceInd=int(Bxshape[1]/2) # currently for y=0
Bxslice = Bx[:,sliceInd,:]
Byslice = By[:,sliceInd,:]
Exslice = Ex[:,sliceInd,:]
Eyslice = Ey[:,sliceInd,:]

# F/e computed for an electron
Fx=-(Exslice-cms*Byslice)
Fy=-(Eyslice+cms*Bxslice)
	       
imField=ax1.imshow(np.sqrt(Fx**2+Fy**2)*kV, cmap='hot',norm=LogNorm(vmin=0.01,vmax=1e4), extent=[info_Ex.xmin*mm,info_Ex.xmax*mm,info_Ex.zmin*mm,info_Ex.zmax*mm],
           aspect='auto', origin='lower') 
figSliceX.colorbar(imField).set_label(label=r'$F_{\perp}$ (kV/m)',size=16)
ax1.set_ylabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$x$ (mm)', fontsize=FontSizeLabelAxis) 
ax1.set_title (r'$y=0$', fontsize=FontSizeLabelAxis) 
plt.savefig(save_dir+fig9994)
plt.savefig(save_dir+fig9994p,dpi=600)

# transverse slice at a given x in (y,z)
figSliceY=plt.figure(9993)  # transverse slice at a given x in (y,z)
gs1 = gridspec.GridSpec(11, 11)
ax1 = figSliceY.add_subplot(gs1[1:10,1:10])

Bxshape = np.shape(Bx)
sliceInd=int(Bxshape[2]/2) # currently for x=0
Bxslice = Bx[:,:,sliceInd]
Byslice = By[:,:,sliceInd]
Exslice = Ex[:,:,sliceInd]
Eyslice = Ey[:,:,sliceInd]

# F/e computed for an electron
Fx=-(Exslice-cms*Byslice)
Fy=-(Eyslice+cms*Bxslice)
	       
imField=ax1.imshow(np.sqrt(Fx**2+Fy**2)*kV, cmap='hot', norm=LogNorm(vmin=0.01,vmax=1e4), extent=[info_Ex.ymin*mm,info_Ex.ymax*mm,info_Ex.zmin*mm,info_Ex.zmax*mm],
           aspect='auto', origin='lower') 
figSliceY.colorbar(imField).set_label(label=r'$F_{\perp}$ (kV/m)',size=16)
ax1.set_ylabel(r'$z$ (mm)', fontsize=FontSizeLabelAxis)
ax1.set_xlabel(r'$y$ (mm)', fontsize=FontSizeLabelAxis) 
ax1.set_title (r'$x=0$', fontsize=FontSizeLabelAxis) 
plt.savefig(save_dir+fig9993)
plt.savefig(save_dir+fig9993p,dpi=600)









