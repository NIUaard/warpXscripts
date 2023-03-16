import numpy as np
import matplotlib.pyplot as plt 
import plotly.graph_objects as go

from scipy import constants

epsilon0= constants.epsilon_0
mu0     = constants.mu_0

# fName='./slab_v3.3d'
fName='./slab_channel_v3.3d'

is3d = False

fh = open(fName,'r')

foundItMacro    = False 
foundItEmbedded = False 


for line in fh: 
#   print (line)
#   print ("---")
# load all contants
   if len(line.split('my_constants.'))>1:
      print (line.split('my_constants.')[1].split('='))     
      zz=line.split('my_constants.')[1].split('=')    
      exec(zz[0]+'='+(zz[1].split('\n')[0]))
      print((eval(zz[0])))
      print ('found and set parameters:'+zz[0]+'='+str(float(eval(zz[1]))))
# check computational domain size
   if len(line.split('geometry.prob_lo'))>1:
      index=0
      bblo = np.zeros((3))
      tmp= line.split('=')[-1].split(' ')
      for i in range(len(tmp)):
         if tmp[i]!='':
            bblo[index]=eval(tmp[i])
            index=index+1
      print ('>>>>>>>>>>>>> geometry.prob_lo:: ', bblo)
   if len(line.split('geometry.prob_hi'))>1:
      index=0
      bbhi = np.zeros((3))
      tmp= line.split('=')[-1].split(' ')
      for i in range(len(tmp)):
         if tmp[i]!='':
            bbhi[index]=eval(tmp[i])
            index=index+1
      print ('>>>>>>>>>>>>> geometry.prob_hi::', bbhi)

# build macroscopic function
   if foundItMacro:
      print ((line.split(' = ')))
      if len(line.split(' = '))>1:
         foundItMacro=False 
      if len(line.split(' = '))<2:
         print (line)
         macrofunc = macrofunc+line
   if line.split('(')[0]=='macroscopic.epsilon_function':
      foundItMacro=True
      print(line)
      macrofunc =  line.split(' = ')[1]
      print (macrofunc)     
      print("foundItMacro OK")

# build macroscopic function
   if foundItEmbedded:
      print ((line.split(' = ')))
      if len(line.split(' = '))>1:
         foundItEmbedded=False 
      if len(line.split(' = '))<2:
         print (line)
         embeddedfunc = embeddedfunc+line
   if line.split('(')[0]=='warpx.eb_implicit_function':
      foundItEmbedded=True
      print(line)
      embeddedfunc =  line.split(' = ')[1]
      print (embeddedfunc)      
      print("foundItEmbedded OK")
      
fh.close()

print ('final macrofunction is:')
print (macrofunc)
print ('final embeddedfunction is:')
print (embeddedfunc)

macrofunc=macrofunc.split('\n')[0].split('"')[1]
print(macrofunc)
#embeddedfunc=embeddedfunc.split('\n')[0].split('"')[1]
embeddedfunc=embeddedfunc.split('"')[1]

def epsilon(x,y,z):
    tmp=(eval(macrofunc))
    return(tmp)

def embedded(x,y,z):
    tmp=(eval(embeddedfunc))
    return(tmp)

print ('---------------------------')
print (epsilon(0.01,0+0.00001,0.1))
print (embedded(0.01,0+0.00001,0.1))
print ('---------------------------')

Nx=51
Ny=51
Nz=101

#reset bllo and bbhi along z because file use moving windows
bblo[2]= 0
bbhi[2]= 0.2280000000

print ("Computational domain:")
print ("bblo:",bblo)
print ("bbhi:",bbhi)
print ("---------------------")
X=np.linspace (bblo[0], bbhi[0], Nx)
Y=np.linspace (bblo[1], bbhi[1], Ny)
Z=np.linspace (bblo[2], bbhi[2], Nz)


# unwrapp array for faster plotting
xx,yy,zz = np.meshgrid (X,Y,Z)

x= xx.flatten() #reshape((1,Nx*Ny*Nz))
y= yy.flatten() #reshape((1,Nx*Ny*Nz))
z= zz.flatten() #reshape((1,Nx*Ny*Nz))


print(foundItMacro)
#if foundItMacro:
ee=eval(macrofunc)
eps=ee.reshape((Nx,Ny,Nz))/np.min(ee)
eps1d=ee/np.min(ee)
   
#if foundItEmbedded:
cc=eval(embeddedfunc)
cond=cc.reshape((Nx,Ny,Nz))
condnorm=0.
   
print ("min eps=", np.min(eps1d))
print ("max eps=", np.max(eps1d))


if is3d: 
   fig = go.Figure(data=go.Volume(
       x=xx.flatten(),
       y=yy.flatten(),
       z=zz.flatten(),
       value=eps1d.flatten(),
       isomin=0,
       isomax=np.max(eps1d),
       opacity=0.1, # needs to be small to see through all surfaces
       surface_count=11, # needs to be a large number for good volume rendering
       ))
   fig.show()
   
if not(is3d):
   plt.figure()
   plt.subplot (2,1,1)
   plt.imshow(eps[:,50,:], aspect='auto')
   plt.subplot (2,1,2)
   plt.imshow(eps[50,:,:], aspect='auto')
   
   plt.figure()
   plt.subplot (2,1,1)
   plt.imshow(cond[:,50,:], aspect='auto')
   plt.subplot (2,1,2)
   plt.imshow(cond[50,:,:], aspect='auto')
   plt.colorbar()
   plt.show()

