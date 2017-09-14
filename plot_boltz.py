#!/usr/bin/python

import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches

case=os.popen('ls *.struct').read().split('.')[0]
filename=case+'.condtens'
f=open(filename)

nline=os.popen('wc -l '+filename).read().split()[0]
nline=int(nline)
nedos=nline-1
f.readline()
data=np.fromfile(f,sep=' ',count=nedos*30,dtype=float).reshape(nedos,30)


energy=data[:,0]
S_xx=data[:,12]
S_yy=data[:,16]
temperature=data[:,1]

itemp=0


color_list=['red','green','cyan','orange','blue','violet','black','magenta']
temperatures=np.linspace(50,750,8)
print temperatures
fig, ax = plot.subplots()
output='Seebeck.png'
color_patch=[]
for index,temp in enumerate(temperatures):
    itemp=index*2
    #ax.plot(energy[itemp::16], S_xx[itemp::16], color=color_list[index], ls='-',lw=2)
    ax.plot(energy[itemp::16], S_yy[itemp::16], color=color_list[index], ls='-',lw=2)
    color_patch.append(mpatches.Patch(color=color_list[index], label='T='+str(temp)+'K'))

plot.legend(bbox_to_anchor=(0.7,0.9),handles=[item for item in color_patch],loc=0,prop={'size':8})

axes=plot.gca()
ymin, ymax = axes.get_ylim()
plot.ylim(ymin,ymax)
ax.set_xlabel('$E-E_f$ (Ry)')
ax.set_ylabel('$Seebeck$ coefficient ($\mu V/K$)')
plot.savefig(output, dpi=300, bbox_inches='tight')

