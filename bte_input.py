#!/usr/bin/python
import os

fw=open("input","w")
f = open(r'POSCAR','r')
f.readline()
f.readline()
#read the lattice constant from the CONTCAR file and calculate the base vectors in the reciprocal space
a1 = f.readline().split()
a2 = f.readline().split()
a3 = f.readline().split()
a1 = [float(i) for i in a1]
a2 = [float(i) for i in a2]
a3 = [float(i) for i in a3]

print >>fw,"       lfactor =  0.1,"
print >>fw,"  lattvec(:,1) =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(a1[0]),float(a1[1]),float(a1[2]))
print >>fw,"  lattvec(:,2) =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(a2[0]),float(a2[1]),float(a2[2]))
print >>fw,"  lattvec(:,3) =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(a3[0]),float(a3[1]),float(a3[2]))

try:
	line=f.readline().split()
	natom=[int(a) for a in line]
	elements=raw_input("please input the element name")
except:
	elements=[a for a in line]
        line=f.readline().split()
        natom=[int(a) for a in line]
total_atom=sum(natom)

for i in range(len(elements)):
	print>>fw,'      elements = "',elements[i],'"',
print>>fw,""
print>>fw,"         types = ",
for i in range(len(elements)):
	for j in range(natom[i]):
		print>>fw,i+1,
print>>fw,""

f.readline()
for i in range(total_atom):
	line=f.readline().split()
	pos=[float(a) for a in line]
	print>>fw,"positions(:,"+str(i+1)+") =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(pos[0]),float(pos[1]),float(pos[2]))

f.close()


result=os.popen('grep -5 DIELECTRIC OUTCAR').readlines()
result=[result[i] for i in range(len(result)-4,len(result)-1)]
for i in range(3):
	dielec=result[i].split()
	print>>fw,"  epsilon(:,"+str(i+1)+") =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(dielec[0]),float(dielec[1]),float(dielec[2]))

natom=2
nline=4*natom+2
result=os.popen('grep -'+str(nline)+' "BORN EFFECTIVE CHARGES (in e, cummulative output)" OUTCAR').readlines()
result=[result[i] for i in range(nline+2,2*nline)]

for i in range(natom):
	born=result[i*4+1].split()
	print>>fw,"   born(:,"+str(i+1)+",1) =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(born[1]),float(born[2]),float(born[3]))
        born=result[i*4+2].split()
        print>>fw,"   born(:,"+str(i+1)+",2) =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(born[1]),float(born[2]),float(born[3]))
        born=result[i*4+3].split()
        print>>fw,"   born(:,"+str(i+1)+",3) =","{0:15.10f} {1:15.10f} {2:15.10f}".format(float(born[1]),float(born[2]),float(born[3]))

print>>fw,"      scell(:) = "

fw.close()
	
	

