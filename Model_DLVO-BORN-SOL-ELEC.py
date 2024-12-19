import math
import numpy as np

pi      =   3.141592653
kb      =   1.38e-23
ee      =   1.602e-19
TT      =   298.15
aa      =   17.5e-9
zeta    =   10e-3
print(str(kb*TT))
eps_cb  =   2.78
eps_nmp =   32.2
hp      =   6.626e-34
v0      =   4e15
n_cb    =   2.0
n_nmp   =   1.47
eps0    =   8.854e-12
epsr    =   32.2

sigma   =   2e-9
nn      =   12
lmd     =   1

As      =   1.5
As      *=  1e-3
ld      =   1e-9

list_x  =   []
list_y  =   []
dx      =   1

ah1 = (0.75*kb*TT) * (((eps_cb-eps_nmp)/(eps_cb+eps_nmp))**2)
ah2 = (0.1325825*hp*v0) * ((((n_cb**2)-(n_nmp**2))**2)/(((n_cb**2)+(n_nmp**2))**1.5))
ah = ah1 + ah2

print('******A******')
print(ah)

print('')
print('*****kBT*****')
print('hh  vdw  born  sol  ele  total')
for i in range(5000):
    hh      =   (2e-12)*(i+1)
    list_x.append(hh)
    dx      =   2e-12
    rc      =   (hh+aa*2)/aa

    vdw1    =   (2*(aa**2))/(hh*((4*aa)+hh))
    vdw2    =   (2*(aa**2))/(((2*aa)+hh)**2)
    vdw3    =   math.log((hh*((4*aa)+hh))/(((2*aa)+hh)**2))
    v_vdw   =   -(ah/6)*(vdw1+vdw2+vdw3)

    born1   =   4*ah*((sigma/aa)**(nn-6))
    born2   =   (math.factorial(nn-6))/(math.factorial(nn-2))
    born3   =   lmd/(lmd+1)
    born4   =   1/((rc-1-lmd)**(nn-5))
    v_born  =   born1*born2*born3*born4

    v_sol   =   pi*As*ld*aa*math.exp(-hh/ld)

    v_ele   =   (4*pi*eps0*eps_nmp*(zeta**2)*(aa**2))/(hh+2*aa)

    hh_p    =   aa
    line    =   str(hh)+' '+str('%.3g' % v_vdw)+' '+str('%.3g' % v_born) +' '\
                +str('%.3g' % v_sol)+' '+str('%.3g' % v_ele)+' '\
                +str('%.3g' % (v_vdw+v_born+v_sol+v_ele)+' *')

    list_y.append(v_vdw+v_born+v_sol+v_ele)
    print(line)


print('**********')
delta1      =   np.diff(list_y, n=1)
list_y1     =   []

for i in delta1:
    j   =   i/dx
    list_y1.append(j)
    #print(j)

delta2      =   np.diff(list_y1, n=1)
list_y2     =   []

for i in delta2:
    j   =   i/dx
    list_y2.append(j)
    print(j)

print('---------------------------------------------------')
Uci         =   0
UU_min      =   1
Fci         =   0
FF_max      =   0
D2U         =   0

y0          =   0
iff         =   -1

for i in range(len(list_y2)):
    y0      =   list_y[i]
    yi      =   list_y[i+1]
    if yi   >   0:
        if yi < y0:
            iff = 1

    if iff == 1:
        if list_y[i + 1] < UU_min:
            Uci = i + 1
            UU_min = list_y[i + 1]
            D2U = list_y2[i]
        if list_y1[i + 1] > FF_max:
            Fci = i + 1
            FF_max = list_y1[i + 1]

    line    =   str(i)+' '+str(list_x[i+1])+' '+str(list_y[i+1])+' '+str(list_y1[i+1])\
                +' '+str(list_y2[i])
    print(line)
k0          =   D2U*list_x[Uci]/FF_max
#k0          =   D2U*(list_x[Uci]**2)/(-UU_min)
print('')
print('hc(e-9):    Uc(kBt):    hc1(e-9):    Fc(e-12):    k0: ')
print(str('%.3g' % (list_x[Uci]/1e-9))+str('        %.3g' % (-UU_min/kb/TT))
      +str('        %.3g' % (list_x[Fci]/1e-9))+str('         %.3g' % (FF_max/1e-12))
      +str('        %.3g' % k0))
