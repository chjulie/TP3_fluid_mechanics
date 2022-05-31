#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 12:41:27 2022

@author: Julie
"""
import matplotlib.pyplot as plt 
import numpy as np
import scipy
from scipy import optimize


data_content1 = np.genfromtxt('Group_1_1.txt', delimiter = ',', skip_header=1, dtype ='float')
data_content2 = np.genfromtxt('Group_1_2.txt', delimiter = ',', skip_header=1, dtype ='float')

zList1 = []
uList1 = []
siguList1= []
zList2 = []
uList2 = []
siguList2 = []

for i in range(0, len(data_content2)):
    zList2.append(data_content2[i][0])
    uList2.append(data_content2[i][1])
    siguList2.append(data_content2[i][2])
    
for i in range(0, len(data_content1)):
    zList1.append(data_content1[i][0])
    uList1.append(data_content1[i][1])
    siguList1.append(data_content1[i][2])
    

'''
1. Mean streamwise velocity for the two boundary layers as a function of the z-coordinate
and report de free stream velocity
'''
fig, ax = plt.subplots()
ax.grid(True)
ax.set_xlim(0,0.5) 

ax.plot(zList1, uList1, label='Dataset 1', color='dimgrey')
ax.plot(zList1, uList2, label='Dataset 2', color='silver')

plt.text(zList1[0], uList1[0], ' 5.0273 m/s', va='bottom')
plt.text(zList2[0], uList2[0], ' 5.0379 m/s', va='top')

plt.xlabel('Depth z [m]')
plt.ylabel('Mean stream velocity U [m/s]')
#plt.title('Mean streamwise velocity')
plt.legend(loc='best')
plt.savefig("plots/mean_streamwise_velocity.png", bbox_inches='tight')

#Finding the free stream velocity (mean of 2 highest value from both datasets)
UFS1 = uList1[0]
UFS2 = uList2[0]


'''
2. Height of the boundary layer in each profile
'''
def interpolation(point1, point2):
    deltax = float(point2[0]) - float(point1[0])
    deltay = float(point2[1]) -float(point1[1])

    m = deltay / deltax
    
    h = point1[1] - (m * point1[0])
    
    return m, h



def boundaryH(UFS, uList, zList):
    
    ubound = 0.99 * UFS

    
    for i, vel in enumerate(uList):
        if vel < ubound:
            lowerIndex = i
            break
    
    p1 = [zList[lowerIndex-1], uList[lowerIndex-1]]
    p2 = [zList[lowerIndex], uList[lowerIndex]]
      
    m, h = interpolation(p2, p1)
    
    BLH = (UFS - h) / (m)
            
    return ubound, BLH


hbl1 = boundaryH(UFS1, uList1, zList1)
hbl2 = boundaryH(UFS2, uList2, zList2)


'''
3. Fit a power law to each velocity profile: U(z) = Uref * (z / zref) ** alpha
Solve to find alpha
'''

def powerLaw(x, alpha):
    #return 5.0273 * np.power((x / 0.5), alpha)
    return 5.03 * np.power((x / 0.5), alpha)

def powerLawforlists(xList, alpha):
    y = []
    '''
    for value in xList:
        y.append(5.0273 * np.power((value / 0.5), alpha))
    '''
    for value in xList:
        y.append(5.03 * np.power((value / 0.5), alpha))
    
    return y


def powerLawFitting(Uref, Zref, zList, uList, siguList, title):

    pars, cov = scipy.optimize.curve_fit(f=powerLaw, xdata=zList, ydata=uList, sigma=siguList, p0 = [0])
    
    alpha = pars[0]
    
    fig2, ax2 = plt.subplots()
    ax2.grid(True)
    ax2.set_ylim(2,6)
    ax2.set_xlim(0,0.5)
    
    ax2.scatter(zList, uList, color='dimgrey',label='Data')
    ax2.plot(zList, powerLawforlists(zList, alpha), color='black', label='Fitted line')
    
    plt.xlabel('Velocity u [m/s]')
    plt.ylabel('Depth z [m]')
    #plt.title('Power Law fitting for dataset ' + title)
    plt.legend(loc='best')
    plt.savefig('plots/powerLawFitting'+title+'.png', bbox_inches='tight')
    
    return alpha
    
print('alpha 1= ', powerLawFitting(UFS1, 0.5, zList1, uList1, siguList1, '1'))
#alpha = 0.1721

print('alpha2= ', powerLawFitting(UFS2, 0.5, zList2, uList2, siguList2, '2'))
#alpha = 0.2102

'''
#4. Fit a logarithmic profile for the lowest 15% of the boundary layer
'''


def logLaw(x, ustar, z_0,):
    return (ustar / 0.41) * np.log(x / z_0)

def logLawFitting(zList, uList, siguList, title):
    z_split = 0.15 * zList[0]
    
    z_values = []
    u_values = []
    sigu_values = []
    for i, val in enumerate(zList):
        if val < z_split:
            z_values.append(val)
            u_values.append(uList[i])
            sigu_values.append(siguList[i])
  
    
    pars, cov = scipy.optimize.curve_fit(f=logLaw, xdata=z_values, ydata=u_values, sigma=sigu_values) #, p0 = [1, 1])
    ustar = pars[0]
    z_0 = pars[1]
    
    
    fig4, ax4 = plt.subplots()
    ax4.grid(True)
    ax4.set_yscale('log')
    ax4.set_ylim(2.3,3.5)
    ax4.set_xlim(0,0.1)
    
    ax4.scatter(z_values, u_values, color='dimgrey',label='Data')
    ax4.plot(z_values, logLaw(z_values, ustar, z_0), color='black', label='Fitted line')
    
    plt.xlabel('Velocity u [m/s]')
    plt.ylabel('Depth z [m] (logarithmic scale)')
    #plt.title('Log law fitting for the lowest 15% part of dataset ' + title)
    plt.legend(loc='best')
    plt.savefig('plots/logLawFitting'+title+'.png', bbox_inches='tight')
    return ustar, z_0

print(logLawFitting(zList1, uList1, siguList1, '1'))
#ustar = 0.1999, z_0 = 5.198e-05

print(logLawFitting(zList2, uList2, siguList2, '2'))
#ustar = 0.2565, z_0 = 3.722e-04

'''
#5. Compute the streamwise turbulence intensity 
'''

def turbIntensity(UFS, sigulist):
    
    Iu = []
    for sig in sigulist:
        Iu.append(sig / UFS)
        
    return Iu

Iu1 = turbIntensity(UFS1, siguList1)
Iu2 = turbIntensity(UFS2, siguList2)

#print('Iu1= ', Iu1)
#print('Iu2= ', Iu2)
fig5, ax5 = plt.subplots()
ax5.grid(True)
ax5.set_xlim(0,0.5)

ax5.plot(zList1, Iu1, label='Dataset 1', color='dimgrey')
ax5.plot(zList2, Iu2, label='Dataset 2', color='silver')

plt.xlabel('Turbulence intesity Iu')
plt.ylabel('Depth z [m]')
#plt.title('Streamwise turbulence intensity')
plt.legend(loc='best')
plt.savefig("plots/Streamwise_turbulence_intesity.png", bbox_inches='tight')


