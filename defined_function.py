###### DEFINED FUNCTIONS ######

import numpy as np

def xy_intersect(a1,a2,b1,b2): #Find intersection of two lines
    s=np.vstack([a1,a2,b1,b2])        
    h=np.hstack((s, np.ones((4, 1)))) 
    l1=np.cross(h[0], h[1]) #first line
    l2=np.cross(h[2], h[3]) #second line
    x,y,z=np.cross(l1, l2) #intersection
    if z == 0:                          
        return (float('inf'), float('inf'))
    return (x/z, y/z)

def NaNlist(n): #Create list of NaN values
    listofzeros = [None] * n
    return listofzeros

def zerolist(n): #Create list of zero values
    listofzeros = [0] * n
    return listofzeros


def find_D(D,a,dt,m,ShorelineChange,i,Qs,Qd):
    D=D+a*dt #Calculate new depth, accounting for sea level rise
    ShorelineChange[i]=-(Qs/D*dt)+(Qd/D*dt) #Store shoreline change values per dt ignoring acccomodation
    D=D+m*(abs(ShorelineChange[i]))
    return D, ShorelineChange[i]

def D_deplenished(D,a,dt,m,ShorelineChange,i,Qs,Qd,Shoreline):
    xs = abs(ShorelineChange[i])
    # shoot the shoreine seaward until we find a depth not equal to zero
    while D <= 0.1:
        xs-=1
        D = D + (a*dt-m*xs)
        
    ShorelineChange[i] = xs
    Shoreline = Shoreline+ShorelineChange[i]
    return D,ShorelineChange[i],Shoreline



