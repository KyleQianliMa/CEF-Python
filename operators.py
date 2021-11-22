import numpy as np
from numpy import sqrt

class Operator():   
 def Jplus(self,J):
    Jplus = np.zeros((int(2*J+1),int(2*J+1)))
    for i in range(int(2*J)):
        Jplus[i+1,i] = sqrt((J-(-J+i))*(J+(-J+i)+1))
    return Jplus

 def Jminus(self,J):
    Jminus =np.zeros((int(2*J+1),int(2*J+1)))
    for i in range(int(2*J)):
        Jminus[i,i+1] = sqrt((J+(-(J-1)+i))*(J-(-(J-1)+i)+1))       
    return Jminus

 def Jz(self,J):    
    Jz = np.zeros((int(2*J+1),int(2*J+1)))
    for i in range(int(2*J+1)):
     Jz[i,i] = (-J+i)
    return Jz

def H(D20,D40,D60,D44,D64):
 H=D20*O20+D40*O40+D60*O60+D44*O44+D64*O64
 Eigenvalues, Eigenvectors = (np.linalg.eigh(H))
 Energy=Eigenvalues-Eigenvalues[0]
 return Energy, Eigenvectors

def scattering(i,j,jx,jy,jz):
    #k=8.6173324*10^(-2)
    x=np.dot(np.dot(i,jx),j)
    y=np.dot(np.dot(i,jy),j)
    z=np.dot(np.dot(i,jz),j)
    S=np.dot(np.conj(x),x)+np.dot(np.conj(y),y)+np.dot(np.conj(z),z)
    return S

#start calculating Nd3+ matrces
Nd=Operator()
J=9/2
jplus=np.matrix(Nd.Jplus(J))
jminus=np.matrix(Nd.Jminus(J))
jz=np.matrix(Nd.Jz(J))
jx=0.5*(jplus+jminus)
jy=(jminus-jplus)/(2*1j)
jsquare=np.dot(jx,np.transpose(np.conj(jx)))+np.dot(jy,np.transpose(np.conj(jy)))+np.dot(jz,np.transpose(np.conj(jz)))
jplussquare = jplus*jplus
jminussquare = jminus*jminus
O20 = 3*(jz**2) - jsquare
O40 = 35*(jz**4) - 30*(jsquare*(jz**2)) + 25*(jz**2) -6*(jsquare) + 3*(jsquare**2)
O60 = 231*(jz**6) - 315*(jsquare*(jz**4)) + 735*(jz**4) + 105*((jsquare**2)*(jz**2)) - 525*(jsquare*(jz**2)) + 294*(jz**2) - 5*(jsquare**3) + 40*(jsquare**2) - 60*(jsquare)
O44 = 0.5*(jplus**4 + jminus**4)
O64 = 0.25*( ((11*jz**2)*(jplus**4+jminus**4) -jsquare*(jplus**4+jminus**4) - 38*(jplus**4+jminus**4)) +   ((jplus**4+jminus**4)*(11*jz**2) -(jplus**4+jminus**4)*jsquare - (jplus**4+jminus**4)*38) )

#Stevens operator conversion to [-1000,1000]range
k20=-281/0.912
k40=-344/(1.25e-2)
k60=-88/(2.09e-4)
k44=93/(-2.82e-2)
k64=104/(-2.77e-3)

#B20  =  212
#B40  =  211
#B44  =  -111
#B60  =  -206
#B64  =  79

def solver(B20,B40,B44,B60,B64):
 D20=B20/k20
 D40=B40/k40
 D44=B44/k44
 D60=B60/k60
 D64=B64/k64
 a=H(D20,D40,D60,D44,D64)
 Energy=a[0]
 Eigenvectors=np.matrix(a[1]).round(10)
 scattering1=scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,2],jx,jy,jz)+scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,3],jx,jy,jz)
 scattering2=scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,4],jx,jy,jz)+scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,5],jx,jy,jz)
 scattering3=scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,6],jx,jy,jz)+scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,7],jx,jy,jz)
 scattering4=scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,8],jx,jy,jz)+scattering(np.transpose(np.conj(Eigenvectors[:,0])),Eigenvectors[:,9],jx,jy,jz)
 N=scattering2
 s1=scattering1/N
 s2=scattering2/N
 s3=scattering3/N
 s4=scattering4/N
 return s1, s2, s3, s4, Energy, Eigenvectors,
