import kylecef as kp
import numpy as np
B20  =  -0.64297937
B40  =  -0.00229734
B44  =  0.02317989
B60  =  -7.056e-05
B64  =  -0.00113058
h    =  0

k20=-281/0.912
k40=-344/(1.25e-2)
k44=93/(-2.82e-2)
k60=-88/(2.09e-4)
k64=104/(-2.77e-3)

B20 = B20*k20
B40 = B40*k40
B44 = B44*k44
B60 = B60*k60
B64 = B64*k64

sol=kp.solver(B20,B40,B44,B60,B64,h) #sol[0:3]=intensities, sol[4]=energy, sol[5] eigenvectors
CalcEnergy=[sol[4][2],sol[4][4],sol[4][6],sol[4][8]]
calcscattering=sol[0:4]
#vertical rolls are eigenvectors
print('Excitation=',sol[4].round(2),'\n','Relative Intensity=\n', np.array(sol[0:4]).round(2).real)
Eigenvectors=sol[5]