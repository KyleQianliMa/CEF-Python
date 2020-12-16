import kylecef as kp
#B20  =  212
#B40  =  211
#B44  =  -111
#B60  =  -206
#B64  =  79
Expenergy        = [11.97, 25.72, 46.88, 85.59]
uncEnergy        = [0.09 ,  0.06,  0.24,  1.26]
Expscattering    = [0.69 ,  1.00,  0.45, 0.040]
UncExpscattering = [0.02 ,  0.00,  0.04, 0.003]

counter=0
bestchi=1e99
step = 5

for B20 in range(150,250,step):
    for B40 in range(0,100,step):
        for B44 in range(-100,0,step):
            for B60 in range(0,100,step):
                for B64 in range(0,100,step):
                  counter=counter+1
                  sol=kp.solver(B20,B40,B44,B60,B64)
                  CalcEnergy=[sol[4][2],sol[4][4],sol[4][6],sol[4][8]]
                  calcscattering=sol[0:4]
                  chi2=0
                  i=0
                  while i <=3:
                      if UncExpscattering[i]==0:
                       i=i+1
                      else:
                       chi2=chi2+((abs(calcscattering[i]-Expscattering[i]))**2)/(UncExpscattering[i]**2)
                       i=i+1
                       
                  i=0
                  while i <=3:
                      if uncEnergy[i]==0:
                       i=i+1
                      else:
                       chi2=chi2+((abs(CalcEnergy[i]-Expenergy[i]))**2)/(uncEnergy[i]**2)
                       i=i+1
                  if chi2 < bestchi:
                   bestchi=chi2
                   solution=[sol,B20,B40,B44,B60,B64]
                   #vertical rolls are eigenvectors
                  if chi2 <5:
                    break 
             
 


                
                  
                      