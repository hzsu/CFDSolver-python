#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0


import numpy
import math

import allClass
import flux
import dissipation
import timeStep
import outTec



rkCoef = [1/4, 1/3, 1/2, 1]


def toPhy(gridIn,far,status):
    
    for i in range(gridIn.num[2]):
        status.Phy[i,3] = status.W[i,0]                                 #Rho
        status.Phy[i,0] = status.W[i,1] / status.W[i,0]                 #U
        status.Phy[i,1] = status.W[i,2] / status.W[i,0]                 #V
        status.Phy[i,4] = status.W[i,3] / status.W[i,0]                 #E
        status.Phy[i,2] = (far.gammaInf - 1) * status.Phy[i,3] * (status.Phy[i,4] - (status.Phy[i,0]**2 + status.Phy[i,1]**2) / 2)
        status.Phy[i,5] = status.Phy[i,4] + status.Phy[i,2] / status.Phy[i,3]  #H

    return status


def rungeKutta(gridIn,far,status,param):

    
    
    # wLast = numpy.zeros((gridIn.num[2],4))
    # rLast = numpy.zeros((gridIn.num[2],4))

    wZero = status.W.copy()
    # status = flux.fluxGo(gridIn,far,status)

    fiderr = open('err.txt','w')

    for i in range(param.step):
        
        # status = flux.rGo(gridIn,far,status,param)
        status = timeStep.timeLocal(gridIn,far,status,param)

        for j in range(4):
            
            if j == 0:
                status = dissipation.dissipationGo(gridIn,far,status,param)
            
            for k in range(gridIn.num[2]):
                status.W[k,:] = wZero[k,:] + rkCoef[j] * status.Time[k] * status.R[k,:]
                
            status = flux.fluxGo(gridIn,far,status)
            status = flux.rGo(gridIn,far,status,param)
            status = toPhy(gridIn,far,status)
            pass

        
        err = 0.0
        for k in range(gridIn.num[2]):
            err += ((status.W[k,0] - wZero[k,0]) / status.Time[k]) ** 2
        if i == 0:
            errZero = err
        err = err / errZero

        wZero = status.W.copy()

        fiderr = open('err.txt','a')
        fiderr.write(str(i)+'\t'+str(err)+'\n')
        print(str(err))  

        if i % 5  == 0:
            fileName = 'output'+str(i)+'.plt'
            outTec.outTecplot(gridIn,far,status,fileName)

        fiderr.close()

    return status



