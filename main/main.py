#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0



import numpy
import math


import allClass
import readGrid
import flux
import dissipation
import timeStep
import scheme
import outTec


def initialization(gridIn,far):
    status = allClass.cell(gridIn)
    
    for i in range(gridIn.num[2]):
        status.Phy[i,:] = [far.U,far.V,far.P,far.Rho,far.E,far.H]
        status.W[i,:] = [far.Rho, far.Rho*far.U, far.Rho*far.V, far.Rho*far.E]

    return status

far = allClass.farField(101330,1.225,1.4,0.8,3) #pInf,rhoInf,gammaInf,maInf,alphaInf
param = allClass.calParam(1,1/32,0.5,10000)       #kTwo,kFour,CFL,step
gridIn = readGrid.readG('naca0012.grd')

status = initialization(gridIn,far)

fileName = 'init.plt'
outTec.outTecplot(gridIn,far,status,fileName)

status = scheme.rungeKutta(gridIn,far,status,param)
