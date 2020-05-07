#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0


import numpy
import math

class grid:
    def __init__(self,num):
        self.num = num.copy()
        self.point = numpy.zeros((num[0],2))
        self.edge = numpy.zeros((num[1],4))     #a       b      cellLeft    cellRight
        self.cell = numpy.zeros((num[2],3))     #edge0   edge1  edge2       edge3
        self.volumn = numpy.zeros((num[2],1))   #P0      P1      P2

class farField:
    def __init__(self,pInf,rhoInf,gammaInf,maInf,alphaInf):
        self.Alpha = alphaInf
        self.c = math.sqrt(gammaInf * pInf / rhoInf)
        self.U = maInf * self.c * math.cos(alphaInf/57.3)
        self.V = maInf * self.c * math.sin(alphaInf/57.3)
        self.P = pInf
        self.Rho = rhoInf
        self.E = pInf / (rhoInf * (gammaInf - 1)) + (self.U ** 2 + self.V ** 2) / 2
        self.H = self.E + self.P / self.Rho
        self.gammaInf = gammaInf
        self.s = self.P / self.Rho ** self.gammaInf
        

class calParam:
    def __init__(self,kTwo,kFour,CFL,step):
        self.kTwo = kTwo
        self.kFour = kFour
        self.CFL = CFL
        self.step = step



class cell:
    def __init__(self,gridIn):
        self.Phy = numpy.zeros((gridIn.num[2],6))   #physical status  U V P Rho E H
        self.Q = numpy.zeros((gridIn.num[2],4))     #flux
        self.D = numpy.zeros((gridIn.num[2],4))     #artificial dissipation
        self.W = numpy.zeros((gridIn.num[2],4))     #rho rhoU rhoV rhoE
        self.R = numpy.zeros((gridIn.num[2],4))     #
        self.Time = numpy.zeros(gridIn.num[2])
        self.alphaI = numpy.zeros(gridIn.num[1])


