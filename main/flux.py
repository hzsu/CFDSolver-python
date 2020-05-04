#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0



import numpy
import math
import allClass

def fluxGo(gridIn,far,status):
    status.Q = numpy.zeros((gridIn.num[2],4))
    
    flux = [0] * 4

    for i in range(gridIn.num[1]):
        edgeP0 = int(gridIn.edge[i,0])
        edgeP1 = int(gridIn.edge[i,1])
        cellK = int(gridIn.edge[i,2])
        cellP = int(gridIn.edge[i,3])
        deltaX = gridIn.point[edgeP1,0] - gridIn.point[edgeP0,0]
        deltaY = gridIn.point[edgeP1,1] - gridIn.point[edgeP0,1]
        deltaS = math.sqrt(deltaX ** 2 + deltaY ** 2)
        
        nX = deltaY / deltaS
        nY = - deltaX / deltaS

        if cellP == -1:         #wall
            U = status.Phy[cellK,0]
            V = status.Phy[cellK,1]
            P = status.Phy[cellK,2]
            Rho = status.Phy[cellK,3]

            c = math.sqrt(far.gammaInf * P / Rho)
            z = U * deltaY - V * deltaX
            status.alphaI[i] = math.fabs(z) + c * deltaS        #artificial dissipation

            flux[0] = 0                     #PPT 8
            flux[1] = deltaY * P
            flux[2] = - deltaX * P
            flux[3] = 0
        
        elif cellP == -2:       #farfield
            U = status.Phy[cellK,0]
            V = status.Phy[cellK,1]
            direction = nX * far.U + nY * far.V
            qNInf = far.U * nX + far.V * nY
            qNE = U * nX + V * nY
            qTInf = - far.U * nY + far.V * nX
            qTE = - U * nY + V * nX
            cInf = far.c
            cE = math.sqrt(far.gammaInf * status.Phy[cellK,2] / status.Phy[cellK,3])

            rEPlus = qNE + 2*cE / (far.gammaInf - 1)                     #PPT 18
            rEMinus = qNE - 2*cE / (far.gammaInf - 1)  
            rInfPlus =  qNInf + 2*cInf / (far.gammaInf - 1)                  #PPT 18
            rInfMinus = qNInf - 2*cInf / (far.gammaInf - 1)

            if ((qNInf / far.c <= 1) and (qNInf / far.c >= -1)):        #subSonic
                if direction <= 0:                                      #inflow
                    qN = (rInfPlus + rEMinus ) / 2
                    qT = qTInf
                    c = (far.gammaInf - 1) * (rInfPlus - rEMinus) / 4
                    s = far.s
                    Rho = (c ** 2 / far.gammaInf / s) ** (1 / (far.gammaInf-1))
                    P = Rho * c ** 2 / far.gammaInf
                else:                                                   #outflow
                    qN = (rEPlus + rInfMinus) / 2
                    qT = qTInf
                    c = (far.gammaInf - 1) * (rEPlus - rInfMinus) / 4
                    s = status.Phy[cellK,2] / status.Phy[cellK,3] ** far.gammaInf
                    Rho = (c ** 2 / far.gammaInf / s) ** (1 / (far.gammaInf-1))
                    P = Rho * c ** 2 / far.gammaInf
            else:                                                       #superSonic
                if direction <= 0:                                      #inflow
                    qN = (rInfPlus + rInfMinus) / 2
                    qT = qTInf
                    c = (far.gammaInf - 1) * (rInfPlus - rInfMinus) / 4
                    s = far.s
                    Rho = (c ** 2 / far.gammaInf / s) ** (1 / (far.gammaInf-1))
                    P = Rho * c ** 2 / far.gammaInf
                else:                                                   #outflow
                    qN = (rEPlus + rEMinus) / 2
                    qT = qTE
                    c = (far.gammaInf - 1) * (rEPlus - rEMinus) / 4
                    s = status.Phy[cellK,2] / status.Phy[cellK,3] ** far.gammaInf
                    Rho = (c ** 2 / far.gammaInf / s) ** (1 / (far.gammaInf-1))
                    P = Rho * c ** 2 / far.gammaInf
            
            U = qN * nX - qT * nY
            V = qN * nY + qT * nX
            E = P / Rho / (far.gammaInf - 1) + (U ** 2 + V ** 2) / 2
            H = E + P / Rho
            
            z = U * deltaY - V * deltaX
            status.alphaI[i] = math.fabs(z) + c * deltaS

            flux[0] = z * Rho
            flux[1] = z * Rho * U + deltaY * P
            flux[2] = z * Rho * V - deltaX * P
            flux[3] = z * Rho * H
        
        else:                       #inflow
            U = (status.Phy[cellK,0] + status.Phy[cellP,0]) / 2
            V = (status.Phy[cellK,1] + status.Phy[cellP,1]) / 2
            P = (status.Phy[cellK,2] + status.Phy[cellP,2]) / 2
            Rho = (status.Phy[cellK,3] + status.Phy[cellP,3]) / 2
            H = (status.Phy[cellK,5] + status.Phy[cellP,5]) / 2

            c = math.sqrt(far.gammaInf * P / Rho)

            z = U * deltaY - V * deltaX
            status.alphaI[i] = math.fabs(z) + c * deltaS
            
            flux[0] = z * Rho
            flux[1] = z * Rho * U + deltaY * P
            flux[2] = z * Rho * V - deltaX * P
            flux[3] = z * Rho * H

        status.Q[cellK,:] += flux[:]
        if ((cellP != -1) and (cellP != -2)):
            status.Q[cellP,:] -= flux[:]

    return status

def rGo(gridIn,far,status,param):

    for i in range(gridIn.num[2]):
        status.R[i,:] = - (status.Q[i,:] - status.D[i,:]) / gridIn.volumn[i]
        
    return status