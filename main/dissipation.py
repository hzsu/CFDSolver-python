#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0


import numpy
import math
import allClass

def dissipationGo(gridIn,far,status,param):
    
    laplacianW = numpy.zeros((gridIn.num[2],4))
    status.D = numpy.zeros((gridIn.num[2],4))

    for i in range(gridIn.num[1]):
        cellK = int(gridIn.edge[i,2])
        cellP = int(gridIn.edge[i,3])

        if ((cellP != -1) and (cellP != -2)):    
            laplacianW[cellK,:]  += (status.W[cellP,:] - status.W[cellK,:]) / 3        #按边，应该是加了三次吧，参考有误。
            laplacianW[cellP,:] += (status.W[cellK,:] - status.W[cellP,:]) / 3

    for i in range(gridIn.num[1]):
        cellK = int(gridIn.edge[i,2])
        cellP = int(gridIn.edge[i,3])

        dTwo = numpy.zeros((gridIn.num[1],4))
        dFour = numpy.zeros((gridIn.num[1],4))

        if ((cellP != -1) and (cellP != -2)):  
            vI = math.fabs((status.Phy[cellP,2] - status.Phy[cellK,2]) / (status.Phy[cellP,2] + status.Phy[cellK,2]))
            epsilonTwo = param.kTwo * vI
            dTwo[i,:] = status.alphaI[i] * epsilonTwo * (status.W[cellP,:] - status.W[cellK,:])


            epsilonFour = param.kFour - epsilonTwo
            if epsilonFour < 0:
                epsilonFour = 0
            dFour[i,:] = - status.alphaI[i] * epsilonFour * (laplacianW[cellP,:] - laplacianW[cellK,:])
        
            status.D[cellK,:] += (dTwo[i,:] + dFour[i,:])
            status.D[cellP,:] -= (dTwo[i,:] + dFour[i,:])           #PPT P12 有误吧，各边相加，实际右侧的cell还是减法

    return status