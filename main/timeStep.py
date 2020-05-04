#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0


import numpy
import math
import allClass

def timeLocal(gridIn,far,status,param):

    status.Time = numpy.zeros(gridIn.num[2])

    for i in range(gridIn.num[1]):
        cellK = int(gridIn.edge[i,2])
        cellP = int(gridIn.edge[i,3])

        status.Time[cellK] += status.alphaI[i]
        if ((cellP != -1) and (cellP != -2)): 
            status.Time[cellP] += status.alphaI[i]

    for i in range(gridIn.num[2]):
        if status.Time[i] <= 0:         #虽然10000随意，但是由于之前时间初始化为0
            status.Time[i] = 10000
        status.Time[i] = param.CFL * gridIn.volumn[i] / status.Time[i]

    return status