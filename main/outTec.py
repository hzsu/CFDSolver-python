#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0


import numpy
import allClass

def outTecplot(gridInIn,far,status,fileName):
    outPoint = numpy.zeros((gridInIn.num[0],7))


    for i in range(gridInIn.num[2]):
        for j in range(3):
            point = int(gridInIn.status[i,j])
            outPoint[point,2] += gridInIn.volumn[i]
            outPoint[point,3] += gridInIn.volumn[i] * status.Phy[i,0]     #U
            outPoint[point,4] += gridInIn.volumn[i] * status.Phy[i,1]     #V
            outPoint[point,5] += gridInIn.volumn[i] * status.Phy[i,2]     #P
            outPoint[point,6] += gridInIn.volumn[i] * status.Phy[i,3]     #Rho
    
    for k in range(gridInIn.num[0]):
        outPoint[k,3] = outPoint[k,3] / outPoint[k,2]
        outPoint[k,4] = outPoint[k,4] / outPoint[k,2]
        outPoint[k,5] = outPoint[k,5] / outPoint[k,2]
        outPoint[k,6] = outPoint[k,6] / outPoint[k,2]
    
    outPoint[:,0] = gridInIn.point[:,0]
    outPoint[:,1] = gridInIn.point[:,1]

    fid = open(fileName,'w')

    raw = f"""
Title = "20200308 0.8Ma 3deg"
Variables= "X","Y","volumn","U","V","P","Density"
Zone N={gridInIn.num[0]},E={gridInIn.num[2]},F=FEPOINT, ET=TRIANGLE
    """

    fid.write(raw)

    for i in range(gridInIn.num[0]):
        for j in range(7):
            fid.write("%.6f" %outPoint[i,j])
            fid.write('\t')
        fid.write('\n')
    
    for i in range(gridInIn.num[2]):
        for j in range(3):
            fid.write(str(gridInIn.status[i,j] + 1))
            fid.write('\t')
        fid.write('\n')

    fid.close() 
    
    pass


def ClCdCal(gridIn,far,status):
    
    cx = 0
    cy = 0
    p = 0.5 * far.Rho * (far.U ** 2 + far.V ** 2)
    Coef = [0,0]
        
    for i in range(gridIn.num[1]):
        
        cellK = int(gridIn.edge[i,2])
        cellP = int(gridIn.edge[i,3])

        if cellP == -1:
            
            edgeP0 = int(gridIn.edge[i,0])
            edgeP1 = int(gridIn.edge[i,1])           
            deltaX = gridIn.point[edgeP1,0] - gridIn.point[edgeP0,0]
            deltaY = gridIn.point[edgeP1,1] - gridIn.point[edgeP0,1]

            cx += status.Phy[cellK,2] / p * deltaY
            cy -= status.Phy[cellK,2] / p * deltaX
        
    Coef[0] = cy * math.cos(far.Alpha / 57.3) - cx * math.sin(far.Alpha / 57.3)
    Coef[1] = cx * math.cos(far.Alpha / 57.3) + cy * math.sin(far.Alpha / 57.3)

    return Coef

