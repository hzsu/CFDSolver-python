#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0


import numpy
import allClass

def outTecplot(gridIn,far,status,fileName):
    outPoint = numpy.zeros((gridIn.num[0],7))


    for i in range(gridIn.num[2]):
        for j in range(3):
            point = int(gridIn.cell[i,j])
            outPoint[point,2] += gridIn.volumn[i]
            outPoint[point,3] += gridIn.volumn[i] * status.Phy[i,0]     #U
            outPoint[point,4] += gridIn.volumn[i] * status.Phy[i,1]     #V
            outPoint[point,5] += gridIn.volumn[i] * status.Phy[i,2]     #P
            outPoint[point,6] += gridIn.volumn[i] * status.Phy[i,3]     #Rho
    
    for k in range(gridIn.num[0]):
        outPoint[k,3] = outPoint[k,3] / outPoint[k,2]
        outPoint[k,4] = outPoint[k,4] / outPoint[k,2]
        outPoint[k,5] = outPoint[k,5] / outPoint[k,2]
        outPoint[k,6] = outPoint[k,6] / outPoint[k,2]
    
    outPoint[:,0] = gridIn.point[:,0]
    outPoint[:,1] = gridIn.point[:,1]

    fid = open(fileName,'w')

    raw = f"""
Title = "20200308 0.8Ma 3deg"
Variables= "X","Y","volumn","U","V","P","Density"
Zone N={gridIn.num[0]},E={gridIn.num[2]},F=FEPOINT, ET=TRIANGLE
    """

    fid.write(raw)

    for i in range(gridIn.num[0]):
        for j in range(7):
            fid.write("%.6f" %outPoint[i,j])
            fid.write('\t')
        fid.write('\n')
    
    for i in range(gridIn.num[2]):
        for j in range(3):
            fid.write(str(gridIn.cell[i,j] + 1))
            fid.write('\t')
        fid.write('\n')

    fid.close() 
    
    pass

