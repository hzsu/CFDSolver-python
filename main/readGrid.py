#南京航空航天大学 GYS课题组
#CFD python计算程序
#By 黄紫, hznt@nuaa.edu.cn
#V1.0



import allClass
import numpy

def readG(fileName):
    
    
    # fileName = input('请输入')
    print(fileName)
    fid = open(fileName,'r')
    num = [0]*3

    basicNum = fid.readline().strip().split("  ")
    num[0] = int(basicNum[0])
    num[1] = int(basicNum[1])
    num[2] = int(basicNum[2])

    gridIn = allClass.grid(num)
    gridIn.point = numpy.zeros((num[0],2))
    gridIn.edge =numpy.zeros((num[1],4))
    gridIn.cell =numpy.zeros((num[2],3))
    gridIn.volumn = numpy.zeros((num[2],1))

    print(str(num[0]),str(num[1]),str(num[2]))
    
    for i in range(num[0]):
        basicNum = fid.readline().strip().split()
        for j in range(2):
            gridIn.point[i][j] = float(basicNum[j])

    for i in range(num[1]):
        basicNum = fid.readline().strip().split()
        for j in range(4):
            gridIn.edge[i][j] = int(basicNum[j])-1          #网格文件中应该是从1开始编号的
        if (gridIn.edge[i][3] == -2) or (gridIn.edge[i][3] == -3):
            gridIn.edge[i][3] += 1


    for i in range(num[2]):
        basicNum = fid.readline().strip().split()
        for j in range(3):
            gridIn.cell[i][j] = int(basicNum[j]) - 1


    for i in range(num[2]):
        basicNum = fid.readline().strip().split()
        for j in range(1):
            gridIn.volumn[i][j] = float(basicNum[j])

    fid.close()

    return gridIn





    


