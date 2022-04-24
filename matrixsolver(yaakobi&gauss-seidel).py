#Orel Dandeker 318554102


##get a matrix and returns new matrix with same values
def copy_matrix(M,rowsCol):
    zerosMat=[]
    for q in range(rowsCol):
        zerosMat.append([])
        for m in range(rowsCol):
            zerosMat[q].append(0)
    for i in range(rowsCol):
        for j in range(rowsCol):
            zerosMat[i][j] = M[i][j]
    return zerosMat

##get a vector and returns new matrix with same values
def copy_vector(V,rowsCol):
    zerosVec=[]
    for q in range(rowsCol):
        zerosVec.append([])
    for m in range(rowsCol):
        zerosVec[m].append(0)
    for i in range(rowsCol):
        zerosVec[i][0] = V[i][0]
    return zerosVec

##pivoting first column
def Pivoting(matrix1,b_vector,rowsCol):

    newMat = copy_matrix(matrix1, rowsCol)
    newVec = copy_vector(b_vector,rowsCol)
    print(newMat)
    print(newVec)
    for i in range(1):
        pivot = matrix1[i][i]
        count = 0
        linetoreplace = 0
        max = pivot
        for j in range(1, rowsCol - i):
            count = count + 1
            if abs(matrix1[i + j][i]) > max:
                max = abs(matrix1[i + j][i])
                linetoreplace = count
        if linetoreplace > 0:
            newMat[i] = matrix1[linetoreplace+i]
            newMat[linetoreplace+i] = matrix1[i]
            newVec[i] = b_vector[linetoreplace+i]
            newVec[linetoreplace + i] = b_vector[i]

            matrix1[i] = newMat[i]
            matrix1[linetoreplace] = newMat[linetoreplace]
            b_vector[i] = newVec[i]
            b_vector[linetoreplace] = newVec[linetoreplace]

    return newMat
#check if matrix has string diagonal and returns boolean value
def maxDiag(mat1 ,rowsCol):
    for i in range(len(mat1)):
        SumAllElements = 0
        for j in range(len(mat1)):
            SumAllElements = SumAllElements + abs(mat1[i][j])
        if abs(mat1[i][i]) < SumAllElements - abs((mat1[i][i])):
            return False
    return True

#yaakobi method
def Yaakobi(mat1,b_Vector,rowsCol):
    x=0
    y=0
    z=0
    iterationNum=1
    afterPivot = Pivoting(mat1,b_Vector,rowsCol)
    if maxDiag(afterPivot, rowsCol):
        print("THE MATRIX IS DOMINATE DIAG:")
    else:
        print("THE MATRIX ISNT HAVE DOMINATE DIAG:")
    while(True):
        x1 = (b_Vector[0][0] - afterPivot[0][1] * y - afterPivot[0][2] * z) / afterPivot[0][0]
        y1 = (b_Vector[1][0] - afterPivot[1][0] * x - afterPivot[1][2] * z) / afterPivot[1][1]
        z1 = (b_Vector[2][0] - afterPivot[2][0] * x - afterPivot[2][1] * y) / afterPivot[2][2]
        if(abs(x1-x) < 0.00001 ) :
            break
        else:
            x = x1
            y = y1
            z = z1
        print(f"iteration #{iterationNum}: x={x},y={y},z={z}")
        iterationNum=iterationNum+1

##guess ziedel method
def GaussSeidel(mat1,b_Vector,rowsCol):
    x=0
    y=0
    z=0
    iterationNum=1
    afterPivot = Pivoting(mat1, b_Vector, rowsCol)
    if maxDiag(afterPivot, rowsCol):
        print("THE MATRIX IS DOMINATE DIAG:")
    else:
        print("THE MATRIX ISNT DOMINATE DIAG:")
    while(True):
        x1 = (b_Vector[0][0] - afterPivot[0][1] * y - afterPivot[0][2] * z) / afterPivot[0][0]
        y1 = (b_Vector[1][0] - afterPivot[1][0] * x1 - afterPivot[1][2] * z) / afterPivot[1][1]
        z1 = (b_Vector[2][0] - afterPivot[2][0] * x1 - afterPivot[2][1] * y1) / afterPivot[2][2]
        if((abs(x1-x) < 0.00001 ) ) :
            break
        else:
            x = x1
            y = y1
            z = z1
        print(f"iteration #{iterationNum} : x={x},y={y},z={z}")
        iterationNum=iterationNum+1



#matrix and vector that you asked to resolve
mat1=[[4,2,0] ,[2,10,4], [0,4,5]]
b_Vector=[[2],[6],[5]]
rowsCol =len(mat1)


option = int(input('press 1 for yaakobi\npress 2 for GaussSeidel\n'))

if option == 1:
    print("YOU USING YAAKOBI METHOD")
    Yaakobi(mat1,b_Vector, rowsCol)
else:
    print("YOU USING GAUSS SEIDEL METHOD")
    GaussSeidel(mat1,b_Vector, rowsCol)