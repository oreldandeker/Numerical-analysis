import math

##Orel Dandeker 318554102


def Linear_interpolation(Table, pointToFindVal):
    IncPoint = math.ceil(pointToFindVal)  ##ROUND UP
    DecPoint = IncPoint - 1  ##ROUND DOWN

    FxDecPoint = Table[DecPoint][1]
    FxIncPoint = Table[IncPoint][1]
    print(FxLinearFormula(FxDecPoint, DecPoint, FxIncPoint, IncPoint, pointToFindVal))  ##send values to the formula


def FxLinearFormula(y1, x1, y2, x2, pointToFindVal):
    result = (((y1 - y2) / (x1 - x2)) * pointToFindVal) + (((y2 * x1) - (y1 * x2)) / (x1 - x2))  ## calculate by formula
    return result


def Polynomial_interpolation(Table, pointToFindVal):
    mat = []  ##set matrix by values
    for i in range(len(Table)):
        row = []
        row.append(1)
        mat.append(row)

    for i in range(len(Table)):
        for j in range(len(Table) - 1):
            x = Table[i][0] ** (j + 1)
            mat[i].append(x)

    RESvector = []
    for i in range(len(Table)):
        RESvector.append([])
        RESvector[i].append(Table[i][1])

    Inversedmat = GeussianElimination(mat, RESvector, len(mat))  ##calculate result vector
    finalAnsVec = MulMatrixVector(Inversedmat, RESvector, len(mat))

    x1 = pointToFindVal

    finalPolAns = finalAnsVec[0][0]  ##calculate polynom in for loop

    for i in range(1, len(finalAnsVec)):
        finalPolAns = finalPolAns + (finalAnsVec[i][0] * x1 ** i)
    print(finalPolAns)


def Spline(Table, pointToFindVal):
    mat = []  ##set matrix by values
    for i in range(len(Table)):
        row = []
        row.append(1)
        mat.append(row)

    for i in range(len(Table)):
        for j in range(len(Table) - 1):
            x = Table[i][0] ** (j + 1)
            mat[i].append(x)

    RESvector = []
    for i in range(len(Table)):
        RESvector.append([])
        RESvector[i].append(Table[i][1])

    Inversedmat = GeussianElimination(mat, RESvector, len(mat))  ##calculate result vector
    finalAnsVec = MulMatrixVector(Inversedmat, RESvector, len(mat))

    x1 = pointToFindVal

    finalPolAns = finalAnsVec[0][0]  ##calculate polynom in for loop

    for i in range(1, len(finalAnsVec)):
        finalPolAns = finalPolAns + (finalAnsVec[i][0] * x1 ** i)
    print(finalPolAns)


def Lagrange_interpolation(Table, pointToFindVal):  ##calculate by lagrange
    yp = 0

    for i in range(len(Table)):
        p = 1
        for j in range(len(Table)):
            if i != j:
                p = p * (pointToFindVal - Table[j][0]) / (Table[i][0] - Table[j][0])
        yp = yp + p * Table[i][1]

    print(yp)


def NevAlgorithm(Table2, pointToFindVal):  ##calcuate Nevil algorithm by the nested formula
    n = len(Table2)
    total = 0
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            Table2[i][1] = ((pointToFindVal - Table2[i - j][0]) * Table2[i][1] - (pointToFindVal - Table2[i][0]) *
                            Table2[i - 1][1]) / (Table2[i][0] - Table2[i - j][0])
    total = Table2[n - 1][1]
    print(total)


# function for matrix multiply (allowed only for matrix from same size) , return a matrix
def multiplyMatrix(matrix1, matrix2, rowsCol):
    result = []
    for q in range(rowsCol):
        result.append([])
        for m in range(rowsCol):  ## creating new deafult matrix
            result[q].append(0)
    for i in range(len(matrix1)):
        for j in range(len(matrix2[0])):
            for k in range(len(matrix2)):  ##multiply rowXcol
                result[i][j] += matrix1[i][k] * matrix2[k][j]
    return result


##This function is used for multiply matrix with 1D vector , return a vector
def MulMatrixVector(InversedMat, b_vector, rows):
    result = []
    for q in range(rows):
        result.append([])  ## creating new deafult matrix
        result[q].append(0)
    for i in range(len(InversedMat)):
        for j in range(len(b_vector[0])):  ##mul inversed matrix with vector b
            for k in range(len(b_vector)):
                result[i][j] += InversedMat[i][k] * b_vector[k][j]
    return result


##This function get a matrix and returns NEW matrix with same values
def copy_matrix(M, rowsCol):
    zerosMat = []
    for q in range(rowsCol):
        zerosMat.append([])
        for m in range(rowsCol):  ## creating new deafult temporary matrix
            zerosMat[q].append(0)
    for i in range(rowsCol):
        for j in range(rowsCol):
            zerosMat[i][j] = M[i][j]  ##copy the elements
    return zerosMat


##This function is calculate The inverse matrix by GeussianElimination method ,the function returns the inversed matrix
def GeussianElimination(matrix1, b_vector, rowsCol):
    elementaricMatrix = []  ##create new elementaric matrix
    for q in range(rowsCol):
        elementaricMatrix.append([])
        for m in range(rowsCol):
            elementaricMatrix[q].append(0)
        elementaricMatrix[q][q] = 1

    allElementericSMultiply = copy_matrix(elementaricMatrix, rowsCol)
    for i in range(
            rowsCol):  ##finding for each column the max index val and if needed replace with row with the pivot row
        ElementaricM = copy_matrix(elementaricMatrix, rowsCol)
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
            ElementaricM[i][i] = 0
            ElementaricM[i + linetoreplace][i + linetoreplace] = 0
            ElementaricM[i + linetoreplace][i] = 1
            ElementaricM[i][i + linetoreplace] = 1
            allElementericSMultiply = multiplyMatrix(ElementaricM, allElementericSMultiply, rowsCol)
            matrix1 = multiplyMatrix(ElementaricM, matrix1, rowsCol)
            ElementaricM = copy_matrix(elementaricMatrix, rowsCol)

        ElementaricM[i][i] = 1 / matrix1[i][i]
        allElementericSMultiply = multiplyMatrix(ElementaricM, allElementericSMultiply, rowsCol)
        matrix1 = multiplyMatrix(ElementaricM, matrix1, rowsCol)

        for j in range(i + 1,
                       rowsCol):  ##multiply by the matched elementeric to make all indexed under the pivot to zeros
            ElementaricM = copy_matrix(elementaricMatrix, rowsCol)
            ElementaricM[j][i] = -matrix1[j][i]
            allElementericSMultiply = multiplyMatrix(ElementaricM, allElementericSMultiply, rowsCol)
            matrix1 = multiplyMatrix(ElementaricM, matrix1, rowsCol)
    for i in range(rowsCol - 1, 0, -1):  ## above all pivots makes them zeros
        for j in range(i - 1, -1, -1):
            ElementaricM = copy_matrix(elementaricMatrix, rowsCol)
            ElementaricM[j][i] = -matrix1[j][i]
            allElementericSMultiply = multiplyMatrix(ElementaricM, allElementericSMultiply, rowsCol)
            matrix1 = multiplyMatrix(ElementaricM, matrix1, rowsCol)

    return allElementericSMultiply  ##return the multiply of the elemnteric whith is that the inverse matrix


# ------------------------main-------------------------------------
Table = [[0.2, 13.7241], [0.35, 13.9776], [0.45, 14.0625], [0.6, 13.9776], [0.75, 13.7241], [0.85, 13.3056], [0.9, 12.7281]]
pointToFindVal = 2.5
pointToFindVal1 = 3
method = int(input("Choose your method : 1-Linear   2-Polynomial    3-Lagrange  4-Nevil  5-Spline\n"))
while 1:
    if method == 1:
        Linear_interpolation(Table, 0.65)
    elif method == 2:
        Polynomial_interpolation(Table, 0.65)
    elif method == 3:
        Lagrange_interpolation(Table, 0.65)
    elif method == 4:
        NevAlgorithm(Table, 0.65)
    elif method == 5:
        Spline(Table, 0.65)
    else:
        print("Invalid choice")
    method = int(input("Choose your method : 1-Linear   2-Polynomial    3-Lagrange  4-Nevil  5-Spline\n"))
