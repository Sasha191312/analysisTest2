import decimal
import math
from sympy.utilities.lambdify import lambdify
import sympy as sp
from sympy import*
from mpmath import ln

# auxiliary functions

# x = sp.symbols('x')


def float_range(a,b,section):
    while a<b:
        yield float(a)
        a += decimal.Decimal(section)


def matrix_multiply(a, b):  # A function that calculates the multiplication of 2 matrices and returns the new matrix
    rows_a = len(a)
    cols_a = len(a[0])
    rows_b = len(b)
    cols_b = len(b[0])
    if cols_a != rows_b:
        print('Number of A columns must equal number of B rows.')
    new_matrix = []
    while len(new_matrix) < rows_a:  # while len small the len rows
        new_matrix.append([])  # add place
        while len(new_matrix[-1]) < cols_b:
            new_matrix[-1].append(0.0)  # add value
    for i in range(rows_a):
        for j in range(cols_b):
            total = 0
            for k in range(cols_a):
                total += a[i][k] * b[k][j]  # mul mat
            new_matrix[i][j] = total
    return new_matrix  # return the A*B=new matrix


def create_i(matrix):  # A function that creates and returns the unit matrix
    i_mtrx = list(range(len(matrix)))  # make it list
    for i in range(len(i_mtrx)):
        i_mtrx[i] = list(range(len(i_mtrx)))

    for i in range(len(i_mtrx)):
        for j in range(len(i_mtrx[i])):
            i_mtrx[i][j] = 0.0  # put the zero

    for i in range(len(i_mtrx)):
        i_mtrx[i][i] = 1.0  # put the pivot
    return i_mtrx  # unit matrix


def inverse(matrix):  # A function that creates and returns the inverse matrix to matrix A
    new_matrix = create_i(matrix)  # Creating the unit matrix
    count = 0
    check = False  # flag
    while count <= len(matrix) and check == False:
        if matrix[count][0] != 0:  # if the val in place not 0
            check = True  # flag
        count = count + 1  # ++
    if not check:
        print("ERROR")
    else:
        temp = matrix[count - 1]
        matrix[count - 1] = matrix[0]  # put zero
        matrix[0] = temp
        temp = new_matrix[count - 1]
        new_matrix[count - 1] = new_matrix[0]
        new_matrix[0] = temp

        for x in range(len(matrix)):
            divider = matrix[x][x]  # find the div val
            if divider == 0:
                divider = 1
            for i in range(len(matrix)):
                matrix[x][i] = matrix[x][i] / divider  # find the new index
                new_matrix[x][i] = new_matrix[x][i] / divider
            for row in range(len(matrix)):
                if row != x:
                    divider = matrix[row][x]
                    for i in range(len(matrix)):
                        matrix[row][i] = matrix[row][i] - divider * matrix[x][i]
                        new_matrix[row][i] = new_matrix[row][i] - divider * new_matrix[x][i]
    return new_matrix  # Return of the inverse matrix

# our methods


def polynomial_interpolation(points, requested_p):
    # creating a new matrix
    mat = list(range(len(points)))
    for i in range(len(mat)):
        mat[i] = list(range(len(mat)))
    for row in range(len(points)):
        mat[row][0] = 1
    for row in range(len(points)):
        for col in range(1, len(points)):
            mat[row][col] = pow(points[row][0], col)
    res_mat = list(range(len(points)))
    for i in range(len(res_mat)):
        res_mat[i] = list(range(1))
    for row in range(len(res_mat)):
        res_mat[row][0] = points[row][1]
    vector_a = matrix_multiply(inverse(mat), res_mat)
    print('a[0]->a[%.f] =' % (len(points)-1), vector_a)
    sum1 = 0
    for i in range(len(vector_a)):
        if i == 0:
            sum1 = vector_a[i][0]
        else:
            sum1 += vector_a[i][0] * requested_p ** i
    print('P%.f(%.2f) = %.10f' % (len(points) - 1, requested_p, sum1))


def neville_interpolation(requested_p):
    n = len(xVectors)
    p = n*[0]
    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[i] = yVectors[i]
            else:
                p[i] = ((requested_p - xVectors[i + k]) * p[i] +
                        (xVectors[i] - requested_p) * p[i + 1]) / \
                       (xVectors[i] - xVectors[i + k])
    print("p(x=%.2f) =" % requested_p, p[0])


def newtonRaphson(f,a,b):
    eps = 0.00000001
    itr = 0
    xr = (a+b)/2
    next_xr = b
    error = (-(ln((eps) / (b-a)) / ln(2)))
    while abs(next_xr-xr)>eps and itr < error:
        if f(xr) == 0:
            print("cant divided by zero")
            return
        temp = next_xr
        next_xr = xr - (f(xr)/f_prime(xr))
        xr = temp
        itr = itr + 1
        print("itr:",itr,"xr:",xr,"f(x):",f(xr),"f'(x):",f_prime(xr))
    print("root:",xr)
    itr = 0


def secantMethod(f,a,b):
    eps = 0.0000000001
    itr = 0
    prev_xr = (a+b)/2
    xr = a + 0.1
    next_xr = xr
    error = (-(ln((eps) / (b-a)) / ln(2)))
    while abs(xr - prev_xr) > eps and itr < error:
        if f(xr) - f(prev_xr) == 0:
            print("cant divided by zero")
            return
        temp = next_xr
        next_xr = (prev_xr*f(xr)-xr*f(prev_xr)) / (f(xr)-f(prev_xr))
        xr = next_xr
        prev_xr = temp
        itr+=1
        print("itr:",itr,"xr:",next_xr)
    print("root:",next_xr)
    itr = 0



def simpson_method(func, sections):
    a = 0
    b = 1

    if sections % 2 != 0:
        return None
    h = (b-a)/sections
    first = func(a)
    last = func(b)
    x = a
    sum = 0
    for i in range(sections-1):
        x += h
        value = func(x)
        if i % 2 == 0:
            sum += 4 * value
        else:
            sum += 2 * value
    total = (h/3)*(first+sum+last)
    print('solution =', total)

    # f(1) is the max value for this function
    error = 1/180 * h**4 * (b-a) * func(1)
    print('error =', error)


def program_for_roots(f,a,b):
    iterationRange = 0.1  # sections
    iterationRange = float_range(a, b, iterationRange)
    for i in iterationRange:
        if f(a) * f(a + 0.1) < 0:
            newtonRaphson(f, i, i + 0.1)
        a = a + 0.1



x = symbols('x')
f = sin(2*x**3 + 5*x**2 - 6) / 2 * math.e ** (-2*x)
f_prime = f.diff(x)
f = lambdify(x,f)
f_prime = lambdify(x,f_prime)
program_for_roots(f,-1, 1.5)
simpson_method(f, 4)
# PART B


points = [[1.2, 3.5095], [1.3,3.6984], [1.4,3.9043], [1.5,4.1293], [1.6,4.3756]]
xVectors = [1.2, 1.3, 1.4, 1.5, 1.6]
yVectors = [3.5095, 3.6984, 3.9043, 4.1293, 4.3756]

print("\n\n\n\ndrill for part b\n my id 321223075\nlast digit = 5\n question 5")

print("###neville interpolation###")
neville_interpolation(1.37)
print("###polynmial interpolation###")
polynomial_interpolation(points, 1.37)
