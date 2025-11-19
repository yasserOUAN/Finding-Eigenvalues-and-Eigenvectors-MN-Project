import cmath
import random

# function to create and fill a 2*2 matrix
def fillmatrix ():
    A = [[0.0, 0.0], [0.0, 0.0]]
    for i in range(2):
        for j in range(2):
            A[i][j]=float (input("A[%d][%d]=" % (i , j )))
    return A

# function to choice the correct method to calculate lamba
def eigenvalue (A):
    """
    For 2x2 matrix: solve λ² - (b)λ + det = 0
    """
    #A = Matrix

    b=A[0][0] * (-1) + A[1][1] * (-1)
    det=(A[0][0] * A[1][1]) - (A[0][1] * A[1][0])
    discriminant = b**2 - 4 * det

    if discriminant>0:
        t1 = (- b + (discriminant ** 0.5)) / 2
        t2 = (- b - (discriminant ** 0.5)) / 2
    elif discriminant == 0 :
        t1 =- b / 2
        t2 = t1
    elif discriminant < 0:
        print("Complex eigenvalues - not supported")
        t1 = (- b - cmath.sqrt(discriminant)) / 2
        t2 = (- b + cmath.sqrt(discriminant)) / 2
    return t1, t2


#function to find eigenvector(s) depending on the nature of lamba
def eigenvector (A,t1,t2): 
    
    """
    Calculates eigenvectors by solving (A - λI)v = 0
    
    For a 2x2 matrix and eigenvalue λ, we get:
    [A[0][0]-λ    A[0][1]  ] [x]   [0]
    [A[1][0]      A[1][1]-λ] [y] = [0] 
    
    """
    
    #A = Matrix
    #t1 eignevector 1
    #t2 eignevector 2

    v2=[0,0]
    a = int( random.randint(1, 10) * random.choice([-1, 1]))
    c = int( random.randint(1, 10) * random.choice([-1, 1]))

    x=(A[0][0] - t1) / A[0][1]

    # lamba 1 is complex
    if type(t1) is complex:
        if (A[0][0] - t1).real > 0 and A[0][1] > 0 or (A[0][0] - t1).real < 0 and (A[0][1]).real < 0:
            x = x * (- 1)
        v1 = [a, complex(a) * x]


    # lamba 1 is int
    else:
        if (A[0][0] - t1) > 0 and A[0][1] > 0 or (A[0][0] - t1) < 0 and A[0][1] < 0 :
            x = x * (- 1)
        v1 = [a, a * x]
    print("from the eigenvalues we deduce that the vector can be written as")
    print(f"-----> x = ( {x} ) * y")
    print(f"Let x = ({a}) so y = {a * x} \n")

    # there is lamba 2
    if t1 != t2:
        # lamba 2 is complex
        if type(t1) is complex:
            y = complex(A[1][0]) / (complex(A[1][1]) - t2)

            if (A[0][0] - t1).real > 0 and A[0][1] > 0 or (A[0][0] - t1).real < 0 and (A[0][1]).real < 0:
                y = y * (- 1)
            v2 = [c, complex(c) * y]

        else: # lamba 2 is int
            y = A[1][0] / (A[1][1] - t2)

            if (A[1][1] - t1) > 0 and A[1][0] > 0 or (A[1][1] - t1) < 0 and A[1][0] < 0:
                y = y * (- 1)
            v2 = [c, c * y]
        print("and the second vector can be written as")
        print(f"-----> X = ( {y} ) * Y")
        print(f"Let X = ({a}) so Y = {a * y} \n")

    return v1 , v2

# function to print the lambas
def values():
    if t1 == t2:
        print(f"The matrix have one eigenvalue and it is: \nλ1 = {t1} ")
    else:
        print(f"The matrix have two eigenvalues and they are: \nλ1 = {t1}  \nλ2 = {t2} ")
    print("\n")


# function to print the vectors
def vectors():
    if v2 == [0, 0]:
        print(f"The matrix have one eigenvector and it is: \n{v1} ")
    else:
        print(f"The matrix have two eigenvectors and they are:\nfor λ1 = {t1} -> {v1}  \nfor λ2 = {t2} -> {v2} ")


# main
# matrix creation
A = [[0.0, 0.0], [0.0, 0.0]]
print("Fill a 2x2 matrix :\n")
A= fillmatrix()
print("\nYour matrix :")
print(f"{A[0]}\n{A[1]}")


# finding and printing eigenvalue
print("\n" + "=" * 50)
print("EIGENVALUES")
print("=" * 50)
t1,t2= eigenvalue (A)
values()



# finding and printing eigenvectors
print("\n" + "=" * 50)
print("EIGENVECTORS")
print("=" * 50)
v1,v2 = eigenvector(A,t1,t2)
vectors()
