import cmath
import random

#                   function to create and fill a 2*2 matrix
def fillmatrix ():
    A = [[0.0, 0.0], [0.0, 0.0]]
    for i in range(2):
        for j in range(2):
            A[i][j]=float (input("A[%d][%d]=" % (i + 1, j + 1)))
    return A

#                   function to choice the correct method to calculate lamba
def eigenvalue (A):

    b=A[0][0] * (-1) + A[1][1] * (-1)
    c=(A[0][0] * A[1][1]) - (A[0][1] * A[1][0])
    U=b ** 2 - 4 * c

    if U>0:
        t1 = (- b + (U ** 0.5)) / 2
        t2 = (- b - (U ** 0.5)) / 2
    elif U == 0 :
        t1 =- b / 2
        t2 = t1
    else:
        t1 = (- b - cmath.sqrt(U)) / 2
        t2 = (- b + cmath.sqrt(U)) / 2
    return t1, t2


#                   function to find eigenvector(s) depending on the nature of lamba
def eigenvector (A,t1,t2):
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


    return v1 , v2

#                   function to print the lambas
def values():
    if t1 == t2:
        print(f"the matrix have one eigenvalue and it is {t1} ")
    else:
        print(f"the matrix have two eigenvalues and they are {t1} and {t2} ")
    print("\n")


#                   function to print the vectors
def vectors():
    if v2 == [0, 0]:
        print(f"the matrix have one eigenvector and it is {v1} ")
    else:
        print(f"the matrix have two eigenvectors and they are {v1} and {v2} ")

     # main
#                   matrix creation
A = [[0.0, 0.0], [0.0, 0.0]]
print("fill a 2x2 matrix \n")
A= fillmatrix()
print(f"{A[0]}\n{A[1]}")

print("===============================================================/n")
#                   finding and printing eigenvalue
t1,t2= eigenvalue (A)
values()

print("===============================================================/n")
#                   finding and printing eigenvectors
print(f"")
v1,v2 = eigenvector(A,t1,t2)
vectors()
