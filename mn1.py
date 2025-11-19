def calculate_eigenvalues(matrix):
    """
    Calculate eigenvalues using characteristic polynomial
    Works for 2x2 and 3x3 matrices
    """
    n = len(matrix)

    if n == 2:
        return calculate_eigenvalues_2x2(matrix)
    elif n == 3:
        return calculate_eigenvalues_3x3(matrix)
    else:
        print("Error: Only 2x2 and 3x3 matrices supported")
        return []


def calculate_eigenvalues_2x2(matrix):
    """
    For 2x2 matrix: solve λ² - (trace)λ + det = 0
    """
    a = matrix[0][0]
    b = matrix[0][1]
    c = matrix[1][0]
    d = matrix[1][1]

    trace = a + d
    det = a * d - b * c

    discriminant = trace**2 - 4 * det

    if discriminant < 0:
        print("Complex eigenvalues - not supported")
        return [trace / 2, trace / 2]

    sqrt_disc = discriminant**0.5
    lam1 = (trace + sqrt_disc) / 2
    lam2 = (trace - sqrt_disc) / 2
    return [lam1, lam2]


def calculate_eigenvalues_3x3(matrix):
    """
    For 3x3 matrix: use the invariant formula for the characteristic polynomial.
    Characteristic polynomial: λ³ - c2 λ² + c1 λ + c0 = 0
    where:
      c2 = trace(A)
      c1 = sum of principal 2x2 minors
      c0 = -det(A)
    """
    a = matrix[0][0]
    b = matrix[0][1]
    c = matrix[0][2]
    d = matrix[1][0]
    e = matrix[1][1]
    f = matrix[1][2]
    g = matrix[2][0]
    h = matrix[2][1]
    i = matrix[2][2]

    # trace
    c2 = a + e + i

    # sum of principal 2x2 minors
    c1 = (e * i - f * h) + (a * i - c * g) + (a * e - b * d)

    # determinant of A
    detA = a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h

    # characteristic polynomial: λ³ - c2 λ² + c1 λ - detA = 0
    eigenvalues = solve_cubic(1, -c2, c1, -detA)
    return eigenvalues


def solve_cubic(a, b, c, d):
    """
    Solve cubic equation: a x³ + b x² + c x + d = 0
    Using simple Newton iteration from a few starting points.
    """
    roots = []
    starts = [-10, 0, 10]

    for start in starts:
        x = start
        for _ in range(100):
            f = a * x**3 + b * x**2 + c * x + d
            fp = 3 * a * x**2 + 2 * b * x + c

            if fp == 0 or (fp < 0.000001 and fp > -0.000001):
                break

            x_new = x - f / fp
            diff = x_new - x

            if diff < 0.000001 and diff > -0.000001:
                # check if this root is new
                new_root = True
                for r in roots:
                    delta = x_new - r
                    if delta < 0.01 and delta > -0.01:
                        new_root = False
                        break
                if new_root:
                    roots.append(x_new)
                break

            x = x_new

    return sorted(roots, reverse=True)


def calculate_eigenvector(matrix, eigenvalue):
    """
    Solve (A - λI)v = 0 for v
    Supports 2x2 and 3x3
    """
    n = len(matrix)
    A_shifted = [[matrix[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        A_shifted[i][i] -= eigenvalue

    if n == 2:
        a = A_shifted[0][0]
        b = A_shifted[0][1]

        if b != 0:
            v = [1, -a / b]
        elif a != 0:
            v = [-b / a, 1]
        else:
            c = A_shifted[1][0]
            d = A_shifted[1][1]
            if d != 0:
                v = [1, -c / d]
            else:
                v = [0, 1]
        return v

    if n == 3:
        v1 = 1
        a11, a12, a13 = A_shifted[0]
        a21, a22, a23 = A_shifted[1]

        if a13 != 0 and a23 != 0:
            denom = a13 * a22 - a23 * a12
            if denom != 0:
                v2 = (a23 * a11 - a13 * a21) * v1 / denom
            else:
                v2 = 0
            v3 = -(a11 * v1 + a12 * v2) / a13

        elif a12 != 0 and a22 != 0:
            v3 = 1
            denom = a11 * a22 - a21 * a12
            if denom != 0:
                v2 = (a21 * a13 - a11 * a23) * v3 / denom
            else:
                v2 = 0
            v1 = -(a12 * v2 + a13 * v3) / a11 if a11 != 0 else 1

        else:
            v2 = 1
            v3 = 1
            if a11 != 0:
                v1 = -(a12 * v2 + a13 * v3) / a11
            else:
                v1 = 1

        return [v1, v2, v3]

    return [1] * n


def matrix_vector_mult(A, v):
    n = len(A)
    return [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]


def print_matrix(matrix):
    for row in matrix:
        print([int(round(x)) for x in row])


def print_vector(vector):
    print([int(round(x)) for x in vector])


def main():
    n = int(input("Enter matrix size (2 or 3): "))
    if n not in [2, 3]:
        print("Error: Only 2x2 and 3x3 matrices supported!")
        return

    matrix = []
    print(f"\nEnter {n}x{n} matrix (separate values with spaces):")
    for i in range(n):
        row = list(map(float, input(f"Row {i+1}: ").split()))
        if len(row) != n:
            print(f"Error: Need exactly {n} values!")
            return
        matrix.append(row)

    print("\nOriginal Matrix:")
    print_matrix(matrix)

    print("\n" + "=" * 50)
    print("EIGENVALUES")
    print("=" * 50)
    eigenvalues = calculate_eigenvalues(matrix)

    print("\nEigenvalues:")
    for i, lam in enumerate(eigenvalues):
        print(f"λ{i+1} = {int(round(lam))}")

    print("\n" + "=" * 50)
    print("EIGENVECTORS")
    print("=" * 50)
    for i, lam in enumerate(eigenvalues):
        print(f"\nEigenvector {i+1} (for λ = {int(round(lam))}):")
        v = calculate_eigenvector(matrix, lam)
        print_vector(v)
        Av = matrix_vector_mult(matrix, v)
        lamv = [lam * x for x in v]
        print(f"  A*v = {[int(round(x)) for x in Av]}")
        print(f"  λ*v = {[int(round(x)) for x in lamv]}")
        


if __name__ == "__main__":
    main()
    #eeeeeee
