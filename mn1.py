def calculate_eigenvalues(matrix, max_iter=1000, tolerance=1e-10):
    """
    Calculate all eigenvalues using QR algorithm
    Returns: list of eigenvalues
    """
    n = len(matrix)
    A_k = [row[:] for row in matrix]
    
    for iteration in range(max_iter):
        Q, R = qr_decomposition(A_k)
        A_new = matrix_multiply(R, Q)
        
        max_off_diag = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    max_off_diag = max(max_off_diag, abs(A_new[i][j]))
        
        if max_off_diag < tolerance:
            break
        
        A_k = A_new
    
    eigenvalues = [A_k[i][i] for i in range(n)]
    return eigenvalues


def calculate_eigenvector(matrix, eigenvalue):
    """
    Calculate eigenvector by solving (A - λI)v = 0
    Returns: eigenvector as a list
    """
    n = len(matrix)
    
    # Create (A - λI)
    A_shifted = [[matrix[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        A_shifted[i][i] -= eigenvalue
    
    # For 2x2 matrix
    if n == 2:
        a = A_shifted[0][0]
        b = A_shifted[0][1]
        
        if abs(b) > 0:
            v = [1, -a/b]
        elif abs(a) > 0:
            v = [-b/a, 1]
        else:
            c = A_shifted[1][0]
            d = A_shifted[1][1]
            if abs(d) > 0:
                v = [1, -c/d]
            else:
                v = [0, 1]
        return v
    
    # For 3x3 matrix
    if n == 3:
        v1 = 1
        
        a11, a12, a13 = A_shifted[0]
        a21, a22, a23 = A_shifted[1]
        
        if abs(a13) > 0 and abs(a23) > 0:
            denom = a13 * a22 - a23 * a12
            if abs(denom) > 0:
                v2 = (a23 * a11 - a13 * a21) * v1 / denom
            else:
                v2 = 0
            v3 = -(a11 * v1 + a12 * v2) / a13
        elif abs(a12) > 0 and abs(a22) > 0:
            v3 = 1
            denom = a11 * a22 - a21 * a12
            if abs(denom) > 0:
                v2 = (a21 * a13 - a11 * a23) * v3 / denom
            else:
                v2 = 0
            v1 = -(a12 * v2 + a13 * v3) / a11 if abs(a11) > 0 else 1
        else:
            v2 = 1
            v3 = 1
            if abs(a11) > 0:
                v1 = -(a12 * v2 + a13 * v3) / a11
            else:
                v1 = 1
        
        return [v1, v2, v3]
    
    # For larger matrices
    v = [1] * n
    for row_idx in range(n):
        non_zero_col = -1
        for col_idx in range(n):
            if abs(A_shifted[row_idx][col_idx]) > 0:
                non_zero_col = col_idx
                break
        
        if non_zero_col >= 0:
            sum_val = 0
            for j in range(n):
                if j != non_zero_col:
                    sum_val += A_shifted[row_idx][j] * v[j]
            
            v[non_zero_col] = -sum_val / A_shifted[row_idx][non_zero_col]
    
    return v

def qr_decomposition(A):
    """QR decomposition using Gram-Schmidt"""
    n = len(A)
    Q = [[0.0] * n for _ in range(n)]
    R = [[0.0] * n for _ in range(n)]
    
    columns = [[A[i][j] for i in range(n)] for j in range(n)]
    
    for j in range(n):
        v = columns[j][:]
        
        for i in range(j):
            R[i][j] = sum(Q[k][i] * columns[j][k] for k in range(n))
            v = [v[k] - R[i][j] * Q[k][i] for k in range(n)]
        
        norm = sum(x**2 for x in v)**0.5
        R[j][j] = norm
        
        if norm > 1e-10:
            for k in range(n):
                Q[k][j] = v[k] / norm
        else:
            Q[j][j] = 1.0
    
    return Q, R


def matrix_multiply(A, B):
    """Multiply two matrices"""
    n = len(A)
    result = [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return result


def matrix_vector_mult(A, v):
    """Multiply matrix by vector"""
    n = len(A)
    result = [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]
    return result


def print_matrix(matrix):
    """Print matrix"""
    for row in matrix:
        print([f"{val:10.4f}" for val in row])


def print_vector(vector):
    """Print vector"""
    print([f"{val:10.4f}" for val in vector])


def main():
    # Get matrix input
    n = int(input("Enter matrix size (n x n): "))
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
    
    # Calculate all eigenvalues
    print("\n" + "="*50)
    print("EIGENVALUES")
    print("="*50)
    eigenvalues = calculate_eigenvalues(matrix)
    
    print("\nEigenvalues:")
    for i, val in enumerate(eigenvalues):
        print(f"λ{i+1} = {int(round(val))}")
    
    # Calculate eigenvectors
    print("\n" + "="*50)
    print("EIGENVECTORS")
    print("="*50)
    
    for i, eigenval in enumerate(eigenvalues):
        print(f"\nEigenvector {i+1} (for λ = {int(round(eigenval))}):")
        eigenvec = calculate_eigenvector(matrix, eigenval)
        print_vector(eigenvec)
        
        # Verify A*v = λ*v
        Av = matrix_vector_mult(matrix, eigenvec)
        lambdav = [eigenval * v for v in eigenvec]
        
        print(f"  A*v = {[int(round(x)) for x in Av]}")
        print(f"  λ*v = {[int(round(x)) for x in lambdav]}")
        


if __name__ == "__main__":
    main()
