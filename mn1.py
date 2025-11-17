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
    Uses row reduction to find null space
    Returns: eigenvector as a list
    """
    n = len(matrix)
    
    # Create (A - λI)
    A_shifted = [[matrix[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        A_shifted[i][i] -= eigenvalue
    
    # For 2x2 matrix, use direct formula
    if n == 2:
        # From equation: a*v1 + b*v2 = 0
        # If b != 0: v1 = 1, v2 = -a/b
        # If b == 0: v1 = 0, v2 = 1
        a = A_shifted[0][0]
        b = A_shifted[0][1]
        
        if abs(b) > 1e-10:
            v = [1.0, -a/b]
        elif abs(a) > 1e-10:
            v = [-b/a, 1.0]
        else:
            # Use second row
            c = A_shifted[1][0]
            d = A_shifted[1][1]
            if abs(d) > 1e-10:
                v = [1.0, -c/d]
            else:
                v = [0.0, 1.0]
        
        return v
    
    # For larger matrices, use general approach
    # Find row with non-zero element
    v = [0.0] * n
    
    # Try each row to find valid null space vector
    for row_idx in range(n):
        # Find largest non-zero coefficient in this row
        max_coef = 0
        max_col = 0
        for col_idx in range(n):
            if abs(A_shifted[row_idx][col_idx]) > max_coef:
                max_coef = abs(A_shifted[row_idx][col_idx])
                max_col = col_idx
        
        if max_coef > 1e-10:
            # Set all variables to 1 except the one with max coefficient
            for i in range(n):
                v[i] = 1.0
            
            # Solve for the variable with max coefficient
            sum_val = 0
            for j in range(n):
                if j != max_col:
                    sum_val += A_shifted[row_idx][j] * v[j]
            
            v[max_col] = -sum_val / A_shifted[row_idx][max_col]
            break
    
    return v


def calculate_largest_eigenvalue(matrix, max_iter=100):
    """
    Calculate largest eigenvalue using power method
    Returns: (eigenvalue, eigenvector)
    """
    n = len(matrix)
    v = [1.0] * n
    
    for iteration in range(max_iter):
        v_new = matrix_vector_mult(matrix, v)
        
        max_val = max(abs(x) for x in v_new)
        if max_val > 1e-10:
            v = [x / max_val for x in v_new]
        else:
            v = v_new
    
    Av = matrix_vector_mult(matrix, v)
    eigenvalue = sum(v[i] * Av[i] for i in range(n)) / sum(v[i]**2 for i in range(n))
    
    return eigenvalue, v


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
        print(f"λ{i+1} = {val:.6f}")
    
    # Calculate eigenvectors
    print("\n" + "="*50)
    print("EIGENVECTORS")
    print("="*50)
    
    for i, eigenval in enumerate(eigenvalues):
        print(f"\nEigenvector {i+1} (for λ = {eigenval:.6f}):")
        eigenvec = calculate_eigenvector(matrix, eigenval)
        print_vector(eigenvec)
        
        # Verify A*v = λ*v
        Av = matrix_vector_mult(matrix, eigenvec)
        lambdav = [eigenval * v for v in eigenvec]
        
        print(f"  A*v = {[f'{x:.4f}' for x in Av]}")
        print(f"  λ*v = {[f'{x:.4f}' for x in lambdav]}")
        


if __name__ == "__main__":
    main()
