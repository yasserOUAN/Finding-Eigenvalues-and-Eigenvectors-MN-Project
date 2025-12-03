#include <iostream>
#include <complex>
#include <cmath>
#include <iomanip>
#include <random>
using namespace std;

// Type aliases for cleaner code
typedef complex<double> Complex;
typedef double Matrix2x2[2][2];

// Function to create and fill a 2x2 matrix
void fillMatrix(Matrix2x2 A) {
    cout << "Fill a 2x2 matrix:\n\n";
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cout << "A[" << i << "][" << j << "] = ";
            cin >> A[i][j];
        }
    }
}

// Function to calculate eigenvalues
void eigenvalue(Matrix2x2 A, Complex& t1, Complex& t2) {
    // For 2x2 matrix: solve λ² - trace(λ) + det = 0
    double trace = A[0][0] + A[1][1];
    double det = (A[0][0] * A[1][1]) - (A[0][1] * A[1][0]);
    double discriminant = trace * trace - 4 * det;

    if (discriminant > 0) {
        // Two real eigenvalues
        t1 = -((-trace + sqrt(discriminant)) / 2.0);
        t2 = -((-trace - sqrt(discriminant)) / 2.0);
    } else if (abs(discriminant) < 1e-10) {
        // One repeated eigenvalue
        t1 = -trace / 2.0;
        t2 = t1;
    } else {
        // Complex eigenvalues
        cout << "Complex eigenvalues detected\n\n";
        Complex sqrtDisc = sqrt(Complex(discriminant, 0));
        t1 = -((-trace + sqrtDisc) / 2.0);
        t2 = -((-trace - sqrtDisc) / 2.0);
    }
}

// Function to calculate eigenvectors
void eigenvector(Matrix2x2 A, Complex t1, Complex t2,
                 Complex v1[2], Complex v2[2], bool& hasTwoVectors) {

    // Random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(1, 10);
    uniform_int_distribution<> sign(0, 1);

    double a = dis(gen) * (sign(gen) ? 1 : -1);
    double c = dis(gen) * (sign(gen) ? 1 : -1);

    // Calculate first eigenvector
    // Solve (A - λI)v = 0
    if (abs(A[0][1]) > 1e-10) {
        Complex x = (A[0][0] - t1) / A[0][1];
        if ((real(A[0][0] - t1) > 0 && A[0][1] > 0) ||
            (real(A[0][0] - t1) < 0 && A[0][1] < 0)) {
            x = -x;
        }
        v1[0] = a;
        v1[1] = Complex(a) * x;

        cout << "From the eigenvalues we deduce that the vector can be written as\n";
        cout << "-----> x = (" << x << ") * y\n";
        cout << "Let x = (" << a << ") so y = " << a * x << "\n\n";
    } else if (abs(A[1][0]) > 1e-10) {
        // Use second row if first row has zero
        Complex x = -(A[1][1] - t1) / A[1][0];
        v1[0] = Complex(a) * x;
        v1[1] = a;
    } else {
        v1[0] = 1;
        v1[1] = 0;
    }

    // Calculate second eigenvector if eigenvalues are different
    if (abs(t1 - t2) > 1e-10) {
        hasTwoVectors = true;

        if (abs(A[1][0]) > 1e-10) {
            Complex y = A[1][0] / (A[1][1] - t2);
            if ((real(A[1][1] - t2) > 0 && A[1][0] > 0) ||
                (real(A[1][1] - t2) < 0 && A[1][0] < 0)) {
                y = -y;
            }
            v2[0] = c;
            v2[1] = Complex(c) * y;

            cout << "And the second vector can be written as\n";
            cout << "-----> X = (" << y << ") * Y\n";
            cout << "Let X = (" << c << ") so Y = " << c * y << "\n\n";
        } else {
            v2[0] = 0;
            v2[1] = 1;
        }
    } else {
        hasTwoVectors = false;
        v2[0] = 0;
        v2[1] = 0;
    }
}

// Function to print eigenvalues
void printValues(Complex t1, Complex t2) {
    cout << fixed << setprecision(4);

    if (abs(t1 - t2) < 1e-10) {
        cout << "The matrix has one eigenvalue:\n";
        cout << "λ1 = " << t1 << "\n";
    } else {
        cout << "The matrix has two eigenvalues:\n";
        cout << "λ1 = " << setprecision(2) << t1 << "\n";
        cout << "λ2 = " << setprecision(2) << t2 << "\n";
    }
    cout << "\n";
}

// Function to print eigenvectors
void printVectors(Complex t1, Complex t2, Complex v1[2], Complex v2[2],
                  bool hasTwoVectors) {
    cout << fixed << setprecision(2);

    if (!hasTwoVectors) {
        cout << "The matrix has one eigenvector:\n";
        cout << "[" << v1[0] << ", " << v1[1] << "]\n";
    } else {
        cout << "The matrix has two eigenvectors:\n";
        cout << "for λ1 = " << setprecision(2) << t1 << " -> [" << setprecision(2) << v1[0] << ", " << setprecision(2) << v1[1] << "]\n";
        cout << "for λ2 = " << setprecision(2) << t2 << " -> [" << setprecision(2) << v2[0] << ", " << setprecision(2) << v2[1] << "]\n";
    }
}

int main() {
    Matrix2x2 A = {{0.0, 0.0}, {0.0, 0.0}};
    Complex t1, t2;
    Complex v1[2], v2[2];
    bool hasTwoVectors;

    // Matrix creation
    fillMatrix(A);

    cout << "\nYour matrix:\n";
    cout << "[" << A[0][0] << ", " << A[0][1] << "]\n";
    cout << "[" << A[1][0] << ", " << A[1][1] << "]\n";

    // Finding and printing eigenvalues
    cout << "\n" << string(50, '=') << "\n";
    cout << "EIGENVALUES\n";
    cout << string(50, '=') << "\n";
    eigenvalue(A, t1, t2);
    printValues(t1, t2);

    // Finding and printing eigenvectors
    cout << "\n" << string(50, '=') << "\n";
    cout << "EIGENVECTORS\n";
    cout << string(50, '=') << "\n";
    eigenvector(A, t1, t2, v1, v2, hasTwoVectors);
    printVectors(t1, t2, v1, v2, hasTwoVectors);

    return 0;
}
