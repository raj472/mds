#include <iostream>
#include <cmath>

using namespace std;

const int MAX_SIZE = 10;

// Function to perform Gaussian elimination
void gaussianElimination(double mat[MAX_SIZE][MAX_SIZE + 1], int n) {
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(mat[j][i]) > fabs(mat[maxRow][i]))
                maxRow = j;
        }
        for (int k = i; k <= n; ++k) {
            double temp = mat[i][k];
            mat[i][k] = mat[maxRow][k];
            mat[maxRow][k] = temp;
        }
        for (int j = i + 1; j < n; ++j) {
            double factor = mat[j][i] / mat[i][i];
            for (int k = i; k <= n; ++k) {
                mat[j][k] -= factor * mat[i][k];
            }
        }
    }
}

// Function to back substitute to find eigenvalues
void backSubstitute(double mat[MAX_SIZE][MAX_SIZE + 1], int n, double eigenvalues[]) {
    for (int i = n - 1; i >= 0; --i) {
        eigenvalues[i] = mat[i][n] / mat[i][i];
        for (int k = i - 1; k >= 0; --k) {
            mat[k][n] -= mat[k][i] * eigenvalues[i];
        }
    }
}

// Function to find eigenvectors
void findEigenvectors(double mat[MAX_SIZE][MAX_SIZE + 1], int n, double eigenvalues[], double eigenvectors[MAX_SIZE][MAX_SIZE]) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            eigenvectors[j][i] = mat[j][i] / (mat[i][i] * eigenvalues[i]);
        }
    }
}

int main() {
    int n;
    double mat[MAX_SIZE][MAX_SIZE + 1];
    double eigenvalues[MAX_SIZE];
    double eigenvectors[MAX_SIZE][MAX_SIZE];

    cout << "Enter the size of square matrix: ";
    cin >> n;

    cout << "Enter the elements of the matrix row-wise: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> mat[i][j];
        }
    }

    // Augmenting matrix with identity matrix
    for (int i = 0; i < n; ++i) {
        for (int j = n; j < n + 1; ++j) {
            mat[i][j] = (i == j - n) ? 1.0 : 0.0;
        }
    }

    gaussianElimination(mat, n);
    backSubstitute(mat, n, eigenvalues);

    cout << "Eigenvalues: ";
    for (int i = 0; i < n; ++i) {
        cout << eigenvalues[i] << " ";
    }
    cout << endl;

    findEigenvectors(mat, n, eigenvalues, eigenvectors);
    cout << "Eigenvectors: " << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << eigenvectors[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
