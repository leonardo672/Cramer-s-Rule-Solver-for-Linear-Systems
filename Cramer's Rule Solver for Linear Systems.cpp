#include <iostream>
#include <vector>

using namespace std;

// Function to calculate the determinant of a 3x3 matrix
double calculateDeterminant(const vector<vector<double>>& matrix) {
    double determinant =
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
        matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
        matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

    return determinant;
}

// Function to solve a system of linear equations using Cramer's method
vector<double> solveByCramer(const vector<vector<double>>& coefficientMatrix, const vector<double>& constantVector) {
    int n = constantVector.size();

    // Calculate the determinant of the coefficient matrix
    double determinant = calculateDeterminant(coefficientMatrix);

    vector<double> solution;
    if (determinant == 0.0) {
        cout << "The system of equations has no unique solution." << endl;
        return solution;
    }

    for (int i = 0; i < n; ++i) {
        vector<vector<double>> tempMatrix = coefficientMatrix;  // Create a copy of the coefficient matrix

        // Replace the i-th column of the tempMatrix with the constant vector
        for (int j = 0; j < n; ++j) {
            tempMatrix[j][i] = constantVector[j];
        }

        // Calculate the solution for the i-th variable
        double variableSolution = calculateDeterminant(tempMatrix) / determinant;
        solution.push_back(variableSolution);
    }

    return solution;
}

int main() {
    int n;

    // Get the size of the system of equations
    cout << "Enter the number of equations: ";
    cin >> n;

    // Coefficient matrix
    vector<vector<double>> coefficientMatrix(n, vector<double>(n));

    // Constant vector
    vector<double> constantVector(n);

    // Get the coefficient matrix
    cout << "Enter the coefficient matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "Equation " << i + 1 << ": ";
        for (int j = 0; j < n; ++j) {
            cin >> coefficientMatrix[i][j];
        }
    }

    // Get the constant vector
    cout << "Enter the constant vector:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "b" << i + 1 << ": ";
        cin >> constantVector[i];
    }

    // Solve the system of linear equations using Cramer's method
    vector<double> solution = solveByCramer(coefficientMatrix, constantVector);

    int numVariables = solution.size();
    if (numVariables > 0) {
        cout << "Solution: ";
        for (int i = 0; i < numVariables; ++i) {
            cout << (char)('x' + i) << " = " << solution[i] << " ";
        }
        cout << endl;
    }

    return 0;
}
