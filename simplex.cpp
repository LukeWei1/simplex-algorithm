// Simplest tableau simplex: maximize c.x subject to A x <= b, x >= 0, b >= 0.
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

typedef vector<vector<double>> Matrix;

// Fills x with the optimal values and returns z.
double simplex(vector<double> c, Matrix A, vector<double> b, vector<double>& x) {
    int m = A.size(), n = c.size();

    // Build tableau: [A | I | b], bottom row [-c | 0 | 0]
    Matrix T(m + 1, vector<double>(n + m + 1, 0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) T[i][j] = A[i][j];
        T[i][n + i] = 1;            // slack variable
        T[i][n + m] = b[i];         // RHS
    }
    for (int j = 0; j < n; j++) T[m][j] = -c[j];

    vector<int> basis(m);           // which variable labels each row
    for (int i = 0; i < m; i++) basis[i] = n + i;

    while (true) {
        // Entering variable: most negative entry in bottom row
        int col = -1;
        double most = -1e-9;
        for (int j = 0; j < n + m; j++)
            if (T[m][j] < most) { most = T[m][j]; col = j; }
        if (col == -1) break;       // no negatives: optimal

        // Ratio test: smallest RHS / entry, over positive entries
        int row = -1;
        double best = 1e18;
        for (int i = 0; i < m; i++)
            if (T[i][col] > 1e-9 && T[i][n + m] / T[i][col] < best) {
                best = T[i][n + m] / T[i][col];
                row = i;
            }
        if (row == -1) { cout << "unbounded\n"; return 1e18; }

        // Pivot: scale pivot row, then zero out the column in other rows
        double p = T[row][col];
        for (double& v : T[row]) v /= p;
        for (int i = 0; i <= m; i++) {
            if (i == row) continue;
            double f = T[i][col];
            for (int j = 0; j <= n + m; j++) T[i][j] -= f * T[row][j];
        }
        basis[row] = col;
    }

    // Basic variables = RHS, nonbasic = 0
    x.assign(n, 0);
    for (int i = 0; i < m; i++)
        if (basis[i] < n) x[basis[i]] = T[i][n + m];
    return T[m][n + m];
}

bool test(const char* name, vector<double> c, Matrix A, vector<double> b,
          vector<double> expectX, double expectZ) {
    vector<double> x;
    double z = simplex(c, A, b, x);
    bool ok = fabs(z - expectZ) < 1e-6;
    for (size_t i = 0; i < x.size(); i++) ok = ok && fabs(x[i] - expectX[i]) < 1e-6;
    cout << name << ": x=" << x[0] << " y=" << x[1] << " z=" << z
         << (ok ? "  PASS" : "  FAIL") << "\n";
    return ok;
}

int main() {
    // Example 1: max 3x + 2y, x + y <= 4, 2x + y <= 6  -> (2,2), z=10
    test("Example 1", {3, 2}, {{1, 1}, {2, 1}}, {4, 6}, {2, 2}, 10);
    // Example 2: max 2x + y, x + y <= 4, x + 3y <= 6   -> (4,0), z=8
    test("Example 2", {2, 1}, {{1, 1}, {1, 3}}, {4, 6}, {4, 0}, 8);
}