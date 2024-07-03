#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col) {
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }

    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}

Matrix sub_matrix(Matrix a, Matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}

Matrix mul_matrix(Matrix a, Matrix b) {
    if (a.cols != b.rows) {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}

Matrix scale_matrix(Matrix a, double k) {
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a) {
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

double det_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    int n = a.rows;

    if (n == 1) {
        return a.data[0][0];
    }

    if (n == 2) {
        return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];
    }

    double det = 0;
    for (int j = 0; j < n; j++) {
        Matrix sub_matrix = create_matrix(n - 1, n - 1);
        for (int row = 1; row < n; row++) {
            int sub_col = 0;
            for (int col = 0; col < n; col++) {
                if (col == j) continue;
                sub_matrix.data[row - 1][sub_col] = a.data[row][col];
                sub_col++;
            }
        }

        double sub_det = det_matrix(sub_matrix); // 递归

        det += (j % 2 == 0 ? 1 : -1) * a.data[0][j] * sub_det;
    }

    return det;
}

Matrix inv_matrix(Matrix a) {
    int n = a.rows;
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }

    double det = det_matrix(a);
    if (fabs(det) < 1e-10) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    Matrix adj = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Matrix sub_matrix = create_matrix(n - 1, n - 1);
            for (int row = 0; row < n; row++) {
                if (row == i) continue;
                int sub_row = row < i ? row : row - 1;
                for (int col = 0; col < n; col++) {
                    if (col == j) continue;
                    int sub_col = col < j ? col : col - 1;
                    sub_matrix.data[sub_row][sub_col] = a.data[row][col];
                }
            }
            adj.data[j][i] = (i + j) % 2 == 0 ? det_matrix(sub_matrix) : -det_matrix(sub_matrix);
        }
    }

    Matrix inv = scale_matrix(adj, 1.0 / det);
    return inv;
}

int rank_matrix(Matrix a) {
    int m = a.rows;
    int n = a.cols;
    int rank = n;

    for (int row = 0; row < rank; row++) {
        if (a.data[row][row]) {
            for (int col = 0; col < m; col++) {
                if (col != row) {
                    double mult = a.data[col][row] / a.data[row][row];
                    for (int i = 0; i < rank; i++) {
                        a.data[col][i] -= mult * a.data[row][i];
                    }
                }
            }
        } else {
            int reduce = 1;
            for (int i = row + 1; i < m; i++) {
                if (a.data[i][row]) {
                    for (int j = 0; j < n; j++) {
                        double temp = a.data[row][j];
                        a.data[row][j] = a.data[i][j];
                        a.data[i][j] = temp;
                    }
                    reduce = 0;
                    break;
                }
            }

            if (reduce) {
                rank--;
                for (int i = 0; i < m; i++) {
                    a.data[i][row] = a.data[i][rank];
                }
            }
            row--;
        }
    }
    return rank;
}

double trace_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    double trace = 0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a) {
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}