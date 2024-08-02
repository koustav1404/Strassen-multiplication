// Problem Statement : Multiplication of two matrices using recursive versions of the Strassen's algorithm

#include <stdio.h>
#define MAX 20 // maximum size of the matrix

// function to add two matrices of size s
void add(int a1[MAX][MAX], int a2[MAX][MAX], int mult[MAX][MAX], int s)
{
    int i, j;
    for (i = 0; i < s; ++i)
    {
        for (j = 0; j < s; ++j)
            mult[i][j] = a1[i][j] + a2[i][j];
    }
}

// function to subtract two matrices of size s
void sub(int a1[MAX][MAX], int a2[MAX][MAX], int mult[MAX][MAX], int s)
{
    int i, j;
    for (i = 0; i < s; ++i)
    {
        for (j = 0; j < s; ++j)
            mult[i][j] = a1[i][j] - a2[i][j];
    }
}

// recursive Strassen's multiplication function
void mul(int a1[MAX][MAX], int a2[MAX][MAX], int res[MAX][MAX], int s1)
{

    // if size is less the equal to 2 then we will implent the Strassen's algorithm of individual elements

    if (s1 <= 2) // if size is less than equal to two
    {
        // Strassen's multiplication formulas
        int p, q, r, s, t, u, v;
        p = (a1[0][0] + a1[1][1]) * (a2[0][0] + a2[1][1]);
        q = (a1[1][0] + a1[1][1]) * a2[0][0];
        r = a1[0][0] * (a2[0][1] - a2[1][1]);
        s = a1[1][1] * (a2[1][0] - a2[0][0]);
        t = (a1[0][0] + a1[0][1]) * a2[1][1];
        u = (a1[1][0] - a1[0][0]) * (a2[0][0] + a2[0][1]);
        v = (a1[0][1] - a1[1][1]) * (a2[1][0] + a2[1][1]);

        // results
        res[0][0] = p + s - t + v;
        res[0][1] = r + t;
        res[1][0] = q + s;
        res[1][1] = p + r - q + u;

        return;
    }

    // if size is greater the to 2 then we will implent the Strassen's algorithm of individual submatrices

    int i, j; // loop variables

    int s2 = s1 / 2; // size of matrix after dividing

    // declaration of submatrices

    // sub matrices of first matrix
    int a11[MAX][MAX], a12[MAX][MAX], a21[MAX][MAX], a22[MAX][MAX];

    // sub matrices of second matrix
    int b11[MAX][MAX], b12[MAX][MAX], b21[MAX][MAX], b22[MAX][MAX];

    // sub matrices of resultant matrix
    int c11[MAX][MAX], c12[MAX][MAX], c21[MAX][MAX], c22[MAX][MAX];

    // matrices to store Strassen's algorithm variables
    int p[MAX][MAX], q[MAX][MAX], r[MAX][MAX], s[MAX][MAX], t[MAX][MAX], u[MAX][MAX], v[MAX][MAX];

    for (i = 0; i < s2; i++)
    {
        for (j = 0; j < s2; j++)
        {

            // Splitting of first matrix

            a11[i][j] = a1[i][j];
            a12[i][j] = a1[i][j + s2];
            a21[i][j] = a1[i + s2][j];
            a22[i][j] = a1[i + s2][j + s2];

            // Splitting of second matrix

            b11[i][j] = a2[i][j];
            b12[i][j] = a2[i][j + s2];
            b21[i][j] = a2[i + s2][j];
            b22[i][j] = a2[i + s2][j + s2];
        }
    }
    int temp1[MAX][MAX], temp2[MAX][MAX]; // temporary matrix to store data

    // Strassen algorithmn implementation on individual submatrices

    // p = (a11 +a22) * (b11 + b22)
    add(a11, a22, temp1, s2);
    add(b11, b22, temp2, s2);
    mul(temp1, temp2, p, s2);

    // q = (a21 +a22) * b11
    add(a21, a22, temp1, s2);
    mul(temp1, b11, q, s2);

    // r = a11 * (b12 - b22)
    sub(b12, b22, temp1, s2);
    mul(a11, temp1, r, s2);

    // s = a22 * (b21 - b11)
    sub(b21, b11, temp1, s2);
    mul(a22, temp1, s, s2);

    // t = (a11 + a12) * b22
    add(a11, a12, temp1, s2);
    mul(temp1, b22, t, s2);

    // u = (a21 - a11) * (b11 + b12)
    sub(a21, a11, temp1, s2);
    add(b11, b12, temp2, s2);
    mul(temp1, temp2, u, s2);

    // v = (a12 - a22) * (b21 + b22)
    sub(a12, a22, temp1, s2);
    add(b21, b22, temp2, s2);
    mul(temp1, temp2, v, s2);

    // calcuations of the resultant submatrices from Strassen's formulas

    // c11 = p + s - t + v
    add(p, s, temp1, s2);
    sub(temp1, t, temp2, s2);
    add(temp2, v, c11, s2);

    //  c12 = r + t
    add(r, t, c12, s2);

    // c21 = q + s
    add(q, s, c21, s2);

    // c22 = p + r - q + u
    add(p, r, temp1, s2);
    sub(temp1, q, temp2, s2);
    add(temp2, u, c22, s2);

    // recombining the results
    for (i = 0; i < s2; i++)
    {
        for (j = 0; j < s2; j++)
        {

            res[i][j] = c11[i][j];
            res[i][j + s2] = c12[i][j];
            res[i + s2][j] = c21[i][j];
            res[i + s2][j + s2] = c22[i][j];
        }
    }
    return;
}
int main()
{
    int a1[MAX][MAX], a2[MAX][MAX], ans[MAX][MAX]; // initilization of matrix 1, matrix 2 ans result matrix
    int s, i, j;                                   // s id size of matrix and i&j are loop variables

    // input of matrix size
    do
    {
        printf("Enter the dimension of the matrix (must be a power of 2) : ");
        scanf("%d", &s);
    } while (s % 2 != 0);

    // Input of first  matrix
    printf("\nEnter elements of matrix 1:\n");
    for (i = 0; i < s; ++i)
    {
        for (j = 0; j < s; ++j)
        {
            scanf("%d", &a1[i][j]);
        }
    }

    // input of second matrix
    printf("\nEnter elements of matrix 2:\n");
    for (i = 0; i < s; ++i)
    {
        for (j = 0; j < s; ++j)
        {
            scanf("%d", &a2[i][j]);
        }
    }
    mul(a1, a2, ans, s); // Strassen's algorithmn calling

    // displaying the result
    printf("\nOutput Matrix:\n");
    for (i = 0; i < s; ++i)
    {
        for (j = 0; j < s; ++j)

            printf("%d  ", ans[i][j]);
        printf("\n");
    }

    return 0;
}