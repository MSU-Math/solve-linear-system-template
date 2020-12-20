#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define DECIMAL 10
#define GIGA_MODIFIER 1e9
#define NANO_MODIFIER 1e-9
#define FIVE 5

double f(int k, int n, int i, int j)
{
    switch (k) {
    case 1: {
        double res;
        double max;
        if (i > j) {
            max = i;
        } else {
            max = j;
        }
        res = n - max + 1;
        return res;
    }

    case 2: {
        if (i > j) {
            return i;
        }
        return j;
    }

    case 3: {
        return abs(i - j);
    }

    case 4: {
        return (1.0 / (i + j - 1.0));
    }
    }
    return 0;
}

struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec - start.tv_nsec) < 0) {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = GIGA_MODIFIER + end.tv_nsec - start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}



int slau(int n) {
	int m, i, j, k;
	m = n + 1;     //расширение матрицы
	FILE *f = fopen("mtr.txt", "r");
	double **A = NULL, **U = NULL, **uu = NULL, *d = NULL, *y = NULL, *x = NULL, sum = 0, deter = 0;
	A = (double**)malloc(m * sizeof(double*));
	U = (double**)malloc(m * sizeof(double*));
	uu = (double**)malloc(m * sizeof(double*));
	d = (double*)malloc(m * sizeof(double*)); //вектор при матрице
	y = (double*)malloc(m * sizeof(double*)); //промежуточный вектор
	x = (double*)malloc(m * sizeof(double*)); //ответ
	for (i = 0; i < n; i++)
	{
		A[i] = (double*)malloc(m * sizeof(double)); //под столбцы
		U[i] = (double*)malloc(m * sizeof(double)); //под столбцы
		uu[i] = (double*)malloc(m * sizeof(double)); //под столбцы
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n + 1; j++)
		{
			fscanf(f, "%lf", &A[i][j]);
		}
	}
	fclose(f);
	for (i = 0; i < n; i++)
	{
		for (j = n; j < n + 1; j++)
		{
			d[i] = A[i][j];
		}
	}
//   определитель для ошибки минус 4
	for (k = 0; k <  - 1; k++)
	{
		for (i = k + 1; i < n; i++)
		{
			double tmp = A[i][k] / A[k][k];
			A[i][k] = 0;
			for (j = k + 1; j < n; j++)
			{
				A[i][j] = A[i][j] - tmp * A[k][j];
			}
		}
	}
	for (i = 0; i < n; i++)
	{
		deter *= A[i][i];
	}
	if (deter <= 0)
	{
		return -4;
	}
//

				// вывод матрицы и значений
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			if(j == n)
			{
				printf(" | ");
			}
			printf("%lf ", A[i][j]);
		}
		printf("\n");
	}

	for(i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (i < j)
			{
				U[i][j] = 0.0;
			}
			if (i == j)
			{
				for (k = 0; k < i; k++)
				{
					sum = sum + U[i][k] * U[i][k];
				}
				U[i][j] = sqrt(A[i][i] - sum);
				sum = 0;
			}
			if (i > j)
			{
				for (k = 0; k < j; k++)
				{
					sum = sum + U[i][k] * U[j][k];
				}
				U[i][j] = (A[i][j] - sum) / U[j][j];
				sum = 0;
			}
		}
	}
/*	printf("Нижнетреугольная матрица:\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%lf ", U[i][j]);
			uu[i][j] = U[j][i];
		}
		printf("\n");
	}
	printf("Верхнетреугольная матрица:\n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%lf ", uu[i][j]);
		}
		printf("\n");
	}
*/
	y[0] = d[0] / U[0][0];
	for (i = 1; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < i; k++)
			{
				sum = sum + U[i][k] * y[k];
			}
			y[i] = (d[i] - sum) / U[i][i];
			sum = 0;
		}
	}
/*	printf("Получаем промежуточный вектор\n");
	for (i = 0; i < n; i++)
	{
		printf("%lf ", y[i]);
	}
	printf("\n");
*/

	x[n - 1] = y[n - 1] / uu[n - 1][n - 1];
	for (i = n - 2; i >= 0; i--)
	{
		for (j = n - 1; j > 0; j--)
		{
			for (k = n - 1; k > i; k--)
			{
				sum = sum + uu[i][k] * x[k];
			}
			x[i] = (y[i] - sum) / uu[i][i];
			sum = 0;
		}
	}
	printf("Solution\n");
	for (i = 0; i < n; i++)
	{
		printf("%lf ", x[i]);
	}
	return 0;
}

int main(int argc, char *argv[])
{
    int n, m, k, check = 0;
    char *endptr;
    struct timespec time_start;
    struct timespec time_end;
    if (argc != 4 && argc != FIVE) {
        return -1;
    }
    n = (int)strtol(argv[1], &endptr, DECIMAL);
    m = (int)strtol(argv[2], &endptr, DECIMAL);
    k = (int)strtol(argv[3], &endptr, DECIMAL);
    if (endptr == argv[1] || endptr == argv[2] || endptr == argv[3]) {
        return -1;
    }
    if ((n <= 0) || (m <= 0) || (m > n) || (k < 0) || (k > 4)) {
        return -1;
    }
    if ((k == 0 && argc == 4) || (argc == FIVE && k != 0)) {
        return -1;
    }
    double *matr;
    matr = (double *)malloc(n * n * sizeof(double));
    if (matr == NULL) {
        return -2;
    }
    if (argc == FIVE) {
        FILE *f;
        f = fopen(argv[4], "r");
        if (f == NULL) {
            free(matr);
            return -3;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int scan = fscanf(f, "%lf", &matr[i * n + j]);
                if (scan == 0 || scan == EOF) {
                    free(matr);
                    fclose(f);
                    return -3;
                }
            }
        }
        fclose(f);
    } else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matr[i * n + j] = f(k, n, i + 1, j + 1);
            }
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &time_start);
    check = slau(n);
    if (check < 0)
    {
	return -4;
    }
    clock_gettime(CLOCK_MONOTONIC, &time_end);
    free(matr);
    time_end = diff(time_start, time_end);
    printf("Time: %lf s\n", (double)(time_end.tv_sec + time_end.tv_nsec * NANO_MODIFIER));
    return 0;
}
