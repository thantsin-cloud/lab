#include <mpi.h>
#include <math.h>
#include <stdio.h>

float fct(float x)
{
    return cos(x);
}

float integral(float a, int n, float h);

void main(argc, argv)
int argc;
char *argv[];
{
    int n, p, i, j, ierr, num;
    float h, result, a, b, pi;
    float my_a, my_range;

    int myid, source, dest, tag;
    MPI_Status status;
    float my_result;

    pi = acos(-1.0);
    a = 0.;
    b = pi * 1. / 2.;
    n = 100000;

    dest = 0;
    tag = 123;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    h = (b - a) / n;
    num = n / p;
    my_range = (b - a) / p;
    my_a = a + myid * my_range;
    my_result = integral(my_a, num, h);

    printf("Process %d has the partial result of %f\n", myid, my_result);

    MPI_Reduce(&my_result, &result, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Finalize();
}

float integral(float a, int n, float h)
{
    int j;
    float h2, aij, integ;

    integ = 0.0;
    h2 = h / 2.;
    for (j = 0; j < n; j++) {
        aij = a + j * h;
        integ += fct(aij + h2) * h;
    }
    return (integ);
}
