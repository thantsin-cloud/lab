#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

int compare_dbls(const void* arg1, const void* arg2){
    double a1 = *(double *) arg1;
    double a2 = *(double *) arg2;
    if (a1 < a2) return -1;
    else if (a1 == a2) return 0;
    else return 1;
}

void qsort_dbls(double *array, int array_len){
    qsort(array, (size_t)array_len, sizeof(double), compare_dbls);
}

int find_bucket(double num, int p_num){
    int x;
    for(x=1; x < p_num+1; x++){
        double bucket_range =(double) x / (double)p_num;
        if(num <= bucket_range){
            return x - 1;
        }
    }
}

int main(int argc, char *argv[]){
    int myrank, P;
    int sub_count[1];

    if(argc != 2){
        printf("\nPlease include N, problem size\n");
        return 0;
    }

    int N = strtol(argv[1], NULL, 10);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    double t1 = MPI_Wtime();
    
    int numpp = N/P;
    double *list = (double*)malloc(numpp*sizeof(double));
    int *count = (int*)calloc(P, sizeof(int));
    int i, bucket;
    double r;
    for(i = 0; i < numpp; i++){

        r = (double)rand() / (double) RAND_MAX;

        list[i] = r;
        bucket = find_bucket(r, P);
        count[bucket]++;
    }

    int *bucket_count = (int*)malloc(P*sizeof(int));
    MPI_Alltoall(count, 1, MPI_INT, bucket_count, 1, MPI_INT, MPI_COMM_WORLD);

    int loc_bcount = 0;
    for(i = 0; i < P; i++){
        loc_bcount+= bucket_count[i]; 
    }

    double *bucket_list = (double*)malloc(loc_bcount*sizeof(double));

    int *displs = (int*)malloc(P*sizeof(int));
    double *dist_list = (double*)malloc(numpp*sizeof(double));
    int *index = (int*)calloc(P,sizeof(int));

    displs[0] = 0;
    for(i = 1; i < P; i++){
        displs[i] = count[i-1] + displs[i-1];
    }

    int *rdispls = (int*)malloc(P*sizeof(int));
    rdispls[0] = 0;
    for(i = 1; i < P; i++){
        rdispls[i] = bucket_count[i-1] + rdispls[i-1];
    }

    for(i = 0; i < numpp; i++){
        bucket = find_bucket(list[i], P);
        dist_list[displs[bucket] + index[bucket]] = list[i];
        index[bucket]++;
    }
    free(list);

    MPI_Alltoallv(dist_list, count, displs, MPI_DOUBLE, bucket_list, bucket_count, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);    

    qsort_dbls(bucket_list, loc_bcount);

    int gathercounts[1];
    gathercounts[0] = loc_bcount;
    MPI_Gather(gathercounts, 1, MPI_INT, count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(myrank==0){
        displs[0] = 0;
        for(i = 1; i < P; i++){
            displs[i] = count[i-1] + displs[i-1];
        }
    }
    
    double* final_list = (double*)malloc(N*sizeof(double));
    MPI_Gatherv(bucket_list,loc_bcount, MPI_DOUBLE, final_list, count, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    double t2 = MPI_Wtime();

    if(myrank == 0){
        int sorted = 1;
        int k;
        for(k = 0; k < N - 2; k++){
            if(final_list[k] > final_list[k+1]){
                sorted = 0;
            }
        }

        if(sorted == 1){
            printf("\nSORTING CORRECT\n");
        }else{
            printf("\nSORTING NOT CORRECT\n");
        }

        printf("\nV2- N: %d P: %d  Execution Time: %f\n",N, P, t2-t1);
    }

    free(index);
    free(displs);
    free(rdispls);
    free(count);
    free(bucket_count);
    free(bucket_list);
    free(final_list);
    MPI_Finalize();
}
