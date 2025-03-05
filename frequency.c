#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ARRAY_SIZE 1000000

int main(int argc, char** argv)
{
    int num = atoi(argv[1]);
    static int list_of_numbers[ARRAY_SIZE];
    int i;
    time_t t;

    srand((unsigned) time(&t));

    for( i = 0 ; i < ARRAY_SIZE ; ++i )
        list_of_numbers[i] = rand() % 100;

    MPI_Status status;

    MPI_Init(NULL, NULL);

    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    int number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    if (pid == 0) {
        int index, i;
        int elements_per_process;
        unsigned long frequency;

        elements_per_process = ARRAY_SIZE / number_of_processes;

        if (number_of_processes > 1) {
            for (i = 1; i < number_of_processes - 1; i++) {
                index = i * elements_per_process;

                MPI_Send(&elements_per_process,
                         1, MPI_INT, i, 0,
                         MPI_COMM_WORLD);
                MPI_Send(&list_of_numbers[index],
                         elements_per_process,
                         MPI_INT, i, 0,
                         MPI_COMM_WORLD);
            }

            index = i * elements_per_process;
            int elements_left = ARRAY_SIZE - index;

            MPI_Send(&elements_left,
                     1, MPI_INT,
                     i, 0,
                     MPI_COMM_WORLD);
            MPI_Send(&list_of_numbers[index],
                     elements_left,
                     MPI_INT, i, 0,
                     MPI_COMM_WORLD);
        }

        frequency = 0;
        for(i = 0; i < elements_per_process; ++i)
            if(list_of_numbers[i] == num)
                frequency += 1;

        unsigned long buffer = 0;
        for (i = 1; i < number_of_processes; i++) {
            MPI_Recv(&buffer, 1, MPI_INT,
                     MPI_ANY_SOURCE, 0,
                     MPI_COMM_WORLD,
                     &status);
            frequency += buffer;
        }

        printf("The frequency of %d in the list of numbers is %ld\n", num, frequency);
    }
    else {
        static int buffer[ARRAY_SIZE];
        int num_of_elements_recieved = 0;
        unsigned long frequency = 0;

        MPI_Recv(&num_of_elements_recieved,
                 1, MPI_INT, 0, 0,
                 MPI_COMM_WORLD,
                 &status);

        MPI_Recv(&buffer, num_of_elements_recieved,
                 MPI_INT, 0, 0,
                 MPI_COMM_WORLD,
                 &status);

        for(i = 0; i < num_of_elements_recieved; ++i)
            if(buffer[i] == num)
                frequency += 1;

        MPI_Send(&frequency, 1, MPI_INT,
                 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
