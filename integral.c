#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<mpi.h>
#include<math.h>

#define INTERVAL_MIN 0.0
#define INTERVAL_MAX 1.0
#define NUM_INTERVALS 1
#define ERROR_THRESHOLD 0.00001

typedef struct
{
    int intervals,rank,size;
    double delta,region;
    double error;
    double min,max;
    double (*func)(double start, double end, double delta);
} Params;

double f(double x)
{
    return sqrt(1-x*x); 
}

double rectangle_rule(double start, double end, double delta)
{
    double area = 0.0;
    double x;
    for(x = start; x <= end; x+= delta)
        area += f(x)*delta;
    return area;
}

double trapezoid_rule(double start, double end, double delta)
{
    double area = 0.0;
    area = 0.5 * ( f(start) + f(end) );
    double x  = 0.0;
    for(x = start + delta; x <= end; x+= delta)
    {
        area += f(x);
    }
    area = area * delta;
    return area;
}

void master(Params * p)
{   
    double last = 0.0;
    double area = p->error+1,local;
    double start, end;
    while(fabs(last - area) > p->error)
    {
        last = area;
        start = MPI_Wtime();
        MPI_Bcast(&p->intervals,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Reduce(&local,&area,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        end = MPI_Wtime();
        p->delta = (p->max-p->min)/p->intervals;
        fprintf(stderr,"Intervals = %8d\tDelta = %f\tResult = %f\tError=%f\tTime = %f\n",p->intervals,p->delta,area,fabs(last-area),end-start);
        p->intervals = p->intervals*2;
    } 
    p->intervals = -1;
    MPI_Bcast(&p->intervals,1,MPI_INT,0,MPI_COMM_WORLD);
    return;
}

void slave(Params *p)
{
    double start, end, area;
    while(1)
    {
        MPI_Bcast(&p->intervals,1,MPI_INT,0,MPI_COMM_WORLD);
        if(p->intervals == -1)
            return;
        p->region = (p->max-p->min)/(double)p->size;
        start = p->min + p->region * (p->rank-1);
        end = start + p->region;
        p->delta = (p->max-p->min)/p->intervals;
        area = p->func(start,end,p->delta);
        MPI_Reduce(&area,&area,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
}

void usage()
{
    printf("integration\n\
            MPI Program to approximate an integral.\n\
            Jharrod LaFon 2011\n\
            Usage: integration [args]\n\
            a <start>\t\tStarting point of interval to be integrated\n\
            b <end>\t\tEnding point of interval to be integrated\n\
            e <error>\t\tMaximum error threshold\n\
            t\t\tUse the trapezoid rule\n\
            r\t\tUse the rectangle rule\n\
            h\t\tPrint this message\n");
}

void parse_args(int argc, char ** argv, Params * p)
{
    int c = 0;
    while((c = getopt(argc,argv,"a:b:e:trh")) != -1)
    {
        switch(c)
        {
            case 'a':
                p->min = atof(optarg);
                break;
            case 'b':
                p->max = atof(optarg);
                break;
            case 'e':
                p->error = atof(optarg);
                break;
            case 't':
                p->func = trapezoid_rule;
                break;
            case 'r':
                p->func = rectangle_rule;
                break;
            case 'h':
                if(p->rank == 0) usage();
                exit(0);
            default:
                break;
        }
    }
    return;
}

int main(int argc, char *argv[])
{
    int rank;
    int size;
    if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        fprintf(stderr, "Unable to initialize MPI!\n");
        return -1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    Params p;
    p.min = INTERVAL_MIN;
    p.max = INTERVAL_MAX;
    p.error = ERROR_THRESHOLD;
    p.intervals = NUM_INTERVALS;
    p.rank = rank;
    p.size = size-1;
    p.func = rectangle_rule;
    
    parse_args(argc,argv,&p);
    
    if(rank == 0)
        master(&p);
    else
        slave(&p);
    MPI_Finalize();
    return 0;
}
