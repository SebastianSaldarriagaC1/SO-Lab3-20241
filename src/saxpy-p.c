/**
 * @defgroup   SAXPY saxpy
 *
 * @brief      This file implements an iterative saxpy operation
 * 
 * @param[in] <-p> {vector size} 
 * @param[in] <-s> {seed}
 * @param[in] <-n> {number of threads to create} 
 * @param[in] <-i> {maximum itertions} 
 *
 * @author     Danny Munera
 * @date       2020
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <pthread.h>

struct saxpy_args {
    double *X; 
    double *Y;
    double *Y_avgs;
    double a;
    int p;
    int max_iters;
    int start_pos;
    int end_pos;
};

void saxpy_setup(double*, double*, double*, double , int , int, int);

void* saxpy_thread(void*);

int main(int argc, char* argv[]){
	// Variables to obtain command line parameters
	unsigned int seed = 1;
  	int p = 10000000;
  	int n_threads = 2;
  	int max_iters = 1000;
  	// Variables to perform SAXPY operation
	double* X;
	double a;
	double* Y;
	double* Y_avgs;
	int i, it;
	// Variables to get execution time
	struct timeval t_start, t_end;
	double exec_time;

	// Getting input values
	int opt;
	while((opt = getopt(argc, argv, ":p:s:n:i:")) != -1){  
		switch(opt){  
			case 'p':  
			printf("vector size: %s\n", optarg);
			p = strtol(optarg, NULL, 10);
			assert(p > 0 && p <= 2147483647);
			break;  
			case 's':  
			printf("seed: %s\n", optarg);
			seed = strtol(optarg, NULL, 10);
			break;
			case 'n':  
			printf("threads number: %s\n", optarg);
			n_threads = strtol(optarg, NULL, 10);
			break;  
			case 'i':  
			printf("max. iterations: %s\n", optarg);
			max_iters = strtol(optarg, NULL, 10);
			break;  
			case ':':  
			printf("option -%c needs a value\n", optopt);  
			break;  
			case '?':  
			fprintf(stderr, "Usage: %s [-p <vector size>] [-s <seed>] [-n <threads number>] [-i <maximum itertions>]\n", argv[0]);
			exit(EXIT_FAILURE);
		}  
	}  
	srand(seed);

	printf("p = %d, seed = %d, n_threads = %d, max_iters = %d\n", \
	 p, seed, n_threads, max_iters);	

	// initializing data
	X = (double*) malloc(sizeof(double) * p);
	Y = (double*) malloc(sizeof(double) * p);
	Y_avgs = (double*) malloc(sizeof(double) * max_iters);

	for(i = 0; i < p; i++){
		X[i] = (double)rand() / RAND_MAX;
		Y[i] = (double)rand() / RAND_MAX;
	}
	for(i = 0; i < max_iters; i++){
		Y_avgs[i] = 0.0;
	}
	a = (double)rand() / RAND_MAX;

#ifdef DEBUG
	printf("vector X= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ",X[i]);
	}
	printf("%f ]\n",X[p-1]);

	printf("vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);

	printf("a= %f \n", a);	
#endif

	/*
	 *	Function to parallelize 
	 */
	gettimeofday(&t_start, NULL);
	//SAXPY iterative SAXPY mfunction
	
    saxpy_setup(X, Y, Y_avgs, a, p, max_iters, n_threads);

	gettimeofday(&t_end, NULL);

#ifdef DEBUG
	printf("RES: final vector Y= [ ");
	for(i = 0; i < p-1; i++){
		printf("%f, ", Y[i]);
	}
	printf("%f ]\n", Y[p-1]);
#endif
	
	// Computing execution time
	exec_time = (t_end.tv_sec - t_start.tv_sec) * 1000.0;  // sec to ms
	exec_time += (t_end.tv_usec - t_start.tv_usec) / 1000.0; // us to ms
	printf("Execution time: %f ms \n", exec_time);
	printf("Last 3 values of Y: %f, %f, %f \n", Y[p-3], Y[p-2], Y[p-1]);
	printf("Last 3 values of Y_avgs: %f, %f, %f \n", Y_avgs[max_iters-3], Y_avgs[max_iters-2], Y_avgs[max_iters-1]);
	return 0;
}	

void saxpy_setup(double *X, double *Y, double *Y_avgs, double a, int p, int max_iters, int n_threads){
    pthread_t threads[n_threads];
    struct saxpy_args args[2];

    args[0].X = X;
    args[0].Y = Y;
    args[0].Y_avgs = Y_avgs;
    args[0].a = a;
    args[0].p = p;
    args[0].max_iters = max_iters;
    args[0].start_pos = 0;
    args[0].end_pos = p/2;

	args[1].X = X;
    args[1].Y = Y;
    args[1].Y_avgs = Y_avgs;
    args[1].a = a;
    args[1].p = p;
    args[1].max_iters = max_iters;
    args[1].start_pos = p/2 + 1;
    args[1].end_pos = p;

    pthread_create(&threads[0], NULL, saxpy_thread, (void *) &args[0]);
	pthread_create(&threads[1], NULL, saxpy_thread, (void *) &args[1]);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);
    
}

void *saxpy_thread(void *arg){
    struct saxpy_args* args = (struct saxpy_args*) arg;
    double *X = args->X;
    double *Y = args->Y;
    double *Y_avgs = args->Y_avgs;
    double a = args->a;
    int p = args->p;
    int max_iters = args->max_iters;
    int start_pos = args->start_pos;
    int end_pos = args->end_pos;

    for(int it = 0; it < max_iters; it++){
		for(int i = start_pos; i < end_pos; i++){
			Y[i] = Y[i] + a * X[i];
			Y_avgs[it] += Y[i];
		}
		Y_avgs[it] = Y_avgs[it] / end_pos;
	}
	printf("Hola desde el hilo con start=%d y end=%d\n", start_pos, end_pos);
}