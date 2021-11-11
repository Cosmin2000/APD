#include <stdlib.h>
#include "genetic_algorithm.h"





struct thread_arg {
	int id;
	int P;
    const sack_object *objects;
    int object_count;
    int generations_count;
    int sack_capacity;
    individual *new;
    individual **current_generation ;
	individual **next_generation;
    pthread_barrier_t *barrier;
    pthread_barrier_t *barrier1;
};

void *f(void *arg) {

    struct thread_arg* th_arg = (struct thread_arg*) arg;
    int id = th_arg->id;
    int P = th_arg->P;  
    int count, cursor;
    int object_count = th_arg->object_count;
    int start = id * (double)object_count / P;
    int end = min(((id + 1) * (double)object_count / P), object_count);

    for (int i = start; i < end; ++i) {
		(*th_arg->current_generation)[i].fitness = 0;
		(*th_arg->current_generation)[i].chromosomes = (int*) calloc(object_count, sizeof(int));
		(*th_arg->current_generation)[i].chromosomes[i] = 1;
		(*th_arg->current_generation)[i].index = i;
		(*th_arg->current_generation)[i].chromosome_length = object_count;
		(*th_arg->next_generation)[i].fitness = 0;
		(*th_arg->next_generation)[i].chromosomes = (int*) calloc(object_count, sizeof(int));
		(*th_arg->next_generation)[i].index = i;
	    (*th_arg->next_generation)[i].chromosome_length = object_count;
	}

    pthread_barrier_wait(th_arg->barrier);

    for (int k = 0; k < th_arg->generations_count; ++k) {
        cursor = 0;
     
        parallel_compute_fitness_function(th_arg->objects, &(*th_arg->current_generation), object_count, th_arg->sack_capacity,start,end);
        pthread_barrier_wait(th_arg->barrier);
    
        mergesort((*th_arg->current_generation),id,object_count,P,th_arg->barrier1, th_arg->new);
        pthread_barrier_wait(th_arg->barrier);
        
        count = object_count * 3 / 10;

        // keep first 30% children (elite children selection)
        start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);
		for (int i = start; i < end; i++) {
			copy_individual((*th_arg->current_generation) + i, (*th_arg->next_generation) + i);
		}
        pthread_barrier_wait(th_arg->barrier);

        cursor = count;

		// mutate first 20% children with the first version of bit string mutation
		count = object_count * 2 / 10;

        // mutate first 20% children with the first version of bit string mutation
        start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);

        for (int i = start; i < end; i++) {
			copy_individual((*th_arg->current_generation) + i, (*th_arg->next_generation) + cursor + i);
			mutate_bit_string_1((*th_arg->next_generation) + cursor + i, k);
		}

        pthread_barrier_wait(th_arg->barrier);

        cursor += count;
        count = object_count * 2 / 10;

        // mutate next 20% children with the second version of bit string mutation
        start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);

        for (int i = start; i < end; i++) {
			copy_individual((*th_arg->current_generation) + i + count, (*th_arg->next_generation) + cursor + i);
			mutate_bit_string_2((*th_arg->next_generation) + cursor + i, k);
		}
        pthread_barrier_wait(th_arg->barrier);

        cursor += count;
        count = object_count * 3 / 10;

        if (count % 2 == 1) {
			 copy_individual((*th_arg->current_generation) + object_count - 1, (*th_arg->next_generation) + cursor + count - 1);
             count--;
		}

		// crossover first 30% parents with one-point crossover
		start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);
        if (start % 2 == 1) {
            end--;
        }

        for (int i = start; i  < end ; i += 2) {
			crossover((*th_arg->current_generation) + i, (*th_arg->next_generation) + cursor + i, k);
		}

        pthread_barrier_wait(th_arg->barrier);

        if (id == 0) {
            
            // switch to new generation
            individual *tmp = NULL;
		    tmp = (*th_arg->current_generation);
		    (*th_arg->current_generation) = (*th_arg->next_generation);
		    (*th_arg->next_generation) = tmp;

        }
        pthread_barrier_wait(th_arg->barrier);

        start = id * (double)object_count / P;
        end = min(((id + 1) * (double)object_count / P), object_count);

        for (int i = start; i < end; i++) {
			(*th_arg->current_generation)[i].index = i;
		}
        pthread_barrier_wait(th_arg->barrier);

        if (id == 0) {
            if (k % 5 == 0) {
			    print_best_fitness((*th_arg->current_generation));
		    }
        }
        pthread_barrier_wait(th_arg->barrier);
    }
   
    start = th_arg->id * (double)object_count / P;
    end = min(((th_arg->id + 1) * (double)object_count / P), object_count);

    parallel_compute_fitness_function(th_arg->objects,&(*th_arg->current_generation), object_count, th_arg->sack_capacity, start, end);
    pthread_barrier_wait(th_arg->barrier);

    mergesort((*th_arg->current_generation),id,object_count,P,th_arg->barrier1, th_arg->new);
    pthread_barrier_wait(th_arg->barrier);

    if (id == 0) {
	    print_best_fitness((*th_arg->current_generation));
    }



    return NULL;

}



int main(int argc, char *argv[]) {

    struct thread_arg *arguments;
    pthread_t *threads;

	// array with all the objects that can be placed in the sack
	sack_object *objects = NULL;

	// number of objects
	int object_count = 0;

	// maximum weight that can be carried in the sack
	int sack_capacity = 0;

	// number of generations
	int generations_count = 0;

    //Number of Threads
    int P = (int)atoi(argv[3]);

    int  r;
    pthread_barrier_t barrier, barrier1;

    pthread_barrier_init(&barrier,NULL, P);

    pthread_barrier_init(&barrier1,NULL, P);


	if (!read_input(&objects, &object_count, &sack_capacity, &generations_count, argc, argv)) {
		return 0;
	}
    

    threads = (pthread_t*) malloc(P * sizeof(pthread_t));
	arguments = (struct thread_arg*) malloc(P * sizeof(struct thread_arg));
    individual *current_generation = (individual*) calloc(object_count, sizeof(individual));
	individual *next_generation = (individual*) calloc(object_count, sizeof(individual));
    individual *new = (individual*) calloc(object_count, sizeof(individual));
    
   

    for (int i = 0; i < P; i++) {
		arguments[i].id = i;
		arguments[i].P = P;
        arguments[i].object_count = object_count;
        arguments[i].sack_capacity = sack_capacity;
		arguments[i].generations_count = generations_count;
        arguments[i].objects = objects;
        arguments[i].current_generation = &current_generation;
        arguments[i].next_generation = &next_generation;
        arguments[i].barrier = &barrier;
        arguments[i].barrier1 = &barrier1;
        arguments[i].new  = new;
        

		r = pthread_create(&threads[i], NULL, f, &arguments[i]);

		if (r) {
			printf("Eroare la crearea thread-ului %d\n", i);
			exit(-1);
		}
	}

    for (int i = 0; i < P; i++) {
		r = pthread_join(threads[i], NULL);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}

	//run_genetic_algorithm(objects, object_count, generations_count, sack_capacity);


    pthread_barrier_destroy(&barrier);

    pthread_barrier_destroy(&barrier1);
	
    // free resources for old generation
	free_generation(current_generation);
	free_generation(next_generation);

    free(current_generation);
    free(next_generation);
    free(objects);
    free(new);
    free(threads);
    free(arguments);

	return 0;
}
