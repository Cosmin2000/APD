#include <stdlib.h>
#include "genetic_algorithm.h"


pthread_barrier_t barrier;

struct my_arg {
	int id;
	int P;
    const sack_object *objects;
    int object_count;
    int generations_count;
    int sack_capacity;
    int *count;
    int *cursor;
    individual *current_generation ;
	individual *next_generation;
};


void *f(void *arg) {
    //int cursor,count;
    struct my_arg* data = (struct my_arg*) arg;
    //int i;
    int start = data->id * (double)data->object_count / data->P;
    int end = min(((data->id + 1) * (double)data->object_count / data->P), data->object_count);

    for (int i = start; i < end; ++i) {
		data->current_generation[i].fitness = 0;
		data->current_generation[i].chromosomes = (int*) calloc(data->object_count, sizeof(int));
		data->current_generation[i].chromosomes[i] = 1;
		data->current_generation[i].index = i;
		data->current_generation[i].chromosome_length = data->object_count;

		data->next_generation[i].fitness = 0;
		data->next_generation[i].chromosomes = (int*) calloc(data->object_count, sizeof(int));
		data->next_generation[i].index = i;
		data->next_generation[i].chromosome_length = data->object_count;
	}

    pthread_barrier_wait(&barrier);

    for (int k = 0; k < data->generations_count; ++k) {
        if(data->id == 0) {
            *data->cursor = 0;
        }
        // compute fitness and sort by it
        parallel_compute_fitness_function(data->objects, data->current_generation, data->object_count, data->sack_capacity,start,end);
        pthread_barrier_wait(&barrier);
        if (data->id == 0) {
            qsort(data->current_generation, data->object_count, sizeof(individual), cmpfunc);    
            //Calculez 30% din numarul de indivizi pentru a-i pastra
            *data->count = data->object_count * 3 / 10;
        }
        pthread_barrier_wait(&barrier);

        // keep first 30% children (elite children selection)
        start = data->id * (double)*data->count / data->P;
        end = min(((data->id + 1) * (double)*data->count / data->P), data->object_count);
		for (int i = start; i < end; ++i) {
			copy_individual(data->current_generation + i, data->next_generation + i);
		}

        pthread_barrier_wait(&barrier);
        if(data->id == 0) {
		    *data->cursor = *data->count;
            // Calculez primii 20% din indivizi
            *data->count = data->object_count * 2 / 10;
        }
        pthread_barrier_wait(&barrier);

        // mutate first 20% children with the first version of bit string mutation
        start = data->id * (double)*data->count / data->P;
        end = min(((data->id + 1) * (double)*data->count / data->P), data->object_count);

        for (int i = start; i < end; ++i) {
			copy_individual(data->current_generation + i, data->next_generation + *data->cursor + i);
			mutate_bit_string_1(data->next_generation + *data->cursor + i, k);
		}

        pthread_barrier_wait(&barrier);
        if (data->id == 0) {
            *data->cursor += *data->count;
            //calculez count pentru urmatorul task
            *data->count = data->object_count * 2 / 10;
        }
        pthread_barrier_wait(&barrier);


        // mutate next 20% children with the second version of bit string mutation
        start = data->id * (double)*data->count / data->P;
        end = min(((data->id + 1) * (double)*data->count / data->P), data->object_count);

        for (int i = start; i < end; ++i) {
			copy_individual(data->current_generation + i + *data->count, data->next_generation + *data->cursor + i);
			mutate_bit_string_2(data->next_generation + *data->cursor + i, k);
		}
        pthread_barrier_wait(&barrier);

        if (data->id == 0) {
            *data->cursor += *data->count;
            //calculez count pentru urmatorul task
            *data->count = data->object_count * 3 / 10;
            // (if there is an odd number of parents, the last one is kept as such)
            if (*data->count % 2 == 1) {
			    copy_individual(data->current_generation + data->object_count - 1, data->next_generation + *data->cursor + *data->count - 1);
			    (*data->count)--;
		    }

        }
        pthread_barrier_wait(&barrier);


		// crossover first 30% parents with one-point crossover
		start = data->id * (double)*data->count / data->P;
        end = min(((data->id + 1) * (double)*data->count / data->P), *data->count);

        //Aici apare seg fault la Tread (id 3) La inputs/in4 5 4 .    ULtimul i a fost 1497
        for (int i = start; i  < end - 1; i += 2) {
			crossover(data->current_generation + i, data->next_generation + (*data->cursor) + i, k);
		}

        pthread_barrier_wait(&barrier);

        if (data->id == 0) {
            // switch to new generation
            individual *tmp = NULL;
		    tmp = data->current_generation;
		    data->current_generation = data->next_generation;
		    data->next_generation = tmp;
        }
        pthread_barrier_wait(&barrier);

        start = data->id * (double)data->object_count / data->P;
        end = min(((data->id + 1) * (double)data->object_count / data->P), data->object_count);

        for (int i = start; i < end; ++i) {
			data->current_generation[i].index = i;
		}
        pthread_barrier_wait(&barrier);

        if (data->id == 0) {
            if (k % 5 == 0) {
			    print_best_fitness(data->current_generation);
		    }
        }
        pthread_barrier_wait(&barrier);
    }

    start = data->id * (double)data->object_count / data->P;
    end = min(((data->id + 1) * (double)data->object_count / data->P), data->object_count);
    parallel_compute_fitness_function(data->objects,data->current_generation, data->object_count, data->sack_capacity,start,end);
    
    pthread_barrier_wait(&barrier);

    if (data->id == 0) {
        qsort(data->current_generation, data->object_count, sizeof(individual), cmpfunc);
	    print_best_fitness(data->current_generation);
    }



    return NULL;

}


int main(int argc, char *argv[]) {

    struct my_arg *arguments;
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

    int  r, count, cursor;

    pthread_barrier_init(&barrier,NULL, P);


	if (!read_input(&objects, &object_count, &sack_capacity, &generations_count, argc, argv)) {
		return 0;
	}
    

    threads = (pthread_t*) malloc(P * sizeof(pthread_t));
	arguments = (struct my_arg*) malloc(P * sizeof(struct my_arg));
    individual *current_generation = (individual*) calloc(object_count, sizeof(individual));
	individual *next_generation = (individual*) calloc(object_count, sizeof(individual));

    for (int i = 0; i < P; i++) {
		arguments[i].id = i;
		arguments[i].P = P;
        arguments[i].count = &count;
        arguments[i].cursor = &cursor;
        arguments[i].object_count = object_count;
        arguments[i].sack_capacity = sack_capacity;
		arguments[i].generations_count = generations_count;
        arguments[i].objects = objects;
        arguments[i].current_generation = current_generation;
        arguments[i].next_generation = next_generation;

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
	
    // free resources for old generation
	free_generation(current_generation);
	free_generation(next_generation);

    free(current_generation);
    free(next_generation);
    free(objects);
    free(threads);
    free(arguments);

	return 0;
}
