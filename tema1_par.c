#include <stdlib.h>
#include "parallel_genetic_algorithm.h"


void *thread_func(void *arg) {


    // Extrag datele pentru fiecare thread
    struct thread_arg* th_arg = (struct thread_arg*) arg;
    int id = th_arg->id;
    int P = th_arg->P;  
    int count, cursor;
    int object_count = th_arg->object_count;

    // Le obtin indicii pentru paralelizarea obtinerii generatiei initiale
    int start = id * (double)object_count / P;
    int end = min(((id + 1) * (double)object_count / P), object_count);

    // paralelizez crearea generatiei initiale
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
        
        // paralelizez calcularea fitness-ului fiecarui individ din generatie
        parallel_compute_fitness_function(th_arg->objects, &(*th_arg->current_generation), object_count, th_arg->sack_capacity,start,end);
        pthread_barrier_wait(th_arg->barrier);
    
        // sortez paralel indivizii descrescator dupa fitness sau crescator dupa numarul de indivizi/index
        parallel_mergesort((*th_arg->current_generation),id,object_count,P,th_arg->barrier_merge, th_arg->new);
        pthread_barrier_wait(th_arg->barrier);
        
        count = object_count * 3 / 10;

        // pastrez primii 30 % din indivizi (selectia elitei). Paralelizez operatia calculand indicii de start si end pentru fiecare thread
        start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);

		for (int i = start; i < end; i++) {
			copy_individual((*th_arg->current_generation) + i, (*th_arg->next_generation) + i);
		}

        cursor = count;

		// Aplic prima varianta de mutatie pe primii 20% din indivizi din generatia curenta (functia mutate_bit_string_1). 
        // Paralelizez operatia calculand indicii de start si end pentru fiecare thread
		count = object_count * 2 / 10;

        start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);

        for (int i = start; i < end; i++) {
			copy_individual((*th_arg->current_generation) + i, (*th_arg->next_generation) + cursor + i);
			mutate_bit_string_1((*th_arg->next_generation) + cursor + i, k);
		}

        cursor += count;
        count = object_count * 2 / 10;

        // Aplic a doua varianta de mutatie pe primii 20% din indivizi din generatia curenta  (functia mutate_bit_string_2). 
        // Paralelizez operatia calculand indicii de start si end pentru fiecare thread
        start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);

        for (int i = start; i < end; i++) {
			copy_individual((*th_arg->current_generation) + i + count, (*th_arg->next_generation) + cursor + i);
			mutate_bit_string_2((*th_arg->next_generation) + cursor + i, k);
		}

        cursor += count;
        count = object_count * 3 / 10;

        // Inainte de a aplica crossover, daca numarul de parinti este impar, pastrez ultimul individ din generatie
        if (count % 2 == 1 && id == 0) {
			 copy_individual((*th_arg->current_generation) + object_count - 1, (*th_arg->next_generation) + cursor + count - 1);
             count--;
		}

		// Paralelizez operatia calculand indicii de start si end pentru fiecare thread.
        // Daca indicele de start pentru un thread este impar, scad indicele de final pentru ca altfel depaseste array-ul.
		start = id * (double)count / P;
        end = min(((id + 1) * (double)count / P), count);
        if (start % 2 == 1) {
            end--;
        }

        // Aplic crossover intr-un punct pe primii 30 % indivizi din generatia curenta (functia crossover).
        for (int i = start; i  < end ; i += 2) {
			crossover((*th_arg->current_generation) + i, (*th_arg->next_generation) + cursor + i, k);
		}

        pthread_barrier_wait(th_arg->barrier);

        // Interschimb generatiile.
        if (id == 0) {    
            
            individual *tmp = NULL;
		    tmp = (*th_arg->current_generation);
		    (*th_arg->current_generation) = (*th_arg->next_generation);
		    (*th_arg->next_generation) = tmp;

        }
        pthread_barrier_wait(th_arg->barrier);

        start = id * (double)object_count / P;
        end = min(((id + 1) * (double)object_count / P), object_count);

        //Setez corect index-ul pentru noua generatie curenta
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

    // paralelizez calcularea fitness-ului fiecarui individ din noua generatie curenta
    parallel_compute_fitness_function(th_arg->objects,&(*th_arg->current_generation), object_count, th_arg->sack_capacity, start, end);
    pthread_barrier_wait(th_arg->barrier);

    // sortez paralel indivizii din noua generatie curenta, descrescator dupa fitness sau crescator dupa numarul de indivizii/index
    parallel_mergesort((*th_arg->current_generation),id,object_count,P,th_arg->barrier_merge, th_arg->new);
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

    pthread_barrier_t barrier, barrier_merge;

    pthread_barrier_init(&barrier,NULL, P);

    pthread_barrier_init(&barrier_merge,NULL, P);


	if (!read_input(&objects, &object_count, &sack_capacity, &generations_count, argc, argv)) {
		return 0;
	}
    

    threads = (pthread_t*) malloc(P * sizeof(pthread_t));
	arguments = (struct thread_arg*) malloc(P * sizeof(struct thread_arg));
    individual *current_generation = (individual*) calloc(object_count, sizeof(individual));
	individual *next_generation = (individual*) calloc(object_count, sizeof(individual));

    // Array folosit pentru mergesort
    individual *new = (individual*) calloc(object_count, sizeof(individual));
    
   
    //Asignez argumentele pentru fiecare thread 
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
        arguments[i].barrier_merge = &barrier_merge;
        arguments[i].new  = new;
        

		r = pthread_create(&threads[i], NULL, thread_func, &arguments[i]);

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

    pthread_barrier_destroy(&barrier);

    pthread_barrier_destroy(&barrier_merge);
	
    // Dezaloc memoria
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
