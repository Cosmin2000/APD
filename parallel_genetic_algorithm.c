#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "parallel_genetic_algorithm.h"


int min (int a, int b) {
	return a < b ? a : b;
}

int read_input(sack_object **objects, int *object_count, int *sack_capacity, int *generations_count, int argc, char *argv[])
{
	FILE *fp;

	if (argc < 4) {
		fprintf(stderr, "Usage:\n\t./tema1 in_file generations_count P\n");
		return 0;
	}

	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		return 0;
	}

	if (fscanf(fp, "%d %d", object_count, sack_capacity) < 2) {
		fclose(fp);
		return 0;
	}

	if (*object_count % 10) {
		fclose(fp);
		return 0;
	}

	sack_object *tmp_objects = (sack_object *) calloc(*object_count, sizeof(sack_object));

	for (int i = 0; i < *object_count; ++i) {
		if (fscanf(fp, "%d %d", &tmp_objects[i].profit, &tmp_objects[i].weight) < 2) {
			free(objects);
			fclose(fp);
			return 0;
		}
	}

	fclose(fp);

	*generations_count = (int) strtol(argv[2], NULL, 10);
	
	if (*generations_count == 0) {
		free(tmp_objects);

		return 0;
	}

	*objects = tmp_objects;

	return 1;
}

void print_objects(const sack_object *objects, int object_count)
{
	for (int i = 0; i < object_count; ++i) {
		printf("%d %d\n", objects[i].weight, objects[i].profit);
	}
}

void print_generation(const individual *generation, int limit)
{
	for (int i = 0; i < limit; ++i) {
		for (int j = 0; j < generation[i].chromosome_length; ++j) {
			printf("%d ", generation[i].chromosomes[j]);
		}

		printf("\n%d - %d\n", i, generation[i].fitness);
	}
}

void print_best_fitness(const individual *generation)
{
	printf("%d\n", generation[0].fitness);
}

void parallel_compute_fitness_function(const sack_object *objects, individual **generation, int object_count, int sack_capacity,int start, int end)
{
	int weight;
	int profit;
	int i;
	int j;
	individual *gener = *generation;
	for (i = start; i < end; i++) {
		weight = 0;
		profit = 0;

		for ( j = 0; j < (gener[i]).chromosome_length; ++j) {
			if ((gener[i]).chromosomes[j]) {
				weight += objects[j].weight;
				profit += objects[j].profit;
			}
		}

		(gener[i]).fitness = (weight <= sack_capacity) ? profit : 0;
	}
}

int compare(individual a, individual b) {
	int k;
	if (a.fitness > b.fitness) {
		return 1;
	} else if (a.fitness == b.fitness) {
		int first_count = 0, second_count = 0;

		for (k = 0; k < a.chromosome_length && k < b.chromosome_length; ++k) {
			first_count += a.chromosomes[k];
			second_count += b.chromosomes[k];
		}

		if ( first_count < second_count) {
			return 1;
		} else if (first_count == second_count) {
			if (a.index < b.index) {
				return 1;
			}
		}
	}
	return -1;
}

void merge(individual *source, int start, int mid, int end, individual *destination) {
	int iA = start; // start-ul primului sub-array
	int iB = mid; // start-ul celui de-al doilea sub-array
	int i;

	for (i = start; i < end; i++) {
		// Daca am terminat al doilea sub-array sau trebuie ales elementul din primul sub-array
		if (end == iB || (iA < mid && compare(source[iA], source[iB]) == 1 )) {
			destination[i] = source[iA];
			iA++;
		} else {
			// este ales elementul din al doilea sub-array
			destination[i] = source[iB];
			iB++;
		}
		
	}
}

void parallel_mergesort(individual *current,int id,int N, int P, pthread_barrier_t *barrier, individual *new) {
	int i, merges;
	int start, end;
	int width;
	
	// folosit pentru interschimbare
	individual *aux;

	if (id == 0) {
		aux = calloc(N , sizeof(individual));
	}
	
	pthread_barrier_wait(barrier);

	
	for (width = 1; width < N; width = 2 * width) {

		// Pentru ca mergesort-ul sa mearga pentru orice dimensiune, daca este cazul, fac un merge in plus. 
		// Folosesc functia ceil pentru asta.
		merges = ceil((double)N / (2*width));
		//Paralelizez merge-urile.
		start = id * merges / P * 2 *width;
    	end = min(((id + 1) * merges/ P * 2 * width), N);
	
		for (i = start; i < end; i = i + 2 * width) {
			merge(current, i,min(i + width,N),min(i + 2 * width, N), new);
		}
		
		// Astept sa termine toate thread-urile
		pthread_barrier_wait(barrier);
		
		if (id == 0) {	
			// Updatez generatia. Folosesc memmove fiindca cu interschimbare simpla cu pointeri, 
			// current si new  pointeaza acceasi zona de memorie.	
			memmove(aux,current, N * sizeof(individual));
			memmove(current,new, N * sizeof(individual));
			memmove(new,aux, N * sizeof(individual));			
		}
		
		pthread_barrier_wait(barrier);
		// Dupa ce toate thread-urile au terminat pot porni urmatorul nivel de merge-uri.
		
	}

	if (id == 0) {
		// Dezaloc vectorul alocat pentru interschimbare.
		free(aux);
	}
}



void mutate_bit_string_1(const individual *ind, int generation_index)
{
	int i, mutation_size;
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	if (ind->index % 2 == 0) {
		// for even-indexed individuals, mutate the first 40% chromosomes by a given step
		mutation_size = ind->chromosome_length * 4 / 10;
		for (i = 0; i < mutation_size; i += step) {
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
		}
	} else {
		// for even-indexed individuals, mutate the last 80% chromosomes by a given step
		mutation_size = ind->chromosome_length * 8 / 10;
		for (i = ind->chromosome_length - mutation_size; i < ind->chromosome_length; i += step) {
			ind->chromosomes[i] = 1 - ind->chromosomes[i];
		}
	}
}

void mutate_bit_string_2(const individual *ind, int generation_index)
{
	int step = 1 + generation_index % (ind->chromosome_length - 2);

	// mutate all chromosomes by a given step
	for (int i = 0; i < ind->chromosome_length; i += step) {
		ind->chromosomes[i] = 1 - ind->chromosomes[i];
	}
}

void crossover(individual *parent1, individual *child1, int generation_index)
{
	individual *parent2 = parent1 + 1;
	individual *child2 = child1 + 1;
	int count = 1 + generation_index % parent1->chromosome_length;

	memcpy(child1->chromosomes, parent1->chromosomes, count * sizeof(int));
	memcpy(child1->chromosomes + count, parent2->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));

	memcpy(child2->chromosomes, parent2->chromosomes, count * sizeof(int));
	memcpy(child2->chromosomes + count, parent1->chromosomes + count, (parent1->chromosome_length - count) * sizeof(int));
}

void copy_individual(const individual *from, const individual *to)
{
	memcpy(to->chromosomes, from->chromosomes, from->chromosome_length * sizeof(int));
}

void free_generation(individual *generation)
{
	int i;

	for (i = 0; i < generation->chromosome_length; ++i) {
		free(generation[i].chromosomes);
		generation[i].chromosomes = NULL;
		generation[i].fitness = 0;
	}
}

