#ifndef HMM_H
#define HMM_H

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

struct hmm_t {
  char *states;
  int states_size;

  char *observables;
  int observables_size;

  int *emission_count;
  int emission_count_max;
  int *emission_offsets;

  double *initial_probs;
  double **transition_probs;
  double **emission_probs;
};

void hmm_free(struct hmm_t *hmm);

struct hmm_t *hmm_read_path(char *path, bool logspace);
struct hmm_t *hmm_read(FILE *file, bool logspace);

void hmm_write(FILE *file, struct hmm_t *hmm);
void hmm_write_path(char *path, struct hmm_t *hmm);

int *hmm_translate_observations_to_indexes(struct hmm_t *hmm, const char *observations, int length);

inline static int cmp(const void* pa, const void* pb);

#endif
