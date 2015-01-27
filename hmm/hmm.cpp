#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#include "hmm.h"

struct hmm_t *hmm_read_path(char *path, bool logspace) {
  FILE *file = fopen(path, "r");
  struct hmm_t *hmm = hmm_read(file, logspace);
  fclose(file);
  return hmm;
}

struct hmm_t *hmm_read(FILE *file, bool logspace) {
  struct hmm_t *hmm = (hmm_t*) malloc(sizeof(struct hmm_t));

  /** STATES **/
  fscanf(file, "states");
  fscanf(file, "%d ", &hmm->states_size);

  hmm->states = (char*) malloc(hmm->states_size * sizeof(int));
  for (int i = 0; i < hmm->states_size; i++) {
    fscanf(file, "%c ", &hmm->states[i]);
  }
  fscanf(file, "\n");

  /** EMISSION COUNT **/
  fscanf(file, "emissionCount");
  hmm->emission_count = (int*) malloc(hmm->states_size * sizeof(int));
  hmm->emission_count_max = INT_MIN;
  for (int i = 0; i < hmm->states_size; i++) {
    fscanf(file, "%i ", &hmm->emission_count[i]);
    if (hmm->emission_count[i] > hmm->emission_count_max) {
      hmm->emission_count_max = hmm->emission_count[i];
    }
  }
  fscanf(file, "\n");

  /** OBSERVABLES **/
  fscanf(file, "observables");
  fscanf(file, "%d ", &hmm->observables_size);

  hmm->observables = (char*) malloc(hmm->observables_size * sizeof(char));
  for (int i = 0; i < hmm->observables_size; i++) {
    fscanf(file, "%s ", &hmm->observables[i]); // eat whitespace.
  }

  /** INITIAL PROBS **/
  fscanf(file, "initProbs");
  hmm->initial_probs = (double*) malloc(hmm->states_size * sizeof(double));
  for (int i = 0; i < hmm->states_size; i++) {
    fscanf(file, "%le ", &hmm->initial_probs[i]);
    if (logspace) {
      hmm->initial_probs[i] = logl(hmm->initial_probs[i]);
    }
  }

  /** TRANSITION PROBS **/
  fscanf(file, "transProbs");
  hmm->transition_probs = (double**) malloc(hmm->states_size * sizeof(double *));
  for (int i = 0; i < hmm->states_size; i++) {
    hmm->transition_probs[i] = (double*) malloc(hmm->states_size * sizeof(double));
    for (int j = 0; j < hmm->states_size; j++) {
      fscanf(file, "%le ", &hmm->transition_probs[i][j]);
      if (logspace) {
        hmm->transition_probs[i][j] = logl(hmm->transition_probs[i][j]);
      }
    }
  }

  /** EMISSION PROBS **/

  // Find the unique emission counts for all states.
  int *emission_count_unique = (int*) malloc(hmm->states_size * sizeof(int));
  memset(emission_count_unique, 0, hmm->states_size * sizeof(int));

  for (int i = 0; i < hmm->states_size; i++) {
    int count = hmm->emission_count[i];
    for (int j = 0; j < hmm->states_size; j++) {
      if (emission_count_unique[j] == count) {
        break;
      } else if (emission_count_unique[j] == 0) {
        emission_count_unique[j] = count;
        break;
      }
    }
  }

  // Find the number of columns we should read.
  int cols = 0;
  for (int i = 0; i < hmm->states_size; i++) {
    if (emission_count_unique[i] == 0) {
      break;
    }
    cols += pow(hmm->observables_size, emission_count_unique[i]);
  }

  qsort(emission_count_unique, hmm->states_size, sizeof(int), &cmp);

  // Compute offsets for emission probabilities for different
  // emission counts in the emission probability matrix.
  hmm->emission_offsets = (int*) malloc(hmm->states_size * sizeof(int));
  for (int i = 0; i < hmm->states_size; i++) {
    if (emission_count_unique[i] == 0) {
      continue;
    } else if (emission_count_unique[i - 1] == 0) {
      hmm->emission_offsets[emission_count_unique[i]] = 0;
    } else {
      hmm->emission_offsets[emission_count_unique[i]] =
        pow(hmm->observables_size, emission_count_unique[i - 1]) +
        hmm->emission_offsets[emission_count_unique[i - 1]];
    }
  }

  free(emission_count_unique);

  fscanf(file, "emProbs");
  hmm->emission_probs = (double**) malloc(hmm->states_size * sizeof(double *));
  for (int i = 0; i < hmm->states_size; i++) {
    hmm->emission_probs[i] = (double*) malloc(cols * sizeof(double));
    for (int j = 0; j < cols; j++) {
      fscanf(file, "%le ", &hmm->emission_probs[i][j]);
      if (logspace) {
        hmm->emission_probs[i][j] = logl(hmm->emission_probs[i][j]);
      }
    }
  }

  return hmm;
}

void hmm_free(struct hmm_t *hmm) {
  if (hmm->initial_probs != NULL) {
    free(hmm->initial_probs);
  }

  for (int i = 0; i < hmm->states_size; i++) {
    if (hmm->transition_probs != NULL && hmm->transition_probs[i] != NULL) {
      free(hmm->transition_probs[i]);
    }
  }
  if (hmm->transition_probs != NULL) {
    free(hmm->transition_probs);
  }

  for (int i = 0; i < hmm->observables_size; i++) {
    if (hmm->emission_probs != NULL && hmm->emission_probs[i] != NULL) {
      free(hmm->emission_probs[i]);
    }
  }
  if (hmm->emission_probs != NULL) {
    free(hmm->emission_probs);
  }

  if (hmm != NULL) {
    free(hmm);
  }
}

void hmm_write_path(char *path, struct hmm_t *hmm) {
  FILE *file = fopen(path, "w");
  hmm_write(file, hmm);
  fclose(file);
}

void hmm_write(FILE *file, struct hmm_t *hmm) {
  /** STATES **/
  fprintf(file, "states\n");
  fprintf(file, "%d\n", hmm->states_size);

  for (int i = 0; i < hmm->states_size; i++) {
    fprintf(file, "%c ", hmm->states[i]);
  }
  fprintf(file, "\n");

  /** EMISSION COUNT **/
  fprintf(file, "emissionCount\n");
  for (int i = 0; i < hmm->states_size; i++) {
    fprintf(file, "%d ", hmm->emission_count[i]);
  }
  fprintf(file, "\n");

  /** OBSERVABLES **/
  fprintf(file, "observables\n");
  fprintf(file, "%d\n", hmm->observables_size);

  for (int i = 0; i < hmm->observables_size; i++) {
    fprintf(file, "%c ", hmm->observables[i]);
  }
  fprintf(file, "\n");

  /** INITIAL PROBS **/
  fprintf(file, "initProbs\n");
  for (int i = 0; i < hmm->states_size; i++) {
    fprintf(file, "%le ", hmm->initial_probs[i]);
  }
  fprintf(file, "\n");

  /** TRANSITION PROBS **/
  fprintf(file, "transProbs\n");
  for (int i = 0; i < hmm->states_size; i++) {
    for (int j = 0; j < hmm->states_size; j++) {
      fprintf(file, "%le ", hmm->transition_probs[i][j]);
    }
    fprintf(file, "\n");
  }

  /** EMISSION PROBS **/
  int emission_cols = hmm->emission_offsets[hmm->emission_count_max]
                    + pow(hmm->observables_size, hmm->emission_count_max);

  fprintf(file, "emProbs\n");
  for (int i = 0; i < hmm->states_size; i++) {
    for (int j = 0; j < emission_cols; j++) {
      fprintf(file, "%le ", hmm->emission_probs[i][j]);
    }
    fprintf(file, "\n");
  }

}

int *hmm_translate_observations_to_indexes(struct hmm_t *hmm, const char *observations, int length) {
  int *observations_idx = (int*) malloc(length * sizeof(int));
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < hmm->observables_size; j++) {
      if (observations[i] == hmm->observables[j]) {
        observations_idx[i] = j;
      }
    }
  }
  return observations_idx;
}

inline static int cmp(const void* pa, const void* pb) {
  int a = *(const int*) pa;
  int b = *(const int*) pb;
  return a - b;
}
