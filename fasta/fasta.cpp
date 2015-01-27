#include <stdlib.h>
#include <ctype.h>

#include "fasta.h"

#define BUFFER_SIZE 4096

struct fasta_t *fasta_new() {
  struct fasta_t *fasta = (fasta_t*) malloc(sizeof(struct fasta_t));
  fasta->entries = NULL;
  fasta->last = NULL;
  fasta->length = 0;
  return fasta;
}

void fasta_free(struct fasta_t *fasta) {
  struct fasta_entry_t *curr_entry = fasta->entries;
  while (curr_entry != NULL) {
    struct fasta_entry_t *next_entry = curr_entry->next;
    free(curr_entry->header);
    free(curr_entry->content);

    curr_entry = next_entry;
  }
}

struct fasta_t *fasta_read_path(char *path) {
  FILE *file = fopen(path, "r");
  if (file == NULL) {
    return NULL;
  }

  struct fasta_t *fasta = fasta_read_file(file);

  fclose(file);
  return fasta;
}

int fasta_line_length(char *line) {
  int i;
  for (i = 0; line[i] != '\n'; i++) { }
  return i;
}

struct fasta_t *fasta_read_file(FILE *file) {
  struct fasta_t *fasta = fasta_new();

  char line[BUFFER_SIZE];

  struct fasta_entry_t *curr_entry = NULL;
  int curr_entry_sequence_pos = 0;
  while (fgets(line, BUFFER_SIZE, file) != NULL) {
    if (line[0] == '>' && curr_entry != NULL) {
      fasta_append_entry(fasta, curr_entry);
      curr_entry = NULL;
    } else if (line[0] == '>') {
      curr_entry = (fasta_entry_t*) malloc(sizeof(struct fasta_entry_t));
      curr_entry->header = (char*) malloc((fasta_line_length(line) - 1) * sizeof(char));
      for (int i = 1; i < fasta_line_length(line) - 1; i++) {
        curr_entry->header[i - 1] = line[i];
      }
    } else if (!isspace(line[0])) {
      int length = fasta_line_length(line);
      if (curr_entry->content == NULL) {
        curr_entry->content = (char*) malloc(length * sizeof(char));
      } else {
        curr_entry->content =
          (char*) realloc(curr_entry->content, (curr_entry_sequence_pos + length) * sizeof(char));
      }
      for (int i = 0; i < length; i++) {
        curr_entry->content[curr_entry_sequence_pos + i] = line[i];
      }
      curr_entry_sequence_pos += length;
    }
  }

  if (curr_entry != NULL) {
    fasta_append_entry(fasta, curr_entry);
  }

  return fasta;
}

void fasta_append_entry(struct fasta_t *fasta, struct fasta_entry_t *entry) {
  if (fasta->last == NULL) {
    fasta->entries = entry;
  } else {
    fasta->last->next = entry;
  }
  fasta->last = entry;
  fasta->length++;
}
