#ifndef FASTA_H
#define FASTA_H

#include <stdio.h>

struct fasta_entry_t {
  char *header;
  char *content;
  struct fasta_entry_t *next;
};

struct fasta_t {
  unsigned int length;
  struct fasta_entry_t *entries;
  struct fasta_entry_t *last;
};

struct fasta_t *fasta_new();
void fasta_free(struct fasta_t *fasta);

struct fasta_t *fasta_read_path(char *path);
struct fasta_t *fasta_read_file(FILE *file);

void fasta_append_entry(struct fasta_t *fasta, struct fasta_entry_t *entry);

#endif
