cmake_minimum_required (VERSION 2.6)
project (fasta)

set(CMAKE_C_FLAGS "-pedantic -Wall -std=c99 -O3")

add_library(fasta fasta.c)

install (TARGETS fasta DESTINATION bin)
install (FILES fasta.h DESTINATION include)