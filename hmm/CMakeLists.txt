cmake_minimum_required (VERSION 2.6)
project (hmm)

set(CMAKE_C_FLAGS "-pedantic -Wall -std=c99 -O3")

add_library(hmm hmm.c)
target_link_libraries(hmm m)

install (TARGETS hmm DESTINATION bin)
install (FILES hmm.h DESTINATION include)