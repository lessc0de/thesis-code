include_directories("${PROJECT_SOURCE_DIR}/fasta" "${PROJECT_SOURCE_DIR}/hmm" "${PROJECT_SOURCE_DIR}")

add_executable(hmm_experiments hmm_experiments.cpp)
target_link_libraries(hmm_experiments fasta_reader hmm_reader)
