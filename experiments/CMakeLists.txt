include_directories("${PROJECT_SOURCE_DIR}/fasta" "${PROJECT_SOURCE_DIR}/hmm" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/zipHMM/zipHMM")

add_executable(hmm_experiments hmm_experiments.cpp)
target_link_libraries(hmm_experiments fasta_reader hmm_reader)

add_executable(hmm_experiments_new hmm_experiments_new.cpp)
target_link_libraries(hmm_experiments_new zipHMM fasta_reader)

add_executable(complexity complexity.cpp)
target_link_libraries(complexity zipHMM)

add_executable(model model.cpp)
target_link_libraries(model zipHMM)

add_executable(sequence sequence.cpp)
target_link_libraries(sequence zipHMM)
