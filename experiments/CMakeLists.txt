include_directories("${PROJECT_SOURCE_DIR}/fasta" "${PROJECT_SOURCE_DIR}/hmm" "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/zipHMM/zipHMM")

add_executable(hmm_experiments hmm_experiments.cpp)
target_link_libraries(hmm_experiments fasta_reader hmm_reader)

add_executable(hmm_experiments_new hmm_experiments_new.cpp)
target_link_libraries(hmm_experiments_new zipHMM fasta_reader)

add_executable(complexity complexity.cpp)
target_link_libraries(complexity zipHMM)

add_executable(bench_viterbi bench_viterbi.cpp)
target_link_libraries(bench_viterbi zipHMM)

add_executable(bench_posterior bench_posterior.cpp)
target_link_libraries(bench_posterior zipHMM)

add_executable(bench_subpath_posterior bench_subpath_posterior.cpp)
target_link_libraries(bench_subpath_posterior zipHMM)

add_executable(backward backward.cpp)
target_link_libraries(backward zipHMM)

add_executable(forward_experiment forward.cpp)
target_link_libraries(forward_experiment zipHMM)

add_executable(posterior posterior.cpp)
target_link_libraries(posterior zipHMM)

add_executable(gen_fibonacci gen_fibonacci.cpp)

add_executable(profile_viterbi profile_viterbi.cpp)
target_link_libraries(profile_viterbi zipHMM profiler)

add_executable(profile_posterior profile_posterior.cpp)
target_link_libraries(profile_posterior zipHMM profiler)

add_executable(profile_subpath_posterior profile_subpath_posterior.cpp)
target_link_libraries(profile_subpath_posterior zipHMM profiler)

add_executable(subpath_posterior subpath_posterior.cpp)
target_link_libraries(subpath_posterior zipHMM)

add_executable(viterbi_example viterbi_example.cpp)
target_link_libraries(viterbi_example zipHMM)
