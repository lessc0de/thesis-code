#include "../zipHMM/hmm_io.hpp"
#include "../zipHMM/hmm_suite.hpp"

int main() {
    // Read input.
    zipHMM::HMMSuite h;
    size_t alphabet_size = 2;
    size_t min_num_of_evals = 500;
    size_t no_states = 2;

    h.read_seq("../zipHMM/test.seq",
               alphabet_size, no_states, min_num_of_evals);

    // Save for future runs.
    h.write_to_directory("some_dir");

    // Read HMM.
    zipHMM::Matrix pi, A, B;
    zipHMM::read_HMM(pi, A, B, "../zipHMM/test.hmm");

    // Run viterbi.
    bool memory_save = false;
    std::vector<unsigned> viterbi_path;
    double likelihood = h.viterbi(pi, A, B, memory_save, viterbi_path);

    // Print result to stdout.
    std::cout << "Likelihood: " << likelihood << std::endl
              << "Path: ";
    for (std::vector<unsigned>::const_iterator it = viterbi_path.begin();
         it != viterbi_path.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    exit(0);
}
