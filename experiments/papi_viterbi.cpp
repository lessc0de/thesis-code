//=========================================================================
// This file is part of HMMlib.
//
// HMMlib is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.
//
// HMMlib is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with HMMlib. If not, see
// <http://www.gnu.org/licenses/>.
//
// Copyright (C) 2010  Bioinformatics Research Centre, Aarhus University.
// Author: Andreas Sand (asand@birc.au.dk)
//=========================================================================

#include "../zipHMM/viterbi.hpp"
#include "../zipHMM/hmm_io.hpp"
#include "../zipHMM/timer.hpp"
#include "../zipHMM/seq_io.hpp"
#include "../zipHMM/hmm_suite.hpp"
#include "/users/muldvang/.local/include/papi.h"

#include <string>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

#ifdef WITH_OMP
#include<omp.h>

#endif

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: sequence model N" << std::endl;
        exit(1);
    }

    std::string sequence_filename = argv[1];
    std::string hmm_filename = argv[2];
    int N = atoi(argv[3]);

    // Read HMM.
    zipHMM::Matrix init_probs;
    zipHMM::Matrix trans_probs;
    zipHMM::Matrix em_probs;
    zipHMM::read_HMM(init_probs, trans_probs, em_probs, hmm_filename);

    std::cout << N << " ";

    // Measure running time
    {
        zipHMM::Timer timer;
        timer.start();

        std::vector<unsigned> viterbi_path;
        zipHMM::HMMSuite h;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 500;
        h.read_seq(sequence_filename, alphabet_size, init_probs.get_height(), min_num_of_evals);
        h.viterbi(init_probs, trans_probs, em_probs, true, viterbi_path);

        timer.stop();
        std::cout << timer.timeElapsed() << " ";
    }

    // Measure cache misses. This is ugly.
    int VPAPI_L1_TCM [1] = {PAPI_L1_TCM};
    int VPAPI_L2_TCM [1] = {PAPI_L2_TCM};
    int VPAPI_L3_TCM [1] = {PAPI_L3_TCM};
    int VPAPI_L1_TCA [1] = {PAPI_L1_TCA};
    int VPAPI_L2_TCA [1] = {PAPI_L2_TCA};
    int VPAPI_L3_TCA [1] = {PAPI_L3_TCA};

    std::vector<int*> events_list;
    events_list.push_back(VPAPI_L1_TCM); //  0x80000006  Yes  Level 1 cache misses
    events_list.push_back(VPAPI_L2_TCM); //  0x80000007  No   Level 2 cache misses
    events_list.push_back(VPAPI_L3_TCM); //  0x80000008  No   Level 3 cache misses
    events_list.push_back(VPAPI_L1_TCA); //  0x80000058  Yes  Level 1 total cache accesses
    events_list.push_back(VPAPI_L2_TCA); //  0x80000059  No   Level 2 total cache accesses
    events_list.push_back(VPAPI_L3_TCA); //  0x8000005a  No   Level 3 total cache accesses

    std::vector<long long> result;

    for (std::vector<int*>::iterator it = events_list.begin();
         it != events_list.end(); ++it) {
        // Start PAPI
        PAPI_start_counters(*it, 1);

        // Run experiment.
        std::vector<unsigned> viterbi_path;
        zipHMM::HMMSuite h;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 500;
        h.read_seq(sequence_filename, alphabet_size, init_probs.get_height(), min_num_of_evals);
        h.viterbi(init_probs, trans_probs, em_probs, true, viterbi_path);

        // Stop PAPI
        std::vector<long long> res(1);
        int errorcode = PAPI_stop_counters(res.data(), res.size());

        // Check for errors
        if (errorcode != PAPI_OK) {
            char errorstring[PAPI_MAX_STR_LEN+1];
            PAPI_perror(errorstring);
            std::cout << "PAPI error " << errorcode << ": " << errorstring
                      << std::endl;
        }

        // Store events
        result.push_back(res[0]);
    }

    for (std::vector<long long>::iterator it = result.begin();
         it != result.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << std::endl;

    exit(0);
}
