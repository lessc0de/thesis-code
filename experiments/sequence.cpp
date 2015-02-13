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

#include <string>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cmath>

#ifdef WITH_OMP
#include<omp.h>

#endif

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: model sequence n." << std::endl;
        exit(1);
    }

    std::string hmm_filename = argv[1];

    std::string sequence_filename = argv[2];

    int n = atof(argv[3]);

    // Read HMM.
    zipHMM::Matrix init_probs;
    zipHMM::Matrix trans_probs;
    zipHMM::Matrix em_probs;
    zipHMM::read_HMM(init_probs, trans_probs, em_probs, hmm_filename);

    // Run experiment.
    std::vector<unsigned> viterbi_path;
    zipHMM::Timer simple_pre_timer;
    zipHMM::Timer zipHMMlib_uncompressed_pre_timer;
    zipHMM::Timer zipHMMlib_uncompressed_path_pre_timer;
    zipHMM::Timer zipHMMlib_pre_timer;
    zipHMM::Timer zipHMMlib_path_pre_timer;
    zipHMM::Timer simple_running_timer;
    zipHMM::Timer zipHMMlib_uncompressed_path_running_timer;
    zipHMM::Timer zipHMMlib_uncompressed_running_timer;
    zipHMM::Timer zipHMMlib_running_timer;
    zipHMM::Timer zipHMMlib_path_running_timer;

    // Simple
    {
        std::cout << n << " ";
        simple_pre_timer.start();
        simple_pre_timer.stop();
        std::cout << simple_pre_timer.timeElapsed() << " ";

        simple_running_timer.start();
        zipHMM::viterbi(sequence_filename, init_probs, trans_probs, em_probs, viterbi_path);
        simple_running_timer.stop();
        std::cout << simple_running_timer.timeElapsed() << " ";
    }

    // zipHMMlib uncompressed
    {
        zipHMMlib_uncompressed_pre_timer.start();
        zipHMM::Viterbi v1;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 0;
        v1.read_seq(sequence_filename, alphabet_size, min_num_of_evals);
        zipHMMlib_uncompressed_pre_timer.stop();
        std::cout << zipHMMlib_uncompressed_pre_timer.timeElapsed() << " ";

        zipHMMlib_uncompressed_running_timer.start();
        v1.viterbi(init_probs, trans_probs, em_probs);
        zipHMMlib_uncompressed_running_timer.stop();
        std::cout << zipHMMlib_uncompressed_running_timer.timeElapsed() << " ";
    }

    // zipHMMlib uncompressed path
    {
        zipHMMlib_uncompressed_path_pre_timer.start();
        zipHMM::Viterbi v1;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 0;
        v1.read_seq(sequence_filename, alphabet_size, min_num_of_evals);
        zipHMMlib_uncompressed_path_pre_timer.stop();
        std::cout << zipHMMlib_uncompressed_path_pre_timer.timeElapsed() << " ";

        zipHMMlib_uncompressed_path_running_timer.start();
        v1.viterbi(init_probs, trans_probs, em_probs, viterbi_path);
        zipHMMlib_uncompressed_path_running_timer.stop();
        std::cout << zipHMMlib_uncompressed_path_running_timer.timeElapsed() << " ";
    }

    // zipHMMlib
    {
        zipHMMlib_pre_timer.start();
        zipHMM::Viterbi v;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 500;
        v.read_seq(sequence_filename, alphabet_size, min_num_of_evals);
        zipHMMlib_pre_timer.stop();
        std::cout << zipHMMlib_pre_timer.timeElapsed() << " ";

        zipHMMlib_running_timer.start();
        v.viterbi(init_probs, trans_probs, em_probs);
        zipHMMlib_running_timer.stop();
        std::cout << zipHMMlib_running_timer.timeElapsed() << " ";
    }

    // zipHMMlib_path
    {
        zipHMMlib_path_pre_timer.start();
        zipHMM::Viterbi v2;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 500;
        v2.read_seq(sequence_filename, alphabet_size, min_num_of_evals);
        zipHMMlib_path_pre_timer.stop();
        std::cout << zipHMMlib_path_pre_timer.timeElapsed() << " ";

        zipHMMlib_path_running_timer.start();
        v2.viterbi(init_probs, trans_probs, em_probs, viterbi_path);
        zipHMMlib_path_running_timer.stop();
        std::cout << zipHMMlib_path_running_timer.timeElapsed() << " ";

        std::cout << std::endl;
    }
    std::cout.flush();
    exit(0);
}
