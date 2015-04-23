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

#include "../zipHMM/hmm_io.hpp"
#include "../zipHMM/timer.hpp"
#include "../zipHMM/seq_io.hpp"
#include "../zipHMM/hmm_suite.hpp"
#include "../zipHMM/posterior_decoding.hpp"

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
    if (argc != 7) {
        std::cerr << "Usage: sequence T model N filename subseq_length" << std::endl;
        exit(1);
    }

    std::string sequence_filename = argv[1];
    int T = atoi(argv[2]);
    std::string hmm_filename = argv[3];
    int N = atof(argv[4]);
    std::ofstream output(argv[5], std::ios::out | std::ios::app);
    int subseq_length = atof(argv[6]);

    if (!output.is_open()) {
        std::cout << "Unable to open output file." << std::endl;
        output.close();
        exit(1);
    }

    if (subseq_length > T) {
        std::cout << "subseq_length cannot be larger than sequence!" << std::endl;
        output.close();
        exit(2);
    }

    std::srand(time(0));
    unsigned i = std::rand() % (T - subseq_length);
    unsigned j = i + subseq_length;

    // Read HMM.
    zipHMM::Matrix pi;
    zipHMM::Matrix A;
    zipHMM::Matrix B;
    zipHMM::read_HMM(pi, A, B, hmm_filename);

    // Run experiment.
    zipHMM::Timer simple_pre_timer;
    zipHMM::Timer one_pre_timer;
    zipHMM::Timer many_pre_timer;

    zipHMM::Timer simple_running_timer;
    zipHMM::Timer one_running_timer;
    zipHMM::Timer many_running_timer;

    // Simple
    {
        simple_pre_timer.start();
        zipHMM::Matrix pd_table;
        std::vector<double> scales;
        std::vector<unsigned> pd_path;
        simple_pre_timer.stop();

        simple_running_timer.start();
        zipHMM::posterior_decoding(sequence_filename, pi, A, B, pd_path, pd_table);
        std::vector<unsigned> subseq_path;
        subseq_path.insert(subseq_path.begin(), pd_path.begin() + i, pd_path.begin() + j);
        simple_running_timer.stop();
    }

    // One
    {
        one_pre_timer.start();

        zipHMM::HMMSuite f;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 1;
        f.read_seq(sequence_filename, alphabet_size, pi.get_height(), min_num_of_evals);

        one_pre_timer.stop();

        one_running_timer.start();
        std::vector<unsigned> pd_path;
        f.indexed_posterior_decoding(pi, A, B, i, j, pd_path);

        one_running_timer.stop();
    }

    // Many
    {
        many_pre_timer.start();

        zipHMM::HMMSuite f;
        size_t alphabet_size = 4;
        size_t min_num_of_evals = 500;
        f.read_seq(sequence_filename, alphabet_size, pi.get_height(), min_num_of_evals);

        many_pre_timer.stop();

        many_running_timer.start();
        std::vector<unsigned> pd_path;
        f.indexed_posterior_decoding(pi, A, B, i, j, pd_path);

        many_running_timer.stop();
    }

    output << N << " "
           << T << " "
           << subseq_length << " "
           << simple_pre_timer.timeElapsed() << " "
           << simple_running_timer.timeElapsed() << " "
           << one_pre_timer.timeElapsed() << " "
           << one_running_timer.timeElapsed() << " "
           << many_pre_timer.timeElapsed() << " "
           << many_running_timer.timeElapsed() << std::endl;
    output.close();
    exit(0);
}
