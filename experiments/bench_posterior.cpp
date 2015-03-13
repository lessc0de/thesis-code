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
    if (argc != 6) {
        std::cerr << "Usage: sequence T model N filename" << std::endl;
        exit(1);
    }

    std::string sequence_filename = argv[1];
    int T = atoi(argv[2]);
    std::string hmm_filename = argv[3];
    int N = atof(argv[4]);
    std::ofstream output(argv[5], std::ios::out | std::ios::app);

    if (!output.is_open()) {
        std::cout << "Unable to open output file." << std::endl;
    }

    // Read HMM.
    zipHMM::Matrix pi;
    zipHMM::Matrix A;
    zipHMM::Matrix B;
    zipHMM::read_HMM(pi, A, B, hmm_filename);

    // Run experiment.
    zipHMM::Timer simple_timer;
    zipHMM::Timer zipHMMlib_timer;

    // Simple
    {
        simple_timer.start();
        zipHMM::Matrix pd_table;
        std::vector<double> scales;
        std::vector<unsigned> pd_path;
        zipHMM::posterior_decoding(sequence_filename, pi, A, B, pd_path, pd_table);
        simple_timer.stop();
    }

    // zipHMMlib
    {
        zipHMMlib_timer.start();

        zipHMM::HMMSuite f;
        int alphabet_size = 4;
        int min_num_of_evals = 0;
        f.read_seq(sequence_filename, alphabet_size, min_num_of_evals);

        std::vector<unsigned> pd_path;
        f.posterior_decoding(pi, A, B, pd_path);

        zipHMMlib_timer.stop();
    }

    output << N << " "
           << T << " "
           << simple_timer.timeElapsed() << " "
           << zipHMMlib_timer.timeElapsed() << std::endl;
    output.close();
    exit(0);
}
