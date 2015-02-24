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
#include "../zipHMM/posterior_decoding.hpp"

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
        std::cerr << "Usage: sequence model_folder k." << std::endl;
        exit(1);
    }

    std::string seq_filename = argv[1];

    std::string data_folder = argv[2];
    std::stringstream s;
    s << data_folder << argv[3] << ".hmm";
    std::string hmm_filename = s.str();

    int k = atof(argv[3]);

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
        zipHMM::posterior_decoding(seq_filename, pi, A, B, pd_path, pd_table);

        simple_timer.stop();
    }

    // zipHMMlib
    {
        zipHMMlib_timer.start();

        // Compute forward table.
        zipHMM::HMMSuite f;
        int alphabet_size = 4;
        int min_num_of_evals = 0;
        f.read_seq(seq_filename, alphabet_size, min_num_of_evals);

        std::vector<unsigned> pd_path;
        f.posterior_decoding(pi, A, B, pd_path);

        zipHMMlib_timer.stop();
    }

    std::cout << k << " "
              << simple_timer.timeElapsed() << " "
              << zipHMMlib_timer.timeElapsed() << std::endl;
    exit(0);
}
