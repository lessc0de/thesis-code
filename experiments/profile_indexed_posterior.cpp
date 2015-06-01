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
#include "../fasta/fasta_reader.hpp"
#include "../zipHMM/forwarder.hpp"
#include "../zipHMM/timer.hpp"
#include "../zipHMM/seq_io.hpp"
#include "../zipHMM/hmm_suite.hpp"

#include <string>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <sys/time.h>

#ifdef WITH_OMP
#include<omp.h>
#endif

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Please provide HMM model and fasta file." << std::endl;
        exit(1);
    }

    // Read HMM.
    zipHMM::Matrix pi;
    zipHMM::Matrix A;
    zipHMM::Matrix B;
    zipHMM::read_HMM(pi, A, B, std::string(argv[1]));

    std::string sequence_filename = argv[2];

    std::vector<unsigned> pd_path;

    zipHMM::HMMSuite h;
    size_t alphabet_size = 4;
    size_t min_num_of_evals = 500;

    h.read_seq(sequence_filename, alphabet_size, pi.get_height(), min_num_of_evals);

    for (int k = 0; k < min_num_of_evals; ++k) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        std::srand(tv.tv_sec * 1000000 + tv.tv_usec);
        int i = std::rand() % (h.get_orig_seq_length() - (h.get_orig_seq_length() / 100));
        int j = i + h.get_orig_seq_length() / 100;

        assert(i < h.get_orig_seq_length());
        assert(j < h.get_orig_seq_length());
        assert(i < j);

        h.indexed_posterior_decoding(pi, A, B, i, j, pd_path);
    }

    exit(0);
}
