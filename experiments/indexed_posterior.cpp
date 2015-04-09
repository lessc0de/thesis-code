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
#include "../zipHMM/timer.hpp"
#include "../zipHMM/seq_io.hpp"
#include "../zipHMM/posterior_decoding.hpp"
#include "../zipHMM/matrix.hpp"
#include "../zipHMM/hmm_suite.hpp"

#include <string>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

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

    std::string seq_filename = argv[2];


    // Reference implementation.
    std::vector<unsigned> ref_pd_path;
    {
        zipHMM::Matrix pd_table;
        std::vector<double> scales;

        zipHMM::posterior_decoding(seq_filename, pi, A, B, ref_pd_path, pd_table);
    }

    // zipHMM implementation
    std::vector<unsigned> zip_pd_path;
    {
        zipHMM::HMMSuite f;
        int alphabet_size = 4;
        int min_num_of_evals = 500;
        f.read_seq(seq_filename, alphabet_size, min_num_of_evals);

        std::srand(std::time(0));
        int i = std::rand() % f.get_orig_seq_length();
        int j = std::rand() % f.get_orig_seq_length();

        if (i > j) {
            int tmp = i;
            i = j;
            j = tmp;
        }
        assert(i <= j);

        f.indexed_posterior_decoding(pi, A, B, i, j, zip_pd_path);
        zip_pd_path.insert(zip_pd_path.begin(), ref_pd_path.begin(), ref_pd_path.begin() + i);
        zip_pd_path.insert(zip_pd_path.end(), ref_pd_path.begin() + j, ref_pd_path.end());
    }

    if (ref_pd_path.size() != zip_pd_path.size()) {
        std::cout << "Paths not identical!" << std::endl;
        std::cout << "Length: " << ref_pd_path.size() << "  ";
        for (size_t i = 0; i < 20 && i < ref_pd_path.size(); ++i) {
            std::cout << ref_pd_path[i] << " ";
        }
        std::cout << "..." << std::endl;
        std::cout << "Length: " << zip_pd_path.size() << "  ";
        for (size_t i = 0; i < 20 && i < zip_pd_path.size(); ++i) {
            std::cout << zip_pd_path[i] << " ";
        }
        std::cout << "..." << std::endl;
        exit(1);
    } else {
        std::cout << "Paths identical!" << std::endl;
        for (size_t i = 0; i < 10; ++i) {
            std::cout << ref_pd_path[i] << " ";
        }
        std::cout << "..." << std::endl;
        for (size_t i = 0; i < 10; ++i) {
            std::cout << zip_pd_path[i] << " ";
        }
        std::cout << "..." << std::endl;
    }

    exit(0);
}
