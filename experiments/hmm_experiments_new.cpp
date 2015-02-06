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

#include <string>
#include <iomanip>
#include <vector>
#include <iostream>
#include <cstdlib>

#ifdef WITH_OMP
#include<omp.h>
#endif

std::vector<unsigned> read_sequence(char *path) {
    struct fasta_t *f = fasta_read_path(path);
    char *seq = (*f->entries).content;

    int seq_length = (*f->entries).length;
    std::vector<unsigned> obs_array(seq_length);

    for (int i = 0; i < seq_length; i++) {
        if (seq[i] == 'A')
            obs_array[i] = 0;
        else if (seq[i] == 'C')
            obs_array[i] = 1;
        else if (seq[i] == 'G')
            obs_array[i] = 2;
        else if (seq[i] == 'T')
            obs_array[i] = 3;
        else {
            std::cout << "Unknown symbol " << seq[i] << std::endl;
            exit(1);
        }
    }
    fasta_free(f);
    return obs_array;
}

bool valid_path(std::vector<unsigned> &sequence, std::vector<unsigned> &path,
                double likelihood, zipHMM::Matrix &init_probs,
                zipHMM::Matrix &trans_probs, zipHMM::Matrix &em_probs) {
    // Init
    double l = std::log(init_probs(path[0], 0) * em_probs(path[0], sequence[0]));

    // Recursion
    for (size_t c = 1; c < sequence.size(); ++c) {
        l += std::log(trans_probs(path[c - 1], path[c]) * em_probs(path[c], sequence[c]));
    }

    return std::abs(l - likelihood) < 1e-3;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Please provide HMM model and fasta file." << std::endl;
        exit(1);
    }

    // Read HMM.
    zipHMM::Matrix init_probs;
    zipHMM::Matrix trans_probs;
    zipHMM::Matrix em_probs;
    zipHMM::read_HMM(init_probs, trans_probs, em_probs, std::string(argv[1]));

    std::string sequence = std::string(argv[2]);

    std::vector<unsigned> ref_viterbi_path;
    zipHMM::Timer ref_timer;
    ref_timer.start();
    double res = zipHMM::viterbi(sequence, init_probs, trans_probs, em_probs, ref_viterbi_path);
    ref_timer.stop();
    std::cout << "Reference:\t" << res << "\t" << ref_timer.timeElapsed() << std::endl;


    std::vector<unsigned> my_viterbi_path;
    zipHMM::Viterbi v;
    size_t alphabet_size = 4;
    size_t min_num_of_evals = 500;
    zipHMM::Timer comp_timer;
    comp_timer.start();
    v.read_seq(sequence, alphabet_size, min_num_of_evals);
    double my_new_res = v.viterbi(init_probs, trans_probs, em_probs, my_viterbi_path);
    comp_timer.stop();
    std::cout << "My new impl.:\t" << my_new_res << "\t" << comp_timer.timeElapsed() << std::endl;

    if (std::abs(res - my_new_res) > 1e-3) {
        std::cout << "Likelihoods not identical!" << std::endl;
        exit(1);
    }
    std::vector<unsigned> seq;
    zipHMM::readSeq(seq, sequence);
    if (!valid_path(seq, my_viterbi_path, my_new_res, init_probs, trans_probs, em_probs)) {
        std::cout << "Viterbi paths not identical!" << std::endl;
        exit(2);
    }

    exit(0);
}
