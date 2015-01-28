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

#include "HMMlib/hmm_vector.hpp"
#include "HMMlib/hmm_table.hpp"
#include "HMMlib/hmm.hpp"

#include "fasta.hpp"
#include "hmm.hpp"

#include <pmmintrin.h>
#include <string>


#ifdef WITH_OMP
#include<omp.h>
#endif

template <typename float_type, typename sse_float_type>
float_type viterbi(hmm_t *h, string seq) {
    using namespace hmmlib;
    int no_states = h->states_size;
    int alphabet_size = h->observables_size;
    int seq_length = seq.size();

    boost::shared_ptr<HMMVector<float_type, sse_float_type> > pi_ptr(new HMMVector<float_type, sse_float_type>(no_states));
    boost::shared_ptr<HMMMatrix<float_type, sse_float_type> > T_ptr(new HMMMatrix<float_type, sse_float_type>(no_states,no_states));
    boost::shared_ptr<HMMMatrix<float_type, sse_float_type> > E_ptr(new HMMMatrix<float_type, sse_float_type>(alphabet_size, no_states));

    HMMVector<float_type, sse_float_type> &pi = *pi_ptr;
    HMMMatrix<float_type, sse_float_type> &T = *T_ptr;
    HMMMatrix<float_type, sse_float_type> &E = *E_ptr;

    for (int i = 0; i < no_states; i++) {
        pi(i) = h->initial_probs[i];
    }

    for (int i = 0; i < no_states; i++) {
        for (int j = 0; j < no_states; j++) {
            T(i, j) = h->transition_probs[i][j];
        }
    }

    for (int i = 0; i < no_states; i++) {
        for (int j = 0; j < alphabet_size; j++) {
            E(j, i) = h->emission_probs[i][j];
        }
    }

    HMM<float_type, sse_float_type> hmm(pi_ptr, T_ptr, E_ptr);

    unsigned obs_array[seq_length];

    for (int i = 0; i < seq_length; i++) {
        if (seq[i] == 'A')
            obs_array[i] = 0;
        else if (seq[i] == 'C')
            obs_array[i] = 1;
        else if (seq[i] == 'G')
            obs_array[i] = 2;
        else if (seq[i] == 'T')
            obs_array[i] = 3;
        else
            exit(1);
    }

    sequence obsseq(obs_array, obs_array + seq_length);
    sequence hiddenseq(seq_length);

    double res  = hmm.viterbi(obsseq, hiddenseq);

    return res;
}

int main(int argc, char **argv) {
    struct hmm_t *hmm = hmm_read_path(argv[1], false);
    struct fasta_t *f = fasta_read_path(argv[2]);
    char *seq = (*f->entries).content;
    string s = string(seq, sizeof(seq));

    std::cout << viterbi<float, float>(hmm, s) << std::endl;

    fasta_free(f);
    hmm_free(hmm);
}
