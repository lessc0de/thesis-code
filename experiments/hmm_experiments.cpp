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

#include "fasta/fasta.h"
#include "hmm/hmm.h"

#include <pmmintrin.h>
using namespace hmmlib;


#ifdef WITH_OMP
#include<omp.h>
#endif

template <typename float_type, typename sse_float_type>
void viterbi_2x2_with_transitions() {
    int no_states = 2;
    int alphabet_size = 2;
    int seq_length = 4;

    boost::shared_ptr<HMMVector<float_type, sse_float_type> > pi_ptr(new HMMVector<float_type, sse_float_type>(no_states));
    boost::shared_ptr<HMMMatrix<float_type, sse_float_type> > T_ptr(new HMMMatrix<float_type, sse_float_type>(no_states,no_states));
    boost::shared_ptr<HMMMatrix<float_type, sse_float_type> > E_ptr(new HMMMatrix<float_type, sse_float_type>(alphabet_size, no_states));

    HMMVector<float_type, sse_float_type> &pi = *pi_ptr;
    HMMMatrix<float_type, sse_float_type> &T = *T_ptr;
    HMMMatrix<float_type, sse_float_type> &E = *E_ptr;

    pi(0)  = 0.2; pi(1) = 0.8;

    // transitions from state 0
    T(0,0) = 0.1; T(0,1) = 0.9;
    // transitions from state 1
    T(1,0) = 0.9; T(1,1) = 0.1;

    // emissions from state 0
    E(0,0) = 0.25; E(1,0) = 0.75;
    // emissions from state 1
    E(0,1) = 0.25; E(1,1) = 0.75;

    HMM<float_type, sse_float_type> hmm(pi_ptr, T_ptr, E_ptr);

    unsigned obs_array[] = { 0, 1, 0, 1 };
    sequence obsseq(obs_array, obs_array + seq_length);
    sequence hiddenseq(seq_length);

    hmm.viterbi(obsseq, hiddenseq);
}

int main() {
    struct hmm_t *hmm = hmm_read_path(argv[1], true);
    struct fasta_t *fasta = fasta_read_path(argv[2]);

    viterbi_2x2_with_transitions<double, double>();
}
