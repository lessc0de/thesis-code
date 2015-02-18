#include "backwarder.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "io_utils.hpp"
#include "hmm_utils.hpp"
#include "PThreadProcessingDevice.hpp"
#include "performance_description.hpp"
#include "Stage2JobControl.hpp"
#include "debug.hpp"
#include "ds.hpp"

#include <utility>
#include <algorithm>
#include <map>
#include <vector>
#include <deque>
#include <stack>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <climits>
#include <dirent.h>


namespace zipHMM {


  void Backwarder::backward(const Matrix &pi, const Matrix &A, const Matrix &B, const std::vector<double> &scales, Matrix &backward_table) const {
    if(pi.get_width() != 1 || pi.get_height() != A.get_width() || A.get_height() != A.get_width() ||
       B.get_height() != A.get_width() || B.get_width() != ds.get_orig_alphabet_size()) {
      std::cerr << "Dimensions of input matrices do not match:" << std::endl;
      std::cerr << "\t" << "pi width:  " << pi.get_width()  << std::endl;
      std::cerr << "\t" << "pi height: " << pi.get_height() << std::endl;
      std::cerr << "\t" << "A width:  "  << A.get_width()   << std::endl;
      std::cerr << "\t" << "A height: "  << A.get_height()  << std::endl;
      std::cerr << "\t" << "B width:  "  << B.get_width()   << std::endl;
      std::cerr << "\t" << "B height: "  << B.get_height()  << std::endl;
      std::exit(-1);
    }

    // find alphabet and seqs for given number of states
    size_t no_states = A.get_width();
    size_t alphabet_size = 0;
    std::vector<std::vector< unsigned> > sequences;
    for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = ds.nStates2seqs.begin(); it != ds.nStates2seqs.end(); ++it) {
      if(it->first >= no_states) {
	sequences = it->second;
	alphabet_size = ds.get_nStates2alphabet_size().find(it->first)->second;
	break;
      }
    }
    if(sequences.empty()) {
      sequences = ds.get_nStates2seqs().rbegin()->second;
      alphabet_size = ds.get_nStates2alphabet_size().rbegin()->second;
    }

    double *symbol2scale = new double[alphabet_size];
    Matrix *symbol2matrix = new Matrix[alphabet_size];

    compute_symbol2scale_and_symbol2matrix(symbol2matrix, symbol2scale, A, B, alphabet_size);

    for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
      const std::vector<unsigned> &sequence = (*it);
      backward_seq(pi, A, B, sequence, scales, symbol2matrix, backward_table);
    }

    delete[] symbol2scale;
    delete[] symbol2matrix;
  }

  void Backwarder::backward_seq(const Matrix &pi, const Matrix &A, const Matrix &B, const std::vector<unsigned> &sequence, const std::vector<double> &scales, const Matrix *symbol2matrix, Matrix &backward_table) const {
    Matrix res;
    Matrix tmp;
    size_t no_states = A.get_height();

    size_t length = sequence.size();
    backward_table.reset(no_states, length);

    // Initialize
    res.reset(1, A.get_width());
    for (size_t i = 0; i < res.get_width(); ++i) {
      res(0, i) = 1.0;
    }
    for (size_t j = 0; j < no_states; ++j) {
      backward_table(j, length - 1) = res(j, 0);
    }

    // multiply matrices across the sequence
    for(size_t c = length - 2; c <= length - 2; --c) { // weird because of unsigned wrap around when negative.
      Matrix::blas_mult(res, symbol2matrix[sequence[c + 1]], tmp);
      Matrix::copy(tmp, res);

      for (size_t j = 0; j < no_states; ++j) {
        backward_table(j, c) = res(j, 0);
      }

      // Scale the values in res using the scales vector.
      for (size_t j = 0; j < res.get_width(); ++j) {
        res(0, j) /= std::exp(scales[c+1]);
      }

    }
  }

  void Backwarder::compute_symbol2scale_and_symbol2matrix(Matrix *symbol2matrix, double *symbol2scale, const Matrix &A, const Matrix &B, const size_t alphabet_size) const{
    // compute C matrices for each symbol in the original alphabet
    make_em_trans_probs_array(symbol2matrix, A, B);

    // compute C matrices for each symbol in the extended alphabet
    for(size_t i = ds.get_orig_alphabet_size(); i < alphabet_size; ++i) {
      const s_pair symbol_pair = ds.get_symbol2pair().find(unsigned(i))->second;
      const unsigned left_symbol = symbol_pair.second; // the multiplication is done in the reverse direction of the sequence
      const unsigned right_symbol  = symbol_pair.first;
      Matrix &left_matrix  = symbol2matrix[left_symbol];
      Matrix &right_matrix = symbol2matrix[right_symbol];

      Matrix::blas_mult(left_matrix, right_matrix, symbol2matrix[i]);
      // symbol2scale[i] = std::log( symbol2matrix[i].normalize() ) + symbol2scale[left_symbol] + symbol2scale[right_symbol];
    }
  }

} // namespace
