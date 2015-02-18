#include "forwarder.hpp"
#include "matrix.hpp"
#include "seq_io.hpp"
#include "io_utils.hpp"
#include "hmm_utils.hpp"
#include "PThreadProcessingDevice.hpp"
#include "performance_description.hpp"
#include "Stage2JobControl.hpp"
#include "debug.hpp"
#include "ds.hpp"

#include <algorithm>
#include <map>
#include <vector>
#include <stack>
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <climits>
#include <dirent.h>


namespace zipHMM {


  double Forwarder::forward_helper(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales, bool compute_path, Matrix &forward_table) const {
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

    double ll = 0.0;
    for(std::vector<std::vector<unsigned> >::const_iterator it = sequences.begin(); it != sequences.end(); ++it) {
      const std::vector<unsigned> &sequence = (*it);
      ll += forward_seq(pi, A, B, sequence, symbol2scale, symbol2matrix, scales, compute_path, forward_table);
    }

    delete[] symbol2scale;
    delete[] symbol2matrix;

    return ll;
  }

  double Forwarder::forward(const Matrix &pi, const Matrix &A, const Matrix &B) const {
    std::vector<double> scales;
    Matrix forward_table;
    return Forwarder::forward_helper(pi, A, B, scales, false, forward_table);
  }

  double Forwarder::forward(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales) const {
    Matrix forward_table;
    return Forwarder::forward_helper(pi, A, B, scales, false, forward_table);
  }

  double Forwarder::forward(const Matrix &pi, const Matrix &A, const Matrix &B, std::vector<double> &scales, Matrix &forward_table) const {
    return Forwarder::forward_helper(pi, A, B, scales, true, forward_table);
  }

  double Forwarder::forward_seq(const Matrix &pi, const Matrix &A, const Matrix &B, const std::vector<unsigned> &sequence, const double *symbol2scale, const Matrix *symbol2matrix, std::vector<double> &scales, bool compute_table, Matrix &forward_table) const {
    Matrix res;
    Matrix tmp;
    double loglikelihood = 0;
    unsigned no_states = A.get_width();

    if (compute_table) {
      forward_table.reset(no_states, sequence.size());
    }

    // compute C_1 and push corresponding scale
    scales.push_back( std::log(init_apply_em_prob(res, pi, B, sequence[0])) );
    if (compute_table) {
      for (size_t j = 0; j < no_states; ++j) {
        forward_table(j, 0) = res(j, 0);
      }
    }

    // multiply matrices across the sequence
    for(size_t i = 1; i < sequence.size(); ++i) {
      Matrix::blas_matrix_vector_mult(symbol2matrix[sequence[i]], res, tmp);
      Matrix::copy(tmp, res);

      if (compute_table) {
        for (size_t j = 0; j < no_states; ++j) {
          forward_table(j, i) = res(j, 0);
        }
      }

      scales.push_back( std::log(res.normalize()) + symbol2scale[sequence[i]] );
    }

    // compute loglikelihood by summing log of scales
    for(std::vector<double>::iterator it = scales.begin(); it != scales.end(); ++it) {
      loglikelihood += (*it);
    }

    return loglikelihood;
  }

  double Forwarder::pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, const DeviceDescriptor &device_descriptor) const {
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

    double *symbol2scale;
    Matrix *symbol2matrix;
    std::vector<ProcessingDevice*> devices;
    size_t numberOfDevices;
    Matrix *result;
    Matrix *temp;
    double loglikelihood;
    size_t head;

    // find alphabet and seqs for given number of states
    size_t no_states = A.get_width();
    size_t alphabet_size = 0;
    const std::vector<std::vector< unsigned> > *sequences;
    for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = ds.get_nStates2seqs().begin(); it != ds.get_nStates2seqs().end(); ++it) {
      if(it->first >= no_states) {
	sequences = &(it->second);
	alphabet_size = ds.get_nStates2alphabet_size().find(it->first)->second;
	break;
      }
    }
    if(sequences == 0) {
      sequences = &(ds.get_nStates2seqs().rbegin()->second);
      alphabet_size = ds.get_nStates2alphabet_size().rbegin()->second;
    }

    symbol2scale = new double[alphabet_size];
    symbol2matrix = new Matrix[alphabet_size];

    compute_symbol2scale_and_symbol2matrix(symbol2matrix, symbol2scale, A, B, alphabet_size);

    result = new Matrix();
    temp = new Matrix();

    std::map<unsigned, s_pair> symbol2pair = ds.get_symbol2pair();
    numberOfDevices = device_descriptor.getNDevices();
    for(unsigned i = 0; i < numberOfDevices; ++i) {
      devices.push_back(device_descriptor.createDevice(i));
      devices[i]->setParameters(&pi, &A, &B, &symbol2pair, symbol2scale, symbol2matrix);
    }

    loglikelihood = 0.0;
    for(std::vector<std::vector<unsigned> >::const_iterator it = sequences->begin(); it != sequences->end(); ++it) {
      for(unsigned i = 0; i < numberOfDevices; ++i) {
	devices[i]->setSeq(&(*it));
      }

      const size_t length = it->size();
      const size_t nBlocks = size_t(std::sqrt(length));
      Stage2JobControl control(length, nBlocks);

      devices[0]->likelihoodVector(control);
      for(unsigned i = 1; i < numberOfDevices; ++i)
	devices[i]->likelihoodMatrix(control);

      for(unsigned i = 0; i < numberOfDevices; ++i)
	devices[i]->join();

      head = control.headBlock - 1;

      Matrix::copy(*control.resultMatrices[head], *result);
      loglikelihood += control.resultLogLikelihoods[head];

      for(size_t i = head + 1; i < nBlocks; ++i) {
	Matrix::mult(*control.resultMatrices[i], *result, *temp);
	std::swap(result, temp);

	loglikelihood += LinearSpace::toLogSpace(result->normalize()) + control.resultLogLikelihoods[i];
      }
    }

    for(unsigned i = 0; i < devices.size(); ++i)
      delete devices[i];

    delete result;
    delete temp;

    delete[] symbol2scale;
    delete[] symbol2matrix;

    return loglikelihood;
  }

  double Forwarder::pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, const std::string &device_filename) const {
    std::vector<DeviceDescriptor> device_descriptors;
    readDescriptors(device_descriptors, device_filename);
    return pthread_forward(pi, A, B, device_descriptors[0]);
  }


  double Forwarder::mr_pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, const DeviceDescriptor &device_descriptor) const {
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

    double *symbol2scale;
    Matrix *symbol2matrix;
    const std::vector<std::vector< unsigned> > *sequences;
    std::vector<ProcessingDevice*> devices;
    size_t numberOfDevices;
    double loglikelihood;

    // find alphabet and seqs for given number of states
    size_t no_states = A.get_width();
    size_t alphabet_size = 0;
    for(std::map<size_t, std::vector<std::vector<unsigned> > >::const_iterator it = ds.get_nStates2seqs().begin(); it != ds.get_nStates2seqs().end(); ++it) {
      if(it->first >= no_states) {
	sequences = &(it->second);
	alphabet_size = ds.get_nStates2alphabet_size().find(it->first)->second;
	break;
      }
    }
    if(sequences == 0) {
      sequences = &(ds.get_nStates2seqs().rbegin()->second);
      alphabet_size = ds.get_nStates2alphabet_size().rbegin()->second;
    }

    symbol2scale = new double[alphabet_size];
    symbol2matrix = new Matrix[alphabet_size];

    compute_symbol2scale_and_symbol2matrix(symbol2matrix, symbol2scale, A, B, alphabet_size);

    std::map<unsigned, s_pair> symbol2pair = ds.get_symbol2pair();
    numberOfDevices = device_descriptor.getNDevices();
    for(unsigned i = 0; i < numberOfDevices; ++i) {
      devices.push_back(device_descriptor.createDevice(i));
      devices[i]->setParameters(&pi, &A, &B, &symbol2pair, symbol2scale, symbol2matrix);
      devices[i]->setSeqs(sequences);
    }

    const size_t length = sequences->size();
    const size_t nBlocks = size_t(std::sqrt(length));

    MapReduceJobControl control(length, nBlocks);

    for(unsigned i = 0; i < numberOfDevices; ++i)
      devices[i]->mapReduceLoglikelihood(control);

    for(unsigned i = 0; i < numberOfDevices; ++i)
      devices[i]->join();

    loglikelihood = 0.0;
    for(size_t i = 0; i < nBlocks; ++i) {
      loglikelihood += control.resultLogLikelihoods[i];
    }

    for(unsigned i = 0; i < devices.size(); ++i)
      delete devices[i];

    delete[] symbol2scale;
    delete[] symbol2matrix;

    return loglikelihood;
  }

  double Forwarder::mr_pthread_forward(const Matrix &pi, const Matrix &A, const Matrix &B, const std::string &device_filename) const {
    std::vector<DeviceDescriptor> device_descriptors;
    readDescriptors(device_descriptors, device_filename);
    return mr_pthread_forward(pi, A, B, device_descriptors[0]);
  }


  void Forwarder::compute_symbol2scale_and_symbol2matrix(Matrix *symbol2matrix, double *symbol2scale, const Matrix &A, const Matrix &B, const size_t alphabet_size) const{
    // compute C matrices for each symbol in the original alphabet
    make_em_trans_probs_array(symbol2scale, symbol2matrix, A, B);

    // compute C matrices for each symbol in the extended alphabet
    for(size_t i = ds.get_orig_alphabet_size(); i < alphabet_size; ++i) {
      const s_pair symbol_pair = ds.get_symbol2pair().find(unsigned(i))->second;
      const unsigned left_symbol = symbol_pair.second; // the multiplication is done in the reverse direction of the sequence
      const unsigned right_symbol  = symbol_pair.first;
      Matrix &left_matrix  = symbol2matrix[left_symbol];
      Matrix &right_matrix = symbol2matrix[right_symbol];

      Matrix::blas_mult(left_matrix, right_matrix, symbol2matrix[i]);
      symbol2scale[i] = std::log( symbol2matrix[i].normalize() ) + symbol2scale[left_symbol] + symbol2scale[right_symbol];
    }
  }

} // namespace
