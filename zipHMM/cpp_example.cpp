#include "hmm_io.hpp"
#include "forwarder.hpp"
#include "matrix.hpp"
#include "posterior_decoding.hpp"
#include "viterbi.hpp"

#include <iostream>
#include <vector>

using namespace zipHMM;

int main(int argc, char **args) {
  
  // Create forwarder object
  Forwarder f1;
  // and read input
  size_t alphabet_size = 3;
  f1.read_seq("example.seq", alphabet_size);

  // Save for future runs
  f1.write_to_directory("example.out");

  // Read in a concrete HMM
  Matrix pi, A, B;
  read_HMM(pi, A, B, "example.hmm");

  // Compute the log-likelihood
  std::cout << "loglikelihood: "
  	    << f1.forward(pi, A, B)
  	    << std::endl;
  std::cout << "pthread loglikelihood: "
     	    << f1.pthread_forward(pi, A, B)
     	    << std::endl;
  std::cout << "mr_pthread loglikelihood: "
     	    << f1.mr_pthread_forward(pi, A, B)
     	    << std::endl;
  
  Forwarder f2;
  f2.read_from_directory("example.out");

  std::cout << "from dir loglikelihood: "
	    << f2.forward(pi, A, B)
	    << std::endl;
  // std::cout << "loglikelihood: "
  //  	    << f2.pthread_forward(pi, A, B)
  //    	    << std::endl;

  Forwarder f_seqs;
  alphabet_size = 3;
  f_seqs.read_seq_directory("example_seq_directory", alphabet_size);
  std::cout << "seqs loglikelihood: "
	    << f_seqs.forward(pi, A, B)
	    << std::endl;
  std::cout << "pthread seqs loglikelihood: "
	    << f_seqs.pthread_forward(pi, A, B)
	    << std::endl;
  std::cout << "mr_pthread loglikelihood: "
     	    << f_seqs.mr_pthread_forward(pi, A, B)
     	    << std::endl;

  // std::vector<unsigned> viterbi_path;
  // double viterbi_ll = viterbi("example.seq", pi, A, B, viterbi_path);
  // 
  // std::cout << "viterbi loglikelihood: " << viterbi_ll << std::endl;
  // 
  std::vector<unsigned> pd_path;
  Matrix pd_table;
  posterior_decoding("example.seq", pi, A, B, pd_path, pd_table);
  std::cout << "posterior path[0:10]: ";
  for(size_t i = 0; i < 10; ++i)
    std::cout << pd_path[i] << " ";
  std::cout << std::endl;
  
  // for(size_t r = 0; r < 2; ++r) {
  //   for(size_t c = 0; c < 10; ++c) {
  //     std::cout << pd_table(r, c) << "\t";
  //   }
  //   std::cout << std::endl;
  // }

  
  std::vector<unsigned> viterbi_path;
  double viterbi_ll = viterbi("example.seq", pi, A, B, viterbi_path);
  std::cout << "viterbi loglikelihood:\t" << viterbi_ll << std::endl;
  std::cout << "viterbi path[0:10]: ";
  for(size_t i = 0; i < 10; ++i)
    std::cout << viterbi_path[i] << " ";
  std::cout << std::endl;

  return 0;

  // pi, A and B can also be created using
  // Matrix pi(4, 0);
  // pi(0,0) = 0.5;
  // ... and so on
}
