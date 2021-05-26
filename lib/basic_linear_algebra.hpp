#ifndef DIFFERENTIAL_LINEAR_EXPERIMENTS_SBOX_EXPERIMENTS_BASIC_LINEAR_ALGEBRA_HPP_
#define DIFFERENTIAL_LINEAR_EXPERIMENTS_SBOX_EXPERIMENTS_BASIC_LINEAR_ALGEBRA_HPP_

#include <bit>
#include <cassert>
#include <vector>

#include <NTL/mat_GF2.h>

namespace basic_linear_algebra {

// Passing from uint to vec_GF2 and back
NTL::vec_GF2 ConvertToNtl(uint64_t x, size_t length);
uint64_t ConvertToUint(const NTL::vec_GF2 &x);
NTL::vec_GF2 Embed(uint64_t x, const NTL::mat_GF2 &matrix);
NTL::mat_GF2 &AppendMatrix(NTL::mat_GF2 &target, const NTL::mat_GF2 &other);

// Calculation of related spaces
NTL::mat_GF2 ComplementSpace(NTL::mat_GF2 &factor_space,
                             bool already_rre = false);
NTL::mat_GF2 OrthogonalComplement(const NTL::mat_GF2 &matrix);

// Reduced row echelon form and related stuff
long Rre(NTL::mat_GF2 &mat);
bool VecGf2Lt(const NTL::vec_GF2 &x, const NTL::vec_GF2 &y);
bool RreLt(const NTL::mat_GF2 &mat, const NTL::mat_GF2 &other);
std::vector<std::vector<NTL::mat_GF2>> ListOfVectorSpaces(long dim,
                                                          long max_subspace_dim);
template<class VecType>
void walsh_hadamard_inplace(VecType &truth_table, size_t bitlength) {
  size_t window_size = ((size_t) (1u) << bitlength);
  size_t number_of_windows = 1;
  for (size_t i = 0; i < bitlength; ++i) {
    window_size >>= 1u;
    for (size_t j = 0; j < number_of_windows; ++j) {
      for (size_t k = 0; k < window_size; ++k) {
        truth_table[(2 * j) * window_size + k] +=
            truth_table[(2 * j + 1) * window_size + k];
        truth_table[(2 * j + 1) * window_size + k] *= -2;
        truth_table[(2 * j + 1) * window_size + k] +=
            truth_table[(2 * j) * window_size + k];
      }
    }
    number_of_windows <<= 1u;
  }
}

};

#endif //DIFFERENTIAL_LINEAR_EXPERIMENTS_SBOX_EXPERIMENTS_BASIC_LINEAR_ALGEBRA_HPP_
