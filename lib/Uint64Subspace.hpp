#ifndef CONDITIONS_SRC_UINT64SUBSPACE_HPP_
#define CONDITIONS_SRC_UINT64SUBSPACE_HPP_
#include <numeric>

#include <NTL/mat_GF2.h>
#include "basic_linear_algebra.hpp"
class Uint64Subspace {
 private:
  bool parity_check_out_of_date = true;
  NTL::mat_GF2 saved_parity_check_matrix;
 public:
  bool init_ = false;
  ssize_t truncation_ = 0;
  friend std::ostream &operator<<(std::ostream &stream,
                                  const Uint64Subspace &x) {
    auto old_flags = stream.flags();
    stream << "Subspace of dimension " << std::dec << x.space_.NumRows() <<
           " with basis: {" << std::hex;
    const auto rows = x.space_.NumRows();
    for (auto y = 0; y < rows; ++y) {
      stream << NtlToUint64(x.space_[y]) << ",";
    }
    stream << "}";
    stream.flags(old_flags);
    return stream;
  }
  static NTL::vec_GF2 Uint64ToNtl(const uint64_t &x, size_t truncation) {
    NTL::vec_GF2 result;
    result.SetLength(truncation);
    const auto limit = 1u << truncation;
    for (uint64_t i = 0; i < limit; ++i) {
      result[i] = static_cast<long>((x >> i) & 1u);
    }
    return result;
  }

  static uint64_t NtlToUint64(const NTL::vec_GF2 &x) {
    return basic_linear_algebra::ConvertToUint(x);
  }
  Uint64Subspace() {
    truncation_ = 0;
    init_ = true;
  }
  Uint64Subspace(NTL::mat_GF2 ints, ssize_t truncation, bool skip_rre = false)
      : truncation_(truncation), space_(std::move(ints)) {
    auto rg = space_.NumRows();
    if (!skip_rre) {
      rg = basic_linear_algebra::Rre(space_);
    }
    space_.SetDims(rg, truncation_);
    init_ = true;
  }

  template<typename RANGE>
  Uint64Subspace(const RANGE &ints, ssize_t truncation, bool skip_rre = false)
      : truncation_(truncation) {
    space_.SetDims(0, truncation);
    // super inefficient
    for (const uint64_t &x: ints) {
      space_.SetDims(space_.NumRows() + 1, truncation);
      space_[space_.NumRows() - 1] = Uint64ToNtl(x, truncation);
    }
    auto rg = space_.NumRows();
    if (!skip_rre) {
      rg = basic_linear_algebra::Rre(space_);
    }
    space_.SetDims(rg, truncation_);
    init_ = true;
  }

  [[nodiscard]] Uint64Subspace OrthogonalComplement() const {
    auto kernel = basic_linear_algebra::OrthogonalComplement(space_);
    Uint64Subspace result(kernel, truncation_);
    return result;
  }
  [[nodiscard]] NTL::vec_GF2 ProjectOntoCanonicalDirectComplement(NTL::vec_GF2 x_ntl) const {
    // WE NEED space_ to be in RRE here
    const auto rows = space_.NumRows();
    for (long i = 0; i < rows; ++i) {
      // Find leading 1
      long leading = -1;
      while (++leading < truncation_ && NTL::IsZero(space_[i][leading]));
      if (!NTL::IsZero(x_ntl[leading]))
        x_ntl += space_[i];
    }
    return x_ntl;
  }

  [[nodiscard]] uint64_t ProjectOntoCanonicalDirectComplement(uint64_t x) const {
    // WE NEED space_ to be in RRE here
    NTL::vec_GF2 x_ntl = basic_linear_algebra::ConvertToNtl(x, truncation_);
    x_ntl = ProjectOntoCanonicalDirectComplement(x_ntl);
    return basic_linear_algebra::ConvertToUint(x_ntl);
  }

  [[nodiscard]] Uint64Subspace CanonicalDirectComplement() const {
    auto copy(space_);
    auto basis = basic_linear_algebra::ComplementSpace(copy);
    Uint64Subspace result(basis, truncation_);
    return result;
  }

  [[nodiscard]] uint64_t ElementK(uint64_t i) const {
    return basic_linear_algebra::ConvertToUint(basic_linear_algebra::Embed(i,
                                                                           space_));
  }

  [[nodiscard]] std::vector<uint64_t> Elements() const {
    std::vector<uint64_t> result(1u << static_cast<uint64_t>(space_.NumRows()));
    std::iota(result.begin(), result.end(), 0);
    std::transform(result.begin(),
                   result.end(),
                   result.begin(),
                   [&](uint64_t &i) {
                     return this->ElementK(i);
                   });
    return result;
  }

  [[nodiscard]] size_t Dimension() const {
    return space_.NumRows();
  }

  bool operator<=(const Uint64Subspace &other) const {
    //assert(truncation_ == other.truncation_);
    const auto dim = Dimension();
    for (size_t i = 0; i < dim; ++i) {
      // super inefficient
      if (!other.ContainsElement(space_[i])) {
        return false;
      }
    }
    return true;
  }

  bool operator<(const Uint64Subspace &other) const {
    return (*this <= other) && (other.Dimension() > Dimension());
  }

  bool operator==(const Uint64Subspace &other) const {
    return (*this <= other) && (other <= *this);
  }

  [[nodiscard]] bool RreLt(const Uint64Subspace &other) const {
    return basic_linear_algebra::RreLt(space_, other.space_);
  }

  [[nodiscard]] bool ContainsElement(NTL::vec_GF2 elm_ntl) const {
    assert(elm_ntl.length() == truncation_);
/*
    auto copy(space_);
    copy.SetDims(copy.NumRows()+1, copy.NumCols());
    copy[space_.NumRows()] = elm_ntl;
    auto rk = gauss(copy);
    return rk == Dimension();
    */
    elm_ntl = ProjectOntoCanonicalDirectComplement(elm_ntl);
    return IsZero(elm_ntl);
  }

  [[nodiscard]] bool ContainsElement(uint64_t elm) const {
    auto elm_ntl = Uint64ToNtl(elm, truncation_);
    return ContainsElement(elm_ntl);
  }

  Uint64Subspace &operator+=(const Uint64Subspace &other) {
    assert(other.truncation_ == truncation_);
    parity_check_out_of_date = true;
    auto old_dim = space_.NumRows();
    space_.SetDims(space_.NumRows() + other.space_.NumRows(),
                   space_.NumCols());
    const auto rows = other.space_.NumRows();
    for (long i = 0; i < rows; ++i) {
      space_[i + old_dim] = other.space_[i];
    }
    auto rg = basic_linear_algebra::Rre(space_);
    space_.SetDims(rg, truncation_);
    return *this;
  }

  Uint64Subspace &operator+=(const NTL::vec_GF2 &other) {
    assert(other.length() == truncation_);
    parity_check_out_of_date = true;
    auto old_dim = space_.NumRows();
    space_.SetDims(space_.NumRows() + 1,
                   space_.NumCols());
    space_[old_dim] = other;
    auto rg = basic_linear_algebra::Rre(space_);
    space_.SetDims(rg, truncation_);
    return *this;
  }

  Uint64Subspace &operator+=(const uint64_t &other) {
    NTL::vec_GF2
        other_ntl = basic_linear_algebra::ConvertToNtl(other, truncation_);
    return operator+=(other_ntl);
  }

  template<class C>
  Uint64Subspace operator+(const C &other) const {
    auto result(*this);
    result += other;
    return result;
  }

  const NTL::mat_GF2 &parity_check_matrix() {
    if (parity_check_out_of_date) {
      auto pc = basic_linear_algebra::OrthogonalComplement(space_);
      saved_parity_check_matrix.swap(pc);
      parity_check_out_of_date = false;
    }
    return saved_parity_check_matrix;
  }

  //const for all intents and purposes.
  Uint64Subspace operator&(Uint64Subspace &other) {
    Uint64Subspace h1(parity_check_matrix(), truncation_);
    Uint64Subspace h2(other.parity_check_matrix(), truncation_);
    h1 += h2;
    return h1.OrthogonalComplement();
  }

  NTL::mat_GF2 space_;
};

#endif //CONDITIONS_SRC_UINT64SUBSPACE_HPP_
