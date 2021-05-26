#ifndef CONDITIONS_VECTORIAL_BOOLEAN_FUNCTION_HPP_
#define CONDITIONS_VECTORIAL_BOOLEAN_FUNCTION_HPP_
#include <cstdint>
#include <vector>
#include <ostream>

class VectorialBooleanFunction {
 private:
  size_t input_size = 0, output_size = 0;
  std::vector<uint64_t> values;
  bool ddt_calculated = false;
 protected:
  std::vector<uint64_t> &GetValuesMutable();
  void TruncateOutputSize(size_t new_output_size);
 public:
  std::vector<std::vector<uint64_t>> ddt;
  bool operator==(const VectorialBooleanFunction &other) const;
  friend std::ostream &operator<<(std::ostream &os,
                                  const VectorialBooleanFunction &vec);
  VectorialBooleanFunction(const std::initializer_list<uint64_t> &x,
                           size_t input_size,
                           size_t output_size);

  template<class V>
  VectorialBooleanFunction(const V &x, size_t input_size,
                           size_t output_size);
  const std::vector<uint64_t> &GetValues() const;
  [[nodiscard]] size_t InputSize() const;
  [[nodiscard]] size_t OutputSize() const;
  uint64_t operator()(uint64_t x) const;

  void calculate_ddt() {
    if (ddt_calculated)
      return;
    ddt.resize(1u << input_size);
    for (auto &x: ddt) x.resize(1u << output_size);
    size_t beta;
    for (size_t i = 0; i < (1u << input_size); ++i) {
      for (size_t x = 0; x < (1u << input_size); ++x) {
        ++ddt[i][values[x] ^ values[x ^ i]];
      }
    }
    ddt_calculated = true;
  }

  size_t differential_uniformity() {
    calculate_ddt();
    size_t uni = 0;
    for (size_t i = 1; i < (1u << input_size); ++i) {
      for (size_t j = 0; j < (1u << output_size); ++j) {
        if (ddt[i][j] > uni) uni = ddt[i][j];
      }
    }
    return uni;
  }

  bool IsAPN() {
    calculate_ddt();
    bool apn = true;
    for (size_t i = 1; i < (1u << input_size) && apn; ++i) {
      for (size_t j = 0; j < (1u << output_size) && apn; ++j) {
        apn &= ddt[i][j] == 0 || ddt[i][j] == 2;
      }
    }
    return apn;
  }
  size_t &OutputSizeMutable();
};

template<class V>
VectorialBooleanFunction::VectorialBooleanFunction(const V &x,
                                                   size_t input_size,
                                                   size_t output_size) :
    values(std::move(x)), input_size(input_size),
    output_size(output_size) {
}

template<class RNG>
VectorialBooleanFunction RandomVBFWithWeight(uint64_t dim,
                                             uint64_t weight,
                                             RNG &rng) {
  std::vector<uint64_t> base;
  base.resize(1u << dim);
  for (size_t i = 0; i < weight; ++i) {
    base[i] = 1;
  }
  for (size_t i = weight; i < 1u << dim; ++i) {
    base[i] = 0;
  }
  std::shuffle(base.begin(), base.end(), rng);
  VectorialBooleanFunction fun(base, dim, 1);
  return fun;
}
template<class RNG>
VectorialBooleanFunction RandomBalancedVBF(uint64_t dim, RNG &rng) {
  return RandomVBFWithWeight(dim, 1u << (dim - 1u), rng);
}
template<class RNG>
VectorialBooleanFunction RandomVBF(uint64_t dim, RNG &rng) {
  std::vector<uint64_t> base;
  base.resize(1u << dim);
  for (size_t i = 0; i < 1u << dim; ++i) {
    base[i] = rng() % 2;
  }
  VectorialBooleanFunction fun(base, dim, 1);
  return fun;
}
#endif
