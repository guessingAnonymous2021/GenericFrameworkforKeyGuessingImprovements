#ifndef CONDITIONS_SRC_COMPONENTCHOICE_HPP_
#define CONDITIONS_SRC_COMPONENTCHOICE_HPP_

#include <bit>
#include <vector>
#include "VectorialBooleanFunction.hpp"

class ComponentChoice: public VectorialBooleanFunction {
 public:
  ComponentChoice(
      const std::vector<uint64_t> &components,
      const VectorialBooleanFunction &v) : VectorialBooleanFunction(v.GetValues(), v.InputSize(), v.OutputSize()) {
    std::vector<uint64_t> function(1u << v.InputSize());
    size_t i = 0;
    for (auto &x: components) {
      for (uint64_t inp = 0; inp < (1u << v.InputSize()); ++inp) {
        function[inp] ^= (std::popcount(x & v(inp)) & 1) << i;
      }
      i++;
    }
    this->GetValuesMutable() = function;
    this->OutputSizeMutable() = i;
  }
};


#endif //CONDITIONS_SRC_COMPONENTCHOICE_HPP_
