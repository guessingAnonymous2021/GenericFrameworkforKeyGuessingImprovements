#ifndef CONDITIONS_SRC_COMPONENTFUNCTION_HPP_
#define CONDITIONS_SRC_COMPONENTFUNCTION_HPP_
#include "VectorialBooleanFunction.hpp"
#include <bit>

class ComponentFunction : public VectorialBooleanFunction {
 public:
  ComponentFunction(const VectorialBooleanFunction &underlying_function,
                    uint64_t component) :
      VectorialBooleanFunction(underlying_function) {
    for (auto &x: GetValuesMutable()) {
      x = std::popcount(x & component) % 2;
    }
    TruncateOutputSize(1);
  }
};

#endif //CONDITIONS_SRC_COMPONENTFUNCTION_HPP_
