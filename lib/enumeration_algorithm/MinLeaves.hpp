#ifndef CONDITIONS_SRC_ENUMERATION_ALGORITHM_MIN_LEAVES_HPP_
#define CONDITIONS_SRC_ENUMERATION_ALGORITHM_MIN_LEAVES_HPP_

#include <iostream>
#include <memory>
#include <optional>

#include "BruteForce.hpp"
#include "../Uint64Subspace.hpp"
#include "../VectorialBooleanFunction.hpp"

namespace MinLeaves {
typedef BruteForce::BinaryDecisionTreePointer BinaryDecisionTreePointer;
typedef BruteForce::BinaryDecisionTree BinaryDecisionTree;

BinaryDecisionTreePointer TreeSearch(const VectorialBooleanFunction &fun,
                                     const NTL::mat_GF2 &vectors_on_path,
                                     const NTL::vec_GF2 &choices,
                                     int last_choice,
                                     const NTL::vec_GF2 &coset_so_far,
                                     double &bound,
                                     int level, bool print= true);

BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun, bool print=true);
BinaryDecisionTreePointer StartSearch(const VectorialBooleanFunction &fun,
                                      double &bound, bool print=true);
BinaryDecisionTreePointer StartSearchWithFixedRoot(const VectorialBooleanFunction &fun,
                                                   double &bound,
                                                   uint64_t root, bool print=true);
BinaryDecisionTreePointer StartSearchWithFixedStump(const VectorialBooleanFunction &fun,
                                                   const BinaryDecisionTreePointer &stump,
                                                   double &bound,
                                                   uint64_t root, bool print=true);
void RecursiveStump(const VectorialBooleanFunction &fun,
                                         const BinaryDecisionTreePointer &stump,
                                         const NTL::mat_GF2 &vectors_on_path,
                                         const NTL::vec_GF2 &choices,
                                         int last_choice,
                                         const NTL::vec_GF2 &coset_so_far,
                                         double &bound,
                                         int level,
                                         bool print);


}

#endif //CONDITIONS_SRC_ENUMERATION_ALGORITHM_MIN_LEAVES_HPP_
