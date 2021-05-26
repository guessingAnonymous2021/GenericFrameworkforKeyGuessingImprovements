#ifndef CONDITIONS_SRC_ENUMERATION_ALGORITHM_EXHAUSTIVEBF_HPP_
#define CONDITIONS_SRC_ENUMERATION_ALGORITHM_EXHAUSTIVEBF_HPP_

#include "BruteForce.hpp"

namespace ExhaustiveBF {
using BinaryDecisionTree = BruteForce::BinaryDecisionTree;
using BinaryDecisionTreePointer = BruteForce::BinaryDecisionTreePointer;
typedef std::vector<BinaryDecisionTreePointer> BinaryDecisionTreeVector;

BinaryDecisionTreeVector TreeSearch(const VectorialBooleanFunction &fun,
                                    const NTL::mat_GF2 &vectors_on_path,
                                    const NTL::vec_GF2 &choices,
                                    int last_choice,
                                    const NTL::vec_GF2 &coset_so_far,
                                    double &bound,
                                    int level, bool print = true);

BinaryDecisionTreeVector StartSearch(const VectorialBooleanFunction &fun, bool print=true);
BinaryDecisionTreeVector StartSearch(const VectorialBooleanFunction &fun,
                                     double &bound, bool print=true);

};

#endif //CONDITIONS_SRC_ENUMERATION_ALGORITHM_EXHAUSTIVEBF_HPP_
