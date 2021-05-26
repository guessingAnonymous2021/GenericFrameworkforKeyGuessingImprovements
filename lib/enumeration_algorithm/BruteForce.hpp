#ifndef CONDITIONS_SRC_ENUMERATION_ALGORITHM_BRUTEFORCE_HPP_
#define CONDITIONS_SRC_ENUMERATION_ALGORITHM_BRUTEFORCE_HPP_

#include <iostream>
#include <memory>
#include <optional>

#include "../Uint64Subspace.hpp"
#include "../VectorialBooleanFunction.hpp"

namespace BruteForce {

bool IsConstantOnSubspace(const VectorialBooleanFunction &fun,
                          const Uint64Subspace &space,
                          uint64_t coset);

class BinaryDecisionTree {
  typedef std::unique_ptr<BinaryDecisionTree> Subtree;
  static void print_spaces(std::ostream &os, int i);

  void recursive_print(std::ostream &os, int i);
 public:
  Subtree left_ = nullptr;
  Subtree right_ = nullptr;

  friend std::ostream &operator<<(std::ostream &os, BinaryDecisionTree &x) {
    os << "Binary tree: \n";
    x.recursive_print(os, 0);
    return os;
  }

  uint64_t value_ = 0;
  double associated_cost_ = 0;
  explicit BinaryDecisionTree() = default;;
  explicit BinaryDecisionTree(uint64_t value, double associated_cost = 0) {
    left_ = {};
    right_ = {};
    value_ = value;
    associated_cost_ = associated_cost;
  }
  BinaryDecisionTree(uint64_t value, std::unique_ptr<BinaryDecisionTree> left,
                     std::unique_ptr<BinaryDecisionTree> right) :
      value_(value), left_(std::move(left)),
      right_(std::move(right)) {
  }

  // Might drop some information, but not value_, the subtrees, and associated cost.
  [[nodiscard]] BinaryDecisionTree DeepImperfectCopy() const {
    if (IsLeaf()) {
      auto result = BinaryDecisionTree(value_);
      result.associated_cost_ = 0;
      return result;
    }
    else {
      auto result = BinaryDecisionTree(value_,
                                left_ ? std::make_unique<BinaryDecisionTree>(left_->DeepImperfectCopy()) : nullptr,
                                right_ ? std::make_unique<BinaryDecisionTree>(right_->DeepImperfectCopy()) : nullptr);
      result.associated_cost_ = associated_cost_;
      return result;
    }
  }

  [[nodiscard]] bool IsLeaf() const;

  [[nodiscard]] size_t Depth() const;

  [[nodiscard]] double AveragePathLength() const;

  [[nodiscard]] uint64_t EvaluateAt(uint64_t x) const;

  [[nodiscard]] size_t BitSize() const;
  [[nodiscard]] size_t OutputBitSize() const;

  [[nodiscard]] VectorialBooleanFunction UnderlyingFunction(size_t size_hint = 0,
                                                            size_t output_size_hint = 0) const;

  bool HasSameLabelsAs(const BinaryDecisionTree& other) {
    if (IsLeaf() != other.IsLeaf()) {
      return false;
    } else if (value_ == other.value_ && !IsLeaf()){
      return left_->HasSameLabelsAs(*other.left_)
          && right_->HasSameLabelsAs(*other.right_);
    } else if (value_ == other.value_) {
      return true;
    }
    return false;
  }
  int leaves() {
    if (IsLeaf()) return 1;
    int total = 0;
    if (right_) {
      total += right_->leaves();
    }
    if (left_) {
      total += left_->leaves();
    }
    return total;
  }

};

typedef std::unique_ptr<BinaryDecisionTree> BinaryDecisionTreePointer;
NTL::vec_GF2 ComplementIn(const NTL::mat_GF2 &parity_check_1,
                          const NTL::vec_GF2 &new_vector);

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

void analyse_component(const VectorialBooleanFunction &function,
                       uint64_t component);
}

#endif //CONDITIONS_SRC_ENUMERATION_ALGORITHM_BRUTEFORCE_HPP_
