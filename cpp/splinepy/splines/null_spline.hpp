#pragma once

#include <algorithm>
#include <array>
#include <memory>

#include <splinepy/splines/splinepy_base.hpp>

namespace splinepy::splines {

/// @brief Null Spline. Placeholder spline that only implements evaluate to
///        return zeros.
///
/// Only meant to be used in backend, where you want to use existing pipeline
/// but don't want to create memory consuming objects.
/// For example, if you have numerous list of
class NullSpline : public splinepy::splines::SplinepyBase {
public:
  /// parametric and physical dimension needs to be known at creation
  /// with these two, for whatever reason, if they are exposed to python side,
  /// it is still possible to call implemented query functions
  const int para_dim_;
  const int dim_;

  /// @brief ctor with parametric and physical dimension
  /// @param para_dim
  /// @param dim
  NullSpline(const int para_dim, const int dim)
      : para_dim_(para_dim),
        dim_(dim) {}

  /// basic implementations
  virtual int SplinepyParaDim() const { return para_dim_; }
  virtual int SplinepyDim() const { return dim_; }
  virtual std::string SplinepySplineName() const { return "NullSpline"; }
  virtual std::string SplinepyWhatAmI() const {
    return "NullSpline, parametric dimension: " + std::to_string(para_dim_)
           + ", physical dimension: " + std::to_string(dim_);
  }
  /// @brief Shouldn't ask NullSpline for has_knot_vectors, and is_rational
  virtual bool SplinepyHasKnotVectors() const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyHasKnotVectors() - invalid function call for",
        SplinepyWhatAmI());
    return false;
  }
  virtual bool SplinepyIsRational() const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyIsRational() - invalid function call for",
        SplinepyWhatAmI());
    return false;
  }
  virtual int SplinepyNumberOfControlPoints() const { return -1; }
  virtual int SplinepyNumberOfSupports() const { return -1; }
  virtual bool SplinepyIsNull() const { return true; }
  virtual void
  SplinepyCurrentProperties(int* degrees,
                            std::vector<std::vector<double>>* knot_vectors,
                            double* control_points,
                            double* weights) const {
    splinepy::utils::PrintAndThrowError(
        "SplinepyCurrentProperties() - invalid function call for",
        SplinepyWhatAmI());
  }

  /// Spline evaluation - fills zeros.
  virtual void SplinepyEvaluate(const double* para_coord,
                                double* evaluated) const {

    std::fill_n(evaluated, dim_, 0.0);
  }
};

// maxdim to pre-create null splines
#ifdef SPLINEPY_MORE
constexpr static int kMaxLookupDim{3};
#else
constexpr static int kMaxLookupDim{10};
#endif
using LookupArray_ =
    std::array<std::array<std::shared_ptr<SplinepyBase>, kMaxLookupDim>,
               kMaxLookupDim>;
/// pre-create null splines up to the dimension we support
static const LookupArray_ kNullSplineLookup = [] {
  LookupArray_ lookup;
  int i{0};
  for (auto& para_dim_array : lookup) {
    int j{0};
    for (auto& dim_element : para_dim_array) {
      dim_element = std::make_shared<NullSpline>(i + 1, j + 1);
      ++j;
    }
    ++i;
  }
  return lookup;
}();

} // namespace splinepy::splines
