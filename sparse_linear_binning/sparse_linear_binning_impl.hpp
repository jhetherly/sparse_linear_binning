#include <cmath>
#include <limits>
#include <vector>

#include "boost/multiprecision/cpp_int.hpp"
#include "sparsepp/spp.h"

using boost::multiprecision::cpp_int;
using spp::sparse_hash_map;
using std::vector;


typedef unsigned long       FiniteNum_t;
typedef unsigned long long  BigFiniteNum_t;
typedef cpp_int             ArbNum_t;


// NOTE: types defining what precision to use
template<class GT, class CT>
struct finite_global_indices_D_fin {
  typedef CT Corner_t;
  static inline unsigned long convert_corner_to_ul (const Corner_t &num)
  {return static_cast<unsigned long>(num);}
  typedef GT GIndex_t;
  // typedef unordered_map<GIndex_t, double> Map_t;
  typedef sparse_hash_map<GIndex_t, double> Map_t;
  static inline unsigned long convert_to_ul (const GIndex_t &num)
  {return static_cast<unsigned long>(num);}
};

template<class GT>
struct finite_global_indices_D_arb {
  typedef ArbNum_t Corner_t;
  static inline unsigned long convert_corner_to_ul (const Corner_t &num)
  {return num.convert_to<unsigned long>();}
  // {return num.template convert_to<unsigned long>();}
  typedef GT GIndex_t;
  // typedef unordered_map<GIndex_t, double> Map_t;
  typedef sparse_hash_map<GIndex_t, double> Map_t;
  static inline unsigned long convert_to_ul (const GIndex_t &num)
  {return static_cast<unsigned long>(num);}
};

template<class CT>
struct arbitrary_precision_global_indices_D_fin {
  // NOTE: this hasher is really bad in terms of performance - a replacement
  //       should be developed
  struct cpp_int_stupid_hasher {
    std::size_t operator()(const cpp_int &num) const
    {
      std::size_t seed = 0;
      spp::hash_combine(seed, num.str());
      return seed;
    }
  };

  typedef CT Corner_t;
  static inline unsigned long convert_corner_to_ul (const Corner_t &num)
  {return static_cast<unsigned long>(num);}
  typedef ArbNum_t GIndex_t;
  // typedef unordered_map<GIndex_t, double, cpp_int_stupid_hasher> Map_t;
  typedef sparse_hash_map<GIndex_t, double, cpp_int_stupid_hasher> Map_t;
  static inline unsigned long convert_to_ul (const GIndex_t &num)
  {return num.convert_to<unsigned long>();}
  // {return num.template convert_to<unsigned long>();}
};

struct arbitrary_precision_global_indices_D_arb {
  // NOTE: this hasher is really bad in terms of performance - a replacement
  //       should be developed
  struct cpp_int_stupid_hasher {
    std::size_t operator()(const cpp_int &num) const
    {
      std::size_t seed = 0;
      spp::hash_combine(seed, num.str());
      return seed;
    }
  };

  typedef ArbNum_t Corner_t;
  static inline unsigned long convert_corner_to_ul (const Corner_t &num)
  {return num.convert_to<unsigned long>();}
  // {return num.template convert_to<unsigned long>();}
  typedef ArbNum_t GIndex_t;
  // typedef unordered_map<GIndex_t, double, cpp_int_stupid_hasher> Map_t;
  typedef sparse_hash_map<GIndex_t, double, cpp_int_stupid_hasher> Map_t;
  static inline unsigned long convert_to_ul (const GIndex_t &num)
  {return num.convert_to<unsigned long>();}
  // {return num.template convert_to<unsigned long>();}
};



// NOTE: adding unsigned integral types without arithmetic (really slow)
template<class T>
void bitwise_add_store_first (T &a, T b, T &carry)
{
  carry = a & b;
  a ^= b;
  while (carry > 0) {
    b = carry <<= 1;
    carry &= a;
    a ^= b;
  }
}


template<class T>
bool almost_equal(const T &x, const T &y, const int &ulp)
{
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
    || std::abs(x-y) < std::numeric_limits<T>::min();
}


template<class T>
inline bool almost_zero(const T &x)
{
  return (x <  std::numeric_limits<T>::epsilon() &&
          x > -std::numeric_limits<T>::epsilon());
}



template <class NumType>
struct StaticLoop {

  template <NumType N> struct uint_{ };

  template <NumType N, typename Lambda, typename IterT>
  static inline void unroller(const Lambda& f, const IterT& iter, uint_<N>) {
    unroller(f, iter, uint_<N-1>());
    f(iter + N);
  }

  template <typename Lambda, typename IterT>
  static inline void unroller(const Lambda& f, const IterT& iter, uint_<0>) {
    f(iter);
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void unroller_template(const Lambda& f, const IterT& iter, uint_<N>) {
    unroller_template(f, iter, uint_<N-1>());
    f.template operator()<N>(iter + N);
  }

  template <typename Lambda, typename IterT>
  static inline void unroller_template(const Lambda& f, const IterT& iter, uint_<0>) {
    f.template operator()<0>(iter);
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void unroller(Lambda& f, const IterT& iter, uint_<N>) {
    unroller(f, iter, uint_<N-1>());
    f(iter + N);
  }

  template <typename Lambda, typename IterT>
  static inline void unroller(Lambda& f, const IterT& iter, uint_<0>) {
    f(iter);
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void unroller_template(Lambda& f, const IterT& iter, uint_<N>) {
    unroller_template(f, iter, uint_<N-1>());
    f.template operator()<N>(iter + N);
  }

  template <typename Lambda, typename IterT>
  static inline void unroller_template(Lambda& f, const IterT& iter, uint_<0>) {
    f.template operator()<0>(iter);
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void loop (Lambda& f, const IterT& iter)
  {
    unroller(f, iter, uint_<N>());
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void loop (const Lambda& f, const IterT& iter)
  {
    unroller(f, iter, uint_<N>());
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void loop_static_num (Lambda& f, const IterT& iter)
  {
    unroller_template(f, iter, uint_<N>());
  }

  template <NumType N, typename Lambda, typename IterT>
  static inline void loop_static_num (const Lambda& f, const IterT& iter)
  {
    unroller_template(f, iter, uint_<N>());
  }

};




template<class GI>
typename GI::GIndex_t compute_global_bin(const vector<unsigned long> &bins,
                                         const vector<unsigned long> &sizes,
                                         const unsigned long &index)
{
  if (index == 0) return 0;

  typename GI::GIndex_t  sum = 0,
                         prod = 1;
  for (unsigned long i = 0; i < index; ++i) {
    if (i > 0) prod *= sizes[i - 1];
    sum += bins[i]*prod;
  }
  return sum;
}


// NOTE: class to compute the global bin number with looping done by the compiler
template<class GI>
struct ComputeGlobalBin {
  const vector<unsigned long>       *bins,
                                    *sizes;
  mutable typename GI::GIndex_t      sum, prod;

  void operator() (const unsigned long &i) const
  {
    if (i > 0) prod *= (*sizes)[i - 1];
    sum += (*bins)[i]*prod;
  }
};

template<unsigned long I, class GI>
struct GlobalBinComputer {
  static typename GI::GIndex_t result (const ComputeGlobalBin<GI> &cgb)
  {
    StaticLoop<unsigned long>::loop<I - 1>( cgb, 0 );
    return cgb.sum;
  }
};
template<class GI>
struct GlobalBinComputer<0, GI> {
  static typename GI::GIndex_t result (const ComputeGlobalBin<GI> &cgb)
  {
    return 0;
  }
};


template<class GI>
void compute_bin_numbers(const unsigned long &D,
                         const typename GI::GIndex_t &nbin,
                         const vector<unsigned long> &sizes, // D
                         vector<unsigned long> &bins) // D
{
  typename GI::GIndex_t denom = 1;

  for (unsigned long i = 0; i < D; ++i) {
    if (i > 0) denom *= sizes[i - 1];

    bins[i] = GI::convert_to_ul((nbin - compute_global_bin<GI>(bins, sizes, i))/denom)%sizes[i];
  }
}

// NOTE: class to compute the bin numbers with looping done by the compiler
template<class GI>
struct ComputeBinNumbers {
  const typename GI::GIndex_t *nbin;
  const vector<unsigned long> *sizes;
  mutable vector<unsigned long>       *bins;
  mutable typename GI::GIndex_t        denom;
  mutable ComputeGlobalBin<GI>         cgb;

  template<unsigned long I>
  void operator() (const unsigned long &i) const
  {
    if (I > 0) denom *= (*sizes)[I - 1];

    this->cgb.sum = 0;
    this->cgb.prod = 1;

    (*bins)[I] = GI::convert_to_ul(((*nbin) - GlobalBinComputer<I, GI>::result(this->cgb))/ denom)%(*sizes)[I];
  }
};

template<unsigned long D, class GI>
void compute_bin_numbers(const ComputeBinNumbers<GI> &cbn)
{
  StaticLoop<unsigned long>::loop_static_num<D - 1>( cbn, 0 );
}


// NOTE: class with the main logic of the algorithm
//        it's used more as a templated namespace than a genuine class and to
//        avoid unnecessary memory allocation
template<class GI>
class ComputeBinsAndVolumeFractionsDyn {
  const unsigned long   D;
  vector<unsigned long> max_sizes,
                        low_bins,
                        shifted_bins;
  vector<double>        max_lengths,
                        low_reminders;
  const typename GI::Corner_t   n_corners;
  typename GI::Corner_t         bad_corners,
                                max_good_corner;
  double                        basis_vol,
                                offset_sample,
                                fob,
                                temp_vol_frac,
                                temp_opposite_vol_frac;
  typename GI::Corner_t         opposite_corner,
                                min_corner;
  bool                          ccol,
                                occol;

  const vector<double>          &extents; // D x 2
  const vector<double>          &bin_sizes; // D
  const vector<unsigned long>   &sizes; // D
  mutable typename GI::Map_t            *result;

public:
  ComputeBinsAndVolumeFractionsDyn (const unsigned long             &d,
                                    const vector<double>            &Extents, // D x 2
                                    const vector<double>            &Bin_sizes, // D
                                    const vector<unsigned long>     &Sizes, // D
                                    typename GI::Map_t              &Result) :
    D(d),
    max_sizes(D), low_bins(D), shifted_bins(D),
    max_lengths(D), low_reminders(D),
    n_corners(static_cast<typename GI::Corner_t>(1) << D),
    bad_corners(0), max_good_corner(0),
    basis_vol(1.0),
    min_corner(0),
    extents(Extents), bin_sizes(Bin_sizes), sizes(Sizes),
    result(&Result)
  {
    for (unsigned long i = 0; i < D; ++i) {
      max_sizes[i] = sizes[i] - 1;
      max_lengths[i] = bin_sizes[i]*max_sizes[i];
    }
  }

  ComputeBinsAndVolumeFractionsDyn& operator() (const double *sample, // D
                                                const double &weight)
  {
    for (unsigned long i = 0; i < D; ++i) {
      offset_sample = sample[i] - extents[2*i];
      fob = floor(offset_sample/bin_sizes[i]);
      if (fob < 0)
        low_bins[i] = 0;
      else
        low_bins[i] = static_cast<unsigned long>(fob);
      if (low_bins[i] > max_sizes[i]) low_bins[i] = max_sizes[i];
      if (offset_sample < 0.0)
        low_reminders[i] = 0.0;
      else if (offset_sample > static_cast<double>(max_lengths[i]))
        low_reminders[i] = static_cast<double>(max_lengths[i]);
      else
        low_reminders[i] = offset_sample;
      low_reminders[i] -= static_cast<double>(low_bins[i])*bin_sizes[i];
      // NOTE: only need to check high end as the low end is already accounted for
      if (max_sizes[i] != low_bins[i]) {
        basis_vol *= bin_sizes[i];
        if (min_corner == 0)
          min_corner = static_cast<typename GI::Corner_t>(1) << i;
      }
      else
        bad_corners += static_cast<typename GI::Corner_t>(1) << i;
    }
    max_good_corner = (~bad_corners)%n_corners;
    basis_vol = 1.0/basis_vol;

    // NOTE: min_corner only improves iteration if lower dimensions are "bad"
    for (typename GI::Corner_t corner = 0; corner < n_corners/2; ++corner) {
      opposite_corner = (corner ^ max_good_corner)%n_corners;
      // skip if this corner shares any dimension with the bad dimensions or
      // if the opposite_corner is less than the current corner (avoids double
      // counting)
      if ((corner & bad_corners)%n_corners > 0 ||
          opposite_corner < corner) continue;
      temp_vol_frac          = basis_vol;
      temp_opposite_vol_frac = basis_vol;
      for (unsigned long i = 0; i < D; ++i) {
        if (!((bad_corners >> i) & 1)) {
          ccol = (corner >> i) & 1;
          occol = (opposite_corner >> i) & 1;
          // NOTE: "swapping" volume fractions to account for opposite corners
          temp_vol_frac           *= std::fabs(occol*bin_sizes[i] - low_reminders[i]);
          // NOTE: if statement for special case of sample beyond grid in all dimensions
          if (corner != opposite_corner)
            temp_opposite_vol_frac  *= std::fabs(ccol*bin_sizes[i]  - low_reminders[i]);
        }
      }
      if (!almost_zero(temp_vol_frac)) {
        for (unsigned long i = 0; i < D; ++i) {
          shifted_bins[i] = GI::convert_corner_to_ul((corner >> i) & 1) + low_bins[i];
        }
        (*result)[compute_global_bin<GI>(shifted_bins, sizes, D)] += weight*temp_vol_frac;
      }
      if (corner != opposite_corner && !almost_zero(temp_opposite_vol_frac)) {
        for (unsigned long i = 0; i < D; ++i) {
          shifted_bins[i] = GI::convert_corner_to_ul((opposite_corner >> i) & 1) + low_bins[i];
        }
        (*result)[compute_global_bin<GI>(shifted_bins, sizes, D)] += weight*temp_opposite_vol_frac;
      }
    }

    return *this;
  }

  void reset ()
  {
    basis_vol = 1.0;
    min_corner = 0;
    bad_corners = 0;
  }
};


template<class GI>
void make_linear_binning (const unsigned long D, // space dimension
                          const double *X, // n_samples x D
                          const double *weights, // n_samples
                          const unsigned long n_samples,
                          const vector<double> &extents, // D x 2
                          const vector<double> &bin_sizes, // D
                          const vector<unsigned long> &sizes, // D
                          double *&result_X, // grid points x D
                          double *&result_weights,  // grid points
                          unsigned long &n_result) // # of grid points
{
  /*
    Linear binning is used for down-sampling data while retaining much
    higher fidelity (in terms of asymptotic behavior) than nearest-neighbor
    binning (the usual type of binning).
    A-----------------------------------B
     |       |                         |
     |                                 |
     |       |                         |
     |- - - -P- - - - - - - - - - - - -|
     |       |                         |
    D-----------------------------------C
    For example, a 2D point P with weight wP:
    Assign a weight to corner A of the proportion of area (times wP)
        between P and C
    Assign a weight to corner B of the proportion of area (times wP)
        between P and D
    Assign a weight to corner C of the proportion of area (times wP)
        between P and A
    Assign a weight to corner D of the proportion of area (times wP)
        between P and B
   */

  typename GI::Map_t              result;

  ComputeBinsAndVolumeFractionsDyn<GI> cbavf(D, extents, bin_sizes, sizes, result);
  for (unsigned long i = 0; i < n_samples; ++i) {
    const double *sample = &X[i*D];
    const double weight = weights[i];

    cbavf(sample, weight).reset();
  }

  delete[] result_X;
  delete[] result_weights;
  n_result = result.size();
  result_X = new double[result.size()*D];
  result_weights = new double[result.size()];
  vector<unsigned long> _bins(D);
  size_t count = 0;
  for (typename GI::Map_t::iterator it = result.begin(); it != result.end(); ++it) {
    const unsigned long index = count*D;
    compute_bin_numbers<GI>(D, it->first, sizes, _bins);
    for (unsigned long i = 0; i < D; ++i)
      result_X[index + i] = _bins[i]*bin_sizes[i] + extents[i*2];
    result_weights[count] = it->second;
    ++count;
  }
}


// Driver function for the linear binning method
void sparse_linear_binning_impl (const unsigned long D, // space dimension
                                 const double *X, // n_samples x D
                                 const double *weights, // n_samples
                                 const unsigned long n_samples,
                                 const double *extents, // D x 2
                                 const double *bin_sizes, // D
                                 const unsigned long *sizes, // D
                                 double *&result_X, // grid points x D
                                 double *&result_weights,  // grid points
                                 unsigned long &n_result); // # of grid points
