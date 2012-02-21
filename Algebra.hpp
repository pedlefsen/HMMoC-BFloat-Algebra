/*
 *    This file is part of HMMoC-BFloat-Algebra, a tiny C++ library supplying
 *    the Algebra<> template, for polymorphic representation of non-negative
 *    floating-point values in either native (eg double), log-transformed, or
 *    BFloat representations.  BFloat is the Buoyant Float type, allowing
 *    fixed-precision representation at an arbitrary precision range.  The
 *    library is based on code from: Lunter G. HMMoC—a compiler for hidden
 *    Markov models. Bioinformatics (2007) 23(18): 2485-2487. The code was
 *    derived by Ian Holmes and Gerton Lunter from Holmes I, Rubin GM. An
 *    expectation maximization algorithm for training hidden substitution
 *    models.  J Mol Biol. 2002 Apr 12;317(5):753-64.  You may use at will,
 *    subject to the license (Apache v2.0), but *please cite those two papers* in
 *    your documentation and publications associated with uses of this library.
 *    Thank you!
 *
 *
 *    Copyright (C) 2011 by Paul T Edlefsen, Fred Hutchinson Cancer Research
 *    Center.
 *
 *    This file is based on algebras.h in HMMoC 1.3, a hidden Markov model
 *    compiler.  Copyright (C) 2007 by Gerton Lunter, Oxford University.
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *        http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
\*/
//
// Algebra.hpp - extended real types
//
// Paul T Edlefsen, June 1, 2011
// Gerton Lunter, 27/8/04
//

#ifndef HMMOC_BFLOAT_ALGEBRA_HPP
#define HMMOC_BFLOAT_ALGEBRA_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#if defined(_MSC_VER) && (_MSC_VER >= 1310)
#  pragma warning (disable : 4675) // suppress ADL warning
#endif

#include <assert.h>

// Do under- and overflow checking in BFloat? (undef for no)
#ifndef NDEBUG
#define BFLOAT_CHECK_UOFLOW
// Do extra checks in operator- to diagnose negative values? (set to 0 for no)
#define ALGEBRA_CHECK_NEGATIVE
#endif // NDEBUG

#ifndef NBOOST_SERIALIZATION // #define NBOOST_SERIALIZATION if you don't want to compile in the boost serialization library
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/version.hpp>
#endif // NBOOST_SERIALIZATION

// TODO: PUT BACK ?  See below (search for seqan)...
//#include <seqan/basic.h>

#include <iostream>
#include <cmath> // for isinf and isnan
#include <limits>

//TAH 11/01/2011  this defines uint32_t, used below
#include <stdint.h>


// You can control some inlining recommentations to the compiler using these:
#define ALGEBRA_INLINE_ARITHMETIC inline
#define ALGEBRA_INLINE_CAST inline
#define ALGEBRA_INLINE_CAST_COMPARE inline
#define ALGEBRA_INLINE_CAST_ARITHMETIC inline

// typedefs
typedef float BFMantissa;
//typedef double BFMantissa;
const BFMantissa cBFloatRange = 20282409603651670423947251286016.0;  // 2.03e+31; 2^104
const BFMantissa cBFloatRangeInv = 1.0/cBFloatRange;
// Aaron E. Darling 6/7/7: need to typecast to avoid compiler warnings about imprecise FP representations
const BFMantissa cBFloatRangeSqrt    = (BFMantissa)1.0e+18;          // Value between square root of the exponent, and the exponent
const BFMantissa cBFloatRangeInvSqrt = (BFMantissa)1.0e-18;          // Square of this should still be representable, with full mantissa!
const BFMantissa logcBFloatRange     = ::log(cBFloatRange);
const int cBFloatDigits              = 7;                 // Number of significant digits (7 for floats, 16 for doubles?)
const int cBFloatInfinity            = 1000000000;        // Tiniest number representable is cBFloatRangeInv ^ BFloatInfinity
// PTE 15/9/08
const BFMantissa cBFloatEpsilon      = ( BFMantissa )1.193e-07;
const int cBFloatConvTableSize       = 100;               // This includes many zero entries, it makes additions a bit faster
const int cBFloatDoubleConvTableSize = 50;                // Table size for bfloat -> double conversion; cBFloatRange^(-size/2) is double 0
const double cLogspaceEpsilon        = 36.7; // TODO: Test & verify

namespace hmmoc {

// PTE 14/9/08
// To help disambiguate ::log( double ) from hmmoc::log(Algebra)
template <typename other_type>
inline
other_type
log ( other_type const & to_be_logged )
{
  return ::log( to_be_logged );
} // log( other_type const & )

// Forward decls
class Logspace;
class DoubleRealspace;
class LongDoubleRealspace;
class FloatRealspace;
class BFloat;
static inline std::ostream& bfloat_print ( std::ostream& out, const BFloat& x );
static inline void assertBFloatValidity ( BFloat const & bf );
static inline void BFloatNormalise ( BFloat& a );

//
// BFloats: more buoyant floats.
//
// NOTE: For now BFloats can represent only positive values.
//
// struct{ float + int } is 8 bytes; nice size makes noticable speed difference
//
// For BFloats, infinity is represented as ( e == cBFloatInfinity, and we never let f be std::numeric_limits<BFMantissa>::infinity() ).  0 is represented as ( f == e == 0 ), and NaN as ( f != f ).
class BFloat {
// Boost serialization
#ifndef NBOOST_SERIALIZATION
private:
  friend class boost::serialization::access;
  template<class Archive>
  inline
  void serialize ( Archive & ar, const unsigned int /* file_version */ )
  {
    ar & BOOST_SERIALIZATION_NVP( f );
    ar & BOOST_SERIALIZATION_NVP( e );
  } // serialize( Archive &, const unsigned int )
#endif // NBOOST_SERIALIZATION

 public:
  static BFMantissa* aConversionLookup;              // used by addition
  static double* aDoubleConversionLookup;            // used by Value()
  BFMantissa f;
  int e;
 public:
  inline
  BFloat ( BFMantissa const & iF, int const & iE, bool normalise = false ) :
    f( iF ), e( iE )
  {
    if( normalise ) {
      BFloatNormalise( *this );
    } else {
      assertBFloatValidity( *this );
    }
  };

  inline
  BFloat () : f( 0.0 ), e( 0 ) {};

  inline
  BFloat ( double const & prob ) { *this = prob; }

  inline
  BFloat ( BFloat const & copy_from ) : f( copy_from.f ), e( copy_from.e ) {
    assertBFloatValidity( *this );
  };

  // Defined later, after these others are defined.
  inline
  BFloat ( Logspace const & copy_from );
  // Defined later, after these others are defined.
  inline
  BFloat ( DoubleRealspace const & copy_from );
  // Defined later, after these others are defined.
  inline
  BFloat ( LongDoubleRealspace const & copy_from );
  // Defined later, after these others are defined.
  inline
  BFloat ( FloatRealspace const & copy_from );

  inline
  ~BFloat () {};

  inline
  BFloat &
  operator= ( BFloat const & copy_from )
  {
    f = copy_from.f;
    e = copy_from.e;

    return *this;
  } // operator=( BFloat const & )

  // Defined later, after these others are defined.
  inline BFloat &
  operator= ( Logspace const & copy_from );
  // Defined later, after these others are defined.
  inline BFloat &
  operator= ( DoubleRealspace const & copy_from );
  // Defined later, after these others are defined.
  inline BFloat &
  operator= ( LongDoubleRealspace const & copy_from );
  // Defined later, after these others are defined.
  inline BFloat &
  operator= ( FloatRealspace const & copy_from );

  inline BFloat &
  operator= ( double prob )
  {
    if( prob != prob ) {
      // NaN
      // For now, don't let it be set to NaN.
      assert( prob == prob ); // NOT NaN
      f = prob;
      e = 0;
    } else
    if( prob == std::numeric_limits<double>::infinity() ) {
      f = 1.0;
      e = cBFloatInfinity;
    } else
    // Simplistic double-to-BFloat conversion - can be slow if 'standard' numbers get very large/small
    if( prob < std::numeric_limits<double>::min() ) {
#ifdef ALGEBRA_CHECK_NEGATIVE
      if( prob < 0.0 ) {
        std::cerr << "BFloat: Negative number: " << prob << std::endl;
        // TODO: ?
        assert( prob >= 0.0 );
      }
#endif
      f = 0.0;
      e = 0;
    } else {
      f = 0.0;
      e = 0;
      while( prob > cBFloatRangeSqrt ) {
        prob *= cBFloatRangeInv;
        e++;
      }
      while( prob < cBFloatRangeInvSqrt ) {
        prob *= cBFloatRange;
        e--;
      }
      f = prob; 
    }
    assertBFloatValidity( *this );
    return *this;
  } // operator=( double const & )

  inline double Value () const { 
    if( e == 0 ) {
      return f;
    }
    if( e == cBFloatInfinity ) {
      return std::numeric_limits<double>::infinity();
    }
    if( f != f ) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (abs(e) < cBFloatDoubleConvTableSize/2) {
      return (double)f * aDoubleConversionLookup[ e + cBFloatDoubleConvTableSize/2 ];
    } else if (e < cBFloatDoubleConvTableSize/2) {
      return 0.0;
    } else {
      return (double)f * exp((double)e * logcBFloatRange);
    }
  } // double Value()

  // sets to zero
  inline
  void clear () {
    f=0.0; e=0;
  }

  inline
  bool isZero () const
  {
    return ( f == 0.0 );
  } // isZero() const

  inline
  bool isOne () const
  {
    return ( f == 1.0 );
  } // isOne() const

  inline
  bool isInfinity () const
  {
    return ( e == cBFloatInfinity );
  } // isInfinity() const

  inline
  bool isNaN () const
  {
    return ( f != f );
  } // isNaN() const

}; // End class BFloat


//
// dummy class to initialise BFloat lookup table
//
class _BFloatInitialize {
public:
  _BFloatInitialize();
};


//
// implementations of BFloat calculations
//

static inline void assertBFloatValidity ( BFloat const & bf )
{
  // Never negative
  assert( bf.f >= 0.0 );

  // TODO: REMOVE?  NaN should be ok to represent, but we might want to raise a red flag...
  assert( !bf.isNaN() );

  // we represent infinity as ( e == cBFloatInfinity ), not as ( f == std::numeric_limits<BFMantissa>::infinity() ).
  assert( bf.f != std::numeric_limits<BFMantissa>::infinity() );

  if( bf.isZero() ) {
    // 0 is represented as both e and f being 0
    assert( bf.e == 0 );
    //} else if( !bf.isInfinity() ) {
    // f must be in this range (AFTER NORMALIZATION)
    //assert( bf.f <= cBFloatRangeSqrt ); 
    //assert( bf.f >= cBFloatRangeInvSqrt ); 
  }
    
  return;
} // static assertBFloatValidity( BFloat const & )

// Normalization of BFloat result of a single operation
#ifdef BFLOAT_CHECK_UOFLOW
static inline void BFloatNormalise ( BFloat& a )
     //#define BFloatNormalise(a)
{
  if (a.f > cBFloatRangeSqrt) {
    a.f *= cBFloatRangeInv;
    a.e++;
  } else if (a.f < cBFloatRangeInvSqrt) {
    if( a.f == 0.0 ) {
      a.e = 0; 
    } else {
      a.f *= cBFloatRange;
      a.e--;
    }
  }
  if( a.e > cBFloatInfinity ) {
    std::cerr << "BFloat: Overflow" << std::endl;
    a.e = cBFloatInfinity;
  } else if( a.e < -cBFloatInfinity ) {
    std::cerr << "BFloat: Underflow" << std::endl;
    a.e = 0;
    a.f = 0.0;
  }
  // TODO: REMOVE
#ifdef ALGEBRA_CHECK_NEGATIVE
    if( a.f < 0 ) {
      std::cerr << "BFloat: Negative number: ";
      bfloat_print( std::cerr, a );
      std::cerr << std::endl;
      std::cerr << "\tin BFloatNormalise( BFMantissa & ) [BFLOAT_CHECK_UOFLOW is true]" << std::endl;
    }
#endif // ALGEBRA_CHECK_NEGATIVE
#ifndef NDEBUG
  assertBFloatValidity( a );
#endif // NDEBUG
}; // BFloatNormalise ( BFloat & ) // when BFLOAT_CHECK_UOFLOW
#else // if BFLOAT_CHECK_UOFLOW .. else ..
static inline void BFloatNormDown ( BFloat& a ) { 
  a.f *= cBFloatRangeInv;
  a.e++;
  assert( a.f <= cBFloatRangeSqrt );
}
static inline void BFloatNormUp ( BFloat& a ) { 
  if (a.f == 0.0) {
    a.e = 0;
  } else {
    a.f *= cBFloatRange;
    a.e--;
    assert( a.f >= cBFloatRangeInvSqrt );
  }
}
static inline void BFloatNormalise ( BFloat& a )
     //#define BFloatNormalise(a) 
{
  if (a.f > cBFloatRangeSqrt) {
    BFloatNormDown(a);
  } else if (a.f < cBFloatRangeInvSqrt) {
    BFloatNormUp(a);
  }
#ifdef ALGEBRA_CHECK_NEGATIVE
    if( a.f < 0 ) {
      std::cerr << "BFloat: Negative number: ";
      bfloat_print( std::cerr, a );
      std::cerr << std::endl;
      std::cerr << "\tin BFloatNormalise( BFMantissa & ) [BFLOAT_CHECK_UOFLOW is false]" << std::endl;
    }
#endif // ALGEBRA_CHECK_NEGATIVE
#ifndef NDEBUG
  assertBFloatValidity( a );
#endif // NDEBUG
};  // BFloatNormalise ( BFloat & ) // when NOT BFLOAT_CHECK_UOFLOW
#endif // End if BFLOAT_CHECK_UOFLOW .. else ..

static inline void DoubleNormalise ( double& f, int& e )
{
  // comparing to 0.0 here fails, because the comparison is done
  // using higher-precision doubles, but the subsequent while-loop
  // uses true doubles, resulting in an infinite loop. (G.L. 3/9/07)
  if( f < std::numeric_limits<double>::min() ) {
#ifdef ALGEBRA_CHECK_NEGATIVE
    if( f < 0.0 ) {
      std::cerr << "BFloat: Negative number: " << f << std::endl;
      // TODO: ?
      assert( f >= 0.0 );
    }
#endif // ALGEBRA_CHECK_NEGATIVE
    f = 0.0; 
    e = 0;
  } else {
    while (f > (double)cBFloatRangeSqrt) {
      f *= (double)cBFloatRangeInv;
      e++;
    }
    while (f < (double)cBFloatRangeInvSqrt) {
      f *= (double)cBFloatRange;
      e--;
    }
  }
 }; // DoubleNormalise( double &, int & )

// Logarithm of a BFloat
static inline double bfloat_doublelog ( const BFloat& a ) { return a.e*logcBFloatRange+::log(a.f); }

// BFloat exp of a double -- could be a teensy tiny value.
static inline BFloat bfloat_doubleexp ( double iA ) 
{
  if( iA == std::numeric_limits<double>::infinity() ) {
    // exp( inf ) is inf
    return BFloat( 1.0, cBFloatInfinity );
  }
  if( iA == -std::numeric_limits<double>::infinity() ) {
    // exp( -inf ) is 0
    return 0.0;
  }
  int iE = (int)floor( iA / ::log(cBFloatRange) );
  iA -= iE * ::log(cBFloatRange);
  return BFloat( ::exp(iA), iE, true );
} // static bfloat_doubleexp( double )

// Returns a double value - or underflow/overflow if it does not fit.
static inline double bfloat2double ( const BFloat & bfloat) { return bfloat.Value(); }

/// TODO: PUT BACK
//template <typename T>
//inline static double
//toDouble ( T const & v )
//{
//  return seqan::convert<T, double>( v );
//} // toDouble( T const & )

inline static double
toDouble ( BFloat const & v )
{
  return bfloat2double( v );
} // toDouble( BFloat const & )

inline static long double
toLongDouble ( BFloat const & v )
{
  return bfloat2double( v );
} // toLongDouble( BFloat const & )

inline static double
toLogDouble ( BFloat const & v )
{
  return bfloat_doublelog( v );
} // toLogDouble( BFloat const & )

static inline BFloat double2bfloat ( double prob )
{
  return prob;
}

static inline BFloat bfloat_pr_product ( const BFloat& a, const BFloat& b )
{ 
  //assertBFloatValidity( a );
  //assertBFloatValidity( b );
  if( a.isZero() ) {
    return 0.0;
  }
  if( b.isZero() ) {
    return 0.0;
  }
  return BFloat( a.f*b.f, a.e+b.e, true );
} // static bfloat_pr_product( BFloat const &, BFloat const & )

static inline BFloat bfloat_pr_double_product ( const BFloat& a, double const & b )
{ 
  //assertBFloatValidity( a );
  if( a.isZero() ) {
    return 0.0;
  }
  if( b == 0.0 ) {
    return 0.0;
  }
  register double mantisse = a.f*b;
  int exponent = a.e;
  DoubleNormalise(mantisse, exponent);
  return BFloat( mantisse, exponent, false );
} // bfloat_pr_double_product ( BFloat const &, double const & )

static inline void bfloat_pr_product_accum ( BFloat& a, const BFloat& b )
{
  if( a.isZero() ) {
    // It's already 0.
    return;
  }
  if( b.isZero() ) {
    // Multiplying by 0..
    a.clear();
    return;
  }
  a.f *= b.f; a.e += b.e; 
  BFloatNormalise( a );
} // static bfloat_pr_product_accum( BFloat &, BFloat const & )

static inline void bfloat_pr_double_product_accum (  BFloat& a, double const & b )
{ 
  if( a.isZero() ) {
    // It's already 0.
    return;
  }
  if( b == 0.0 ) {
    // Multiplying by 0..
    a.clear();
    return;
  }
  register double mantisse = a.f*b;
  DoubleNormalise( mantisse, a.e );
  a.f = mantisse;
  assertBFloatValidity( a );
  return;
} // static bfloat_pr_double_product_accum( BFloat &, double const & )

// PTE 14/9/08
static inline BFloat bfloat_pr_power (const BFloat& a, const BFloat& b) 
{ 
  return BFloat( ::pow(a.f,b.f), (a.e*b.e), true );
}

static inline BFloat bfloat_pr_double_power ( const BFloat& a, double const & b )
{
  // TODO: test, make sure it's safe.  Probably it's not.
  register double mantisse = ::pow(a.f, b);
  int exponent = a.e;
  DoubleNormalise(mantisse, exponent);
  return BFloat(mantisse, exponent);
}

static inline BFloat bfloat_pr_quotient ( const BFloat& a, const BFloat& b )
{
  if( b.isZero() ) {
    return std::numeric_limits<double>::quiet_NaN();
  } else if( b.isInfinity() ) {
    return 0.0;
  }
  return BFloat( a.f/b.f, a.e-b.e, true );
} // bfloat_pr_product ( BFloat const &, BFloat const & )
  
static inline void bfloat_pr_quotient_accum ( BFloat& a, const BFloat& b ) 
{ 
  if( b.isZero() ) {
    a = std::numeric_limits<double>::quiet_NaN();
  } else if( b.isInfinity() ) {
    a = 0.0;
  } else {
    a.f /= b.f; 
    a.e -= b.e;
  }
  BFloatNormalise( a ); 
} // bfloat_pr_quotient_accum( BFloat &, BFloat const & )

static inline BFloat bfloat_pr_sum ( const BFloat& a, const BFloat& b )
{
  if( b.isZero() ) {
    assertBFloatValidity( a );
    return a;
  }
  if( a.isZero() ) {
    assertBFloatValidity( b );
    return b;
  }
  if( b.isInfinity() || a.isInfinity() ) {
    return std::numeric_limits<double>::infinity();
  }

  BFloat r;
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize) {
      assertBFloatValidity( a );
      return a;
    } else {
      r.f = a.f + b.f * BFloat::aConversionLookup[ a.e - b.e ];
      r.e = a.e;
      assertBFloatValidity( r );
      // TODO: REMOVE?
      //BFloatNormalise( r );
      return r;
    }
  } else {
    if (a.e <= b.e - cBFloatConvTableSize) {
      assertBFloatValidity( b );
      return b;
    } else {
      r.f = b.f + a.f * BFloat::aConversionLookup[ b.e - a.e ];
      r.e = b.e;
      // TODO: REMOVE?
      //BFloatNormalise( r );
      return r;
    }
  }
} // bfloat_pr_sum ( BFloat const &, BFloat const & )
 
static inline void bfloat_pr_sum_accum ( BFloat& a, const BFloat& b) 
{
  //assertBFloatValidity( a );
  //assertBFloatValidity( b );
  if( b.isZero() ) {
    // Adding 0 does nothing.
    return;
  }
  if( a.isInfinity() ) {
    // Adding to infinity does nothing.
    return;
  }
  if( a.isZero() ) {
    // Adding to 0 is an assignment.
    a.e = b.e;
    a.f = b.f;
    assertBFloatValidity( a );
    return;
  }
  if( b.isInfinity() ) {
    // Adding infinity makes a infinite.
    a.e = b.e;
    a.f = b.f;
    return;
  }
  if (a.e >= b.e) {
    if (a.e < b.e + cBFloatConvTableSize)
      a.f += b.f * BFloat::aConversionLookup[ a.e - b.e ];
  } else {
    if (a.e > b.e - cBFloatConvTableSize) {
      a.f = b.f + a.f * BFloat::aConversionLookup[ b.e - a.e ];
      a.e = b.e;
    } else {
      // Assign b to a.
      a.f = b.f;
      a.e = b.e;
    }
  }
  // TODO: REMOVE?
  //BFloatNormalise( a );
} // bfloat_pr_sum_accum( BFloat &, BFloat const & )

static inline void bfloat_print ( const BFloat& x ) 
{
  bfloat_print( std::cout, x );
}
static inline std::ostream& bfloat_print ( std::ostream& out, const BFloat& x ) 
{
  static const double log10 = ::log(10.0);
  static const double maxmantisse = 10.0 * (1.0 - 0.55 * exp(-cBFloatDigits * log10));
  //out.setf(ios::fixed,ios::floatfield);
  out.precision( cBFloatDigits );
  if( x.isInfinity() ) {
    out << 1.0 << "e+Inf";
  } else if( x.isNaN() ) {
    out << "NaN";
  } else if( x.isZero() ) {
    out << "0.0";
  } else if( x.isOne() ) {
    out << "1.0";
  } else if( x.e == 0 ) {
    out << x.f;
  } else {
    double iM = (::log(x.f) + logcBFloatRange*(double)x.e) / log10;
    //long iExp = long(floor(iM));
    // PTE 16/9/08
    double iExp = floor(iM);
    iM = exp((iM - iExp) * log10);
    if (iM > maxmantisse) {
      iExp += 1;
      iM = 1.0;
    }
    out << iM << ( iExp<0 ? "e" : "e+" ) << iExp;
  }
  //out.setf(ios::fixed,ios::floatfield);  // default
  // TODO: MAGIC #
  out.precision( 6 );           // default
  return out;
} // bfloat_print ( ostream &, BFloat const & )

// forward decl.
static inline bool bfloat_equal ( const BFloat& a, const BFloat& b );

// PTE 14/9/08
static inline BFloat bfloat_pr_diff (const BFloat& a, const BFloat& b) 
{
  if( b.isZero() ) {
    // subtracting zero does nothing.
    assertBFloatValidity( a );
    return a;
  } else if( a.isZero() ) {
    // a is 0.  Negative!
#ifdef ALGEBRA_CHECK_NEGATIVE
      std::cerr << "BFloat: Negative number" << std::endl;
      std::cerr << "\tin bfloat_pr_diff( zero, ";
      bfloat_print( std::cerr, b );
      std::cerr << " )" << std::endl;
#endif // ALGEBRA_CHECK_NEGATIVE
    return a; // 0
  } // End if a is 0
  if( ( a.e == b.e ) && ( a.f <= b.f ) ) {
#ifdef ALGEBRA_CHECK_NEGATIVE
    if(
      //( (a.e == b.e) && ( ( a.f - b.f ) <= -cBFloatEpsilon ) )
      ( (a.e == b.e) && ( ( a.f - b.f ) <= -(100.0f*cBFloatEpsilon) ) )
    ) {
      std::cerr << "BFloat: Negative number" << std::endl;
      std::cerr << "\tin bfloat_pr_diff( ";
      bfloat_print( std::cerr, a );
      std::cerr << ", ";
      bfloat_print( std::cerr, b );
      std::cerr << " )" << std::endl;
      std::cerr << "\t\t(a.e - b.e) is zero" << std::endl;
      std::cerr << "\t\t(a.f - b.f) is " << ( a.f - b.f ) << std::endl;
      return 0.0;
    } // End if it's negative ..
#endif // ALGEBRA_CHECK_NEGATIVE
    return 0.0;
  }
  if( (a.e > b.e) || ( (a.e == b.e) && ( ( a.f - b.f ) > -cBFloatEpsilon ) ) ) {
    // a positive result...
    if( a.e >= b.e + cBFloatConvTableSize ) {
      assertBFloatValidity( a );
      return a;
    } else {
      return BFloat( a.f - b.f * BFloat::aConversionLookup[ a.e - b.e ], a.e, true );
    }
  } else {
    // It could still be positive
    if( a.e > ( b.e - cBFloatConvTableSize ) ) {
      BFMantissa tmp_f = ( a.f * BFloat::aConversionLookup[ b.e - a.e ] ) - b.f;
      if( tmp_f > 0 ) {
        return BFloat( tmp_f, b.e, true );
      } // else negative
    }

    if( bfloat_equal( a, b ) ) {
      return 0.0;
    }

    // a negative result..
#ifdef ALGEBRA_CHECK_NEGATIVE
      std::cerr << "BFloat: Negative number" << std::endl;
      std::cerr << "\tin bfloat_pr_diff( ";
      bfloat_print( std::cerr, a );
      std::cerr << ", ";
      bfloat_print( std::cerr, b );
      std::cerr << " )" << std::endl;
      std::cerr << "\t\t(a.e - b.e) is " << ( a.e - b.e ) << std::endl;
      std::cerr << "\t\t(a.f - b.f) is " << ( a.f - b.f ) << std::endl;
#endif // ALGEBRA_CHECK_NEGATIVE
    return 0.0;

    //if (a.e <= b.e - cBFloatConvTableSize)
    //  return b;
    //else
    //  return BFloat( b.f - a.f * BFloat::aConversionLookup[ b.e - a.e ], b.e );
  } // End if positive .. else negative ..
} // bfloat_pr_diff ( const BFloat &, const BFloat & )
 
// PTE 14/9/08
static inline void bfloat_pr_diff_accum ( BFloat& a, const BFloat& b) 
{
  assertBFloatValidity( a );
  assertBFloatValidity( b );
  if( b.isZero() ) {
    // subtracting zero does nothing.
    return;
  }
  if( a.isZero() ) {
    // a is 0.  Negative!
#ifdef ALGEBRA_CHECK_NEGATIVE
      std::cerr << "BFloat: Negative number" << std::endl;
      std::cerr << "\tin bfloat_pr_diff_accum( zero, ";
      bfloat_print( std::cerr, b );
      std::cerr << " )" << std::endl;
#endif // ALGEBRA_CHECK_NEGATIVE
    return;
  } // End if a is 0
  if( ( a.e == b.e ) && ( a.f <= b.f ) ) {
#ifdef ALGEBRA_CHECK_NEGATIVE
    if(
      ( ( a.f - b.f ) <= -cBFloatEpsilon )
    ) {
      std::cerr << "BFloat: Negative number" << std::endl;
      std::cerr << "\tin bfloat_pr_diff_accum( ";
      bfloat_print( std::cerr, a );
      std::cerr << " [ a.f = " << a.f << ", a.e = " << a.e << " ]";
      std::cerr << ", ";
      bfloat_print( std::cerr, b );
      std::cerr << " [ b.f = " << b.f << ", b.e = " << b.e << " ]";
      std::cerr << " )" << std::endl;
      a.f = 0.0;
      a.e = 0;
      return;
    } // End if it's negative.
#endif // ALGEBRA_CHECK_NEGATIVE
    a.f = 0.0;
    a.e = 0;
    return;
  }
  if( (a.e > b.e) || ( (a.e == b.e) && ( ( a.f - b.f ) > -cBFloatEpsilon ) ) ) {
    // a positive result...
    if( a.e < b.e + cBFloatConvTableSize ) {
      a.f -= b.f * BFloat::aConversionLookup[ a.e - b.e ];
      BFloatNormalise( a );
    }
  } else {
    // It could still be positive
    if( a.e > b.e - cBFloatConvTableSize ) {
      if( a.f * BFloat::aConversionLookup[ b.e - a.e ] > b.f ) {
        a.f = a.f * BFloat::aConversionLookup[ b.e - a.e ] - b.f;
        a.e = b.e;
        BFloatNormalise( a );
        return;
      } // else negative..
    }

    // a negative result..
#ifdef ALGEBRA_CHECK_NEGATIVE
      std::cerr << "BFloat: Negative number" << std::endl;
      std::cerr << "\tin bfloat_pr_diff_accum( ";
      bfloat_print( std::cerr, a );
      std::cerr << " [ a.f = " << a.f << ", a.e = " << a.e << " ]";
      std::cerr << ", ";
      bfloat_print( std::cerr, b );
      std::cerr << " [ b.f = " << b.f << ", b.e = " << b.e << " ]";
      std::cerr << " )" << std::endl;
#endif // ALGEBRA_CHECK_NEGATIVE
    a.f = 0.0;
    a.e = 0;
  } // End if positive .. else negative ..
} // bfloat_pr_diff_accum ( BFloat &, const BFloat & )

static inline bool bfloat_less ( const BFloat& a, const BFloat& b) 
{
  // ASSUMPTION: a and b are non-negative.
  assertBFloatValidity( a );
  assertBFloatValidity( b );
  if( a.isZero() ) {
    return ( !b.isZero() );
  }
  if( b.isZero() ) {
    return false;
  }
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return false;
    else
      return a.f < b.f * BFloat::aConversionLookup[ a.e - b.e ];
  }
  if (a.e <= b.e - cBFloatConvTableSize)
    return true;
  else
    return a.f * BFloat::aConversionLookup[ b.e - a.e ] < b.f;
} // static bfloat_less ( BFloat const &, BFloat const & )
  
static inline bool bfloat_equal ( const BFloat& a, const BFloat& b )
{
  assertBFloatValidity( a );
  assertBFloatValidity( b );
  if( a.isZero() ) {
    return ( b.isZero() );
  }
  if( b.isZero() ) {
    return false;
  }
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return false;
    else
      return a.f == b.f * BFloat::aConversionLookup[ a.e - b.e ];
  }
  if (a.e <= b.e - cBFloatConvTableSize)
    return false;
  else
    return a.f * BFloat::aConversionLookup[ b.e - a.e ] == b.f;
} // static bfloat_equal( BFloat const &, BFloat const & )

static inline bool bfloat_lessequal ( const BFloat& a, const BFloat& b )
{
  assertBFloatValidity( a );
  assertBFloatValidity( b );
  if( a.isZero() ) {
    return true;
  }
  if( b.isZero() ) {
    return false;
  }
  if (a.e > b.e) {
    if (a.e >= b.e + cBFloatConvTableSize)
      return false;
    else
      return a.f <= b.f * BFloat::aConversionLookup[ a.e - b.e ];
  }
  if (a.e <= b.e - cBFloatConvTableSize)
    return true;
  else
    return a.f * BFloat::aConversionLookup[ b.e - a.e ] <= b.f;
} // static bfloat_lessequal( BFloat const &, BFloat const & )


//
// Wrapper to allow BFloats to be used by Algebra template
//
struct BFloatMethods
{
  typedef BFloat Value;
  static inline double to_double (BFloat const & iX) { return bfloat2double(iX); }
  static inline BFloat from_double (double iP) { return double2bfloat(iP); }
  static inline BFloat pmul ( BFloat const & iX, BFloat const & iY) { return bfloat_pr_product(iX,iY); }
  static inline BFloat pmuldouble ( BFloat const & iX, double iY) { return bfloat_pr_double_product(iX,iY); }
  // PTE 14/9/08
  static inline BFloat ppow ( BFloat const & iX, BFloat const & iY) { return bfloat_pr_power(iX,iY); }
  // PTE 14/9/08
  static inline BFloat ppowdouble ( BFloat const & iX, double iY) { return bfloat_pr_double_power(iX,iY); }
  static inline BFloat pdiv ( BFloat const & iX, BFloat const & iY) { return bfloat_pr_quotient(iX,iY); }
  static inline BFloat psum ( BFloat const & iX, BFloat const & iY) { return bfloat_pr_sum(iX,iY); }
  static inline BFloat pdiff ( BFloat const & iX, BFloat const & iY) { return bfloat_pr_diff(iX,iY); }
  static inline BFloat doubleexp ( double iX ) { return bfloat_doubleexp(iX); }
  static inline double doublelog ( BFloat const & iX ) { return bfloat_doublelog(iX); }
  static inline void pmulacc ( BFloat& iX, BFloat const & iY ) { bfloat_pr_product_accum(iX,iY); }
  static inline void pmulaccdouble ( BFloat& iX, double iY ) { bfloat_pr_double_product_accum(iX,iY); }
  static inline void pdivacc ( BFloat& iX, BFloat const & iY ) { bfloat_pr_quotient_accum(iX,iY); }
  static inline void psumacc ( BFloat& iX, BFloat const & iY ) { bfloat_pr_sum_accum(iX,iY); }
  static inline void pdiffacc ( BFloat& iX, BFloat const & iY ) { bfloat_pr_diff_accum(iX,iY); }
  static inline bool less ( BFloat const & iX, BFloat const & iY ) { return bfloat_less(iX,iY); }
  static inline bool equal ( BFloat const & iX, BFloat const & iY ) { return bfloat_equal(iX,iY); }
  static inline bool lessequal ( BFloat const & iX, BFloat const & iY ) { return bfloat_lessequal(iX,iY); }
  static inline std::ostream& print ( std::ostream& iOut, BFloat const & iX ) { return bfloat_print( iOut, iX ); }
};


//
// Wrapper to use Algebra with a double or a float.
//
// PTE 15/9/08
template <typename RealspaceType>
struct RealspaceMethods
{
  typedef RealspaceType Value;
  static inline double to_double ( Value const & iX ) { return (double)iX.x; }
  static inline Value from_double ( double const & iP ) { return Value(iP); }
  static inline Value pmul ( Value const & iX, Value const & iY ) { return iX.x*iY.x; }
  static inline Value pmuldouble ( Value const & iX, double const & iY ) { return iX.x*iY; }
  // PTE 14/9/08
  static inline Value ppow ( Value const & iX, Value const & iY ) { return ::pow(iX.x,iY.x); }
  // PTE 14/9/08
  static inline Value ppowdouble ( Value const & iX, double const & iY ) { return ::pow(iX.x,iY); }
  static inline Value pdiv ( Value const & iX, Value const & iY ) { return iX.x/iY.x; }
  static inline Value psum ( Value const & iX, Value const & iY ) { return iX.x+iY.x; }
  static inline Value pdiff ( Value const & iX, Value const & iY ) { return iX.x-iY.x; }
  static inline Value doubleexp ( double const & iX ) { return Value(exp(iX)); }
  static inline double doublelog ( Value const & iX ) { return (double)log(iX.x); }
  static inline void pmulacc ( Value& iX, Value const & iY ) { iX.x*=iY.x; }
  static inline void pmulaccdouble ( Value& iX, double const & iY ) { iX.x*=static_cast<double>(Value(iY)); }
  static inline void pdivacc ( Value& iX, Value const & iY ) { iX.x /= iY.x; }
  static inline void psumacc ( Value& iX, Value const & iY ) { iX.x += iY.x; }
  static inline void pdiffacc ( Value& iX, Value const & iY ) { iX.x -= iY.x; }
  static inline bool less ( Value const & iX, Value const & iY ) { return ( !equal( iX, iY ) && (iX.x<iY.x) ); } //{ return iX.x<iY.x; }
  // TODO: Implement a numeric_limits<DoubleRealspace> that inherits from numeric_limits<double>, and likewise for FloatRealspace, etc.
  //static inline bool equal ( Value const & iX, Value const & iY ) { return ( ( iX.x == iY.x ) || ( ( iX.x > iY.x ) ? ( ( iX.x-iY.x ) <= std::numeric_limits<RealspaceType>::epsilon() ) : ( ( iY.x-iX.x ) <= std::numeric_limits<RealspaceType>::epsilon() ) ) ); } //{ return iX.x==iY.x; }
  //static inline bool equal ( Value const & iX, Value const & iY ) { return ( ( iX.x == iY.x ) || ( ( iX.x > iY.x ) ? ( ( iX.x-iY.x ) <= std::numeric_limits<float>::epsilon() ) : ( ( iY.x-iX.x ) <= std::numeric_limits<float>::epsilon() ) ) ); } //{ return iX.x==iY.x; }
  static inline bool equal ( Value const & iX, Value const & iY ) { return ( ( iX.x == iY.x ) || ( ( iX.x > iY.x ) ? ( ( iX.x-iY.x ) <= std::numeric_limits<typename Value::XType>::epsilon() ) : ( ( iY.x-iX.x ) <= std::numeric_limits<typename Value::XType>::epsilon() ) ) ); } //{ return iX.x==iY.x; }
  static inline bool lessequal ( Value const & iX, Value const & iY ) { return iX.x<=iY.x; }
  static inline std::ostream& print ( std::ostream& iOut, Value const & iX ) { if( iX.x < 0 ) { iOut << '-'; return bfloat_print( iOut, double2bfloat( -iX.x ) ); } return bfloat_print( iOut, double2bfloat(iX.x) ); }
}; // End class RealspaceMethods

// PTE 15/9/08
class DoubleRealspace {
// Boost serialization
#ifndef NBOOST_SERIALIZATION
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize ( Archive & ar, const unsigned int /* file_version */ )
  {
    ar & BOOST_SERIALIZATION_NVP( x );
  } // serialize( Archive &, const unsigned int )
#endif // NBOOST_SERIALIZATION

 public:
  typedef double XType;
  double x;

 public:
  DoubleRealspace () {}
  DoubleRealspace ( double x ) : x (x) {}
  template <typename other_type>
  DoubleRealspace ( other_type const & from_val )
  {
    x = toDouble( from_val );
  }
  //operator double& () { return x; }
  operator double () const { return x; }
  void clear () {
    x = 0.0;
  }
  DoubleRealspace &
  operator= ( double const & d )
  {
    x = d;
    return *this;
  } // operator=( double )
  template <typename other_type>
  DoubleRealspace &
  operator= ( other_type const & other_val )
  {
    x = toDouble( other_val );

    return *this;
  } // operator= ( other_type const & )
  
  inline double Value () const {
    return x;
  } // double Value()

  bool isZero () const
  {
    return ( x == 0.0 );
  } // isZero() const

  bool isOne () const
  {
    return ( x == 1.0 );
  } // isOne() const

  bool isInfinity () const
  {
    return std::isinf( x );
  } // isInfinity() const

  bool isNaN () const
  {
    return ( x != x );
  } // isNaN() const

}; // End class DoubleRealspace

// Added by Paul T Edlefsen 15/9/08
inline static double
toDouble ( DoubleRealspace const & v )
{
  return v.x;
} // toDouble( DoubleRealspace const & )

inline static long double
toLongDouble ( DoubleRealspace const & v )
{
  return v.x;
} // toLongDouble( DoubleRealspace const & )

inline static double
toLogDouble ( DoubleRealspace const & v )
{
  return log( v.x );
} // toLogDouble( DoubleReaspace const & )

///// mark
// PTE 15/9/08
class LongDoubleRealspace {
// Boost serialization
#ifndef NBOOST_SERIALIZATION
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize ( Archive & ar, const unsigned int /* file_version */ )
  {
    ar & BOOST_SERIALIZATION_NVP( x );
  } // serialize( Archive &, const unsigned int )
#endif // NBOOST_SERIALIZATION

 public:
  typedef long double XType;
  long double x;

 public:
  LongDoubleRealspace () {}
  LongDoubleRealspace ( double x ) : x (x) {}
  LongDoubleRealspace ( long double x ) : x (x) {}
  LongDoubleRealspace ( BFloat const & from_val )
  {
    x = toLongDouble( from_val );
  }
  template <typename other_type>
  LongDoubleRealspace ( other_type const & from_val )
  {
    x = toLongDouble( from_val );
  }
  //operator double& () { return x; }
  operator double () const { return static_cast<double>( x ); }
  operator long double () const { return x; }
  void clear () {
    x = 0.0;
  }
  LongDoubleRealspace &
  operator= ( long double const & d )
  {
    x = d;
    return *this;
  } // operator=( long double )
  template <typename other_type>
  LongDoubleRealspace &
  operator= ( other_type const & other_val )
  {
    x = toLongDouble( other_val );

    return *this;
  } // operator= ( other_type const & )
  
  inline long double Value () const {
    return x;
  } // long double Value()

  bool isZero () const
  {
    return ( x == static_cast<long double>( 0.0 ) );
  } // isZero() const

  bool isOne () const
  {
    return ( x == static_cast<long double>( 1.0 ) );
  } // isOne() const

  bool isInfinity () const
  {
    return std::isinf( x );
  } // isInfinity() const

  bool isNaN () const
  {
    return ( x != x );
  } // isNaN() const

}; // End class LongDoubleRealspace

// Added by Paul T Edlefsen 15/9/08
inline static double
toDouble ( LongDoubleRealspace const & v )
{
  return static_cast<double>( v.x );
} // toDouble( LongDoubleRealspace const & )

inline static double
toLongDouble ( LongDoubleRealspace const & v )
{
  return v.x;
} // toDouble( LongDoubleRealspace const & )

inline static double
toLogDouble ( LongDoubleRealspace const & v )
{
  return log( v.x );
} // toLogDouble( DoubleReaspace const & )
///// endmark

///////////////////////
class FloatRealspace {
// Boost serialization
#ifndef NBOOST_SERIALIZATION
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize ( Archive & ar, const unsigned int /* file_version */ )
  {
    ar & BOOST_SERIALIZATION_NVP( x );
  } // serialize( Archive &, const unsigned int )
#endif // NBOOST_SERIALIZATION

 public:
  typedef float XType;
  float x;

 public:
  FloatRealspace () {}
  FloatRealspace ( float x ) : x (x) {}
  FloatRealspace ( double x ) : x (x) {}
  template <typename other_type>
  FloatRealspace ( other_type const & from_val )
  {
    x = toDouble( from_val );
  }
  //operator double& () { return x; }
  operator double () const { return x; }
  void clear () {
    x = 0.0;
  }
  FloatRealspace &
  operator= ( float const & d )
  {
    x = d;
    return *this;
  } // operator=( float )
  FloatRealspace &
  operator= ( double const & d )
  {
    x = d;
    return *this;
  } // operator=( double )
  template <typename other_type>
  FloatRealspace &
  operator= ( other_type const & other_val )
  {
    x = toDouble( other_val );

    return *this;
  } // operator= ( other_type const & )
  
  inline double Value () const {
    return x;
  } // double Value()

  bool isZero () const
  {
    return ( x == 0.0f );
  } // isZero() const

  bool isOne () const
  {
    return ( x == 1.0f );
  } // isOne() const

  bool isInfinity () const
  {
    return std::isinf( x );
  } // isInfinity() const

  bool isNaN () const
  {
    return ( x != x );
  } // isNaN() const

}; // End class FloatRealspace

// Added by Paul T Edlefsen 15/9/08
inline static double
toDouble ( FloatRealspace const & v )
{
  return v.x;
} // toDouble( FloatRealspace const & )

inline static long double
toLongDouble ( FloatRealspace const & v )
{
  return v.x;
} // toLongDouble( FloatRealspace const & )

inline static double
toLogDouble ( FloatRealspace const & v )
{
  return log( v.x );
} // toLogDouble( FloatReaspace const & )

///////////////////////
//
// Simple log-space numbers - don't use, except possibly for Viterbi
//
template <typename LogspaceType>
struct LogspaceMethods
{
  typedef LogspaceType Value;
  static inline double to_double ( Value const & iX ) {
    if( iX.isInfinity() ) {
      return std::numeric_limits<double>::infinity();
    } else if( iX.isNaN() ) {
      return std::numeric_limits<double>::quiet_NaN();
    } else if( iX.isZero() ) {
      // log( 0 ) = -inf
      return 0.0;
    } // else {
    return ::exp(iX.x);
    //}
  } // to_double( Value const & )

  static inline Value from_double ( double const & iP ) {
    if( iP == std::numeric_limits<double>::infinity() ) {
      return Value( iP );
    } else if( iP != iP /*NaN*/ ) {
      // ?
      assert( iP == iP );
      return Value( iP );
    } else if( iP == 0.0 ) {
      // log( 0 ) = -inf
      return Value( -std::numeric_limits<double>::infinity() );
    } else {
      return Value( ::log( iP ) );
    }
  } // from_double( double const & )

  static inline Value pmul ( Value const & iX, Value const & iY ) { return iX.x+iY.x; }
  static inline Value pmuldouble ( Value const & iX, double const & iY) { return iX.x+::log(iY); }
  // PTE 14/9/08
  static inline Value ppow ( Value const & iX, Value const & iY) { return iX.x*iY.x; }
  // PTE 14/9/08
  static inline Value ppowdouble ( Value const & iX, double const & iY) { return iX.x*::log(iY); }
  static inline Value pdiv ( Value const & iX, Value const & iY) { return iX.x-iY.x; }
  static inline Value psum ( Value const & iX, Value const & iY) { return logspace_add(iX.x,iY.x); }
  static inline Value pdiff ( Value const & iX, Value const & iY) { return logspace_diff(iX.x,iY.x); }
  static inline Value doubleexp ( double const & iX ) { return iX; }
  static inline double doublelog ( Value const & iX ) { return iX.x; }
  static inline void pmulacc ( Value& iX, Value const & iY) { iX.x+=iY.x; }
  static inline void pmulaccdouble ( Value& iX, double const & iY) { iX.x+=::log(iY); }
  static inline void pdivacc ( Value& iX, Value const & iY) { iX.x -= iY.x; }
  static inline void psumacc ( Value& iX, Value const & iY) { iX.x = logspace_add(iX.x,iY.x); }
  static inline void pdiffacc ( Value& iX, Value const & iY) { iX.x = logspace_diff(iX.x,iY.x); }
  static inline bool less ( Value const & iX, Value const & iY) { return iX.x<iY.x; }
  static inline bool equal ( Value const & iX, Value const & iY) { return iX.x==iY.x; }
  static inline bool lessequal ( Value const & iX, Value const & iY) { return iX.x<=iY.x; }
  static inline std::ostream& print ( std::ostream& iOut, Value const & iX ) {
    if( ( iX.x != iX.x ) /* NaN */ ) {
      iOut << iX.x;
      return iOut;
    } else if( iX.x == std::numeric_limits<double>::infinity() ) {
      iOut << iX.x;
      return iOut;
    } else if( iX.x == -std::numeric_limits<double>::infinity() ) {
      // -inf is log( 0 )
      iOut << "0.0";
      return iOut;
    }
    //iOut << iX.x;
    //return iOut;
    return bfloat_print( iOut, bfloat_doubleexp(iX.x) );
  }
};

class Logspace {
// Boost serialization
#ifndef NBOOST_SERIALIZATION
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize ( Archive & ar, const unsigned int /* file_version */ )
  {
    ar & BOOST_SERIALIZATION_NVP( x );
  } // serialize( Archive &, const unsigned int )
#endif // NBOOST_SERIALIZATION

  // PTE 14/9/08
 public:
  double x;

 public:
  Logspace () {}
  // TODO: ERE I AM.  Can't do this.  Need to expect arg double is already logged.
  Logspace ( double const & other_val ) {
    //operator=( other_val );
    x = other_val;
  }
  Logspace ( Logspace const & other_logspace )
  {
    x = other_logspace.x;
  }
  Logspace ( LongDoubleRealspace const & other_val )
  {
    operator=( other_val );
  } // <init>( LongDoubleRealspace const & )
  Logspace ( FloatRealspace const & other_val )
  {
    operator=( other_val );
  } // <init>( FloatRealspace const & )
  Logspace ( BFloat const & other_val )
  {
    operator=( other_val );
  } // <init>( BFloat const & )
  //operator double&(){ return x; }
  // PTE 15/9/08
  //operator double () const { return x; } //return ::exp(x); }
  void clear() {
    // Set to log( 0 ).
    x = -std::numeric_limits<double>::infinity();
  }

  Logspace &
  operator= ( Logspace const & other_logspace )
  {
    x = other_logspace.x;

    return *this;
  } // operator= ( Logspace const & )
  Logspace &
  operator= ( LongDoubleRealspace const & other_val )
  {
    return operator=( toDouble( other_val ) );
  } // operator= ( LongDoubleRealspace const & )
  Logspace &
  operator= ( FloatRealspace const & other_val )
  {
    return operator=( toDouble( other_val ) );
  } // operator= ( FloatRealspace const & )
  Logspace &
  operator= ( BFloat const & other_val )
  {
    if( other_val.isZero() ) {
      x = -std::numeric_limits<double>::infinity();
    } else if( other_val.isInfinity() ) {
      x = std::numeric_limits<double>::infinity();
    } else {
      x = toLogDouble( other_val );
    }
    return *this;
  } // operator= ( BFloat const & )
  Logspace &
  operator= ( double const & other_val )
  {
    if( other_val == std::numeric_limits<double>::infinity() ) {
      x = other_val;
    } else if( other_val != other_val /* NaN */ ) {
      // ?
      assert( other_val == other_val );
      x = other_val;
    } else if( other_val == 0.0 ) {
      // log( 0 ) = -inf
      x = -std::numeric_limits<double>::infinity();
    } else {
      x = ::log( other_val );
    }

    return *this;
  } // operator= ( double const & )

  bool isZero () const
  {
    return ( x == -std::numeric_limits<double>::infinity() );
  } // isZero () const

  bool isOne () const
  {
    return ( x == 0 );
  } // isOne () const

  bool isInfinity () const
  {
    return ( x == std::numeric_limits<double>::infinity() );
  } // isInfinity () const

  bool isNaN () const
  {
    return ( x != x );
  } // isNaN () const

  //friend double logspace_add( double iX, double iY );
  //friend double logspace_diff( double iX, double iY );
  //friend bool LogspaceMethods<Logspace>::less( Logspace iX, Logspace iY);
}; // End class Logspace

// OLD.
//inline Logspace logspace_addsmall( Logspace iX, Logspace iY ) {
//  if (iX.x - iY.x > 36.7) return iX;
//  return Logspace( iX.x + ::log(1.0+exp(iY.x-iX.x)) );
//}

// Paul changed from "iX>iY" to "iX.x>iY.x" to work around a g++ "i686-apple-darwin8-g++-4.0.1 (GCC) 4.0.1 (Apple Computer, Inc. build 5367)" complaint of an "ambiguous overload for operator>".
// PTE 15/9/08
inline double
logspace_add ( double lhs, double const & rhs )
{
  static const double log2 = ::log( 2.0 );
  if( lhs == rhs ) {
    lhs += log2;
  } else if( lhs > rhs ) {
    lhs += ::log1p(::exp(rhs-lhs));
  } else if( lhs < rhs ) {
    lhs = rhs + ::log1p(::exp(lhs-rhs));
  }
  return lhs;
} // logspace_add( double, double const & )

inline double
logspace_diff ( double lhs, double const & rhs )
{
  if( lhs <= rhs ) {
#ifdef ALGEBRA_CHECK_NEGATIVE
    if(
      ( ( lhs - rhs ) < -cLogspaceEpsilon )
    ) {
      std::cerr << "Logspace: Negative number" << std::endl;
      std::cerr << "\tin logspace_diff( " << lhs << ", " << rhs << " )" << std::endl;
    }
#endif // ALGEBRA_CHECK_NEGATIVE
    return -std::numeric_limits<double>::infinity();
  }
  // assert( lhs > rhs )
  return ( lhs + ::log1p(-exp(rhs-lhs)) );
} // logspace_diff( double, double const & )

// Added by Paul T Edlefsen 15/9/08
inline static double
toDouble ( Logspace const & v )
{
  return ::exp( v.x );
} // toDouble( Logspace const & )

inline static long double
toLongDouble ( Logspace const & v )
{
  return ::exp( static_cast<long double>( v.x ) );
} // toLongDouble( Logspace const & )

inline static double
toLogDouble ( Logspace const & v )
{
  return v.x;
} // toLogDouble( Logspace const & )

inline
BFloat::
  BFloat ( Logspace const & copy_from ) {
    *this = bfloat_doubleexp( copy_from.x );
  };
inline
BFloat::
  BFloat ( LongDoubleRealspace const & copy_from ) {
    *this = copy_from.x;
  };
inline
BFloat::
  BFloat ( DoubleRealspace const & copy_from ) {
    *this = copy_from.x;
  };
inline
BFloat::
  BFloat ( FloatRealspace const & copy_from ) {
    *this = copy_from.x;
  };

//TAH 11/01/11 Modified next four functions to compile with g++ v>4.5
inline
BFloat &
  BFloat::operator= ( Logspace const & copy_from ) {
    *this = bfloat_doubleexp( copy_from.x );
    return *this;
  };

inline
BFloat &
  BFloat::operator= ( DoubleRealspace const & copy_from ) {
    *this = copy_from.x;
    return *this;
  };
inline
BFloat &
  BFloat::operator= ( LongDoubleRealspace const & copy_from ) {
    *this = static_cast<double>( copy_from.x );
    return *this;
  };
inline
BFloat &
  BFloat::operator= ( FloatRealspace const & copy_from ) {
    *this = copy_from.x;
    return *this;
  };


//
// Algebra class - Wrapper for overloading all arithmetic operators to use a different algebra.
//
// Gerton Lunter, 19/3/03
// Based on logprob.h by by Ian Holmes.
//

template <class AlgebraMethods>
class Algebra {
// Boost serialization
#ifndef NBOOST_SERIALIZATION
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize ( Archive & ar, const unsigned int /* file_version */ )
  {
    ar & BOOST_SERIALIZATION_NVP( val );
  } // serialize( Archive &, const unsigned int )
#endif // NBOOST_SERIALIZATION

public:
  // typedef
  typedef typename AlgebraMethods::Value Value;

  // value
  Value val;

public:
  // constructors
  Algebra() { }  // no initialisation, for speed
  Algebra (double px) : val(from_double(px)) { }
  //Algebra (const Algebra<AlgebraMethods> & lx) : val(lx.val) { }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  Algebra (const Algebra<AnyAlgebraMethods> & lx) : val(lx.val) { }
  Algebra (const BFloat & v) : val(v) { }
  // PTE 14/9/08
  Algebra (const Logspace & v) : val(v) { }
  // PTE 15/9/08
  Algebra (const LongDoubleRealspace & v) : val(v) { }
  // PTE 15/9/08
  Algebra (const FloatRealspace & v) : val(v) { }

  // PTE 13/9/08
  template <class other_type>
  explicit Algebra (const other_type & px) : val(from_double(px)) { }
  // PTE 13/9/08
  //Algebra (int px) : val(from_double((double)px)) { }

  // fast initialization
  void clear() { val.clear(); }

  // assignment operators
  //inline Algebra& operator= (const Algebra& lx) { val = lx.val; return *this; }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline Algebra<AlgebraMethods> & operator= (const Algebra<AnyAlgebraMethods> & lx) { val = lx.val; return *this; }
  inline Algebra& operator= (double px) { val = from_double(px); return *this; }
  // PTE 13/9/08
  //template <class other_type>
  //inline Algebra& operator= (const other_type & px) { val = from_double(px); return *this; }
  // PTE 13/9/08
  inline Algebra& operator= (int px) { val = from_double((double)px); return *this; }
  // PTE 13/9/08
  inline Algebra& operator= (uint32_t px) { val = from_double((double)px); return *this; }

  // arithmetic operators; all combinations of Algebra and double are covered
  //inline friend Algebra operator+ (const Algebra& lx, const Algebra& ly) { return from_value (psum (lx.val, ly.val)); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend Algebra<AlgebraMethods> operator+ (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return from_value (psum (lx.val, ly.val)); }
  inline friend Algebra operator+ (const Algebra& lx, double py) { return from_value (psum (lx.val, from_double(py))); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator+ (const Algebra& lx, const other_type & py) { return from_value (psum (lx.val, from_double(py))); }
  inline friend Algebra operator+ (double px, const Algebra& ly) { return from_value (psum (from_double(px), ly.val)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator+ (const other_type & px, const Algebra& ly) { return from_value (psum (from_double(px), ly.val)); }
  //inline Algebra& operator+= (const Algebra& lx) { psumacc (val, lx.val); return *this; }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline Algebra<AlgebraMethods> & operator+= (const Algebra<AnyAlgebraMethods> & lx) { psumacc (val, lx.val); return *this; }
  inline Algebra& operator+= (double px) { psumacc (val, from_double(px)); return *this; }
  // PTE 13/9/08
  //template <class other_type>
  //inline Algebra& operator+= (const other_type & px) { psumacc (val, from_double(px)); return *this; }

  //inline friend Algebra operator- (const Algebra& lx, const Algebra& ly) { return from_value (pdiff (lx.val, ly.val)); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend Algebra<AlgebraMethods> operator- (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return from_value (pdiff (lx.val, ly.val)); }

  inline friend Algebra operator- (const Algebra& lx, double py) { return from_value (pdiff (lx.val, from_double(py))); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator- (const Algebra& lx, const other_type & py) { return from_value (pdiff (lx.val, from_double(py))); }
  inline friend Algebra operator- (double px, const Algebra& ly) { return from_value (pdiff (from_double(px), ly.val)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator- (const other_type & px, const Algebra& ly) { return from_value (pdiff (from_double(px), ly.val)); }
  //inline Algebra& operator-= (const Algebra& lx) { pdiffacc (val, lx.val); return *this; }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline Algebra<AlgebraMethods> & operator-= (const Algebra<AnyAlgebraMethods> & lx) { pdiffacc (val, lx.val); return *this; }
  inline Algebra& operator-= (double px) { pdiffacc (val, from_double(px)); return *this; }
  // PTE 13/9/08
  //template <class other_type>
  //inline Algebra& operator-= (const other_type & px) { pdiffacc (val, from_double(px)); return *this; }

  //inline friend Algebra operator* (const Algebra& lx, const Algebra& ly) { return from_value (pmul (lx.val, ly.val)); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend Algebra<AlgebraMethods> operator* (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return from_value (pmul (lx.val, ly.val)); }
  inline friend Algebra operator* (const Algebra& lx, double py) { return from_value (pmuldouble (lx.val, py)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator* (const Algebra& lx, const other_type & py) { return from_value (pmuldouble (lx.val, py)); }
  inline friend Algebra operator* (double px, const Algebra& ly) { return from_value (pmuldouble (ly.val, px)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator* (const other_type & px, const Algebra& ly) { return from_value (pmuldouble (ly.val, px)); }
  //inline Algebra& operator*= (const Algebra& lx) { pmulacc (val, lx.val); return *this; }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline Algebra<AlgebraMethods> & operator*= (const Algebra<AnyAlgebraMethods> & lx) { pmulacc(val, lx.val); return *this; }
  inline Algebra& operator*= (double px) { pmulaccdouble (val, px); return *this; }
  // PTE 13/9/08
  //template <class other_type>
  //inline Algebra& operator*= (const other_type & px) { pmulaccdouble (val, px); return *this; }

  //inline friend Algebra operator/ (const Algebra& lx, const Algebra& ly) { return from_value (pdiv (lx.val, ly.val)); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend Algebra<AlgebraMethods> operator/ (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return from_value (pdiv (lx.val, ly.val)); }
  inline friend Algebra operator/ (const Algebra& lx, double py) { return from_value (pdiv (lx.val, from_double(py))); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator/ (const Algebra& lx, const other_type & py) { return from_value (pdiv (lx.val, from_double(py))); }
  inline friend Algebra operator/ (double px, const Algebra& ly) { return from_value (pdiv (from_double(px), ly.val)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend Algebra operator/ (const other_type & px, const Algebra& ly) { return from_value (pdiv (from_double(px), ly.val)); }
  //inline Algebra& operator/= (const Algebra& lx) { pdivacc (val, lx.val); return *this; }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline Algebra<AlgebraMethods> & operator/= (const Algebra<AnyAlgebraMethods> & lx) { pdivacc (val, lx.val); return *this; }
  inline Algebra& operator/= (double px) { pdivacc (val, from_double(px)); return *this; }
  // PTE 13/9/08
  //template <class other_type>
  //inline Algebra& operator/= (const other_type & px) { pdivacc (val, from_double(px)); return *this; }

  // miscellaneous operators
  inline friend double log ( const Algebra& lx ) { return doublelog( lx.val ); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend const other_type & log( const Algebra& lx ) { return doublelog( lx.val ); }
  //static inline Algebra exp ( const Algebra & lx ) { return doubleexp( to_double( lx ) ); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  static inline Algebra<AlgebraMethods> exp ( const Algebra<AnyAlgebraMethods> & lx ) { return doubleexp( to_double( lx ) ); }
  // PTE 14/9/08
  template <class other_type>
  static inline Algebra exp ( const other_type & lx ) { return doubleexp( ( double )lx ); }

  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  //static inline Algebra<AlgebraMethods> pow ( const Algebra<AlgebraMethods> & lhs, const Algebra<AnyAlgebraMethods> & rhs ) { return lhs.ppow( rhs ); }
  inline friend Algebra<AlgebraMethods> pow ( const Algebra<AlgebraMethods> & lhs, const Algebra<AnyAlgebraMethods> & rhs ) { return from_value( ppow( lhs.val, rhs.val ) ); }
  // PTE 14/9/08
  template <typename other_type>
  //static inline Algebra pow ( const Algebra& lhs, const other_type & rhs ) { return lhs.ppowdouble( (double)rhs ); }
  inline friend Algebra pow ( const Algebra& lhs, const other_type & rhs ) { return from_value( ppowdouble( lhs.val, (double)rhs ) ); }

//   // PTE 13/9/08
//   inline friend bool _isinf () {
//     // TODO: Fix.
//     return isinf( prob() );
//   } // isinf( Algebra const & )
// 
//   // PTE 13/9/08
//   inline friend bool _isnan () {
//     // TODO: Fix.
//     return isnan( prob() );
//   } // isnan( Algebra const & )

  // increment & decremement
  Algebra& operator++ () { *this += 1.; return *this; }
  Algebra operator++ (int) { Algebra tmp (*this); ++(*this); return tmp; }

  Algebra& operator--() { *this -= 1.; return *this; }
  Algebra operator-- (int) { Algebra tmp (*this); --(*this); return tmp; }

  // relational operators
  //inline friend int operator== (const Algebra& lx, const Algebra& ly) { return equal(lx.val, ly.val); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend int operator== (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return equal(lx.val, ly.val); }
  inline friend int operator== (const Algebra& lx, const double py) { return equal(lx.val, from_double(py)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator== (const Algebra& lx, const other_type & py) { return equal(lx.val, from_double(py)); }
  inline friend int operator== (const double px, const Algebra& ly) { return equal(from_double(px), ly.val); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator== (const other_type & px, const Algebra& ly) { return equal(from_double(px), ly.val); }

  //inline friend int operator!= (const Algebra& lx, const Algebra& ly) { return !equal(lx.val, ly.val); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend int operator!= (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return !equal(lx.val, ly.val); }
  inline friend int operator!= (const Algebra& lx, const double py) { return !equal(lx.val, from_double(py)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator!= (const Algebra& lx, const other_type & py) { return !equal(lx.val, from_double(py)); }
  inline friend int operator!= (const double px, const Algebra& ly) { return !equal(from_double(px), ly.val); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator!= (const other_type & px, const Algebra& ly) { return !equal(from_double(px), ly.val); }

  //inline friend int operator< (const Algebra& lx, const Algebra& ly) { return less(lx.val, ly.val); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend int operator< (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return less(lx.val, ly.val); }
  inline friend int operator< (const Algebra& lx, const double py) { return less(lx.val, from_double(py)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator< (const Algebra& lx, const other_type & py) { return less(lx.val, from_double(py)); }
  inline friend int operator< (const double px, const Algebra& ly) { return less(from_double(px), ly.val); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator< (const other_type & px, const Algebra& ly) { return less(from_double(px), ly.val); }

  //inline friend int operator> (const Algebra& lx, const Algebra& ly) { return less(ly.val, lx.val); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend int operator> (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return less(ly.val, lx.val); }
  inline friend int operator> (const Algebra& lx, const double py) { return less(from_double(py), lx.val); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator> (const Algebra& lx, const other_type & py) { return less(from_double(py), lx.val); }
  inline friend int operator> (const double px, const Algebra& ly) { return less(ly.val, from_double(px)); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator> (const other_type & px, const Algebra& ly) { return less(ly.val, from_double(px)); }

  //inline friend int operator<= (const Algebra& lx, const Algebra& ly) { return lessequal(lx.val, ly.val); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend int operator<= (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return lessequal(lx.val, ly.val); }
  inline friend int operator<= (const Algebra& lx, const double py) { return lessequal( lx.val, from_double(py) ); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator<= (const Algebra& lx, const other_type & py) { return lessequal( lx.val, from_double(py) ); }
  inline friend int operator<= (const double px, const Algebra& ly) { return lessequal( from_double(px), ly.val); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator<= (const other_type & px, const Algebra& ly) { return lessequal( from_double(px), ly.val); }

  //inline friend int operator>= (const Algebra& lx, const Algebra& ly) { return lessequal( ly.val, lx.val); }
  // PTE 15/9/08
  template <typename AnyAlgebraMethods>
  inline friend int operator>= (const Algebra<AlgebraMethods> & lx, const Algebra<AnyAlgebraMethods> & ly) { return lessequal( ly.val, lx.val); }
  inline friend int operator>= (const Algebra& lx, const double py) { return lessequal( from_double(py), lx.val ); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator>= (const Algebra& lx, const other_type & py) { return lessequal( from_double(py), lx.val ); }
  inline friend int operator>= (const double px, const Algebra& ly) { return lessequal( ly.val, from_double(px) ); }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend int operator>= (const other_type & px, const Algebra& ly) { return lessequal( ly.val, from_double(px) ); }

  // stream operators
  inline friend std::ostream& operator<< (std::ostream& out, const Algebra& lx) { return AlgebraMethods::print(out, lx.val); }
  inline friend std::istream& operator>> (std::istream& in, Algebra& lx) { double px; in >> px; lx.val = px; return in; }
  // PTE 13/9/08
  //template <class other_type>
  //inline friend istream& operator>> (istream& in, const Algebra& lx) { const other_type & px; in >> px; lx.val = px; return in; }

  // cast operators
  inline double prob () const { return to_double( val ); }
  // PTE 13/9/08
  //inline operator double() const { return to_double (val); }

private:
  // private AlgebraMethods method wrappers
  static inline double to_double (Value const & X) { return AlgebraMethods::to_double (X); }
  static inline Value from_double (double P) { return AlgebraMethods::from_double (P); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline Value pmul (Value const & X, AnyValue const & Y) { return AlgebraMethods::pmul (X, Value(Y)); }
  static inline Value pmul (Value const & X, Value const & Y) { return AlgebraMethods::pmul (X, Y); }
  static inline Value pmuldouble (Value const & X, double Y) { return AlgebraMethods::pmuldouble (X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline Value ppow (Value const & X, AnyValue const & Y) { return AlgebraMethods::ppow (X, Value(Y)); }
  // PTE 14/9/08
  static inline Value ppow (Value const & X, Value const & Y) { return AlgebraMethods::ppow (X, Y); }
  // PTE 14/9/08
  static inline Value ppowdouble (Value X, double Y) { return AlgebraMethods::ppowdouble (X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline Value pdiv (Value const & X, AnyValue const & Y) { return AlgebraMethods::pdiv( X, Value(Y)); }
  static inline Value pdiv (Value const & X, Value const & Y) { return AlgebraMethods::pdiv( X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline Value psum (Value const & X, AnyValue const & Y) { return AlgebraMethods::psum (X, Value(Y)); }
  static inline Value psum (Value const & X, Value const & Y) { return AlgebraMethods::psum (X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline Value pdiff (Value const & X, AnyValue const & Y) { return AlgebraMethods::pdiff (X, Value(Y)); }
  static inline Value pdiff (Value const & X, Value const & Y) { return AlgebraMethods::pdiff (X, Y); }
  static inline Value doubleexp (double X) { return AlgebraMethods::doubleexp( X ); }
  static inline double doublelog (Value const & X) { return AlgebraMethods::doublelog( X ); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline void pmulacc (Value& X, AnyValue const & Y) { AlgebraMethods::pmulacc (X, Value(Y)); }
  static inline void pmulacc (Value& X, Value const & Y) { AlgebraMethods::pmulacc (X, Y); }
  static inline void pmulaccdouble (Value& X, double Y) { AlgebraMethods::pmulaccdouble (X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline void pdivacc ( Value& X, AnyValue const & Y ) { AlgebraMethods::pdivacc( X, Value(Y) ); }
  static inline void pdivacc ( Value& X, Value const & Y ) { AlgebraMethods::pdivacc( X, Y ); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline void psumacc (Value& X, AnyValue const & Y) { AlgebraMethods::psumacc(X, Value(Y)); }
  static inline void psumacc (Value& X, Value const & Y) { AlgebraMethods::psumacc(X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline void pdiffacc (Value& X, AnyValue const & Y) { AlgebraMethods::pdiffacc (X, Value(Y)); }
  static inline void pdiffacc (Value& X, Value const & Y) { AlgebraMethods::pdiffacc (X, Y); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline bool less (Value const & X, AnyValue const & Y ) { return AlgebraMethods::less( X, Value(Y) ); }
  static inline bool less (Value const & X, Value const & Y ) { return AlgebraMethods::less( X, Y ); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline bool equal (Value const & X, AnyValue const & Y ) { return AlgebraMethods::equal( X, Value(Y)); }
  static inline bool equal (Value const & X, Value const & Y ) { return AlgebraMethods::equal( X, Y ); }
  // PTE 15/9/08
  template <typename AnyValue>
  static inline bool lessequal( Value const & X, AnyValue const & Y ) { return AlgebraMethods::lessequal( X, Value(Y) ); }
  static inline bool lessequal( Value const & X, Value const & Y ) { return AlgebraMethods::lessequal( X, Y ); }

public:
  // static constructor from Value type value
  // PTE 15/9/08
  template <typename AnyValue>
  static inline Algebra from_value (AnyValue const & X) { Algebra lx; lx.val = Value(X); return lx; }
  static inline Algebra from_value (Value const & X) { Algebra lx; lx.val = X; return lx; }
}; // End class Algebra


  template <typename AlgebraMethods>
  inline static int
  isinf ( Algebra<AlgebraMethods> const & a )
  {
    return a.val.isInfinity(); //std::isinf( toLogDouble( a ) );
  } // isinf( Algebra<AlgebraMethods> const & )

  template <typename AnyOtherType>
  inline static bool
  isinf ( AnyOtherType const & v )
  {
    return std::isinf( v );
  } // isinf( AnyOtherType const & )

  template <typename AlgebraMethods>
  inline static int
  isnan ( Algebra<AlgebraMethods> const & a )
  {
    return a.val.isNaN(); //std::isnan( toLogDouble( a ) );
  } // isnan( Algebra<AlgebraMethods> const & )

  template <typename AnyOtherType>
  inline static int
  isnan ( AnyOtherType const & v )
  {
    return std::isnan( v );
  } // isnan( AnyOtherType const & )


//// TODO: Why won't this work? (ANSWER: it seems that isnan and isinf are macros!)
////template <class AlgebraMethods>
////inline static bool
////::isinf ( Algebra<AlgebraMethods> const & a )
////{
////  // TODO: Fix.
////  return ::isinf( a.prob() );
////} // ::isinf( Algebra const & )
//
//inline static bool
//::isinf ( Algebra<BFloatMethods> const & a )
//{
//  // TODO: Fix.
//  return ::isinf( a.prob() );
//} // ::isinf( Algebra const & )
//
//inline static bool
//::isinf ( Algebra<LogspaceMethods> const & a )
//{
//  // TODO: Fix.
//  return ::isinf( a.prob() );
//} // ::isinf( Algebra const & )
//
//template <class AlgebraMethods>
//inline static bool
//::isnan ( Algebra<AlgebraMethods> const & a )
//{
//  // TODO: Test/fix.
//  return ::isnan( a.prob() );
//} // ::isnan( algebra const & )

} // End namespace hmmoc

// TODO: Finish these

namespace std
{

// Added by Paul T Edlefsen 13/9/08
template <class AlgebraMethods>
inline static double
toDouble ( hmmoc::Algebra<AlgebraMethods> const & a )
{
  return a.prob();
} // toDouble( Algebra const & )

template <class AnyOtherType>
inline static double
toDouble ( AnyOtherType const & v )
{
  return static_cast<double>( v );
} // toDouble( AnyOtherType const & )

template <class AlgebraMethods>
inline static long double
toLongDouble ( hmmoc::Algebra<AlgebraMethods> const & a )
{
  return a.prob();
} // toLongDouble( Algebra const & )

template <class AnyOtherType>
inline static long double
toLongDouble ( AnyOtherType const & v )
{
  return static_cast<long double>( v );
} // toDouble( AnyOtherType const & )

template <class AlgebraMethods>
inline static double
toLogDouble ( hmmoc::Algebra<AlgebraMethods> const & a )
{
  return log( a );
} // toLogDouble( Algebra const & )

/// Paul added for (eg) prob - .1
template <typename AlgebraMethods,
          typename T>
ALGEBRA_INLINE_CAST_ARITHMETIC
typename hmmoc::Algebra<AlgebraMethods>
operator - (
  const hmmoc::Algebra<AlgebraMethods>& lhs,
  const T& rhs
)
{
  hmmoc::Algebra<AlgebraMethods> p(lhs);
  p -= rhs;
  return p;
}

/// Paul added for (eg) 1.0 - prob
template <typename T,
          typename AlgebraMethods>
ALGEBRA_INLINE_CAST_ARITHMETIC
hmmoc::Algebra<AlgebraMethods>
operator - (
  const T& lhs,
  const hmmoc::Algebra<AlgebraMethods> & rhs)
{
  hmmoc::Algebra<AlgebraMethods> p(lhs);
  p -= rhs;
  return p;
}

// NOTE: This forces the result of all multiplications into the log domain!
// TODO: Change!
/// Multiplication of probabilities.
template <typename LAlgebraMethods,
          typename RAlgebraMethods>
ALGEBRA_INLINE_CAST_ARITHMETIC
hmmoc::Algebra<LAlgebraMethods>
operator * (
  const hmmoc::Algebra<LAlgebraMethods> & lhs,
  const hmmoc::Algebra<RAlgebraMethods> & rhs
)
{
  hmmoc::Algebra<LAlgebraMethods> p(lhs);
  p *= rhs;
  // TODO: REMOVE
  //cout << "ProbOp*2*: linear?: " << lhs << ", log?: " << rhs << ", rval: " << p << std::endl;
  return p;
}

// NOTE: This forces the result of all divisions into the log domain!
// TODO: Change!
/// Paul added for (eg) .5 * prob
template <typename T,
          typename AlgebraMethods>
ALGEBRA_INLINE_CAST_ARITHMETIC
typename hmmoc::Algebra<AlgebraMethods>
operator * (
  const T& lhs,
  const hmmoc::Algebra<AlgebraMethods> & rhs
)
{
  hmmoc::Algebra<AlgebraMethods> p(lhs);
  p *= rhs;
  // TODO: REMOVE
  //cout << "ProbOp*3*: linear?: " << lhs << ", log?: " << rhs << ", rval: " << p << std::endl;
  return p;
}

// NOTE: This forces the result of all divisions into the log domain!
// TODO: Change!
/// Paul added for (eg) prob * .5
template <typename AlgebraMethods,
          typename T>
ALGEBRA_INLINE_CAST_ARITHMETIC
hmmoc::Algebra<AlgebraMethods>
operator * (
  const hmmoc::Algebra<AlgebraMethods>& lhs,
  const T& rhs
)
{
  hmmoc::Algebra<AlgebraMethods> p(lhs);
  p *= rhs;
  // TODO: REMOVE
  //cout << "ProbOp*4*: linear?: " << lhs << ", log?: " << rhs << ", rval: " << p << std::endl;
  return p;
}

/// Paul added for (eg) p_dbl = .5; p_dbl *= prob;
// Casts rhs to double and multiplies it into lhs.
template <typename AlgebraMethods>
ALGEBRA_INLINE_CAST_ARITHMETIC
double &
operator *= (
  double& lhs,
  const hmmoc::Algebra<AlgebraMethods> & rhs
)
{
  lhs *= toDouble( rhs );
  return lhs;
}

// NOTE: This forces the result of all divisions into the log domain!
// TODO: Change!
/// Division of probabilities.
template <typename LAlgebraMethods,
          typename RAlgebraMethods>
ALGEBRA_INLINE_CAST_ARITHMETIC
hmmoc::Algebra<LAlgebraMethods>
operator / (
  const hmmoc::Algebra<LAlgebraMethods>& lhs,
  const hmmoc::Algebra<RAlgebraMethods> & rhs
)
{
  hmmoc::Algebra<LAlgebraMethods> p(lhs);
  p /= rhs;
  return p;
}

// NOTE: This forces the result of all divisions into the log domain!
// TODO: Change!
/// Paul added for (eg) 1.0 / prob
template <typename T,
          typename AlgebraMethods>
ALGEBRA_INLINE_CAST_ARITHMETIC
hmmoc::Algebra<AlgebraMethods>
operator / (
  const T& lhs,
  const hmmoc::Algebra<AlgebraMethods> & rhs
)
{
  hmmoc::Algebra<AlgebraMethods> p(lhs);
  p /= rhs;
  return p;
}

// NOTE: This forces the result of all divisions into the log domain!
// TODO: Change!
/// Paul added for (eg) prob / 1.0
template <typename AlgebraMethods,
          typename T>
ALGEBRA_INLINE_CAST_ARITHMETIC
hmmoc::Algebra<AlgebraMethods>
operator / (
  const hmmoc::Algebra<AlgebraMethods>& lhs,
  const T& rhs
)
{
  hmmoc::Algebra<AlgebraMethods> p(lhs);
  p /= rhs;
  return p;
}
//@}

//   /**
//    * @brief Specialization of @c std::numeric_limits for probabilities.
//    *
//    * In almost all instances, the functions and values encapsulated
//    * for @c std::numeric_limits<hmmoc::probabilities::probability> are
//    * simply delegated to the corresponding functions and values for
//    * the underlying @c value_type.  The following, however, are
//    * notable exceptions.
//    *
//    * @li @c is_signed is false regardless of the @c value_type as
//    * probabilities cannot be negative.
//    *
//    * @li @c min() and @c max() return the extreme values for the
//    * internal value of a probability in the relevant domain.  The
//    * range corresponds to @f$ [0,1] @f$ and @f$ [-\infty, 0] @f$ for
//    * the hmmoc::probabilities::linear_domain and
//    * hmmoc::probabilities::log_domain, respectively.
//    *
//    * @li @c min() returns the value of @c
//    * std::numeric_limits<Value>::min() for any domain.  For types
//    * without denormalization, this corresponds to the minimum finite
//    * value representable by the internal type; for types with
//    * denormalization this corresponds to the minimum positive
//    * normalized value representable by the internal type.  This seems
//    * to be the most appropriate interpretation of the intent of @c
//    * std:numeric_limits<>::min().
//    *
//    * @li @c max() returns a value of unity (1) for eny domain.  This
//    * corresponds to the maximum value in the linear domain, which
//    * seems to be the most appropriate interpretation of the intent of
//    * @c std:numeric_limits<>::max().
//    */
//   template <typename Domain, typename Value, typename Validator>
//   struct std::numeric_limits< hmmoc::probabilities::
//     probability<Domain,Value,Validator> >
//   {
//     static const bool is_specialized = true;

//     static Value min() throw() { return std::numeric_limits<Value>::min(); }
//     static Value max() throw() { return Value(1); }

//     static const int digits = std::numeric_limits<Value>::digits;
//     static const int digits10 = std::numeric_limits<Value>::digits10;
//     static const bool is_signed = false;
//     static const bool is_integer = std::numeric_limits<Value>::is_integer;
//     static const bool is_exact = std::numeric_limits<Value>::is_exact;
//     static const int radix = std::numeric_limits<Value>::radix;
//     static Value epsilon() throw() { return std::numeric_limits<Value>::epsilon(); }
//     static Value round_error() throw()
//     { return std::numeric_limits<Value>::round_error(); }

//     static const int min_exponent = std::numeric_limits<Value>::min_exponent;
//     static const int min_exponent10 = std::numeric_limits<Value>::min_exponent10;
//     static const int max_exponent = std::numeric_limits<Value>::max_exponent;
//     static const int max_exponent10 = std::numeric_limits<Value>::max_exponent10;

//     static const bool has_infinity = std::numeric_limits<Value>::has_infinity;
//     static const bool has_quiet_NaN = std::numeric_limits<Value>::has_quiet_NaN;
//     static const bool has_signaling_NaN
//     = std::numeric_limits<Value>::has_signaling_NaN;
//     static const float_denorm_style has_denorm
//     = std::numeric_limits<Value>::has_denorm;
//     static const bool has_denorm_loss = std::numeric_limits<Value>::has_denorm_loss;

//     static Value infinity() throw()
//     { return std::numeric_limits<Value>::infinity(); }
//     static Value quiet_NaN() throw()
//     { return std::numeric_limits<Value>::quiet_NaN(); }
//     static Value signaling_NaN() throw()
//     { return std::numeric_limits<Value>::signaling_NaN(); }
//     static Value denorm_min() throw()
//     { return std::numeric_limits<Value>::denorm_min(); }

//     static const bool is_iec559 = std::numeric_limits<Value>::is_iec559;
//     static const bool is_bounded = std::numeric_limits<Value>::is_bounded;
//     static const bool is_modulo = std::numeric_limits<Value>::is_modulo;

//     static const bool traps = std::numeric_limits<Value>::traps;
//     static const bool tinyness_before = std::numeric_limits<Value>::tinyness_before;
//     static const float_round_style round_style
//     = std::numeric_limits<Value>::round_style;
//   };

//   // Paul added for min( prob a, prob b )
//   template <typename Domain, typename Value, typename Validator>
//   hmmoc::probabilities::probability<Domain,Value,Validator>
//   min (const hmmoc::probabilities::probability<Domain,Value,Validator>& a, const hmmoc::probabilities::probability<Domain,Value,Validator>& b)
//   {
//     if( a > b ) {
//       return b;
//     }
//     return a;
//   }
  
//   // Paul added for (eg) min( prob, .5 )
//   template <typename Domain, typename Value, typename Validator,
//             typename T>
//   Value
//   min (const hmmoc::probabilities::probability<Domain,Value,Validator>& a, const T& b)
//   {
//     return std::min( a.value_cast(hmmoc::probabilities::linear_domain()), b );
//   }

//   // Paul added for (eg) min( .5, prob )
//   template <typename T,
//             typename Domain, typename Value, typename Validator>
//   Value
//   min (const T&a, const hmmoc::probabilities::probability<Domain,Value,Validator>& b)
//   {
//     return min( a, b.value_cast(hmmoc::probabilities::linear_domain()) );
//   }
  
} // End namespace std

//
// and these are the the things that we'll use:
//

#define bfloat hmmoc::Algebra<hmmoc::BFloatMethods>

#define logspace hmmoc::Algebra<hmmoc::LogspaceMethods<hmmoc::Logspace> >

#define doublerealspace hmmoc::Algebra<hmmoc::RealspaceMethods<hmmoc::DoubleRealspace> >
#define longdoublerealspace hmmoc::Algebra<hmmoc::RealspaceMethods<hmmoc::LongDoubleRealspace> >
#define floatrealspace hmmoc::Algebra<hmmoc::RealspaceMethods<hmmoc::FloatRealspace> >
#define realspace doublerealspace

#endif // HMMOC_BFLOAT_ALGEBRA_HPP
