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
 *    This file is based on algebras.cc in HMMoC 1.3, a hidden Markov model
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
// Algebra.cpp - extended real types
//
// Paul T Edlefsen, June 1, 2011
// Gerton Lunter, 27/8/04
//
//

#include "Algebra.hpp"

using namespace hmmoc;

BFMantissa *BFloat::aConversionLookup;           // Actual location of the static members of BFloat class
double *BFloat::aDoubleConversionLookup;


_BFloatInitialize _dummyInitializer;             // This initializes aConversionLookup and aDoubleConversionLookup


_BFloatInitialize::_BFloatInitialize () {

  BFloat::aConversionLookup = new BFMantissa[cBFloatConvTableSize];
  BFloat::aDoubleConversionLookup = new double[cBFloatDoubleConvTableSize];

  BFMantissa iBFM = 1.0;
  for (int i = 0; i < cBFloatConvTableSize; i++) {
    BFloat::aConversionLookup[ i ] = iBFM;
    iBFM *= cBFloatRangeInv;
  }

  for (int i = 0; i < cBFloatDoubleConvTableSize; i++) {
    BFloat::aDoubleConversionLookup[ i ] = exp( (i-cBFloatDoubleConvTableSize/2) * logcBFloatRange );
  }

} // _BFloatInitialize(..)


