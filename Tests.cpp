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
#include "Algebra.hpp"

#include <limits>
using std::numeric_limits;

#include <sstream> // For stringstream.
#include <ctime> // for std::time

#include <iostream>
using std::cout;
using std::endl;

using namespace hmmoc;

int
test_algebras ( int argc, char **argv )
{
  clock_t start = clock();

  bfloat bfloat0 = 1.0;
  cout << "[0] bfloat0: " << bfloat0 << endl;
  bfloat0 *= 1E-10;
  cout << "[1] bfloat0 after *= 1E-10: " << bfloat0 << endl;
  cout << "[1] log(bfloat0) is now: " << log( bfloat0 ) << endl;
  for( int i = 2; i <= 30; i++ ) {
    bfloat0 *= bfloat0;
    cout << "[" << i << "] bfloat0 after *= bfloat0: " << bfloat0 << endl;
    cout << "[" << i << "] log(bfloat0) is now: " << log( bfloat0 ) << endl;
    cout << "[" << i << "] 1.0 - bfloat0 is now: " << ( 1.0 - bfloat0 ) << endl;
  }
  
  bfloat bfloat1 = 1.0;
  cout << "[0] bfloat1: " << bfloat1 << endl;
  bfloat1 *= 1E-1;
  cout << "[1] bfloat1 after *= 1E-1: " << bfloat1 << endl;
  cout << "[1] log(bfloat1) is now: " << log( bfloat1 ) << endl;
  for( int i = 2; i <= 10; i++ ) {
    bfloat1 *= 1E-1;
    cout << "[" << i << "] bfloat1 after *= 1E-1: " << bfloat1 << endl;
    cout << "[" << i << "] log(bfloat1) is now: " << log( bfloat1 ) << endl;
    cout << "[" << i << "] 1.0 - bfloat1 is now: " << ( 1.0 - bfloat1 ) << endl;
  }

  cout << "largest float: " << numeric_limits<float>::max() << endl;
  cout << "smallest float: " << numeric_limits<float>::min() << endl;

  bfloat bfloat2 = ( bfloat1 * 1E-5 );
  cout << " bfloat1 is " << bfloat1 << endl;
  cout << " bfloat2 is " << bfloat2 << endl;
  //cout << " ( bfloat2 - bfloat1 ) = " << ( bfloat2 - bfloat1 ) << endl;
  cout << " ( bfloat1 - bfloat2 ) = " << ( bfloat1 - bfloat2 ) << endl;
  cout << " bfloat0 is " << bfloat0 << endl;
  //cout << " ( bfloat0 - bfloat1 ) = " << ( bfloat0 - bfloat1 ) << endl;
  cout << " ( bfloat1 - bfloat0 ) = " << ( bfloat1 - bfloat0 ) << endl;

  clock_t end = clock();
  double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  cout << "Elapsed CPU time: " << cpu_time_used << endl;
  return 0;
} // test_algebras(..)

int
main ( int argc, char **argv )
{
  return
    test_algebras( argc, argv );
} // main( int, char** )


