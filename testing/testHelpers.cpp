// =========================================================================
// File:       testHelpers.cpp
// Purpose:    Helper module unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include <cstdlib>
#include "tledUnitTest.h"
#include "tledHelper.h"

using namespace std;
using namespace tledUnitTest;

int main(void) {
  InitUnitTest();

  {
    const int numEntries = 1000;

    vector<int> testVec(numEntries), outVec(0);
    vector<int>::iterator i_e;
    vector<int>::const_iterator ic_ref;

    for (i_e = testVec.begin(); i_e < testVec.end(); i_e++) *i_e = rand()%(numEntries/10);
    outVec = tledHelper::MakeSortedUnique(testVec);
    assert(outVec.size() < testVec.size());
    
    for (i_e = outVec.begin(); i_e < outVec.end(); i_e++) {
      ic_ref = find(testVec.begin(), testVec.end(), *i_e);
      tledUnitTestAssert(ic_ref < testVec.end());
      if (i_e > outVec.begin()) {
	tledUnitTestAssert(*(i_e - 1) < *i_e);
      }
    }
  }

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
