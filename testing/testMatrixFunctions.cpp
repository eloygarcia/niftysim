// =========================================================================
// File:       testVectorArithmetic.cpp
// Purpose:    CPU vector arithmetic unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2011
//
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
//
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledVectorArithmetic.h"
#include "tledMatrixFunctions.h"

#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <limits>

using namespace tledUnitTest;

/** Decomposes a matrix that is already SDP -> SVD same as eigenvalue decomp. */
static void _TestSVDOnSPD() {
  using namespace tledVectorArithmetic;

  int numFound, numTotal;
  float rotError;

  numFound = numTotal = 0;
  rotError = 0;
  for (int t = 0; t < 100000; t++) {
    float EV[3*3], E[3];
    float A[3*3], tmp[3];
    float U[3*3], V[3*3], S[3];

    for (float *p_v = EV; p_v < EV + 3; p_v++) *p_v = (drand48() - drand48());
    for (float *p_v = EV + 3; p_v < EV + 6; p_v++) *p_v = (drand48() - drand48());
    ScalarDiv(EV, Norm(EV));
    Sub(EV + 3, EV + 3, ScalarMul(tmp, EV, Dot(EV, EV + 3)));
    ScalarDiv(EV + 3, Norm(EV + 3));
    ScalarDiv(EV + 6, Norm(Cross(EV + 6, EV, EV + 3)));

    for (float *p_e = E; p_e < E + 3; p_e++) *p_e = 10*drand48();
    std::sort(E, E + 3, std::greater<float>());

    for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
	A[3*r+c] = 0;
	for (float const *pc_vl = EV + 3*r, *pc_e = E, *pc_vr = EV + 3*c; pc_vl < EV + 3*r + 3; pc_vl++, pc_e++, pc_vr++) {
	    A[r*3+c] += (*pc_vr)*(*pc_e)*(*pc_vl);
	  }
      }

    SVD3x3(U, S, V, A);

    for (int e = 0; e < 3; e++) {
      /* Can afford to miss some close singular values as main application is determining rigid body rotations, where actual s. vals not of any interest */
      if (std::fabs(S[e] - E[e]) < 1e-3*(E[0] + E[1] + E[2])) numFound += 1;
      numTotal += 1;
    }
    for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
	rotError += std::fabs(U[3*r+c] - V[3*r+c]);
	tledUnitTestAssert(!std::isnan(U[3*r+c]) && !std::isnan(V[3*r+c]));
      }
  }

  tledUnitTestAssert(numTotal - numTotal < 1e-3*numTotal);
  tledUnitTestAssert(rotError < 1e-4*numTotal);
  std::cout << "Found " << numFound << "/" << numTotal << std::endl;
}

/** Decomposes a matrix that is already SDP -> SVD same as eigenvalue decomp. */
static void _TestSpecial() {
  static const float A[][9] = {{5.55068159f, 0.000297248363f, -0.00282490253f, 0.000297278166f, 5.55150747f, 0.00143945217f, -0.00282490253f, 0.00143939257f, 5.54389668f},
			       {3.85472012f, -1.09679937f, -3.08530664f, -1.09679937f, 0.751185894f, 1.4400984f, -3.08530664f, 1.44009829f, 5.55591631f},
			       {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
  static const float S[][3] = {{5.5519f, 5.5515f, 5.5426f},
			       {8.3377f, 1.5049f, 0.3193f}};

  float const (*pc_sol)[3] = S;

  for (float const (*pc_input)[9] = A; !std::isnan((*pc_input)[0]); pc_input++, pc_sol++) {
    float U[3*3], V[3*3], S[3];

    SVD3x3(U, S, V, *pc_input);
    for (int i = 0; i < 3; i++) {
      tledUnitTestAssert(std::fabs(S[i] - (*pc_sol)[i]) < 1e-3f*std::fabs((*pc_sol)[i]) || std::fabs((*pc_sol)[i]) < 1e-5f);
    }
  }
}

template <typename TScalar>
static void _ReadQRDTest(std::vector<TScalar> &r_dst, const std::string &pfx, const int numRows, const int numCols) {
  std::ostringstream oss;
  std::ifstream fin;

  oss << "QRD_test_" << pfx << numRows << "x" << numCols << ".txt";
  fin.open(tledUnitTest::GetResourcePath(oss.str()).c_str());
  std::copy(std::istream_iterator<double>(fin), std::istream_iterator<double>(), std::back_inserter(r_dst));
  assert(r_dst.size() > 0 && r_dst.size()%numCols == 0 && r_dst.size()%numRows == 0);
}

/** Mostly a simple regression test against Matlab results. */
static void _TestQRDInverse() {
  const int numRows[] = {10, 15, -1};
  const int numCols[] = {10, 10, -1};

  for (int t = 0; numRows[t] > 0; t++) {
    std::vector<double> A;

    {
      std::vector<double> R, QTest(numRows[t]*numRows[t]), RTest(numRows[t]*numCols[t]);
      std::vector<double> QtA(numRows[t]*numCols[t]);

      _ReadQRDTest(A, "", numRows[t], numCols[t]);
      QRDMxN(&QTest.front(), &RTest.front(), &A.front(), numRows[t], numCols[t]);
      MatMultAB(&QTest.front(), numRows[t], numRows[t], &A.front(), numRows[t], numCols[t], &QtA.front());
      _ReadQRDTest(R, "R_ref_", numRows[t], numCols[t]);
      for (int r = 0; r < numRows[t]; r++) for (int c = r; c < numCols[t]; c++) {
	  tledUnitTestAssert(std::fabs(std::fabs(QtA[r*numCols[t]+c]) - std::fabs(R[r*numCols[t]+c])) < 1e-2*std::fabs(R[r*numCols[t]+c]));
	  tledUnitTestAssert(std::fabs(std::fabs(RTest[r*numCols[t]+c]) - std::fabs(R[r*numCols[t]+c])) < 1e-2*std::fabs(R[r*numCols[t]+c]));
	}
    }

    if (numRows[t] == numCols[t]) {
      std::vector<double> AInv(numRows[t]*numRows[t]), ITest(numRows[t]*numRows[t]);

      std::copy(A.begin(), A.end(), AInv.begin());
      MatInverse(&AInv.front(), numRows[t]);
      MatMultAB(&AInv.front(), numRows[t], numRows[t], &A.front(), numRows[t], numRows[t], &ITest.front());
      for (int r = 0; r < numRows[t]; r++) for (int c = 0; c < numRows[t]; c++) {
	  if (r == c) {
	    tledUnitTestAssert(std::fabs(1 - ITest[r*numRows[t]+c]) < 1e-3);
	  } else {
	    tledUnitTestAssert(std::fabs(ITest[r*numRows[t]+c]) < 1e-3);
	  }
	}
    }
  }
}

int main(void) {
  InitUnitTest();

  _TestSpecial();
  _TestSVDOnSPD();
  _TestQRDInverse();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
