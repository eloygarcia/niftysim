// =========================================================================
// File:       testPolynomial.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2012
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledHelper.h"
#include "tledPolynomial.h"

#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <vector>
#include <limits>
#include <fstream>
#include <cmath>

using namespace tledPolynomial;

static int _TestSolutions(const float test_sols[], const float ref_sols[], const int numSols) {
  int numFoundSolutions;

  numFoundSolutions = 0;
  for (int s = 0, t; s < numSols; s++) {
    for (t = 0; t < numSols && std::fabs(test_sols[t] - ref_sols[s]) > 1e-3*std::fabs(ref_sols[s]); t++) {
      if (std::fabs(ref_sols[s]) < 1e-5 && std::fabs(test_sols[t]) < 1e-5) break;
    }
    if (t < numSols) numFoundSolutions += 1;
  }

  return numFoundSolutions;
}

static void _TestQuadraticSolver() {
  int numFound, numTotal;

  numFound = numTotal = 0; 
  for (int tInd = 0; tInd < 10000; tInd++) {    
    float x[2], a, coeffs[3], tx[2];

    /* Solve polynomials eq's with real solution */
    a = 10*(drand48() - drand48());
    x[0] = 10*(drand48() - drand48());
    x[1] = 10*(drand48() - drand48());

    coeffs[0] = a;
    coeffs[1] = -a*(x[0] + x[1]);
    coeffs[2] = a*x[0]*x[1];

    tledPolynomial::FindQuadraticRoots(tx, coeffs);
    numFound += _TestSolutions(x, tx, 2);
    numTotal += 2;
  }
  tledUnitTestAssert(numTotal - numFound < 1e-3f*numTotal);

  for (int tInd = 0; tInd < 1000; tInd++) {    
    float coeffs[3], tx[2];
    int gt;

    /* Solve polynomials eq's with complex solutions */
    for (gt = 0; gt < 100; gt++) {
      coeffs[1] = 10*(drand48() - drand48());
      if (std::fabs(coeffs[1]) > 1e-3) break;      
    } 
    if (gt == 100) continue;
    
    for (gt = 0; gt < 100; gt++) {
      coeffs[0] = 10*(drand48() - drand48());
      coeffs[2] = 10*(drand48() - drand48());
      if (coeffs[1]*coeffs[1] < 4*coeffs[0]*coeffs[2]) break;
    } 
    if (gt == 100) continue;

    tledPolynomial::FindQuadraticRoots(tx, coeffs);
    for (float const *pc_s = tx; pc_s < tx + 2; pc_s++) {
      tledUnitTestAssert(std::isnan(*pc_s));
    }
  }
}

/*
 * The cubic polynomial solver isn't 100% reliable, but it's good enough for our purposes, hence we only require the error rate to be less than 1 per mille
 */
static void _TestCubicSolver() {
  int numFound, numTotal;
  
  numFound = numTotal = 0;
  for (int tInd = 0; tInd < 3000; tInd++) {
    static const float interval[] = { -10.1, 10.1 };

    float ref_sol[3], test_sol[3];
    float cs[4], s, cs2[4];
    
    ref_sol[0] = 10*(drand48() - drand48());
    ref_sol[1] = 10*(drand48() - drand48());

    if (std::fabs(ref_sol[0] - ref_sol[1]) < 1e-3*std::max(std::fabs(ref_sol[0]), std::fabs(ref_sol[1]))) {
      ref_sol[0] += 1e-2*ref_sol[1];
    }

    cs[1] = 1;
    cs[2] = -(ref_sol[0] + ref_sol[1]);
    cs[3] = ref_sol[0]*ref_sol[1];

    s = (drand48()*5 + 0.001)*(rand()%2 == 0? -1 : 1);
    for (float *p_c = cs + 1; p_c < cs + 4; p_c++) *p_c *= s;

    FindQuadraticRoots(test_sol, cs + 1);
    numTotal += 2;
    numFound += _TestSolutions(test_sol, ref_sol, 2);

    /* a*x^3 + b*x^2 + c*x + d = 0 with |a| << |b|, |c|, |d| -> solution roughly same as b*x^2 + c*x + d = 0 */
    cs[0] = 1e-15*drand48()*std::min(std::fabs(cs[1]), std::min(std::fabs(cs[2]), std::fabs(cs[3])));
    FindCubicRoots(test_sol, cs, interval);
    numTotal += 2;
    numFound += _TestSolutions(test_sol, ref_sol, 2);

    /* a*x^3 + b*x^2 + c*x = 0, one solution is 0, others the same as a*x^2 + b*x + c = 0 */
    for (float *p_c = cs; p_c < cs + 3; p_c++) *p_c = *(p_c + 1);
    cs[3] = 1e-15*drand48()*s;
    ref_sol[2] = 0;
    FindCubicRoots(test_sol, cs, interval);
    numTotal += 3;
    numFound += _TestSolutions(test_sol, ref_sol, 3);

    ref_sol[2] = 5*(drand48() - drand48());

    if (std::fabs(ref_sol[0] - ref_sol[2]) < 1e-3*std::max(std::fabs(ref_sol[0]), std::fabs(ref_sol[2]))) {
      ref_sol[2] += 1e-2*ref_sol[0];
    }

    if (std::fabs(ref_sol[1] - ref_sol[2]) < 1e-3*std::max(std::fabs(ref_sol[1]), std::fabs(ref_sol[2]))) {
      ref_sol[2] += 1e-2*ref_sol[1];
    }

    cs2[0] = cs[0];
    cs2[1] = -ref_sol[2]*cs[0] + cs[1];
    cs2[2] = -ref_sol[2]*cs[1] + cs[2];
    cs2[3] = -ref_sol[2]*cs[2];
    FindCubicRoots(test_sol, cs2, interval);
    numTotal += 3;
    numFound += _TestSolutions(test_sol, ref_sol, 3);
  }
  std::cout << "total = " << numTotal << ", found = " << numFound << std::endl;
  tledUnitTestAssert(1.0f - (float)numFound/numTotal < 1e-3);
}

static void _TestCubicSolverSpecificPolys() {
  static const float specialSols[][4] = {{6.017860413, 4.134849548, -1.064121366, 1.310336828},
					 {-2.20806313, -3.63424349, -2.40320373, 1.22414041},
					 {7.53415831443273, 3.52366239792023, -1.74538379048788, 4.016951561},
					 {-5.07708883, -0.394624919, -0.416434854, 0.354697943},					 
					 {0, 1.62344804e-10, 0.00797788426, -0.00797689334},
					 {4.5232482, 0.0382397212, 0.890572548, -4.52178144},
					 {-6.81551266, -6.84867811, 0.504556358, -2.81646109},
					 {-2.31708241, 1.70305347, 0.114719778, 1.0},					 
					 {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
  static const float interval[] = {-10, 10};

  for (float const *pc_sol = &specialSols[0][0]; !std::isnan(pc_sol[0]); pc_sol += 4) {
    float coeffs[4], test_sol[3];
    
    coeffs[0] = pc_sol[3];
    coeffs[1] = (-(pc_sol[0] + pc_sol[1] + pc_sol[2]))*pc_sol[3];
    coeffs[2] = (pc_sol[0]*pc_sol[1] + pc_sol[0]*pc_sol[2] + pc_sol[1]*pc_sol[2])*pc_sol[3];
    coeffs[3] = (-pc_sol[0]*pc_sol[1]*pc_sol[2])*pc_sol[3];

    FindCubicRoots(test_sol, coeffs, interval);
    tledUnitTestAssert(_TestSolutions(test_sol, pc_sol, 3) == 3);
  }
}

static void _TestCubicSolverSpecificQuadraticPolys() {
  static const float coeffs[][4] = {{3.92767263e-10, 4.75208521, -6.95810509, -9.63133717},
				    {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
  static const float refSols[][2] = {{-0.868747234, 2.33296871},
				     {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
  static const float interval[] = {-10, 10};
  
  for (float const *pc_c = &coeffs[0][0], *pc_rs = &refSols[0][0]; !std::isnan(*pc_c); pc_c += 4, pc_rs += 3) {
    float testSol[3];

    FindCubicRoots(testSol, pc_c, interval);
    tledUnitTestAssert(std::isnan(testSol[2]));
    tledUnitTestAssert(_TestSolutions(testSol, pc_rs, 2) == 2);    
  }
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestQuadraticSolver();

  _TestCubicSolverSpecificPolys();
  _TestCubicSolverSpecificQuadraticPolys();
  _TestCubicSolver();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
