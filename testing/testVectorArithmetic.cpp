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

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <limits>

using namespace std;
using namespace tledUnitTest;
using namespace tledVectorArithmetic;

static void _TestSoln(const float m[3][3], const float r[3], const float s[3]) {
  int cInd, rInd;
  float max_s;
  
  max_s = -numeric_limits<float>::max();
  for (cInd = 0; cInd < 3; cInd++) if (max_s < fabs(s[cInd])) max_s = fabs(s[cInd]);

  for (rInd = 0; rInd < 3; rInd++) {
    float t;

    t = 0;
    for (cInd = 0; cInd < 3; cInd++) t += m[cInd][rInd]*s[cInd];
    tledUnitTestAssert(fabs(t - r[rInd]) < 1e-5*max_s);
  }
} /* _TestSoln */

/*
 * Solves numerous equations with known solution, and checks the relative error.
 * Elements are intently set to 0 to verify the correctness of the maximum column-element strategy of the 
 * linear solver.
 */
static void _Test3x3Solver(void) {
  const int numTests = 1000;

  int tInd;

  for (tInd = 0; tInd < numTests; tInd++) {
    float m[3][3], mr[3][3], m2[3][3], r[3], s[3];
    int cInd, rInd;
    
    for (cInd = 0; cInd < 3; cInd++) for (rInd = 0; rInd < 3; rInd++) {
	m[cInd][rInd] = (drand48() - drand48());
      }

    for (rInd = 0; rInd < 3; rInd++) {
      r[rInd] = (drand48() - drand48());
      m[rInd][rInd] += 1;
    }
    
    memcpy(mr, m, sizeof(m));
    memcpy(m2, m, sizeof(m));
    
    {
      int numZeros, zInd, used_rInd, used_cInd;

      numZeros = rand()%3;
      used_rInd = used_cInd = 1000;
      for (zInd = 0; zInd < numZeros; zInd++) {
	rInd = rand()%3, cInd = rand()%3;
	while (rInd == used_rInd) rInd = rand()%3;
	while (cInd == used_cInd) cInd = rand()%3;
	used_rInd = rInd, used_cInd = cInd;
	mr[cInd][rInd] = 0;
      }
    }

    SolveEquation3x3(s[0], s[1], s[2], m[0], m[1], m[2], r);
    _TestSoln(m, r, s);

    m2[0][0] = 0, m2[0][1] = 0, m2[0][2] += 1;
    SolveEquation3x3(s[0], s[1], s[2], m2[0], m2[1], m2[2], r);
    _TestSoln(m2, r, s);

    SolveEquation3x3(s[0], s[1], s[2], mr[0], mr[1], mr[2], r);
    _TestSoln(mr, r, s);
  }
} /* _Test3x3Solver */

/*
 * Tests:
 * \f$(\vec{x}, \vec{y}) = 0, ||\vec{x}|| > 0 \& ||\vec{y}|| > 0: (\vec{x} \times \vec{y}, \vec{x}) = (\vec{x} \times \vec{y}, \vec{y}) = 0\f$
 * And in doing so the following functions: Dot, Cross, Sub, Norm
 */
static void _TestCross(void) {
  const int numTests = 1000;

  int tInd, cInd;
  float x[3], y[3], z[3];

  for (tInd = 0; tInd < numTests; tInd++) {
    for (cInd = 0; cInd < 3; cInd++) {
      x[cInd] = drand48();
      y[cInd] = drand48();
    }

    ScalarDiv(x, Norm(x));
    Sub(y, y, ScalarMul(z, x, Dot(x, y)));
    ScalarDiv(y, Norm(y));
    Cross(z, x, y);

    tledUnitTestAssert(std::fabs(Dot(x, z)) < 1e-4f && std::fabs(Dot(y, z)) < 1e-4f);
  }
} /* _TestCross */    

static void _TestComputeAngle(void) {
  static const int numTests = 7;
  static const float angles[numTests] = {30.f, 30.f, 30.f, 120.f, 45.f,45.f, 90.f};
  static const float x[numTests][2][3] = {{{1.f, 0.f, 0.f}, {0.866025403784f, 0.5f, 0.f}}, 
					  {{1.f, 0.f, 0.f}, {0.866025403784f, 0.f, 0.5f}},
					  {{0.f, 1.f, 0.f}, {0.f, 0.866025403784f, 0.5f}},
					  {{1.f, 0.f, 0.f}, {-0.5f, 0.866025403784f, 0.f}},
					  {{1.f, 0.f, 0.f}, {0.707106781187f, 0.707106781187f, 0.f}},
					  {{0.f, 1.f, 0.f}, {0.f, 0.707106781187f, 0.707106781187f}},
					  {{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}}};    

  for (int tInd = 0; tInd < numTests; tInd++) {
    float xScaled[2][3];
    
    tledUnitTestAssert(std::fabs(ComputeAngleNormalised(x[tInd][0], x[tInd][1]) - angles[tInd]*tledPi/180.0) < std::fabs(angles[tInd]*tledPi/180.0)*1e-3);

    for (int xInd = 0; xInd < 2; xInd++) ScalarMul(xScaled[xInd], x[tInd][xInd], drand48()*10);
    tledUnitTestAssert(std::fabs(ComputeAngle(xScaled[0], xScaled[1]) - angles[tInd]*tledPi/180.0) < std::fabs(angles[tInd]*tledPi/180.0)*1e-3);
  }
}

static void _TestSingleComponentQuaternionConversion(void) {
  for (int qc = 0; qc < 3; qc++) {
    float q[4], rot[3][3];

    std::fill(q, q + 4, 0.f);
    q[qc] = 1.f/std::sqrt(2.f);
    q[3] = 1.f/std::sqrt(2.f);
    tledVectorArithmetic::ComputeRotationFromQuaternion(&rot[0][0], q);
    for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
	if (r == c && r == qc) {
	  tledUnitTestAssert(std::fabs(rot[r][c] - 1) < 1e-5f);
	} else {
	  if (c != qc && r != qc && r != c) {
	    tledUnitTestAssert(std::fabs(std::fabs(rot[r][c]) - 1) < 1e-5f);
	  } else {
	    tledUnitTestAssert(std::fabs(rot[r][c]) < 1e-5f);
	  }
	}
      }
    tledVectorArithmetic::ComputeQuaternionFromRotation(q, &rot[0][0]);
    for (int c = 0; c < 3; c++) {
      if (c == qc) {
	tledUnitTestAssert(std::fabs(q[c] - 1/std::sqrt(2.f)) < 1e-5f);
      } else {
	tledUnitTestAssert(std::fabs(q[c]) < 1e-5f);
      }
    }
    tledUnitTestAssert(std::fabs(std::fabs(q[3]) - 1.f/std::sqrt(2.f)) < 1e-5f);
  }
}

static void _Test90DegQuaternion(void) {
  for (int qc = 0; qc < 3; qc++) {
    const int nc = (qc + 1)%3, oc = (qc + 2)%3;

    float q[4], rot[3][3];

    std::fill(&rot[0][0], &rot[0][0] + 9, 0.f);
    rot[qc][nc] = 1.f;
    rot[nc][qc] = -1.f;
    tledVectorArithmetic::Cross(rot[oc], rot[qc], rot[nc]);

    tledVectorArithmetic::ComputeQuaternionFromRotation(q, &rot[0][0]);
    tledUnitTestAssert(std::fabs(std::fabs(q[oc]) - 1.f/std::sqrt(2.f)) < 1e-4f);
    tledUnitTestAssert(std::fabs(q[3] - 1.f/std::sqrt(2.f)) < 1e-4f);

    tledVectorArithmetic::ComputeRotationFromQuaternion(&rot[0][0], q);
    tledUnitTestAssert(std::fabs(std::fabs(rot[qc][nc]) - 1.f) < 1e-4f);
    tledUnitTestAssert(std::fabs(std::fabs(rot[nc][qc]) - 1.f) < 1e-4f);
    tledUnitTestAssert(std::fabs(std::fabs(rot[oc][oc]) - 1.f) < 1e-4f);
  }
}

static void _TestRandomAngleQuaternion(void) {
  for (int a = 0; a < 3; a++) {
    const int a0 = (a + 1)%3;
    const int a1 = (a + 2)%3;

    float rot[3][3], q[4];
    float angle = float(drand48()*tledPi/2);
    float ca = std::cos(angle), sa = std::sin(angle);

    std::fill(&rot[0][0], &rot[0][0] + 3*3, 0.f);
    rot[a][a] = 1.f;
    rot[a0][a0] = ca;
    rot[a1][a0] = -sa;
    rot[a0][a1] = sa;
    rot[a1][a1] = ca;

    tledVectorArithmetic::ComputeQuaternionFromRotation(q, &rot[0][0]);
    tledUnitTestAssert(std::fabs(q[3] - std::cos(angle/2)) < 1e-3f);
    tledUnitTestAssert(std::fabs(q[a] - std::sin(angle/2)) < 1e-3f);
    tledUnitTestAssert(std::fabs(q[a0]) + std::fabs(q[a1]) < 1e-4f);
  }
}

static void _TestSpecialQuaternion(void) {
  const float rotMats[][3*3] = {{0.78392005f, -0.0135232108f, 0.620714545f, 0.535165906f, -0.492116153f, -0.686599731f, 0.314748675f, 0.870424569f, -0.378542542f},
				{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};

  for (int t = 0; !std::isnan(rotMats[t][0]); t++) {
    float q[4];

    tledVectorArithmetic::ComputeQuaternionFromRotation(q, rotMats[t]);
    tledUnitTestAssert(std::fabs(Norm(q, 4) - 1) < 1e-3f);
  }
}

static void _TestQuaternionConversion(void) { 
  _TestSpecialQuaternion();
  _TestSingleComponentQuaternionConversion();
  _Test90DegQuaternion();
  _TestRandomAngleQuaternion();
}

int main(void) {
  InitUnitTest();
  
  _Test3x3Solver();
  _TestCross();
  _TestComputeAngle();
  _TestQuaternionConversion();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
