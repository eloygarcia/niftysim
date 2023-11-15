// =========================================================================
// File:       testBVHTraverserCPU.cpp
// Purpose:    tledBVHTraverserCPU unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledVectorArithmetic.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledDeformableDeformableBVHTraverserCPU.h"
#include "tledSelfCollisionBVH.h"

#include <iostream>

using namespace tledVectorArithmetic;

typedef tledDeformableContactSurfaceT3CPU _TestSurface;
typedef tledSelfCollisionBVHImpl<_TestSurface, tledAABB<2> > _TestBVH;

class _TestBVHTraverserCPU : public tledDeformableDeformableBVHTraverserImplCPU<_TestBVH> {
public:
  typedef tledDeformableDeformableBVHTraverserImplCPU<_TestBVH> Superclass;

  /**
   * \name Tests
   * @{
   */
public:
  static void RunDiscreteEdgeEdgeRandomTest(void) {
    const int numTests = 10000;

    int numFail = 0;

    for (int t = 0; t < numTests; t++) {
      float r, q, testR, testQ;
      float u[3], v[3], A[3], B[3], C[3], D[3], d[3], ab[3], cd[3];

      for (int c = 0; c < 3; c++) {
	C[c] = 10.0f*(drand48() - drand48());

	u[c] = 2*(drand48() - drand48());
	v[c] = 2*(drand48() - drand48());
      }
      r = 0.98f*drand48() + 0.01f;	
      q = 0.98f*drand48() + 0.01f;	

      Cross(d, u, v);      
      ScalarMul(d, 0.5f*drand48()/Norm(d));
      assert(std::fabs(Dot(u, d)) < 1e-4f*Dot(u, u));
      assert(std::fabs(Dot(v, d)) < 1e-4f*Dot(v, v));
      
      Add(cd, C, ScalarMul(cd, v, q));
      Add(ab, d, cd);

      Add(D, v, C);      
      Sub(A, ab, ScalarMul(A, u, r));
      Add(B, A, u);           

      if (Superclass::ComputeEdgeEdgeClosestPointParameters(testR, testQ, A, B, C, D)) {
	if (std::fabs(std::fabs(Dot(u, v))/(Norm(u)*Norm(v)) - 1) < 1e-3f) {
	  tledUnitTestAssert(testR >= 0 && testR <= 1 && testQ >= 0 && testQ <= 1);
	} else {
	  numFail += !(std::fabs(testR - r) < 1e-2f*std::fabs(r) || std::fabs(testR - r) < 1e-4f) || !(std::fabs(testQ - q) < 1e-2f*std::fabs(q) || std::fabs(testR - q) < 1e-4f);
	}
      } else numFail += 1;
    } /* for tests */

    std::cout << numFail << " failed out of " << numTests << std::endl;
    tledUnitTestAssert(numFail < numTests/1000); 
  } 

  static void RunDiscreteEdgeEdgeSpecialTests(void) {
    const float As[][3] = {{7.18895f, 8.84292f, 1.47905f},
			   {4.14093399f, 4.791646f, -6.62481356f},
			   {-3.61556f, -1.88081f, 0.645997f},
			   {0.0393397f, -2.62175f, -0.243502f},
			   {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
    const float Bs[][3] = {{6.50241f, 9.05044f, 1.36094f},
			   {3.33034968f, 5.35124493f, -7.12709236f},
			   {-2.03036f, -2.84619f, -0.0644127f},
			   {0.536913f, -2.78745f, -1.89951f},
			   {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
    const float Cs[][3] = {{7.41181f, 8.34768f, 1.52679f},
			   {3.07415318f, 5.135252f, -6.83298635f},
			   {-3.30588f, -2.40231f, 0.470821f},
			   {0.545805f, -2.86354f, -0.99453f},
			   {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
    const float Ds[][3] = {{6.57105f, 8.60311f, 1.42579f},
			   {3.47620273f, 4.85968208f, -6.58247614f},
			   {-3.40658f, -2.31192f, 0.0492441f},
			   {-0.523466f, -3.07054f, -0.761424f},
			   {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}};
    const float rs[] = {0.140762f,
			0.5f,
			0.284597f,
			0.444304f,
			std::numeric_limits<float>::quiet_NaN()};
    const float qs[] = {0.522909f,
			std::numeric_limits<float>::quiet_NaN(),
			0.0991764f,
			0.220797f,
			std::numeric_limits<float>::quiet_NaN()};

    for (int t = 0; !std::isnan(rs[t]); t++) {
      float testR, testQ;

      tledUnitTestAssert(Superclass::ComputeEdgeEdgeClosestPointParameters(testR, testQ, As[t], Bs[t], Cs[t], Ds[t]));
      tledUnitTestAssert(std::fabs(testR - rs[t]) < 1e-2f*std::fabs(rs[t]) || std::fabs(testR - rs[t]) < 1e-4f);
      tledUnitTestAssert((std::isnan(qs[t]) && testQ > 0 && testQ < 1)
			 || std::fabs(testQ - qs[t]) < 1e-2f*std::fabs(qs[t]) || std::fabs(testR - qs[t]) < 1e-4f);
    }
  } 
  /** @} */
};

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestBVHTraverserCPU::RunDiscreteEdgeEdgeSpecialTests();
  _TestBVHTraverserCPU::RunDiscreteEdgeEdgeRandomTest();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
