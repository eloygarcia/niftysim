// =========================================================================
// File:       tledSelfCollisionBVHCreator.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledSelfCollisionBVHCreator.h"

namespace tledSelfCollisionBVHCreator_internal {
  template <class TAABB>
  static void InitialiseChildPrimitiveSetsSpecAABB(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const TAABB &bv, const std::vector<float> &primitiveCentroids) {
    int maxAxisInd;

    {
      float maxAxisLen, axisLen;
      int aInd;

      maxAxisLen = -1;
      for (aInd = 0; aInd < 3; aInd++) {
	assert(bv.Bounds[aInd][1] - bv.Bounds[aInd][0] > 0);
	if ((axisLen = bv.Bounds[aInd][1] - bv.Bounds[aInd][0]) > maxAxisLen) {
	  maxAxisLen = axisLen;
	  maxAxisInd = aInd;
	}
      }

      {
	std::vector<int>::iterator i_minPInd, i_maxPInd, i_cPInd;	
	float maxAxisComp, minAxisComp;

	maxAxisComp = -(minAxisComp = std::numeric_limits<float>::max());
	for (i_cPInd = pIndsBegin; i_cPInd < pIndsEnd; i_cPInd++) {
	  const float axisComp = *(&primitiveCentroids.front() + 3*(*i_cPInd) + maxAxisInd);

	  if (axisComp > maxAxisComp) {
	    maxAxisComp = axisComp;
	    i_maxPInd = i_cPInd;
	  } 
	
	  if (axisComp < minAxisComp) {
	    minAxisComp = axisComp;
	    i_minPInd = i_cPInd;
	  }
	}

	std::iter_swap(i_minPInd, pIndsBegin);
	std::iter_swap(i_maxPInd, pIndsEnd - 1);
      }
    }
  }

  template <class TOBB>
  static void InitialiseChildPrimitiveSetsSpecOBB(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const TOBB &bv, const std::vector<float> &primitiveCentroids) {
    using namespace tledVectorArithmetic;

    std::vector<int>::iterator i_minPInd = pIndsBegin, i_maxPInd = pIndsEnd - 1;	
    float maxAxisComp, minAxisComp;

    assert(bv.Extents[0] >= bv.Extents[1] && bv.Extents[0] >= bv.Extents[2]);
    maxAxisComp = -(minAxisComp = std::numeric_limits<float>::max());
    for (std::vector<int>::iterator i_cPInd = pIndsBegin; i_cPInd < pIndsEnd; i_cPInd++) {
      float axisComp, dCnt[3];

      assert(bv.IsInside(&primitiveCentroids.front() + 3*(*i_cPInd)));
      Sub(dCnt, &primitiveCentroids.front() + 3*(*i_cPInd), bv.Centroid);
      axisComp = Dot(dCnt, bv.Axes[0]);
      if (axisComp > maxAxisComp) {
	maxAxisComp = axisComp;
	i_maxPInd = i_cPInd;
      } 
	
      if (axisComp < minAxisComp) {
	minAxisComp = axisComp;
	i_minPInd = i_cPInd;
      }
    }

    std::iter_swap(i_minPInd, pIndsBegin);
    std::iter_swap(i_maxPInd, pIndsEnd - 1);
  }

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledAABB<2> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledAABB<2> > &bv, const std::vector<float> &primitiveCentroids) {
    InitialiseChildPrimitiveSetsSpecAABB<tledAABB<2> >(pIndsBegin, pIndsEnd, bv, primitiveCentroids);
  }

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledAABB<4> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledAABB<4> > &bv, const std::vector<float> &primitiveCentroids) {
    InitialiseChildPrimitiveSetsSpecAABB<tledAABB<4> >(pIndsBegin, pIndsEnd, bv, primitiveCentroids);
  }

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledOBB<2> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledOBB<2> > &bv, const std::vector<float> &primitiveCentroids) {
    InitialiseChildPrimitiveSetsSpecOBB<tledOBB<2> >(pIndsBegin, pIndsEnd, bv, primitiveCentroids);
  }

  template <>
  void InitialiseChildPrimitiveSets<tledSelfCollisionBV<tledOBB<4> > >(const std::vector<int>::iterator &pIndsBegin, const std::vector<int>::iterator &pIndsEnd, const tledSelfCollisionBV<tledOBB<4> > &bv, const std::vector<float> &primitiveCentroids) {
    InitialiseChildPrimitiveSetsSpecOBB<tledOBB<4> >(pIndsBegin, pIndsEnd, bv, primitiveCentroids);
  }

#ifndef NDEBUG
  template <class TAABB>
  static void _CheckAABB(const TAABB &aabb) {
    assert(aabb.ComputeVolume() > 0);
  }

  template <class TOBB>
  static void _CheckOBB(const TOBB &obb) {
    using namespace tledVectorArithmetic;

    assert(obb.ComputeVolume() > 0);
    assert(Norm(obb.Centroid) >= 0);
    for (int c = 0; c < 3; c++) {
      assert(std::fabs(1 - Norm(obb.Axes[c])) < 1e-3f);
      assert(c == 0 || std::fabs(Dot(obb.Axes[c-1], obb.Axes[c])) < 1e-3f);
    }
  }

  template <>
  void CheckBV<tledAABB<2> >(const tledAABB<2> &bv) {
    _CheckAABB(bv);
  }

  template <>
  void CheckBV<tledAABB<4> >(const tledAABB<4> &bv) {
    _CheckAABB(bv);
  }

  template <>
  void CheckBV<tledOBB<2> >(const tledOBB<2> &bv) {
    _CheckOBB(bv);
  }

  template <>
  void CheckBV<tledOBB<4> >(const tledOBB<4> &bv) {
    _CheckOBB(bv);
  }
#endif
}
