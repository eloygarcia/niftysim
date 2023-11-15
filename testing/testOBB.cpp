// =========================================================================
// File:       testOBB.cpp
// Purpose:    tledOBB unit test
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledUnitTest.h"
#include "tledOBB.h"
#include "tledVectorArithmetic.h"

static void _TestBasicTriangleFit() {
  using namespace tledVectorArithmetic;

  const int nodeInds[][3] = {{0, 1, 2},
			     {0, 1, 2},
			     {2, 0, 1},
			     {1, 0, 2},			     
			     {-1, -1, -1}};
  const float axes[][2][3] = {{{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}},
			      {{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}},
			      {{-0.6651f, -0.7395f, -0.1037f}, {0.1121f, -0.2362f, 0.9652f}},
			      {{-0.6708f, -0.6615f, -0.3354f}, {0.2894f, -0.6499f, 0.7028f}}};
  const float origins[][3] = {{0.f, 0.f, 0.f},
			      {1.f, 2.f, 3.f},
			      {-1.3f, -2.923f, 10.539f},
			      {-3.631f, 4.27f, 89.351f}};
  const float axisWeights[][3] = {{2.f, 1.f, 0.5f},
				  {2.f, 1.f, 0.5f},
				  {2.4f, 0.9f, 0.35f},
				  {3.f, 0.5f, 0.71f}};

  for (int t = 0; nodeInds[t][0] >= 0; t++) {
    float vertices[3][3], tmp[3];
    tledOBB<2> obb;

    std::copy(origins[t], origins[t] + 3, vertices[0]);
    Add(vertices[1], ScalarMul(tmp, axes[t][0], axisWeights[t][0]), origins[t]);
    Add(vertices[2], Add(vertices[2], ScalarMul(tmp, axes[t][1], axisWeights[t][1]), ScalarMul(vertices[2], axes[t][0], axisWeights[t][0]*axisWeights[t][2])), origins[t]);
    obb.ComputeFromNodeList(nodeInds[t], nodeInds[t] + 3, &vertices[0][0]);
    obb.AddMargin(1e-4f);
    
    tledUnitTestAssert(std::fabs(Dot(obb.Axes[0], obb.Axes[1])) < 1e-3f);
    tledUnitTestAssert(std::fabs(Dot(obb.Axes[0], obb.Axes[2])) < 1e-3f);
    tledUnitTestAssert(std::fabs(Dot(obb.Axes[1], obb.Axes[2])) < 1e-3f);

    for (int a = 0; a < 2; a++) {
      tledUnitTestAssert(std::fabs(Dot(obb.Axes[a], axes[t][a]) - 1.f) < 1e-3 || std::fabs(Dot(obb.Axes[a], axes[t][a]) + 1.f) < 1e-3);
    }

    for (int v = 0; v < 3; v++) {
      tledUnitTestAssert(obb.IsInside(vertices[v]));
    }
  }
}

static void _TestMergeContained() {
  using namespace tledVectorArithmetic;

  static const float baseAxes[][2][2][3] = {{{{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}},
					     {{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}}},
					    {{{0.6651f, 0.7395f, 0.1037f}, {-0.1121f, 0.2362f, -0.9652f}},
					     {{0.6651f, 0.7395f, 0.1037f}, {-0.1121f, 0.2362f, -0.9652f}}},
					    {{{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}},
					     {{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}}}};
  static const float extents[][2][3] = {{{1.1f, 0.2f, 0.1f},
					 {0.55f, 0.1f, 0.05f}},
					{{2.1f, 1.2f, 1.1f},
					 {1.55f, 1.1f, 0.55f}}};
  static const float centroids[][2][3] = {{{1.f, 2.f, 3.f},
					   {1.f, 2.f, 3.f}},
					  {{1.f, 2.f, 3.f},
					   {1.1f, 2.1f, 3.1f}}};

  for (int t = 0; !std::isnan(baseAxes[t][0][0][0]); t++) {
    tledOBB<2> obb0, obb1, merged;
    float d[3];

    for (int a = 0; a < 2; a++) {
      ScalarDiv(obb0.Axes[a], baseAxes[t][0][a], Norm(baseAxes[t][0][a]));
      ScalarDiv(obb1.Axes[a], baseAxes[t][1][a], Norm(baseAxes[t][1][a]));
    }
    assert(std::fabs(Dot(obb0.Axes[1], obb0.Axes[0])) < 1e-4f);
    assert(std::fabs(Dot(obb1.Axes[1], obb1.Axes[0])) < 1e-4f);
    
    Cross(obb0.Axes[2], obb0.Axes[0], obb0.Axes[1]);
    Cross(obb1.Axes[2], obb1.Axes[0], obb1.Axes[1]);

    std::copy(centroids[t][0], centroids[t][0] + 3, obb0.Centroid);
    std::copy(centroids[t][1], centroids[t][1] + 3, obb1.Centroid);

    std::copy(extents[t][0], extents[t][0] + 3, obb0.Extents);
    std::copy(extents[t][1], extents[t][1] + 3, obb1.Extents);

    merged = obb0;
    tledOBB<2>::Merge(merged, obb1);
    
    for (int a = 0; a < 3; a++) {
      tledUnitTestAssert(Norm(Sub(d, obb0.Axes[a], merged.Axes[a])) < 1e-3f);
    }
    tledUnitTestAssert(Norm(Sub(d, obb0.Extents, merged.Extents)) < 1e-1f);
    tledUnitTestAssert(Norm(Sub(d, obb0.Centroid, merged.Centroid)) < 1e-3f || Norm(Sub(d, obb0.Centroid, obb1.Centroid)) > 1e-3f);
  }
}

static void _TestLargeOBBFit() {
  using namespace tledVectorArithmetic;

  static const int numPts = 30;
  static const float boxAxes[][2][3] = {{{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}},
					{{0.5699f, 0.2601f, 0.7794f}, {0.5395f, -0.8339f, -0.1162f}},
					{{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}}};
  static const float centres[][3] = {{0.f, 1.f, 2.f},
				     {0.92f, 9.26f, 3.71f}};
  static const float extents[][3] = {{3.f, 2.f, 1.f},
				     {102.3f, 50.32f, 1.432f}};
  static const std::vector<int> inds = tledSequenceGenerator::MakeSequence(0, numPts);
  
  for (int t = 0; !std::isnan(boxAxes[t][0][0]); t++) {
    float pts[numPts][3], axes[3][3];    
    tledOBB<2> obb;

    for (int a = 0; a < 2; a++) {
      assert(a == 0 || std::fabs(Dot(boxAxes[t][a-1], boxAxes[t][a])) < 1e-3f);
      ScalarDiv(axes[a], boxAxes[t][a], Norm(boxAxes[t][a]));
    }
    Cross(axes[2], axes[0], axes[1]);

    for (int p = 0; p < numPts; p++) {
      std::copy(centres[t], centres[t] + 3, pts[p]);

      for (int c = 0; c < 3; c++) {
	float tmp[3];

	Add(pts[p], pts[p], ScalarMul(tmp, axes[c], (1 - 2*(rand()%2))*float((drand48()/4 + 0.75f)*extents[t][c])));	
      }
    }
    obb.ComputeFromNodeList(&inds.front(), &inds.back() + 1, &pts[0][0]);
    
    for (int p = 0; p < numPts; p++) {
      tledUnitTestAssert(obb.IsInside(pts[p]));
    }
    
    for (int a = 0; a < 3; a++) {
      tledUnitTestAssert(std::fabs(Dot(obb.Axes[a], axes[a])) >= .75f);
    }
  }
}

static void _TestConstruction() {
  _TestBasicTriangleFit();
  _TestMergeContained();
  _TestLargeOBBFit();
}

static void _TestParallelAxes() {
  using namespace tledVectorArithmetic;

  static const float boxAxes[][2][3] = {{{1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}},
					{{0.5699f, 0.2601f, 0.7794f}, {0.5395f, -0.8339f, -0.1162f}},
					{{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}, {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}}};
  static const float extents[][3] = {{12.1f, 4.54f, 3.71f},
				     {50.1f, 19.2f, 5.1f}};
  static const float centres[][3] = {{0.f, 0.f, 0.f},
				     {4.98f, 1.38f, 0.97f}};

  for (int t = 0; !std::isnan(boxAxes[t][0][0]); t++) {
    tledOBB<2> obb0, obb1;

    for (int a = 0; a < 2; a++) {
      assert(a == 0 || std::fabs(Dot(boxAxes[t][a], boxAxes[t][a-1])) < 1e-4f);
      ScalarDiv(obb0.Axes[a], boxAxes[t][a], Norm(boxAxes[t][a]));
    }
    Cross(obb0.Axes[2], obb0.Axes[0], obb0.Axes[1]);
    std::copy(extents[t], extents[t] + 3, obb0.Extents);
    std::copy(centres[t], centres[t] + 3, obb0.Centroid);
    obb1 = obb0;

    ScalarMul(obb1.Extents, 0.5f);
    for (int a = 0; a < 3; a++) {
      float d[3];

      ScalarMul(d, obb0.Axes[a], 1.05f*(obb0.Extents[a] + obb1.Extents[a]));
      Add(obb1.Centroid, obb0.Centroid, d);
      tledUnitTestAssert(!obb0.DoesIntersect(obb1));
      
      ScalarMul(d, obb0.Axes[a], 0.95f*(obb0.Extents[a] + obb1.Extents[a]));
      Add(obb1.Centroid, obb0.Centroid, d);
      tledUnitTestAssert(obb0.DoesIntersect(obb1));
    }
  }
}

static void _TestOverlapTest() {
  _TestParallelAxes();
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestConstruction();
  _TestOverlapTest();

  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
}
