// =========================================================================
// File:       testContactSolverCPU.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    November 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================

#include "tledUnitTest.h"
#include "tledSimulator.h"
#include "tledTimeStepperCPU.h"
#include "tledDeformableContactSurfaceCPU.h"
#include "tledCentralDifferenceTimeStepperCPU.h"
#include "tledDeformableDeformableContactSolverCPU.h"
#include "tledVectorArithmetic.h"

#include <cstdlib>
#include <string>

class _TestRateResponseContactSolverCPU : public tledDeformableDeformableContactSolverImplCPU<tledDeformableContactSurfaceT3CPU> {
public:
  typedef NodeFacetConstraintItem Constraint;

public:
  void ComputeTestRateResponses(float *p_f, NodeFacetConstraintItem &r_ci, const float uNexts[], const float uCurrs[], const int nodeIndex, const int facetIndex) {
    float xi[3];

    this->GetMesh().ProjectOntoFacetCached(xi, GetMesh().GetNodeCoordinates(nodeIndex), facetIndex);
    this->GetMesh().ComputeShapeValues(r_ci.ShapeValues, xi[0], xi[1]);    
    r_ci.ContactNodeIndices[0] = nodeIndex;
    for (int i = 0; i < 3; i++) r_ci.ContactNodeIndices[1+i] = GetMesh().GetFacet(facetIndex).NodeIndices[i];
    std::copy(GetMesh().GetFacetNormalCached(facetIndex), GetMesh().GetFacetNormalCached(facetIndex) + 3, r_ci.Normal);
    r_ci.GapValue = 0;
    ComputeNodeFacetRateResponse<true, true>(r_ci, uNexts, uCurrs);
    
    for (int i = 0; i < 4; i++) {
      tledVectorArithmetic::Add(p_f + 3*this->GetMesh().MapSurface2VolumeNode(r_ci.ContactNodeIndices[i]), p_f + 3*this->GetMesh().MapSurface2VolumeNode(r_ci.ContactNodeIndices[i]), r_ci.ContactForces[i]);
    }
  }

public:
  _TestRateResponseContactSolverCPU(tledUnstructuredContactManager &r_contactRes) : tledDeformableDeformableContactSolverImplCPU<tledDeformableContactSurfaceT3CPU>(r_contactRes) {
    if (!this->DoSelfCollision()) this->ToggleDoSelfCollision();
    this->Init();
  }
};

static void _TestTetraRateConstraint(void) {
  using namespace tledVectorArithmetic;

  tledModel model(tledUnitTest::GetResourcePath("contact_2_tetras.xml").c_str());
  tledSimulator sim(&model);
  tledUnstructuredContactManager &r_manager = sim.GetContactManager()->GetUnstructuredContactManager();
  _TestRateResponseContactSolverCPU testSolver(r_manager);
  tledCentralDifferenceTimeStepperCPU &r_stepper = static_cast<tledCentralDifferenceTimeStepperCPU&>(static_cast<tledSolverCPU*>(sim.GetSolver())->GetTimeStepper());
  _TestRateResponseContactSolverCPU::Constraint ci;
  std::vector<float> f(model.GetNumNodes()*3, 0);
  float vRel[3], vRelStart[3], diff[3];

  std::fill(r_stepper.GetNextDisplacements(), r_stepper.GetNextDisplacements() + 3*8, 0.f);
  std::fill(r_stepper.GetCurrentDisplacements(), r_stepper.GetCurrentDisplacements() + 3*8, 0.f);
  for (int i = 0; i < 4; i++) {
    r_stepper.GetPreviousDisplacements()[3*i] = 1.0f;
    r_stepper.GetCurrentDisplacements()[3*i] = 1.0f;
    r_stepper.GetNextDisplacements()[3*i] = 1.0765f;
  }

  testSolver.ComputeTestRateResponses(&f.front(), ci, r_stepper.GetNextDisplacements(), r_stepper.GetCurrentDisplacements(), 1, 6);

  Sub(vRelStart, r_stepper.GetNextDisplacements() + 3*ci.ContactNodeIndices[0], r_stepper.GetCurrentDisplacements() + 3*ci.ContactNodeIndices[0]);
  for (int v = 0; v < 3; v++) {
    float tmp[3];

    Sub(vRelStart, vRelStart, ScalarMul(tmp, Sub(tmp, r_stepper.GetNextDisplacements() + 3*ci.ContactNodeIndices[1+v], r_stepper.GetCurrentDisplacements() + 3*ci.ContactNodeIndices[1+v]), ci.ShapeValues[v]));
  }

  for (std::vector<float>::iterator i_f = f.begin(); i_f < f.end(); i_f++) *i_f *= -1;
  r_stepper.EvolveDisplacements(&f.front());

  Sub(vRel, r_stepper.GetNextDisplacements() + 3*ci.ContactNodeIndices[0], r_stepper.GetCurrentDisplacements() + 3*ci.ContactNodeIndices[0]);
  for (int v = 0; v < 3; v++) {
    float tmp[3];

    Sub(vRel, vRel, ScalarMul(tmp, Sub(tmp, r_stepper.GetNextDisplacements() + 3*ci.ContactNodeIndices[1+v], r_stepper.GetCurrentDisplacements() + 3*ci.ContactNodeIndices[1+v]), ci.ShapeValues[v]));
  }

  tledUnitTestAssert(Norm(Add(diff, vRel, vRelStart)) < 1e-2*Norm(vRelStart));
}

int main(int argc, char *argv[]) {
  tledUnitTest::InitUnitTest();

  _TestTetraRateConstraint();
  tledUnitTestPrintSuccess;

  return EXIT_SUCCESS;
} /* main */
