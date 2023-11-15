// =========================================================================
// File:       tledCentralDifferenceTimeStepperGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledCentralDifferenceTimeStepperGPU_H
#define tledCentralDifferenceTimeStepperGPU_H

#include "tledTimeStepperGPU.h"
#include "tledCentralDifferenceTimeStepper.h"

/**
 * \brief GPU central difference time integration
 * \ingroup solver
 */
class tledCentralDifferenceTimeStepperGPU : public tledCentralDifferenceTimeStepper<tledTimeStepperGPU> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledCentralDifferenceTimeStepper<tledTimeStepperGPU> Superclass;
  /** @} */

  /**
   * \name Solutions
   * @{
   */
public:
  struct SolutionVariables : public Superclass::SolutionVariables {
    float4 *PreviousDisplacements;

  public:
    SolutionVariables(float4 *dp_UCurr, float4 *dp_UNext, float4 *dp_UPrev);
  };

private:
  float4 *mdp_CD;
  float4 *mdp_PreviousDisplacements;

private:
  SolutionVariables*& _GetOnDeviceRepresentation(void) { return reinterpret_cast<SolutionVariables*&>(this->GetOnDeviceRepresentation()); }
  const SolutionVariables* _GetOnDeviceRepresentation(void) const { return reinterpret_cast<const SolutionVariables*>(this->GetOnDeviceRepresentation()); }

public:
  virtual void RetrieveSolutionFromDevice(void);
  /** @} */

  /**
   * \brief Displacement Evolution
   * @{
   */
private:
  void _UpdateDeviceDataStructure(void);

public:
  virtual void FinishTimeStep(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledCentralDifferenceTimeStepperGPU(const int numNodes, const float dt, const float alpha, const float M[]);
  virtual ~tledCentralDifferenceTimeStepperGPU(void);
  /** @} */
};

#endif
