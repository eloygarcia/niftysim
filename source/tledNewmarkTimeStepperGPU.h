// =========================================================================
// File:       tledNewmarkTimeStepperGPU.h
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
#ifndef tledNewmarkTimeStepperGPU_H
#define tledNewmarkTimeStepperGPU_H

#include "tledTimeStepperGPU.h"
#include "tledNewmarkTimeStepper.h"

/**
 * \brief GPU Newmark time integration
 * \ingroup solver
 */
class tledNewmarkTimeStepperGPU : public tledNewmarkTimeStepper<tledTimeStepperGPU> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledNewmarkTimeStepper<tledTimeStepperGPU> Superclass;
  /** @} */

  /**
   * \name Solutions
   * @{
   */
public:
  /** \brief On device solution data structure */      
  struct SolutionVariables : public Superclass::SolutionVariables {
    float4 *CurrentAccelerations;
    float4 *PreviousAccelerations;
    float4 *CurrentVelocities;
    float4 *PreviousVelocities;

  public:
    SolutionVariables(float4 *p_uCurr, float4 *p_uNext, float4 *p_aCurr, float4 *p_aPrev, float4 *p_vCurr, float4 *p_vPrev);
  };

private:
  float4 *mdp_CurrentAccelerations, *mdp_PreviousAccelerations;
  float4 *mdp_CurrentVelocities, *mdp_PreviousVelocities;
  float *mdp_M;

private:
  void _UpdateDeviceDataStructure(void);

  SolutionVariables*& _GetOnDeviceRepresentation(void) { return reinterpret_cast<SolutionVariables*&>(this->GetOnDeviceRepresentation()); }
  const SolutionVariables* _GetOnDeviceRepresentation(void) const { return reinterpret_cast<const SolutionVariables*>(this->GetOnDeviceRepresentation()); }

public:
  virtual void RetrieveSolutionFromDevice(void);  
  /** @} */

  /**
   * \name Displacement Evolution
   * @{
   */
public:
  virtual void FinishTimeStep(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledNewmarkTimeStepperGPU(const int numNodes, const float dt, const float alpha, const float M[]);
  virtual ~tledNewmarkTimeStepperGPU(void);
  /** @} */  

};

#endif
