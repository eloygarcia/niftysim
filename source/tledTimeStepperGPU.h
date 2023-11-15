// =========================================================================
// File:       tledTimeStepperGPU.h
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
#ifndef tledTimeStepperGPU_H
#define tledTimeStepperGPU_H

#include "tledTimeStepper.h"

#include <vector_types.h>

/**
 * \brief GPU Time ODE solver
 * \ingroup solver
 */
class tledTimeStepperGPU : public tledTimeStepper {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledTimeStepper Superclass;
  /** @} */

  /**
   * \name Solutions
   * @{
   */
public:
  /** \brief On device solution data structure */      
  struct SolutionVariables {
    float4 *CurrentDisplacements;
    float4 *NextDisplacements;
  };

private:
  float4 *mp_HostDeviceCopyBuffer;
  float4 *mdp_CurrentDisplacements, *mdp_NextDisplacements;
  SolutionVariables *mdp_OnDeviceRepresentation;

protected:
  /** Retrieves a float4 array from the device and converts it to a standard 3-components-per-node host float array */
  static void ConvertToHostFloatArray(float *p_outBuffer, float4 *hp_tmpBuffer, const float4 *dp_deviceArray, const int numNodes);

  float4*& GetHostDeviceCopyBuffer(void) { return mp_HostDeviceCopyBuffer; }
  const float4* GetHostDeviceCopyBuffer(void) const { return mp_HostDeviceCopyBuffer; }

  SolutionVariables*& GetOnDeviceRepresentation(void) { return mdp_OnDeviceRepresentation; }
  const SolutionVariables* GetOnDeviceRepresentation(void) const { return mdp_OnDeviceRepresentation; }

  float4*& GetOnDeviceNextDisplacements(void) { return mdp_NextDisplacements; }
  float4*& GetOnDeviceCurrentDisplacements(void) { return mdp_CurrentDisplacements; }

public:
  const float4* GetOnDeviceNextDisplacements(void) const { return mdp_NextDisplacements; }
  const float4* GetOnDeviceCurrentDisplacements(void) const { return mdp_CurrentDisplacements; }

  /** Retrieves the on-device solutions and copies them to their respective host buffers */
  virtual void RetrieveSolutionFromDevice(void);  

  SolutionVariables* GetSolutionDeviceBuffer(void) { return mdp_OnDeviceRepresentation; }

  virtual void SetCurrentDisplacements(const float U[]);
  /** @} */

  /**
   * \name Displacement Evolution
   * @{
   */
public:
  virtual void FinishTimeStep(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledTimeStepperGPU(const int numNodes);
  virtual ~tledTimeStepperGPU(void);
  /** @} */  
};
#endif
