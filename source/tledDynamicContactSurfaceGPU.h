// =========================================================================
// File:       tledDynamicContactSurfaceGPU.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   CUDA
// Created:    September 2014
// 
// Copyright (c) 2014, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledDynamicContactSurfaceGPU_H
#define tledDynamicContactSurfaceGPU_H

#ifdef _GPU_

#include "tledContactSurfaceGPU.h"

/**
 * \brief Decorator for GPU dynamic contact surfaces providing some facilities useful for contact surfaces with non-static node positions.
 * \ingroup contact
 */
template <class TBaseSurface>
class tledDynamicContactSurfaceGPU : public TBaseSurface {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseSurface Superclass;
  typedef typename Superclass::Facet Facet;
  /** @} */

  /**
   * \name Device Memory
   * @{
   */
public:
  struct GPUSurface : public Superclass::GPUSurface {
    float3 *OldNodeCoordinates;
    float3 *NodeCoordinates0;
  };

private:
  float3 *mdp_OldCoordinateBuffer;

private:
  GPUSurface& _GetHostSurface(void) { return static_cast<GPUSurface&>(this->GetHostGPUSurface()); }
  const GPUSurface& _GetHostSurface(void) const { return static_cast<const GPUSurface&>(this->GetHostGPUSurface()); }
  GPUSurface* _GetDeviceSurface(void) { return static_cast<GPUSurface*>(this->GetDeviceGPUSurface()); }

protected:
  /** How many coordinate vectors need to be stored? */
  virtual int GetCoordinateHistoryBufferMultiplier(void) const = 0;

  /** Previous node positions (not necessarily last time step) buffer base pointer (R/W) */
  float3* GetOnDeviceOldNodeCoordinateBuffer(void) { return mdp_OldCoordinateBuffer; }

  /** Initial configuration node positions (R/W) */
  float3* GetAllOnDeviceNodeCoordinates0(void) { return _GetHostSurface().NodeCoordinates0; }

  virtual void AllocateDeviceSurface(void);
  virtual void CopySurfaceToDevice(void);
  virtual void InitNodes(void);

public:
  /** Previous node positions (not necessarily last time step) buffer base pointer (R/O) */
  const float3* GetAllOnDeviceOldNodeCoordinates(void) const { return _GetHostSurface().OldNodeCoordinates; }

  /** Initial configuration node positions (R/O) */
  const float3* GetAllOnDeviceNodeCoordinates0(void) const { return _GetHostSurface().NodeCoordinates0; }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void LoadFromXMLPostloadHook(void);

  tledDynamicContactSurfaceGPU(void);
  virtual ~tledDynamicContactSurfaceGPU(void);
  /** @} */
};

#include "tledDynamicContactSurfaceGPU.tpp"

#endif
#endif
