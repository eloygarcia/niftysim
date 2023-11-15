// =========================================================================
// File:       tledBVHGPU.h
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
#ifndef tledBVHGPU_H
#define tledBVHGPU_H

#include "tledBVH.h"
#include "tledCUDAHelpers.h"

class tledBVHGPU {
  /**
   * \name On-Device Representation
   * @{
   */
protected:
  /** Device memory allocation and initialisation hook. */
  virtual void InitDeviceMemory(void) = 0;

public:
  /** Copies the BVs to the GPU */
  virtual void CopyBVsToGPU(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
protected:
  tledBVHGPU(void) {}

public:
  virtual ~tledBVHGPU(void);
  /** @} */
};

/**
 * \brief Adds GPU-specifc member functions to BVH API
 * \ingroup contact
 */
template <class TBaseBVH>
class tledBVHImplGPU : public TBaseBVH {
  /**
   * \name Types
   * @{
   */
public:
  typedef TBaseBVH Superclass;
  typedef typename Superclass::BoundingVolume BoundingVolume;
  typedef typename Superclass::ContactMesh ContactMesh;
  /** @} */

  /**
   * \name Device Representation
   * @{
   */
private:
  typename BoundingVolume::GPUBV *mhp_GPUBVs, *mdp_GPUBVs;

protected:
  virtual void InitDeviceMemory(void);

public:
  /** R/W access to host rep. of GPU BVs, can be used in conjunction with CopyBVsToGPU to modify the GPU BVH on the host. */
  typename BoundingVolume::GPUBV& GetHostGPUBV(const int bvInd) { return mhp_GPUBVs[bvInd]; }

  typename BoundingVolume::GPUBV* GetOnDeviceBVs(void) { return mdp_GPUBVs; }
  const typename BoundingVolume::GPUBV* GetOnDeviceBVs(void) const { return mdp_GPUBVs; }

  virtual void CopyBVsToGPU(void);  

  virtual void Init(tledBVHCreator &r_bvhBuilder);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual void LoadFromXMLPostloadHook(void);

  tledBVHImplGPU(ContactMesh &r_mesh) : Superclass(r_mesh) {}
  virtual ~tledBVHImplGPU(void);
  /** @} */
};

#include "tledBVHGPU.tpp"
#endif
