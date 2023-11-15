// =========================================================================
// File:       tledContactVolumeSurfaceExtractor.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    March 2012
// 
// Copyright (c) 2012, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledContactVolumeSurfaceExtractor_H
#define tledContactVolumeSurfaceExtractor_H

#include "tledVolumeSurfaceExtractor.h"
#include "tledContactSurfaceCreator.h"

#include <algorithm>

/**
 * \brief Surface extractor for contact meshes.
 * \ingroup surface
 */
template <class TSurface> 
class tledContactVolumeSurfaceExtractor : public tledContactSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledContactSurfaceCreator<TSurface> Superclass;
  typedef TSurface Surface;
  /** @} */

  /**
   * \name Main 
   * @{
   */
private:
  class _BaseExtractor : public tledVolumeSurfaceExtractor<Surface> {
  protected:
    virtual void InitNodes();
    
  public:
    _BaseExtractor(const tledMesh &mesh) : tledVolumeSurfaceExtractor<Surface>(mesh) {}
  };

private:
  _BaseExtractor m_BaseExtractor;

public:
  virtual void Create(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledContactVolumeSurfaceExtractor(const tledMesh &mesh) : Superclass(m_BaseExtractor), m_BaseExtractor(mesh) {}

  virtual ~tledContactVolumeSurfaceExtractor(void) {}
  /** @} */
};

#include "tledContactVolumeSurfaceExtractor.tpp"
#endif
