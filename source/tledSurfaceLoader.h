// =========================================================================
// File:       tledSurfaceLoader.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2011
// 
// Copyright (c) 2011, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledSurfaceLoader_H
#define tledSurfaceLoader_H

#include "tledBasicMeshFileReader.h"
#include "tledBasicSurfaceCreator.h"
#include "tledMatrixFunctions.h"
#include "tledVectorArithmetic.h"
#include "tledSurface.h"

/**
 * \brief Surface mesh file-reader base class
 * \ingroup surface
 */
template <class TSurface>
class tledSurfaceLoader : public tledBasicMeshFileReader<TSurface> {
public:
  typedef TSurface Surface;

  /**
   * \name Transforms
   * @{
   */
protected:
  virtual void ApplyTransforms(void);
  /** @} */

  /**
   * \name FS I/O
   * @{
   */
protected:
  /** 
   * Reads the mesh from file has to be implemented by all format-specific subclasses.<br>
   * Application of transforms, etc. happens after mesh loading, in the Read member function and does not have to be implemented in subclasses. 
  */
  virtual void ReadFile(void) = 0;

public:
  virtual void Read(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSurfaceLoader(void) {}
  virtual ~tledSurfaceLoader(void) {}
  /** @} */
};

/**
 * \brief Takes an initialised surface loader object and makes it covariant with tledBasicSurfaceCreator so that it can be used e.g. with tledContactSurfaceCreator
 * \ingroup surface
 *
 * The output buffer must be set on the adapter, not on the loader!
 */
template <class TSurface>
class tledSurfaceLoaderSurfaceCreatorAdapter : public tledBasicSurfaceCreator<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledBasicSurfaceCreator<TSurface> Superclass;
  typedef tledSurfaceLoader<TSurface> Loader;
  /** @} */

  /**
   * \name Loader
   * @{
   */
private:
  Loader *mp_Loader;

public:
  Loader& GetLoader(void) { return *mp_Loader; }
  const Loader& GetLoader(void) const { return *mp_Loader; }
  /** @} */
  
  /**
   * \name Main 
   * @{
   */
public:
  virtual void Create(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledSurfaceLoaderSurfaceCreatorAdapter(Loader &r_loader) : mp_Loader(&r_loader) {}
  virtual ~tledSurfaceLoaderSurfaceCreatorAdapter(void) {}
  /** @} */
};

#include "tledSurfaceLoader.tpp"
#endif
