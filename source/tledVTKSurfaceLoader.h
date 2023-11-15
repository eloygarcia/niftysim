// =========================================================================
// File:       tledVTKSurfaceLoader.h
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
#ifndef tledVTKSurfaceLoader_H
#define tledVTKSurfaceLoader_H

#include "tledHelper.h"
#include "tledSurfaceLoader.h"

#include <iostream>
#include <algorithm>
#include <cassert>

#ifdef _Visualisation_
#include <vtkTriangleFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkCell.h>
#include <vtkSmartPointer.h>
#include <vtkErrorCode.h>
#endif

/**
 * \brief Reader for VTK surface files
 * \ingroup import
 */
template <class TSurface>
class tledVTKSurfaceLoader : public tledSurfaceLoader<TSurface> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledSurfaceLoader<TSurface> Superclass;
  typedef TSurface Surface;
  /** @} */

  /**
   * \name Loading
   * @{
   */
protected:
  virtual void ReadFile(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledVTKSurfaceLoader(void) {}
  /** @} */
};

#include "tledVTKSurfaceLoader.tpp"
#endif
