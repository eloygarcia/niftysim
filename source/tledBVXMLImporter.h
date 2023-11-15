// =========================================================================
// File:       tledBVXMLImporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    May 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledBVXMLImporter_H
#define tledBVXMLImporter_H

#include "tledXMLImporter.h"

#include <algorithm>

template <class TBV>
class tledBVXMLImporter : public tledXMLImporter<TBV> {
  /**
   * \name Types
   * @{
   */
public:
  typedef tledXMLImporter<TBV> Superclass;
  typedef TBV BoundingVolume;
  /** @} */

  /**
   * \name Processing
   * @{
   */
public:
  virtual void Import(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledBVXMLImporter(void) : Superclass() {}
  virtual ~tledBVXMLImporter(void) {}
  /** @} */
};

/**
 * \brief Dummy class used in BVH importers to get the right importer for a given BV type.
 */
template <class TBV>
class tledBVHImporterBVXMLImporterAdapter : public tledBVXMLImporter<TBV> {
  /**
   * \name Processing
   * @{
   */
public:
  virtual void Import(void) { std::abort(); }
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledBVHImporterBVXMLImporterAdapter(void) {}
  /** @} */  
};

#include "tledBVXMLImporter.tpp"
#endif
