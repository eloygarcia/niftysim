// =========================================================================
// File:       tledXMLImporter.h
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
#ifndef tledXMLImporter_H
#define tledXMLImporter_H

#include "tledModel.h"
#include "xmlParser.h"

#include <typeinfo>
#include <limits>
#include <sstream>
#include <iostream>
#include <cstdlib>

/**
 * \brief Base class for XML importers
 * \ingroup model
 * \ingroup fileio
 */
template <class TOuput>
class tledXMLImporter {
  /**
   * \name Ouput
   * @{
   */
private:
  TOuput *mp_Output;

public:
  void SetOuputObject(TOuput &r_ouput) { mp_Output = &r_ouput; }
  TOuput& GetOutput(void) { return *mp_Output; }
  const TOuput& GetOutput(void) const { return *mp_Output; }
  /** @} */

  /**
   * \name Input
   * @{
   */
private:
  const XMLNode *mpc_Root;

public:
  void SetRootNode(const XMLNode &root) { mpc_Root = &root; }
  const XMLNode& GetRootNode(void) const { return *mpc_Root; }
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  template <typename TValue>
  static TValue GetNumericAttribute(const XMLNode &node, const std::string &name, const bool abortOnError = true);

  template <typename TValue>
  static TValue GetNumericElementValue(const XMLNode &node, const bool abortOnError = true);

  /** Gets a unique (by name) child of the root node, by default aborts with an informative error messages if that node isn't found. */
  XMLNode GetUniqueChild(const std::string &name, const bool abortOnError = true) const;

  /** Returns the value of an optional attribute of the output object, or returns a default value if it isn't found. */
  template <typename TValue>
  TValue GetOptionalNumericElementValue(const std::string &name, const TValue defaultValue) const;

public:
  virtual void Import(void) = 0;
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  virtual ~tledXMLImporter(void) {}
  /** @} */
};

#include "tledXMLImporter.tpp"
#endif
