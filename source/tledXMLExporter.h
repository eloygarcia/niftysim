// =========================================================================
// File:       tledXMLExporter.h
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    April 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#ifndef tledXMLExporter_H
#define tledXMLExporter_H

#include "xmlParser.h"

#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <vector>

/**
 * \brief Base class for exporters for pre-computed data for repeated simulations.
 * \ingroup model
 * \ingroup fileio
 */
template <class TInput>
class tledXMLExporter {
  /**
   * \name XML Parser Interface
   * @{
   */
private:
  XMLNode m_XMLNode;

public:
  virtual const char* GetRootElementName(void) const = 0;

  /** 
   * \brief Root node for the input object's XML representation. 
   * 
   * The XML tree is not destroyed with exporter object; it is the client's responsibility to release the memory
   * allocated for the XML representation through XMLNode::deleteNodeContent, when it is done with it.
   */ 
  const XMLNode& GetRootNode(void) const { return m_XMLNode; }
  XMLNode& GetRootNode(void) { return m_XMLNode; }
  /** @} */

  /**
   * \name Input
   * @{
   */
private:
  const TInput *mpc_Input;

public:
  const TInput& GetInput(void) const { return *mpc_Input; }
  void SetInput(const TInput &input) { mpc_Input = &input; }
  /** @} */

  /**
   * \name Processing
   * @{
   */
protected:
  void CreateRoot(void);  
  virtual void WriteBody(void) = 0;
  
  template <typename TValue>
  static XMLNode AddNumericAttribute(XMLNode &r_node, const std::string &attribName, const TValue value);

  /** Creates a textual XML node directly below the root */
  XMLNode CreateTextNode(const std::string &name, const std::string &value);
  
  /** Creates an XML node holding a numeric value (int, float, ..) directly below the root */
  template <typename TValue>
  XMLNode CreateNumericNode(const std::string &name, const TValue value);

  /** Creates an XML node holding a list of numeric value (int, float, ..) directly below the root */
  template <typename TValue>
  XMLNode CreateNumericListNode(const std::string &name, const std::vector<TValue> &values);

  /** Creates an XML node holding a list of numeric value (int, float, ..) directly below the root */
  template <typename TValue>
  XMLNode CreateNumericListNode(const std::string &name, const TValue *values, const int numValues);

public:
  virtual void Export(void);
  /** @} */

  /**
   * \name Construction, Destruction
   * @{
   */
public:
  tledXMLExporter(void) : mpc_Input(NULL) {}
  virtual ~tledXMLExporter(void) {}
  /** @} */
};

#include "tledXMLExporter.tpp"
#endif
