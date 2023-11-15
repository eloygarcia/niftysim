// =========================================================================
// File:       tledXMLExporter.tpp
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

template <class TInput>
void tledXMLExporter<TInput>::CreateRoot(void) {
  m_XMLNode = XMLNode::createXMLTopNode(this->GetRootElementName());
}

template <class TInput>
template <typename TValue>
XMLNode tledXMLExporter<TInput>::AddNumericAttribute(XMLNode &r_node, const std::string &name, const TValue value) {
  r_node.addAttribute(name.c_str(), static_cast<std::ostringstream&>(std::ostringstream() << value).str().c_str());

  return r_node;
}

template <class TInput>
XMLNode tledXMLExporter<TInput>::CreateTextNode(const std::string &name, const std::string &value) {
  XMLNode node;

  (node = this->GetRootNode().addChild(name.c_str())).addText(value.c_str());

  return node;
}

template <class TInput>
template <typename TValue>
XMLNode tledXMLExporter<TInput>::CreateNumericNode(const std::string &name, const TValue value) {
  std::ostringstream oss;

  oss << value;

  return this->CreateTextNode(name, oss.str());
}

template <class TInput>
template <typename TValue>
XMLNode tledXMLExporter<TInput>::CreateNumericListNode(const std::string &name, const TValue *values, const int numValues) {
  std::ostringstream oss;

  std::copy(values, values + numValues, std::ostream_iterator<TValue>(oss, "\n"));

  return this->CreateTextNode(name, oss.str());
}

template <class TInput>
template <typename TValue>
XMLNode tledXMLExporter<TInput>::CreateNumericListNode(const std::string &name, const std::vector<TValue> &values) {
  return this->CreateNumericListNode(name, &values.front(), values.size());
}

template <class TInput>
void tledXMLExporter<TInput>::Export(void) {
  assert(mpc_Input != NULL);
  this->CreateRoot();
  this->WriteBody();
}
