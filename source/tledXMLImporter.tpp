// =========================================================================
// File:       tledXMLImporter.tpp
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

template <class TOutput>
template <typename TValue>
TValue tledXMLImporter<TOutput>::GetNumericAttribute(const XMLNode &node, const std::string &name, const bool abortOnError) {
  if (node.getAttribute(name.c_str()) != NULL) {
    std::istringstream iss(node.getAttribute(name.c_str()));
    TValue value;
    
    iss >> value;
    if (iss.fail()) {
      if ((typeid(TValue) == typeid(float) || typeid(TValue) == typeid(double)) && iss.str() == "nan") {
	value = std::numeric_limits<TValue>::quiet_NaN();
      } else if (abortOnError) {
	std::cerr << "Parsing of value of attribute " << name << " on " << node.getName() << " failed. Value: " << node.getAttribute(name.c_str()) << std::endl;
	std::abort();
      }
    } else return value;
  } else if (abortOnError) {
    std::cerr << "No attribute " << name << " set on " << node.getName() << std::endl;
    std::abort();
  }

  return 0;
}

template <class TOutput>
template <typename TValue>
TValue tledXMLImporter<TOutput>::GetNumericElementValue(const XMLNode &node, const bool abortOnError) {
  std::istringstream iss(node.getText());
  TValue value;
    
  iss >> value;
  if (iss.fail()) {
    if ((typeid(TValue) == typeid(float) || typeid(TValue) == typeid(double)) && iss.str() == "nan") {
      value = std::numeric_limits<TValue>::quiet_NaN();
    } else if (abortOnError) {
      std::cerr << "Parsing of value of element " << node.getName() << " failed. Value: " << node.getText() << std::endl;
      std::abort();
    } 
  }

  return value;  
}

template <class TOutput>
template <typename TValue>
TValue tledXMLImporter<TOutput>::GetOptionalNumericElementValue(const std::string &name, const TValue defaultValue) const {
  TValue value = defaultValue;

  if (this->GetRootNode().nChildNode(name.c_str()) > 0) {
    XMLNode node = this->GetRootNode().getChildNode(name.c_str());

    if (node.getAttribute(name.c_str()) != NULL) {
      std::istringstream iss(node.getAttribute(name.c_str()));

      iss >> value;
    }
  }

  return value;
}

template <class TOutput>
XMLNode tledXMLImporter<TOutput>::GetUniqueChild(const std::string &name, const bool abortOnError) const {
  if (abortOnError && this->GetRootNode().nChildNode(name.c_str()) != 1) {
    std::cerr << "Expected to find exactly 1 \"" << name << "\" element in " << this->GetRootNode().getName() << ", found " << this->GetRootNode().nChildNode(name.c_str()) << std::endl;
    std::abort();
  }

  return this->GetRootNode().getChildNode(name.c_str());
}
