// =========================================================================
// File:       tledNodeRejector.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    October 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledNodeRejector.h"
#include "tledBoxNodeRejector.h"
#include "tledSphereNodeRejector.h"

#include <string>

void tledNodeRejector::ProceedWithNextRejector() { 
  if (this->HasNext()) {
    this->GetNext().SetNodeIndices(this->GetNodeIndices());
    this->GetNext().SetNodes(this->GetNodes());
    this->GetNext().RunRejection();
  }
}

void tledNodeRejector::RunRejection() {
  std::vector<int> retain;

  retain.reserve(this->GetNodeIndices().size());
  for (std::vector<int>::const_iterator ic_n = this->GetNodeIndices().begin(); ic_n < this->GetNodeIndices().end(); ic_n++) {
    if (!this->DoReject(*ic_n)) retain.push_back(*ic_n);
  }
  this->GetNodeIndices().swap(retain);
  this->ProceedWithNextRejector();
}

void tledNodeRejector::AppendRejectorToChain(tledNodeRejector *p_rejector) {
  if (!this->HasNext()) mp_Next = p_rejector;
  else mp_Next->AppendRejectorToChain(p_rejector);
}

tledNodeRejector::~tledNodeRejector() {
  if (this->HasNext()) {
    delete mp_Next;
  }
}

tledNodeRejector* tledNodeRejector::CreateRejector(const XMLNode rootNode) {
  tledNodeRejector *p_rejector = NULL;

  if (rootNode.getAttribute("Type") == NULL) {
    tledFatalError("No node-rejector type specified.");
  } else {
    std::string type = rootNode.getAttribute("Type");

    if (type == "Box") {
      p_rejector = new tledBoxNodeRejector();
    } else if (type == "Sphere") {
      p_rejector = new tledSphereNodeRejector();
    } else {
      tledLogErrorStream(tledHelper::FatalError() << "Node-rejector type " << type << " not recognised.");
    }

    if (p_rejector != NULL) p_rejector->InitialiseFromXMLSpec(rootNode);
  }

  return p_rejector;
}
