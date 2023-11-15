// =========================================================================
// File:       tledHelper.cpp
// Purpose:    
// Package:    NiftySim: Nonlinear finite element program
// Author:     Stian Johnsen
// Language:   C++
// Created:    July 2013
// 
// Copyright (c) 2013, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// stian.johnsen.09@ucl.ac.uk
// =========================================================================
#include "tledHelper.h"

#include <iostream>
#include <cstdlib>

namespace tledHelper {
  Error& Error::operator<<(const std::string &str) {
    this->GetMsgStream() << str;

    return *this;
  }

  Error& Error::operator<<(const int num) {
    this->GetMsgStream() << num;

    return *this;
  }

  Error& Error::operator<<(const double num) {
    this->GetMsgStream() << num;

    return *this;
  }

  Error& Error::operator<<(const char c) {
    this->GetMsgStream() << c;

    return *this;    
  }

  Error& Error::operator<<(const void *ptr) {
    this->GetMsgStream() << ptr;

    return *this;    
  }

  void Error::Log(Error &r_errorStream, const char *whereFile, const int whereLine, const char *whereFnct) {
    r_errorStream.GetOutStream() << r_errorStream.GetPrefix() << " in " << whereFnct << " @ " << whereFile << ":" << whereLine << ": " << r_errorStream.GetMsgStream().str() << std::endl;
    if (r_errorStream.IsFatal()) std::abort();
  }
}
