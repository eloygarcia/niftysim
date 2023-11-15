// =========================================================================
// File:       tledConstraint.cpp
// Purpose:    Constraint base class
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    June 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================


#include "tledConstraint.h"

#include <cmath>
#include <limits>
#include <cstring>

using namespace std;

vector<int> tledConstraint::emptyIntVec;
vector<float> tledConstraint::emptyFloatVec;

float tledConstraint::ComputeAmplitude(const double TR, const enum loadShape ls) const
{
   float Amp;
   if (ls == POLY345)
   {
      double TR2 = TR*TR;
      Amp = (float)( TR2*(10*TR - TR2*(15 - 6*TR)) );
   }
   else if (ls == RAMP)
   {
      Amp = (float)TR;
   }
   else if (ls == STEP)
   {
      Amp = 1;
   }
   else if (ls == HILLY)
   {
      const float cycles = 4;
      const float pi = 3.14159f;
      Amp = (float)( TR*exp(cos(cycles*pi*TR)-1) );
   } 
   else if (ls == RAMPHOLD) {
     if (TR >= 0.1) Amp = 1;
     else Amp = (float)(TR/0.1);
   }
   else // Shouldn't ever be reached
   {
     Amp = std::numeric_limits<float>::quiet_NaN();
   }

   return Amp;
}

enum loadShape atols(const char* str)
{
   if (!strcmp(str,"POLY345"))
      return POLY345;
   else if (!strcmp(str,"RAMP"))
      return RAMP;
   else if (!strcmp(str,"RAMPHOLD"))
      return RAMPHOLD;
   else if (!strcmp(str,"STEP"))
      return STEP;
   else if (!strcmp(str,"HILLY"))
      return HILLY;
   else
      cerr << "!!! Invalid load shape " << str << " --> using POLY345" << endl;
   
   return POLY345;
}
