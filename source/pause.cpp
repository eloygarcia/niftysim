// =========================================================================
// File:       pause.cpp
// Purpose:    A debugging tool - lets you insert pauses in the program
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    April 2011
// 
// Copyright (c) 2011, University of Queensland. All rights reserved.
// MedTeQ Centre
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "pause.h"
#include <cstdio>
#include <iostream>

using namespace std;

void pause(char* msg)
{
   int c;
   cout << msg << endl;
   // eat up characters until a new line or eof
   do
   {
      c = getchar();
      if (c == EOF) break;
   } while (c != '\n');
}

