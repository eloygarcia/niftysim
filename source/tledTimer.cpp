// =========================================================================
// File:       tledTimer.cpp
// Purpose:    Class for timing computations
// Package:    NiftySim: Nonlinear finite element program
// Author:     Zeike Taylor
// Language:   C++
// Created:    August 2008
// 
// Copyright (c) 2010, University College London. All rights reserved.
// Centre for Medical Image Computing (CMIC)
// See the LICENSE.txt file in the root folder
// 
// ztaylor@itee.uq.edu.au
// =========================================================================

#include "tledTimer.h"

tledTimer::tledTimer()
{
#ifndef _USE_BOOST_
   start = std::clock();
   finish = std::clock();
#else
   finish = start = boost::posix_time::microsec_clock::universal_time();
#endif
}

void tledTimer::StartTimer(void)
{
#ifndef _USE_BOOST_
   start = std::clock();
#else
   start = boost::posix_time::microsec_clock::universal_time();
#endif
}

void tledTimer::StopTimer(void)
{
#ifndef _USE_BOOST_
   finish = std::clock();
#else
   finish = boost::posix_time::microsec_clock::universal_time();
#endif
}

double tledTimer::GetDuration(void) const
{
#ifndef _USE_BOOST_
   return (finish - start)*1e3/CLOCKS_PER_SEC;
#else
   return double((finish - start).total_milliseconds());
#endif
}
