// =========================================================================
// File:       tledTimer.h
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


#ifndef tledTimer_H
#define tledTimer_H

#include <iostream>
#ifndef _USE_BOOST_
#include <ctime>
#else
#include <boost/date_time/posix_time/posix_time.hpp>
#endif

class tledTimer
{
public:
   tledTimer();
   ~tledTimer() {;}
   
   void StartTimer(void);
   void StopTimer(void);

   /** Time in milliseconds */
   double GetDuration(void) const;

private:
#ifndef _USE_BOOST_
  std::clock_t start;
  std::clock_t finish;
#else
  boost::posix_time::ptime start;
  boost::posix_time::ptime finish;
#endif
};


#endif // tledTimer_H
