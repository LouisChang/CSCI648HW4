//
//  rvgs.h
//  CSCI648HW4
//
//  Created by Sidi Chang on 3/1/15.
//  Copyright (c) 2015 Sidi Chang. All rights reserved.
//

#ifndef __CSCI648HW4__rvgs__
#define __CSCI648HW4__rvgs__

#include <stdio.h>
#if !defined( _RVGS_ )
#define _RVGS_

long Bernoulli(double p);
long Binomial(long n, double p);
long Equilikely(long a, long b);
long Geometric(double p);
long Pascal(long n, double p);
long Poisson(double m);

double Uniform(double a, double b);
double Exponential(double m);
double Erlang(long n, double b);
double Normal(double m, double s);
double Lognormal(double a, double b);
double Chisquare(long n);
double Student(long n);

#endif
#endif /* defined(__CSCI648HW4__rvgs__) */
