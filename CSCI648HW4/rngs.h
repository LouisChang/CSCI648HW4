//
//  rngs.h
//  CSCI648HW4
//
//  Created by Sidi Chang on 3/1/15.
//  Copyright (c) 2015 Sidi Chang. All rights reserved.
//

#ifndef __CSCI648HW4__rngs__
#define __CSCI648HW4__rngs__
#if !defined( _RNGS_ )
#define _RNGS_

double Random(void);
void   GetSeed(long *x);
void   PutSeed(long x);
void   PlantSeeds(long x);
void   SelectStream(int s);
void   TestRandom(void);

#endif

#endif /* defined(__CSCI648HW4__rngs__) */
