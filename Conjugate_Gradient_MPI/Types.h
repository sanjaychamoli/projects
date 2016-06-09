/*
 * Types and defines
 */

#ifndef TYPES_H_
#define TYPES_H_

#include <cmath>
#include <vector>

#define K   2*M_PI

#define K_2 4*M_PI*M_PI

#define ROOT 0

/* This typedef makes it possible to switch between float and double accuracy */
typedef double real_t;

/* Short hand for STL vector */
template <typename T>
using vec_t = std::vector <T>;

#endif /* TYPES_H_ */
