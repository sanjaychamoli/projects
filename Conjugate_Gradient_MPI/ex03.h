/*
 * ex03.h
 *
 *  Created on: Dec 27, 2014
 *      Author: ken
 */

#ifndef EX03_H_
#define EX03_H_

#include "Types.h"

inline real_t f ( real_t x, real_t y)
{
  return K_2 * std::sin( K*x ) * std::sinh ( K*y );
}

inline real_t U ( real_t x, real_t y)
{
  return std::sin( K*x ) * std::sinh ( K*y );
}


#endif /* EX03_H_ */
