// This is really -*- C++ -*-

#ifndef _coord_h
#define _coord_h

#include <cmath>

typedef std::pair<double, double> coordPair;

/** Tolerance for coordinate comparison . . .

    This could be a bit smaller and still not run into the IEEE 754
    standard for double truncation but this should be good enough.
*/
constexpr static double ctol = 1.0e-14;

size_t coordHash(coordPair P);

//! Hash object for galaxy model coordinates
namespace std
{
  template<>
  struct hash<coordPair> 
  {
    //! The hash function
    size_t operator()(coordPair s) const 
    { 
      return coordHash(s); 
    }
  };
}

//! Comparison object for coordinates
struct eqcoord
{
  //! Operator for comparing coordinate pairs
  bool operator()(coordPair c1, coordPair c2) const
  {
    return ( std::fabs(c1.first  - c2.first ) < ctol &&
	     std::fabs(c1.second - c2.second) < ctol   );
  }
};


#endif
