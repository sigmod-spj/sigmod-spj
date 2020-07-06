#ifndef _HASHER_H
#define _HASHER_H

/**
 * double hasher maps a int/string to a double uniformly distributed in [0,1]
 */
class double_hasher {
  uint64_t seed;
public:
  double_hasher();
  double_hasher(unsigned _seed);
  double operator()(unsigned int x);
  double operator()(std::string x);
};

/**
 * die_hashser maps a int/string to a geometrically distributed int.
 * With probability 0.5 it is 0, with probability 0.25 it is 1 etc..
 */
class die_hasher {
  uint64_t seed;
public:
  die_hasher();
  die_hasher(unsigned _seed);
  uint64_t operator()(unsigned int x);
  uint64_t operator()(std::string x);
};

#endif
