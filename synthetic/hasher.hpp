#ifndef _HASHER_H
#define _HASHER_H

class double_hasher {
  uint64_t seed;
public:
  double_hasher();
  double_hasher(unsigned _seed);
  double operator()(unsigned int x);
};

class die_hasher {
  uint64_t seed;
public:
  die_hasher();
  die_hasher(unsigned _seed);
  uint64_t operator()(unsigned int x);
};

#endif
