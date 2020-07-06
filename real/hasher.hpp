#ifndef _HASHER_H
#define _HASHER_H

class double_hasher {
  uint64_t seed;
public:
  double_hasher();
  double_hasher(unsigned _seed);
  double operator()(unsigned int x);
  double operator()(std::string x);
};

class die_hasher {
  uint64_t seed;
public:
  die_hasher();
  die_hasher(unsigned _seed);
  uint64_t operator()(unsigned int x);
  uint64_t operator()(std::string x);
};

#endif
