#include "MurmurHash3.h"
#include "hasher.hpp"
#include <random>
#include <chrono>

double_hasher::double_hasher() {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
}

double_hasher::double_hasher(unsigned int _seed) {
  seed = _seed;
}

double double_hasher::operator()(unsigned int x){
  uint64_t input[1]= {x};
  uint64_t result[1]= {0};
  MurmurHash3_x64_128(input,sizeof(unsigned int),seed,result);
  return result[0] / static_cast<double>(std::numeric_limits<uint64_t>::max());
}

die_hasher::die_hasher() {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
}

die_hasher::die_hasher(unsigned int _seed) {
  seed = _seed;
}

uint64_t die_hasher::operator()(unsigned int x){
  uint64_t input[1]= {x};
  uint64_t result[1]= {0};
  MurmurHash3_x64_128(input,sizeof(unsigned int),seed,result);
  return __builtin_clz(result[0]);
}
