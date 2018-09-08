#include "./bignum.hpp"

int main() {
  Bignum a((uint64_t)0xaaaaaaaaaaaaaaaa);

  Bignum x = a.mul(a);
  x.print();
  // 151236607520417094856213830793044048100
  // 0x71c71c71c71c71c638e38e38e38e38e4
  return 0;
}
