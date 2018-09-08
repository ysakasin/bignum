#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>

#include <cstdio>
#include <iostream>

// immutable class
class Bignum {
private:
  enum class Sign {
    Positive,
    Negative,
  };

  Bignum(Sign sign, size_t length, uint32_t *digits)
      : sign(sign), length(length), digits(digits) {}

public:
  Bignum(const Bignum &src) : sign(src.sign), length(src.length) {
    digits = new uint32_t[src.length];
    std::memcpy(digits, src.digits, length);
  }

  Bignum(Bignum &&src) : sign(src.sign), length(src.length) {
    digits = src.digits;
    src.digits = nullptr;
  }

  Bignum(long long n) : sign(n < 0 ? Sign::Negative : Sign::Positive) {
    digits = new uint32_t[2];
    digits[0] = (uint32_t)abs(n);
    digits[1] = (uint32_t)(abs(n) >> 32);

    if (digits[1] == 0) {
      length = 1;
    } else {
      length = 2;
    }
  }

  Bignum(uint64_t n) : sign(Sign::Positive) {
    digits = new uint32_t[2];
    digits[0] = (uint32_t)n;
    digits[1] = (uint32_t)(n >> 32);

    if (digits[1] == 0) {
      length = 1;
    } else {
      length = 2;
    }
  }

  ~Bignum() {
    if (digits != nullptr) {
      delete[] digits;
    }
  }

  Bignum &operator=(Bignum &&src) {
    this->sign = src.sign;
    this->length = src.length;
    this->digits = src.digits;
    src.digits = nullptr;
    return *this;
  }

  Bignum add(const Bignum &rhs) const {
    if (this->sign == rhs.sign) {
      auto pair = this->sign_ignored_add(rhs);
      return std::move(Bignum(this->sign, pair.first, pair.second));
    }

    if (this->sign_ignored_less_than(rhs)) {
      auto pair = rhs.sign_ignored_sub(*this);
      return std::move(Bignum(rhs.sign, pair.first, pair.second));
    } else {
      auto pair = this->sign_ignored_sub(rhs);
      return std::move(Bignum(this->sign, pair.first, pair.second));
    }
  }

  Bignum sub(const Bignum &rhs) const {
    Sign negate_sign =
        (rhs.sign == Sign::Positive ? Sign::Negative : Sign::Positive);
    Bignum negate_rhs(negate_sign, rhs.length, rhs.digits);

    return std::move(this->add(negate_rhs));
  }

  Bignum mul(const Bignum &rhs) const {
    auto pair = sign_ignored_mul_simple(rhs);
    Sign new_sign = (this->sign == rhs.sign ? Sign::Positive : Sign::Negative);

    return std::move(Bignum(new_sign, pair.first, pair.second));
  }

  void print() const {
    for (long long i = length - 1; i >= 0; i--) {
      printf("%x", digits[i]);
    }
    printf("\n");
  }

private:
  static bool is_overflow(uint64_t x, uint64_t y) {
    return x <= std::numeric_limits<uint64_t>::max() - y;
  }

  static bool is_overflow(uint32_t x, uint32_t y) {
    return x <= std::numeric_limits<uint32_t>::max() - y;
  }

  uint32_t at(size_t pos) const {
    if (pos < length) {
      return digits[pos];
    } else {
      return 0;
    }
  }

  std::pair<size_t, uint32_t *> sign_ignored_add(const Bignum &rhs) const {
    size_t max_length = std::max(length, rhs.length);
    uint32_t *new_digits = new uint32_t[max_length + 1];

    uint32_t carry = 0; // 0 or 1
    for (size_t i = 0; i < max_length; i++) {
      uint64_t x = this->at(i);
      uint64_t y = rhs.at(i);

      uint64_t t = x + y + carry;
      new_digits[i] = (uint32_t)t;
      carry = (uint32_t)(t >> 32);
    }

    new_digits[max_length] = carry;
    size_t new_length = (carry == 1 ? max_length + 1 : max_length);
    return std::make_pair(new_length, new_digits);
  }

  bool sign_ignored_less_than(const Bignum &rhs) const {
    return sign_ignored_comp(rhs) == -1;
  }

  bool sign_ignored_greater_than(const Bignum &rhs) const {
    return sign_ignored_comp(rhs) == 1;
  }

  // this <  rhs  => -1
  // this == rhs  =>  0
  // this >  rhs  =>  1
  int sign_ignored_comp(const Bignum &rhs) const {
    if (this->length != rhs.length) {
      if (this->length < rhs.length) {
        return -1;
      } else {
        return 1;
      }
    }

    for (size_t i = this->length; i > 0; i--) {
      size_t index = i - 1;
      if (this->digits[index] != rhs.digits[index]) {
        if (this->digits[index] < rhs.digits[index]) {
          return -1;
        } else {
          return 1;
        }
      }
    }

    // this equals rhs
    return 0;
  }

  std::pair<size_t, uint32_t *> sign_ignored_sub(const Bignum &rhs) const {
    size_t max_length = std::max(length, rhs.length);
    uint32_t *new_digits = new uint32_t[max_length];
    size_t new_length = 1;

    uint32_t borrow = 0; // 0 or 1
    for (size_t i = 0; i < max_length; i++) {
      uint32_t x = this->at(i);
      uint32_t y = rhs.at(i);

      new_digits[i] = x - y - borrow;
      if (x < y) {
        borrow = 1;
      } else if (borrow == 1 && x == y) {
        borrow = 1;
      } else {
        borrow = 0;
      }

      if (new_digits[i] != 0) {
        new_length = i + 1;
      }
    }
    // assert(borrow == 0);

    return std::make_pair(new_length, new_digits);
  }

  static size_t sanitize_length(size_t length, const uint32_t *digits) {
    while (digits[length - 1] == 0 && length > 1) {
      length--;
    }
    return length;
  }

  std::pair<size_t, uint32_t *>
  sign_ignored_mul_simple(const Bignum &rhs) const {
    const size_t n = this->length;
    const size_t m = rhs.length;

    size_t new_length = this->length * rhs.length;
    uint32_t *new_digits = new uint32_t[new_length];
    std::fill_n(new_digits, new_length, 0);

    for (size_t j = 0; j < n; j++) {
      if (this->at(j) == 0) {
        new_digits[j + m] = 0;
        continue;
      }

      uint32_t k = 0;
      for (size_t i = 0; i < m; i++) {
        uint64_t x = this->at(j);
        uint64_t y = rhs.at(i);

        uint64_t t = x * y + new_digits[i + j] + k;
        new_digits[i + j] = (uint32_t)t;
        k = (uint32_t)(t >> 32);
      }
      new_digits[j + m] = k;
    }

    new_length = sanitize_length(new_length, new_digits);
    return std::make_pair(new_length, new_digits);
  }

  static void shift_digits(uint32_t *digits, size_t length, size_t shift) {
    if (shift <= 0) {
      return;
    }
    for (size_t i = 0, k = 0; i < length; i++) {
      uint32_t x = digits[i];
      digits[i] = (x << shift) | k;
      k = x >> (32 - shift);
    }
    // assert(k == 0);
  }

  // p.257
  std::pair<size_t, uint32_t *>
  sign_ignored_div_simple(const Bignum &rhs) const {
    size_t new_length;
    uint32_t *new_digits;

    if (sign_ignored_less_than(rhs)) {
      new_length = 1;
      new_digits = new uint32_t[1];
      new_digits[0] = 0;
      return std::make_pair(new_length, new_digits);
    }

    const size_t n = rhs.length;
    const size_t m = this->length - rhs.length;
    const uint64_t b = std::numeric_limits<uint32_t>::max();

    new_length = m + 1;
    new_digits = new uint32_t[new_length];

    // D1
    int shift = 0;
    uint32_t v = rhs.digits[rhs.length - 1];
    while (v < std::numeric_limits<uint32_t>::max() / 2) {
      v <<= 1;
      shift++;
    }

    uint32_t *rhs_digits_copy = new uint32_t[rhs.length + 1];
    std::memcpy(rhs_digits_copy, rhs.digits, rhs.length);
    rhs_digits_copy[rhs.length] = 0;

    uint32_t *this_digits_copy = new uint32_t[this->length + 1];
    std::memcpy(this_digits_copy, this->digits, this->length);
    this_digits_copy[this->length] = 0;

    shift_digits(rhs_digits_copy, rhs.length, shift);
    shift_digits(this_digits_copy, this->length + 1, shift);

    // D2, D7
    for (size_t j = m; j >= 0; j--) {

      // D3
      uint64_t tmp =
          (uint64_t)this_digits_copy[j + n] * b + this_digits_copy[j + n - 1];
      uint64_t q = tmp / rhs_digits_copy[n - 1];
      uint64_t r = tmp % rhs_digits_copy[n - 1];

      while (q >= b ||
             q * rhs_digits_copy[n - 2] > b * r + this_digits_copy[j + n - 2]) {
        q -= 1;
        r += rhs_digits_copy[n - 1];
        if (r < b) {
          break;
        }
      }

      // D4
      uint32_t carry = 0;
      uint32_t borrow = 0;
      for (size_t i = 0; i <= n; i++) {
        uint64_t t = rhs_digits_copy[i] * q + carry;
        carry = (uint32_t)(t >> 32);

        uint32_t x = this_digits_copy[j + i];
        uint32_t y = (uint32_t)t;
        this_digits_copy[j + i] = x - y - borrow;

        if (x < y) {
          borrow = 1;
        } else if (borrow == 1 && x == y) {
          borrow = 1;
        } else {
          borrow = 0;
        }
      }

      // D5
      if (borrow != 0) {
        // D6
        q -= 1;
        uint32_t carry = 0;
        for (size_t i = 0; i <= n; i++) {
          uint64_t x = this_digits_copy[j + i];
          uint64_t y = rhs_digits_copy[i];

          uint64_t t = x + y + carry;
          this_digits_copy[j + i] = (uint32_t)t;
          carry = (uint32_t)(t >> 32);
        }
      }

      new_digits[j] = q;
    }

    delete[] rhs_digits_copy;
    delete[] this_digits_copy;

    return std::make_pair(new_length, new_digits);
  }

private:
  Sign sign;
  size_t length;
  uint32_t *digits;
};
