/*
 * Computation of the n^th decimal digit of pi with constant memory using
 * only 32-bit integer arithmetic.
 * By David McWilliams, 2021.
 *
 * This program is optimized for mapping to Factorio combinators.
 * The first integer overflow occurs when 3*N > sqrt(INT32_MAX).
 * Only the first 17400 digits are guaranteed to be accurate.
 *
 * Based on pi1.c by Fabrice Bellard, 1997.
 * https://bellard.org/pi/
 *
 * Uses the hypergeometric series by Bill Gosper, 1974.
 * pi = sum( (50*n-6)/(binomial(3*n,n)*2^n), n=0..infinity )
 * https://arxiv.org/abs/math/0110238
 *
 * Uses the constant memory algorithm by Simon Plouffe, 1996.
 * https://arxiv.org/abs/0912.0303
 *
 * See also the fastest known n^th decimal digit program by Xavier Gourdon, 2003.
 * http://numbers.computation.free.fr/Constants/Algorithms/pidec.cpp
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Return (a^b) mod m */
int pow_mod(int a, int b, int m) {
  int result = 1;
  while (b > 0) {
    if (b & 1) {
      result = result * a % m;
    }
    a = a * a % m;
    b >>= 1;
  }
  return result;
}

/* Solve for x: (a * x) % m == 1
 * https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Modular_integers
 */
int inv_mod(int a, int m) {
  if (a < 0) {
    a += m;
  }
  int b = m;
  int x = 1;
  int y = 0;
  int q;
  // 11 iterations is enough if m is a prime power less than sqrt(INT32_MAX)
  // Longest test case: a=17711, m=28657
  for (int i = 0; i < 11; i++) {
    q = (a == 0) ? 0 : b / a;
    b -= a * q;
    y -= x * q;
    q = (b == 0) ? 0 : a / b;
    a -= b * q;
    x -= y * q;
  }
  if (b == 0) {
    return x;
  } else {
    return y + m;
  }
}

/* Check if n is prime, 2 <= n < sqrt(INT32_MAX) */
int is_prime(int n) {
  static const int SMALL_PRIMES[47] = {
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
     31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211
  };
  int count = 0;
  for (int i = 0; i < 47; i++) {
    // Prime test fails if it is divisible by another prime
    if (n % SMALL_PRIMES[i] != 0) {
      count++;
    }
    // Exception: a prime is allowed to be divisible by itself
    if (n == SMALL_PRIMES[i]) {
      count++;
    }
  }
  // If we passed all 47 tests, it is prime
  return (count == 47);
}

/* Return the prime number immediately after n */
int next_prime(int n) {
  do {
    n++;
  } while (!is_prime(n));
  return n;
}

int prime_power[10];
int prime_power_count;

/* Remove prime factors from n and count how many were removed */
int factor_count(int *n) {
  for (int i = prime_power_count - 1; i >= 0; i--) {
    if (*n % prime_power[i] == 0) {
      *n /= prime_power[i];
      return i;
    }
  }
}

/* Calculate sum = (sum + n/d) and store the decimal part in fixed-point format
 * with 18 decimal places across two 32-bit integers.
 *
 * This is equivalent to the floating point one-liner:
 * sum = fmod(sum + (double)n / (double)d, 1.0);
 *
 * d must be less than sqrt(INT32_MAX).
 * An exception is made for powers of 2, where d may be up to 8388608.
 */
void fixed_point_sum(int n, int d, int *hi, int *lo) {
  // Avoid overflow for large powers of 2
  int r = 0;
  if (d > 60000) {
    d = d / 256;
    r = n % 256 * 125;
    n = n / 256;
  }

  // Digits 1 to 9
  int a = n * 32000 + r;
  *hi += a / d * 31250;
  int b = a % d * 31250;
  *hi += b / d;

  // Digits 10 to 18
  int c = b % d * 32000;
  *lo += c / d * 31250;
  *lo += c % d * 31250 / d;

  // Carry
  if (*lo > 1000000000) {
    *hi += 1;
  }

  // Discard overflow digits
  *hi = *hi % 1000000000;
  *lo = *lo % 1000000000;
}

/* Return 9 digits of pi */
int pi_digits(int start_digit) {
  int sum = 0;
  int sum_low = 0;
  // N = (start_digit + 19) / log10(13.5)
  // log10(13.5) is approximately equal to 269/238
  int N = (start_digit + 19) * 238 / 269;

  // Factor the Gosper series into fractions over prime powers
  for (int prime = 2; prime <= (3 * N); prime = next_prime(prime)) {
    // Compute the first few prime powers
    // Only 10 powers are needed if start_digit < 17500
    // Only powers up to 50000 are needed if start_digit < 17500
    static const int ROOT_50K[10] = {50000, 50000, 223, 36, 14, 8, 6, 4, 3, 3};
    prime_power_count = 0;
    for (int i = 0; i < 10; i++) {
      if (prime <= ROOT_50K[i]) {
        prime_power[i] = (int)pow(prime, i);
        prime_power_count++;
      }
    }

    // For small primes, use a prime power with exponent greater than 1
    int exponent = -1;
    for (int i = 0; i < prime_power_count; i++) {
      if (prime_power[i] <= (3 * N)) {
        exponent++;
      }
    }
    int m = (int)pow(prime, exponent);

    if (prime == 2) {
      // Add the 2^N term in the denominator.
      exponent += N - 1;
      // We have some more powers of 2 in the 10^start_digit decimal shift
      // in the numerator. Use them to cancel out the 2^N term.
      m = (int)pow(prime, exponent - start_digit);
      // Since start_digit grows faster than N, eventually we will
      // cancel the entire exponent and m will become 0.
      if (m == 0) {
        continue;
      }
    }

    // Multiply by 10^start_digit to move the target digit
    // to the most significant decimal place.
    int decimal = 10;
    if (prime == 2) {
      // We already used those powers of 2
      decimal = 5;
    }
    int decimal_shift = pow_mod(decimal, start_digit, m);

    // Main loop
    int subtotal = 0;
    int numerator = 1;
    int denominator = 1;
    //printf("N=%d\n", N);
    for (int k = 1; k <= N; k++) {
      // Terms for the numerator
      int t1 = 2 * k;
      int t2 = 2 * k - 1;
      exponent += factor_count(&t1);
      exponent += factor_count(&t2);
      int terms = (t1 % m) * (t2 % m) % m;
      numerator = numerator * terms % m;

      // Terms for the denominator
      int t3 = 6 * k - 4;
      int t4 = 9 * k - 3;
      exponent -= factor_count(&t3);
      exponent -= factor_count(&t4);
      terms = (t3 % m) * (t4 % m) % m;
      denominator = denominator * terms % m;
      //printf("m=%d k=%d numerator=%d denominator=%d\n", m, k, numerator, denominator);

      // Multiply all parts together
      int t = (50 * k - 6) % m;
      t = t * (int)pow(prime, exponent) % m;
      t = t * numerator % m;
      t = t * inv_mod(denominator, m);

      subtotal = (subtotal + t) % m;
    }
    subtotal = subtotal * decimal_shift % m;

    // We have a fraction over a prime power, add it to the final sum
    //printf("prime %d = %d/%d\n", prime, subtotal, m);
    fixed_point_sum(subtotal, m, &sum, &sum_low);
  }
  return sum;
}

int main(int argc, char *argv[]) {
  // Display help message
  if (argc < 2) {
    printf("This program computes digits of pi.\n");
    printf("Usage: pifactory <START_DIGIT> [END_DIGIT]\n");
    return 0;
  }

  // Read command line arguments
  int start = atoi(argv[1]);
  int end = start;
  if (argc >= 3) {
    end = atoi(argv[2]);
  }

  // Turn off output buffer
  setvbuf(stdout, NULL, _IONBF, 0);

  // Print digits of pi
  if (start == 0) {
    printf("3.");
    start++;
  }
  for (int i = start - 1; i < end; i += 9) {
    printf("%09d", pi_digits(i));
  }
  printf("\n");

  return 0;
}
