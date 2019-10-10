# Search for 6-term Machin-like formulas for pi
# Algorithm by Jörg Arndt https://www.jjj.de/arctan/

from sympy import factorint
from sympy import prime
from numpy import arctan
import cypari2
pari = cypari2.Pari()

# I am only interested in arctangents between 1/100 and 1/463
MIN_TERM = 100
MAX_TERM = 463

# Generate and factor a^2+1
factors = {}
for a in range(1, MAX_TERM + 1):
    factors[a] = factorint(a**2 + 1)

# Generate primes of the form 4k+1
primes = []
i = 0
while True:
    i += 1
    p = prime(i)
    if p > MAX_TERM:
        break
    if p % 4 == 1:
        primes.append(p)

# Generate all 5-term combinations of primes
for x1 in range(0, len(primes)):
    for x2 in range(x1 + 1, len(primes)):
        for x3 in range(x2 + 1, len(primes)):
            for x4 in range(x3 + 1, len(primes)):
                for x5 in range(x4 + 1, len(primes)):
                    primes_subset = {2: True, primes[x1]: True, primes[x2]: True,
                                     primes[x3]: True, primes[x4]: True, primes[x5]: True}

                    # Filter squares by selected primes
                    squares = []
                    for a in range(MIN_TERM, MAX_TERM + 1):
                        fits = True
                        for i in factors[a]:
                            if i not in primes_subset:
                                fits = False
                                break
                        if fits:
                            squares.append(a)

                    # Restart if we don't have enough squares
                    if len(squares) < len(primes_subset):
                        continue

                    # Remove 2
                    del primes_subset[2]

                    # Generate all 6-term combinations of squares
                    for y1 in range(0, len(squares)):
                        for y2 in range(y1 + 1, len(squares)):
                            for y3 in range(y2 + 1, len(squares)):
                                for y4 in range(y3 + 1, len(squares)):
                                    for y5 in range(y4 + 1, len(squares)):
                                        for y6 in range(y5 + 1, len(squares)):
                                            squares_subset = [squares[y1], squares[y2], squares[y3], squares[y4],
                                                              squares[y5], squares[y6]]

                                            # Build matrix of prime factors
                                            m = "["
                                            for p in primes_subset:
                                                if m != "[":
                                                    m += ";"
                                                for j in range(len(squares_subset)):
                                                    if j > 0:
                                                        m += ","
                                                    a = squares_subset[j]
                                                    n = 0
                                                    if p in factors[a]:
                                                        n = factors[a][p]
                                                    # Choose sign
                                                    if a % p < p / 2:
                                                        n *= -1
                                                    m += str(n)
                                            m += "]"

                                            # Solve matrix
                                            nullspace = pari.Mat(m).matkerint()
                                            coefficients = nullspace[0]

                                            # Test formula
                                            total = 0
                                            for i in range(len(squares_subset)):
                                                total += coefficients[i]*arctan(1/squares_subset[i]) * 4
                                            if abs(total) < 0.01:
                                                continue

                                            # Print formula
                                            print(f"{total} / 4 = ", end="")
                                            for i in range(len(squares_subset)):
                                                if i > 0:
                                                    print(" + ", end="")
                                                print(f"{coefficients[i]}*arctan(1/{squares_subset[i]})", end="")
                                            print()
