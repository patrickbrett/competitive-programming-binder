# Binary exponentiation - compute a^b mod m. Complexity O(log(n))
def binary_exponentiation(a, b, m):
  res = 1 % m
  a = a % m
  while b > 0:
    if b & 1:
      res = res * a % m
      a = a * a % m
    b //= 2
  return res

# Extended Euclidean Algorithm. Finds x,y such that
# ax + by = gcd(a,b). Returns (gcd(a,b), x, y). Compexity: O(log(min(a,b)))
def gcd(a, b, x=0, y=0):
  if b == 0:
    y = 0
    x = -1 if a < 0 else 1
    return (abs(a), x, y)
  else:
    g = gcd(b, a % b, y, x)
    y -= a / b * x
    return (g, x, y)

# Multiplicative inverse of a % m, for a, m coprime.
# Complexity: O(log(a))
def inv(a, m):
  _, x, y = gcd(m, a)
  return ((y % m) + m) % m


# Chinese remainder algorithm. Solves x = a[i] mod m[i] for x mod lcm(m)
# for m[i] pairwise coprime. In general x= x0 + t * lcm(m) for all t.
# a and m are lists of ints.
def chinese_remainder_algorithm(a, m):
  n = len(a)
  u = a[0]
  v = m[0]
  for i in range(1, n):
    r, p, q = gcd(v, m[i])
    t = v
    if (a[i] - u) % r != 0:
      return -1 # No solution
    v = v / r * m[i]
    u = ((a[i] - u) / r * p * t + u) % v
  if u < 0:
    u += v
  return u


# Prime factors in O(log(n)) using precomputed fast_sieve(N >= n)
def fast_factors(n, fac):
  res = []
  while n > 1:
    f = fac[n]
    while n % f == 0:
      n /= f
    res.append(f)
  return res


# Prime factors in O(sqrt(n)) with no precomputation
def slow_factors(n):
  res = []
  i = 2
  while i * i <= n:
    if n % i == 0:
      res.append(i)
      while n % i == 0:
        n /= i
    i += 1
  if n > 1:
    res.append(n)
  return res


