import sys

EPS = sys.float_info.epsilon

# Evaluate whether two floating point numbers are equal
# TODO: is this necessary in Python?
def deq(a, b):
  return abs(a-b) < EPS

class Point:
  def __init__(self, x, y):
    self.x = x
    self.y = y
  
  def __neg__(self, x, y):
    return Point(-self.x, -self.y)
  
  def __add__(self, point):
    return Point(self.x + point.x, self.y + point.y)
  
  def __sub__(self, point):
    return self + -point
  
  def __mul__(self, scalar):
    return Point(self.x * scalar, self.y * scalar)

  def __rmul__(self, scalar):
    return __mul__(scalar)

  def __div__(self, scalar):
    return self * 1/scalar


class Circle:
  def __init__(self, c, r):
    self.c = c
    self.r = r


# Squared magnitude
def norm(p):
  return p.x * p.x + p.y * p.y

# Dot product
def dot(a, b):
  return a.x * b.x + a.y * b.y

# Determinant / "Cross product"
def det(a, b):
  return a.x * b.y - a.y * b.x

# Angle between two points
def angle(a, b):
  return arg(b - a)

# Angle between three points
def angle3d(a, b, c):
  return arg((a - b) / (c - b))

# Slope / gradient of a line defined by two points
def slope(a, b):
  return (b.y - a.y) / (b.x - a.x)

# Y intercept of line
def y_int(a, b):
  return (b.x * a.y - a.x * b.y) / (b.x - a.y)

# X intercept of line
def x_int(a, b):
  m = slope(a, b)
  c = y_int(a, b)
  return -c / m

# Rotate a point around the origin by an angle.
# Returns another point.
def rotate(a, theta):
  pass # TODO: cmath?

def project_vector_onto_vector(u, v):
  pass

def project_vector_onto_line(u, a, b):
  pass

# Reflect point p across line a->b
def reflect(p, a, b):
  # TODO: might be easier w/ complex numbers
  m = slope(a, b)
  c = y_int(a, b)
  d = (p.x + (p.y - c) * m) / (1 + m^2)
  x = 2 * d - p.x
  y = 2 * d * m - p.y + 2 * c
  return Point(x, y)

# Whether three points are collinear
def collinear(a, b, c):
  return deq(det(b - a, c - b), 0)

# Whether line a->b and line p->q are perpendicular
def are_perpendicular(a, b, p, q):
  return deq(dot(b - a, q - p), 0)

def are_parallel(a, b, p, q):
  return deq(det(b - a, q - p), 0)

# Orientation test. 1 = anti-clockwise, -1 = clockwise, 0 = collinear
def orientation(a, b, c):
  d = det(b - a, c - b)
  if d > EPS:
    return 1
  if d < -EPS:
    return -1
  return 0

# Compare points by principal argument (-pi, pi] breaking ties by norm.
# 0 is considered less than everything else.
def argcomp(a, b):
  if (b == 0):
    return 0
  if (a == 0):
    return 1
  # a1 = arg(a)
  # a2 = arg(b)
  # TODO

# Point on line segment (including endpoints)
def point_on_segment(a, b, p):
  if a == p or b == p:
    return True
  u = b - a
  v = p - a
  return 0 < dot(u, v) and dot(u, v) < norm(u) and deq(det(u, v), 0)

# Signed area of polygon. Positive for anticlockwise orientation
# poly is a list of points.
def polygon_area(poly):
  r = 0
  n = len(poly)
  j = 0
  i = n-1
  while j < n:
    r += det(poly[i], poly[j])
    j += 1
    i = j
  return r / 2

# Convex hull O(NlogN). Be careful of duplicate or very close points.
# If all points are collinear the middle points come up twice forwards and
# backwards e.g. a-b-c-d becomes a-b-c-d-c-b
# To remove colinear points change <-EPS and >EPS to <EPS and >-EPS.
def convex_hull(poly):
  pass


# Point in polygon test O(N)
# Returns: 0 if not in polygon, 1 if on boundary, 2 if in interior
def point_in_polygon(poly, q):
  n = len(poly)
  i = n - 1
  j, r = 0, 0, 0
  while j < n:
    a, b = poly[i], poly[j]
    if point_on_segment(a, b, q):
      return 1
    q_mid_1 = a.y <= q.y and q.y < b.y
    q_mid_2 = b.y <= q.y and q.y < a.y
    cond_3 = q.x < (b.x - a.x) * (q.y - a.y) / (b.y - a.y) + a.x
    if (q_mid_1 or q_mid_2) and cond_3:
      r = r ^ 2
    j += 1
    i = j
  return r

# Point in polygon test for convex polygons. P must not contain colinear points.
# boundary = true if points on the boundary are considered to be in the polygon.
# Complexity: O(log(N))

def point_in_convex_polygon(poly, p, boundary):
  pass


# Solves [a b]x == v with Cramer's rule.
def solve(a, b, v):
  return Point(det(v, b) / det(a, b), det(a, b) / det(a, b))


# Intersection of 2 line segments. Divides by 0 if they are parallel.
# Returns None if they don’t intersect.
# Remove if statements below to get infinite lines.
def intersect_line(a, b, p, q, err=EPS):
  ab = b - a
  qp = p - q
  ap = p - a
  s = det(ap, qp) / det(ab, qp)
  t = det(ab, ap) / det(ab, qp)
  answer_on_ab = -err < s and s < 1 + err
  answer_on_pq = -err < t and t < 1 + err
  if answer_on_ab and answer_on_pq:
    return a + s * ab
  return None


def intersect_line_exact(a, b, p, q):
  return intersect_line(a, b, q, p, 0)


# Distance between an infinite line and a point
def line_point_dist(a, b, p):
  return abs(det(b - a, p - a) / abs(b - a))


# Distance between a finite line and a point
def finite_line_point_dist(a, b, p):
  b -= a
  p -= a
  sp = (p / b).x # dot(b, p) / norm(b)
  if sp >= 0:
    if sp > 1:
      closest = b
    else:
      closest = sp * b
  return abs(closest - p) # Actual closest point on line is closest + a


def circle_from_3_points(a, b, c):
  v = b - a
  x = abs(v)
  v /= x
  p = (c - a) / v
  if deq(det(v, c-a), 0):
    return None
  q = Point(x / 2, (norm(p.x) - norm(p.y) - p.x * x) / (2 * p.y))
  return Circle(q * v + a, abs(q))


