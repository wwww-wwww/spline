import numpy, math

op = open("original_spline", "r").readlines()
op = [line.split(" ") for line in op]
op = [(float(p[0]), float(p[1])) for p in op]
op = [(p[0], p[1]) for p in op]


def dist(p1, p2):
  return math.sqrt(math.pow(p1[0] - p2[0], 2) + math.pow(p1[1] - p2[1], 2))


def length(x, y):
  return math.sqrt(math.pow(x, 2) + math.pow(y, 2))


def catmull_rom(points, start, t):
  alpha = 0.5

  if start == 0:
    p0 = points[0]
  else:
    p0 = points[start - 1]

  p1 = points[start + 0]
  p2 = points[start + 1]

  if start == len(points) - 2:
    p3 = p2
  else:
    p3 = points[start + 2]

  a1 = -p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]
  a2 = -p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]
  b1 = 2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]
  b2 = 2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]
  c1 = -p0[0] + p2[0]
  c2 = -p0[1] + p2[1]
  d1 = 2 * p1[0]
  d2 = 2 * p1[1]

  x = (a1 * (t**3) + b1 * (t**2) + c1 * t + d1) * alpha
  y = (a2 * (t**3) + b2 * (t**2) + c2 * t + d2) * alpha

  dx = (3 * (t**2) * a1 + 2 * t * b1 + c1) * alpha
  dy = (3 * (t**2) * a2 + 2 * t * b2 + c2) * alpha

  ddx = (6 * t * a1 + 2 * b1) * alpha
  ddy = (6 * t * a2 + 2 * b2) * alpha

  dl = length(dx, dy)

  dx /= dl
  dy /= dl

  return (x, y, (dx, dy), (ddx, ddy))


def catmull_rom_x(points, start, x):
  alpha = 0.5

  if start == 0:
    p0 = points[0]
  else:
    p0 = points[start - 1]

  p1 = points[start + 0]
  p2 = points[start + 1]

  if start == len(points) - 2:
    p3 = p2
  else:
    p3 = points[start + 2]

  a1 = -p0[0] + 3 * p1[0] - 3 * p2[0] + p3[0]
  b1 = 2 * p0[0] - 5 * p1[0] + 4 * p2[0] - p3[0]
  c1 = -p0[0] + p2[0]
  d1 = 2 * p1[0]

  t = numpy.roots([a1, b1, c1, d1 - x / alpha])
  t = [t.real for t in t if numpy.isreal(t)]
  t = [t for t in t if t >= 0 and t <= 1]

  return t

def catmull_rom_y(points, start, x):
  alpha = 0.5

  if start == 0:
    p0 = points[0]
  else:
    p0 = points[start - 1]

  p1 = points[start + 0]
  p2 = points[start + 1]

  if start == len(points) - 2:
    p3 = p2
  else:
    p3 = points[start + 2]

  a1 = -p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]
  b1 = 2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]
  c1 = -p0[1] + p2[1]
  d1 = 2 * p1[1]

  t = numpy.roots([a1, b1, c1, d1 - x / alpha])
  t = [t.real for t in t if numpy.isreal(t)]
  t = [t for t in t if t >= 0 and t <= 1]

  return t

def find_on_spline(points, point):
  t = []
  for i in range(len(points) - 1):
    ts = catmull_rom_x(points, i, point[0]) + catmull_rom_y(points, i, point[1])
    for t1 in ts:
      print(t1)
      t.append((i, catmull_rom(points, i, t1)))

  if len(t) > 0:
    return min(t, key=lambda y: dist(y[1], point))

  return None

del op[3]

print(find_on_spline(op, (1163, 444)))

print("\n".join([f"{p[0]} {p[1]}" for p in op]))