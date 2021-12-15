import svgpathtools, random, math, numpy, subprocess
import xml.etree.ElementTree as ET
from svgpathtools import CubicBezier, QuadraticBezier, Arc, Line

file = "input.svg"

paths, attributes, svg_attribs = svgpathtools.svg2paths2(file)

tree = ET.parse(file)
root = tree.getroot()
viewBox = root.attrib["viewBox"].split(" ")
viewBox = [float(x) for x in viewBox]
off_x = -viewBox[0]
off_y = -viewBox[1]

if "width" in root.attrib:
  width = float(root.attrib["width"])
else:
  width = viewBox[2]

if "height" in root.attrib:
  height = float(root.attrib["height"])
else:
  height = viewBox[3]

scale_x = width / viewBox[2]
scale_y = height / viewBox[3]

width = round(width)
height = round(height)

bitdepth = 3

scale = 1


def lerp(a, b, t):
  return round(a + (b - a) * t)


def round(x, n=0):
  return int(x * (10**n) + 0.5 * (1 if x > 0 else -1)) / (10**n)


def sim(a, b):
  return round(a[0]) == round(b[0]) and round(a[1]) == round(b[1])


def dist(p1, p2, euclid=False):
  # max manhattan distance to determine if two points can occupy the same space when shifted by 1
  if euclid:
    return math.sqrt(math.pow(p1[0] - p2[0], 2) + math.pow(p1[1] - p2[1], 2))

  return max(abs(p2[0] - p1[0]), abs(p2[1] - p1[1]))


def length(x, y):
  return math.sqrt(math.pow(x, 2) + math.pow(y, 2))


def jitter(x, y):
  angle = random.random() * (2 * math.pi)
  return (round(x + math.cos(angle)), round(y + math.sin(angle)))


def midpoint(a, b):
  return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)


def ang(p1, p2):
  return math.atan2(p2[1] - p1[1], p2[0] - p1[0])


def chebyshev(i, n):
  return (math.cos((n - i) / n * math.pi) + 1) / 2


def bezier(p0, p1, p2, p3, t):
  return (p0[0] * (1 - t)**3 + 3 * p1[0] * t * (1 - t)**2 + 3 * p2[0] * t**2 *
          (1 - t) + p3[0] * t**3, p0[1] * (1 - t)**3 + 3 * p1[1] * t *
          (1 - t)**2 + 3 * p2[1] * t**2 * (1 - t) + p3[1] * t**3)


def bezier_dt(p0, p1, p2, p3, t):
  c1x = p3[0] - (3 * p2[0]) + (3 * p1[0]) - p0[0]
  c2x = (3 * p2[0]) - (6 * p1[0]) + (3 * p0[0])
  c3x = (3 * p1[0]) - (3 * p0[0])

  c1y = p3[1] - (3 * p2[1]) + (3 * p1[1]) - p0[1]
  c2y = (3 * p2[1]) - (6 * p1[1]) + (3 * p0[1])
  c3y = (3 * p1[1]) - (3 * p0[1])

  return ((3 * c1x * t**2) + (2 * c2x * t) + c3x,
          (3 * c1y * t**2) + (2 * c2y * t) + c3y)


def cubic_from_bezier(p0, p1, p2, p3, spline_len):
  start = [bezier(p0, p1, p2, p3, 0)]
  end = [bezier(p0, p1, p2, p3, 1)]

  # equidistant points along the curve
  cp = []
  t = 0
  use_uniform = False
  while t <= 1:
    cp.append(bezier(p0, p1, p2, p3, t))
    dt = bezier_dt(p0, p1, p2, p3, t)
    l = length(dt[0], dt[1])
    if l == 0:
      use_uniform = True
      break
    t += 2 / l

  if use_uniform:
    cp.clear()
    cp.extend([
      bezier(p0, p1, p2, p3, i / int(spline_len))
      for i in range(int(spline_len))
    ])

  while cp and cp[0] == start[0]:
    del cp[0]

  while cp and cp[-1] == end[0]:
    del cp[-1]

  points = start + cp + end

  return points


def cubic_from_arc(center, radius, angle1, angle2, npoints=None):
  if npoints:
    cubic_points = []
    for i in range(npoints + 1):
      angle = (angle1 + i / npoints * angle2)
      point = (center[0] + math.cos(angle) * radius[0],
               center[1] + math.sin(angle) * radius[1])
      point = (round(point[0] * scale_x), round(point[1] * scale_y))
      cubic_points.append(point)
    return cubic_points

  npoints = 128
  points = cubic_from_arc(center, radius, angle1, angle2, npoints)

  return points


def solve_catmull_rom(points, start, x, y):
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

  a2 = -p0[1] + 3 * p1[1] - 3 * p2[1] + p3[1]
  b2 = 2 * p0[1] - 5 * p1[1] + 4 * p2[1] - p3[1]
  c2 = -p0[1] + p2[1]
  d2 = 2 * p1[1]

  t = numpy.roots([a1, b1, c1, d1 - x / alpha])
  t = [t.real for t in t if numpy.isreal(t)]
  t = [t for t in t if t >= 0 and t <= 1]

  t2 = numpy.roots([a2, b2, c2, d2 - y / alpha])
  t2 = [t.real for t in t2 if numpy.isreal(t)]
  t2 = [t for t in t2 if t >= 0 and t <= 1]

  return t + t2


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

  dx /= dl or 1
  dy /= dl or 1

  return (x, y, (dx, dy), (ddx, ddy))


def optimize(points):
  points = prune_invalid_points(points)
  if len(points) < 2: return []
  t = [f"{round(p[0] * scale, 12)} {round(p[1] * scale, 12)}" for p in points]
  t = "\n".join(t)
  with open("original_spline", "w") as f:
    f.write(t)

  subprocess.run("optimize")
  subprocess.run("round")
  t = open("quantized_spline", "r").readlines()
  points = [line.split(" ") for line in t]
  points = [(float(p[0]), float(p[1])) for p in points]
  points = prune_invalid_points(points)
  return points


def prune_invalid_points(points):
  if len(points) < 2: return []
  new_points = []

  for point in points[:-1]:
    if new_points and sim(new_points[-1], point):
      continue
    new_points.append(point)

  new_points.append((points[-1][0], points[-1][1]))
  if len(new_points) < 2: return []

  while len(new_points) > 2 and sim(new_points[0], new_points[1]):
    del new_points[1]

  while len(new_points) > 2 and sim(new_points[-2], new_points[-1]):
    del new_points[-2]

  if len(new_points) < 2: return []
  if len(new_points) == 2 and sim(new_points[0], new_points[1]):
    return []
  return new_points


def create_spline(points, color, do_optimize=True):
  if len(points) > 2 and do_optimize:
    points = optimize(points)

  points = prune_invalid_points(points)
  if len(points) < 2:
    return ""
  if len(points) == 2 and sim(points[0], points[1]):
    return ""

  new_color = [round(c / 255 * bitdepth * 2 - bitdepth, 2) for c in color]

  spline = [
    "Spline",
    f"{new_color[0]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{new_color[1]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{new_color[2]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
  ]

  spline.extend([f"{round(p[0], 2)} {round(p[1], 2)}" for p in points])
  spline.append("EndSpline")

  return spline


splines = []

for path, attribute in zip(paths, attributes):
  color = "#000000"
  # TODO: classes

  if "fill" in attribute:
    color = attribute["fill"]
    if color.startswith("url"):
      continue

  if "stroke" in attribute:
    color = attribute["stroke"]

  if color == "#FFFFFF": continue

  color = color.lstrip("#")
  if len(color) == 6:
    color = [int(color[i:i + 2], 16) for i in (0, 2, 4)]
  elif len(color) == 3:
    color = [int(color[i:i + 1], 16) * 8 for i in (0, 1, 2)]
  else:
    color = [0, 0, 0]

  if color[0] > 200 and color[1] > 200 and color[2] > 200: continue

  if "fill-opacity" in attribute:
    if float(attribute["fill-opacity"]) == 0: continue

  current_spline = []

  isplines = []
  for segment in path:
    if type(segment) == Line:
      start = (segment.start.real + off_x, segment.start.imag + off_y)
      end = (segment.end.real + off_x, segment.end.imag + off_y)
      start = (start[0] * scale_x, start[1] * scale_y)
      end = (end[0] * scale_x, end[1] * scale_y)
      if sim(start, end):  # this is a point
        continue

      d = [end[0] - start[0], end[1] - start[1]]
      l = length(d[0], d[1])
      d[0] /= l
      d[1] /= l

      p1 = (start[0] + d[0], start[1] + d[1])
      p2 = (end[0] - d[0], end[1] - d[1])
      cubic_points = [start, p1, p2, end]

      if current_spline:
        if dist(current_spline[0], cubic_points[0]) < 1 / scale:
          current_spline.reverse()

        if dist(current_spline[-1], cubic_points[-1]) < 1 / scale:
          cubic_points.reverse()

        if dist(current_spline[-1], cubic_points[0]) > 1 / scale:
          isplines.append(create_spline(current_spline, color))
          current_spline.clear()

      if current_spline:
        cubic_points.pop(0)

      current_spline.extend(cubic_points)

    elif type(segment) in [CubicBezier, QuadraticBezier]:
      p0 = (segment.start.real + off_x, segment.start.imag + off_y)
      p1 = (segment.control1.real + off_x, segment.control1.imag + off_y)
      if type(segment) == QuadraticBezier:
        p2 = (segment.control.real + off_x, segment.control.imag + off_y)
      else:
        p2 = (segment.control2.real + off_x, segment.control2.imag + off_y)
      p3 = (segment.end.real + off_x, segment.end.imag + off_y)

      cubic_points = cubic_from_bezier(p0, p1, p2, p3, segment.length())
      cubic_points = [(scale_x * p[0], scale_y * p[1]) for p in cubic_points]

      if current_spline:
        if dist(current_spline[0], cubic_points[0]) < 1 / scale:
          current_spline.reverse()

        if dist(current_spline[-1], cubic_points[-1]) < 1 / scale:
          cubic_points.reverse()

        if dist(current_spline[-1], cubic_points[0]) > 1 / scale:
          isplines.append(create_spline(current_spline, color))
          current_spline.clear()

      if current_spline:
        cubic_points.pop(0)

      current_spline.extend(cubic_points)

    elif type(segment) == Arc:
      print("arc")
      center = (segment.center.real + off_x, segment.center.imag + off_y)
      radius = (segment.radius.real, segment.radius.imag)

  if current_spline:
    isplines.append(create_spline(current_spline, color))

  splines.extend(isplines)

lines = [
  f"Width {int(round(width * scale))}",
  f"Height {int(round(height * scale))}",
  f"Bitdepth {bitdepth}",
]

lines.extend(["\n".join(spline) for spline in splines])
lines.append(f"- Set + {2**bitdepth - 1}")

text = "\n".join(lines)
open("spline", "w+").write(text)
