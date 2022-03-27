import argparse
import math
import random
import subprocess
import svgpathtools
import sys
import xml.etree.ElementTree
from svgpathtools import CubicBezier, QuadraticBezier, Arc, Line
import optimize

# used for noise
gradients = [
  [0, 1],
  [1, 0],
  [0, -1],
  [-1, 0],
  [1, 1],
  [1, -1],
  [-1, 1],
  [-1, -1],
]


def bezier1d(p0, p1, p2, p3, t):
  return p0 * (1 - t)**3 + 3 * p1 * t * (1 - t)**2 + 3 * p2 * t**2 * (
    1 - t) + p3 * t**3


def bezier2d(p0, p1, p2, p3, t):
  return (
    bezier1d(p0[0], p1[0], p2[0], p3[0], t),
    bezier1d(p0[1], p1[1], p2[1], p3[1], t),
  )


def bezier1d_dt(p0, p1, p2, p3, t):
  return 3 * (1 - t)**2 * (p1 - p0) + 6 * t * (1 - t) * (
    p2 - p1) + 3 * t**2 * (p3 - p2)


def bezier2d_dt(p0, p1, p2, p3, t):
  return (
    bezier1d_dt(p0[0], p1[0], p2[0], p3[0], t),
    bezier1d_dt(p0[1], p1[1], p2[1], p3[1], t),
  )


def sample_bezier(p0, p1, p2, p3):
  points = []

  t = 0
  while t < 1:
    if t > 1: break
    points.append(bezier2d(p0, p1, p2, p3, t))
    dt = bezier2d_dt(p0, p1, p2, p3, t)
    l = length(dt[0], dt[1])
    l = l or 10
    t += 1 / l  # 1 is the absolute distance per step

  if t > 1:
    points.append(bezier2d(p0, p1, p2, p3, 1))

  return points


def diff(p1, p2):
  return (p1[0] - p2[0], p1[1] - p2[1])


def norm(p1):
  length = math.sqrt(math.pow(p1[0], 2) + math.pow(p1[1], 2)) or 1
  return (p1[0] / length, p1[1] / length)


def length(x, y):
  return math.sqrt(math.pow(x, 2) + math.pow(y, 2))


def dist(p1, p2):
  return length(p1[0] - p2[0], p1[1] - p2[1])


def sim(p1, p2):
  return round(p1[0]) == round(p2[0]) and round(p1[1]) == round(p2[1])


def sample_line(p0, p1):
  points = []
  d = diff(p1, p0)
  n = norm(d)
  l = dist(p0, p1)
  t = 0
  while t < l:
    points.append((p0[0] + t * n[0], p0[1] + t * n[1]))
    t += 1
  return points or [p0, p1]

def create_spline(points, error, scale):
  if len(points) > 2:
    old_len = len(points)
    points = optimize.optimize(points, error * (scale ** 2))
    print(f"{old_len} -> {len(points)}")

  if len(points) < 2:
    return ""

  if len(points) == 2 and dist(points[0], points[1]) < 1:
    # optionally add noise or ignore this case
    return ""

  points = [(round(p[0]), round(p[1])) for p in points]

  for i in range(len(points) - 1):
    if sim(points[i], points[i + 1]):
      direction = random.choice(gradients)
      points[i +
             1] = [points[i][0] + direction[0], points[i][1] + direction[1]]

  #new_color = [round(c / 255 * bitdepth * 2 - bitdepth, 2) for c in color]
  new_color = [-7, -7, -7]
  spline = [
    "Spline",
    f"{new_color[0]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{new_color[1]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{new_color[2]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    "0.7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
  ]

  spline.extend([f"{round(p[0], 2)} {round(p[1], 2)}" for p in points])
  spline.append("EndSpline")

  return spline


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="Input SVG file")
  parser.add_argument("output", help="Output jxl tree")
  parser.add_argument("-s", "--scale", type=float, default=1)
  parser.add_argument("-u", "--upsample", type=int, default=1)
  parser.add_argument("-e",
                      "--error",
                      help="Allowed error, multiplied by scale^2",
                      type=float,
                      default=1)

  args = parser.parse_args()

  paths, attributes, svg_attribs = svgpathtools.svg2paths2(args.input)

  tree = xml.etree.ElementTree.parse(args.input)
  root = tree.getroot()

  viewBox = [float(x) for x in root.attrib["viewBox"].split(" ")]

  off_x = -viewBox[0]
  off_y = -viewBox[1]

  try:
    width = float(root.attrib["width"])
  except:
    width = viewBox[2]

  try:
    height = float(root.attrib["height"])
  except:
    height = viewBox[3]

  width *= args.scale
  height *= args.scale

  header = f"""Bitdepth 4
  Upsample {args.upsample}
  Width {round(width)}
  Height {round(height)}
  RCT 0
  """

  scale_x = width / viewBox[2]
  scale_y = height / viewBox[3]

  with open(args.output, "w+") as output:
    output.write(header)
    splines = []

    for path, attribute in zip(paths, attributes):
      current_spline = []
      for segment in path:
        if type(segment) in [CubicBezier, QuadraticBezier]:
          p0 = (segment.start.real + off_x, segment.start.imag + off_y)
          if type(segment) == QuadraticBezier:
            p1 = (segment.control.real + off_x, segment.control.imag + off_y)
            p2 = (segment.control.real + off_x, segment.control.imag + off_y)
          else:
            p1 = (segment.control1.real + off_x, segment.control1.imag + off_y)
            p2 = (segment.control2.real + off_x, segment.control2.imag + off_y)
          p3 = (segment.end.real + off_x, segment.end.imag + off_y)

          p0 = (p0[0] * scale_x, p0[1] * scale_y)
          p1 = (p1[0] * scale_x, p1[1] * scale_y)
          p2 = (p2[0] * scale_x, p2[1] * scale_y)
          p3 = (p3[0] * scale_x, p3[1] * scale_y)

          cubic_points = sample_bezier(p0, p1, p2, p3)

          if current_spline:
            if dist(current_spline[0], cubic_points[0]) < 1:
              current_spline.reverse()

            if dist(current_spline[-1], cubic_points[-1]) < 1:
              cubic_points.reverse()

            if dist(current_spline[-1], cubic_points[0]) > 1:
              spline = create_spline(current_spline, args.error, args.scale)
              if spline:
                spline_str = "\n".join(spline)
                if spline_str not in splines:
                  splines.append(spline_str)
                  output.write(spline_str)
                  output.write("\n")
                  output.flush()
              current_spline.clear()

          if current_spline:
            cubic_points.pop(0)

          current_spline.extend(cubic_points)
        elif type(segment) == Line:
          p0 = (segment.start.real + off_x, segment.start.imag + off_y)
          p1 = (segment.end.real + off_x, segment.end.imag + off_y)
          p0 = (p0[0] * scale_x, p0[1] * scale_y)
          p1 = (p1[0] * scale_x, p1[1] * scale_y)

          cubic_points = sample_line(p0, p1)

          if cubic_points:
            if current_spline:
              if dist(current_spline[0], cubic_points[0]) < 1:
                current_spline.reverse()

              if dist(current_spline[-1], cubic_points[-1]) < 1:
                cubic_points.reverse()

              if dist(current_spline[-1], cubic_points[0]) > 1:
                spline = create_spline(current_spline, args.error, args.scale)
                if spline:
                  spline_str = "\n".join(spline)
                  if spline_str not in splines:
                    splines.append(spline_str)
                    output.write(spline_str)
                    output.write("\n")
                    output.flush()
                current_spline.clear()

            if current_spline:
              cubic_points.pop(0)

            current_spline.extend(cubic_points)
        else:
          print(type(segment))

      if current_spline:
        spline = create_spline(current_spline, args.error, args.scale)
        if spline:
          spline_str = "\n".join(spline)
          if spline_str not in splines:
            splines.append(spline_str)
            output.write(spline_str)
            output.write("\n")
            output.flush()

    output.write("- Set 15")
