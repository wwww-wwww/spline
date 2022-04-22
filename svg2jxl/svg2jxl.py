import argparse
import math
import random
import svgpathtools
import sys
import io
import zlib
import base64
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


def length(x, y):
  return math.sqrt(math.pow(x, 2) + math.pow(y, 2))


def norm(p1):
  p1_length = length(p1[0], p1[1]) or 1
  return (p1[0] / p1_length, p1[1] / p1_length)


# euclidean distance
def dist(p1, p2):
  return length(p1[0] - p2[0], p1[1] - p2[1])


# same rounded point
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


def create_spline(points, args, color, thickness):
  if len(points) > 2:
    old_len = len(points)
    points = optimize.optimize(points, args.error * (args.scale**2))
    if args.verbose:
      print(f"Optimized spline points: {old_len} -> {len(points)}")

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

  if color[0] != color[1] or color[0] != color[2]:
    thickness *= args.thicknessc
    sat = args.saturationc
  else:
    thickness *= args.thickness
    sat = args.saturation

  bitdepth = args.bitdepth

  if args.sharpness and thickness > 0.5:
    thickness -= 0.25
    bitdepth += 1

  new_color = [
    round((c / 255 - 1) * (2**(bitdepth - 1) - 1)) * sat for c in color
  ]

  spline = [
    "Spline",
    f"{new_color[0]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{new_color[1]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{new_color[2]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{thickness} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
  ]

  spline.extend([f"{round(p[0], 2)} {round(p[1], 2)}" for p in points])
  spline.append("EndSpline")

  return spline


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="Input SVG file")
  parser.add_argument("output", help="Output JXL tree", nargs='?')
  parser.add_argument("-s",
                      "--scale",
                      type=float,
                      default=1,
                      help="Scale for working resolution (default: 1)")
  parser.add_argument("-u",
                      "--upsample",
                      type=int,
                      default=1,
                      help="JXL tree upsample (default: 1)")
  parser.add_argument("-e",
                      "--error",
                      help="Allowed error, multiplied by scale^2  (default: 1)",
                      type=float,
                      default=1)
  parser.add_argument("-b",
                      "--bitdepth",
                      type=int,
                      default=4,
                      help="Bitdepth for colors and JXL (default: 4)")
  parser.add_argument("-t",
                      "--thickness",
                      type=float,
                      default=0.25,
                      help="Thickness for black lines (default: 0.25)")
  parser.add_argument("-tc",
                      "--thicknessc",
                      type=float,
                      default=0.25,
                      help="Thickness for colored lines (default: 0.25)")
  parser.add_argument("-sat",
                      "--saturation",
                      type=float,
                      default=1,
                      help="Saturation for black lines (default: 1)")
  parser.add_argument("-satc",
                      "--saturationc",
                      type=float,
                      default=1,
                      help="Saturation for colored lines (default: 1)")
  parser.add_argument("-sh",
                      "--sharpness",
                      type=int,
                      default=1,
                      help="Add sharpness to thick lines (default: 1)")
  parser.add_argument("-v",
                      "--verbose",
                      action="store_true",
                      help="Print more information to stdout")

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

  header = f"""Bitdepth {args.bitdepth}
Upsample {args.upsample}
Width {round(width)}
Height {round(height)}
RCT 0
"""

  scale_x = width / viewBox[2]
  scale_y = height / viewBox[3]

  output = io.StringIO() if args.output is None else open(args.output, "w+")

  try:
    output.write(header)
    splines = []

    for path, attribute in zip(paths, attributes):
      color = "#000000"
      thickness = 1
      # Currently does not support classes (url(#id))

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

      if "stroke-width" in attribute:
        thickness = float(attribute["stroke-width"])

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

          if not cubic_points:
            continue

          if current_spline:
            if dist(current_spline[0], cubic_points[0]) < 1:
              current_spline.reverse()

            if dist(current_spline[-1], cubic_points[-1]) < 1:
              cubic_points.reverse()

            if dist(current_spline[-1], cubic_points[0]) > 1:
              spline = create_spline(current_spline, args, color, thickness)
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

          if not cubic_points:
            continue

          if current_spline:
            if dist(current_spline[0], cubic_points[0]) < 1:
              current_spline.reverse()

            if dist(current_spline[-1], cubic_points[-1]) < 1:
              cubic_points.reverse()

            if dist(current_spline[-1], cubic_points[0]) > 1:
              spline = create_spline(current_spline, args, color, thickness)
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
        elif type(segment) == Arc:
          path.extend(segment.as_cubic_curves(curves=2))

      if current_spline:
        spline = create_spline(current_spline, args, color, thickness)
        if spline:
          spline_str = "\n".join(spline)
          if spline_str not in splines:
            splines.append(spline_str)
            output.write(spline_str)
            output.write("\n")
            output.flush()

    output.write(f"- Set {2 ** args.bitdepth - 1}")

    if args.output is None:
      output_bytes = output.getvalue().encode()
      output_gzip_bytes = zlib.compress(output_bytes, level=9)
      output_b64_bytes = base64.b64encode(output_gzip_bytes[2:-4], altchars=b"-_")
      zcode = output_b64_bytes.decode().replace("=", "")
      print(f"https://jxl-art.surma.technology/?zcode={zcode}")
  finally:
    output.close()
