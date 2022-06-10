import cv2
import numpy as np
import optimize
import random, math

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


# same rounded point
def sim(p1, p2):
  return round(p1[0]) == round(p2[0]) and round(p1[1]) == round(p2[1])


def kmeans_color_quantization(image, clusters=8, rounds=1):
  h, w = image.shape[:2]
  samples = np.zeros([h * w, 3], dtype=np.float32)
  count = 0

  for x in range(h):
    for y in range(w):
      samples[count] = image[x][y]
      count += 1

  compactness, labels, centers = cv2.kmeans(
    samples, clusters, None,
    (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10000, 0.0001),
    rounds, cv2.KMEANS_RANDOM_CENTERS)

  centers = np.uint8(centers)
  res = centers[labels.flatten()]
  return res.reshape((image.shape))


im = cv2.imread("PNG_demo_Banana.png")
(h, w) = im.shape[:2]
kernel = np.ones((3, 3), np.uint8)
kernel2 = np.ones((5, 5), np.uint8)

scale = 1

im = cv2.resize(im, (w // scale, h // scale), interpolation=cv2.INTER_AREA)
im = cv2.fastNlMeansDenoisingColored(im, None, 16, 16, 5, 12)
im = kmeans_color_quantization(im, clusters=32)

im = im

grad = [[0, 1], [1, 1], [1, 0], [1, -1], [0, -1], [-1, -1], [-1, 0], [-1, 1]]

layer1 = [
  f"""Bitdepth {8}
Width {round(w // scale)}
Height {round(h // scale)}
Upsample {2}
RCT 0
Animation
Duration 0
NotLast"""
]

layer2 = ["BlendMode kMul\nNotLast"]

layer3 = ["BlendMode kAdd"]


def length(x, y):
  return math.sqrt(math.pow(x, 2) + math.pow(y, 2))


# euclidean distance
def dist(p1, p2):
  return length(p1[0] - p2[0], p1[1] - p2[1])


def add(p1, p2):
  return (p1[0] + p2[0], p1[1] + p2[1])


def sub(p1, p2):
  return (p1[0] - p2[0], p1[1] - p2[1])


def norm(p1):
  p1_length = length(p1[0], p1[1]) or 1
  return (p1[0] / p1_length, p1[1] / p1_length)


splines1 = []
splines2 = []
splines3 = []


def create_spline(splines, points, color, thickness=0.6):
  points = [(round(p[0]), round(p[1])) for p in points]

  spline = [
    "Spline",
    f"{color[2]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{color[1]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{color[0]} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
    f"{thickness} 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
  ]

  i = 1
  while i < len(points):
    point = points[i - 1]
    if point == points[i]:
      del points[i]
    else:
      i += 1

  if len(points) < 2: return [""]

  if len(points) == 2 and dist(points[0], points[1]) < 2: return [""]

  points_str = "".join(["".join([str(p) for p in point]) for point in points])
  if points_str in splines: return [""]

  for point in points:
    point = " ".join([str(p) for p in point])
    spline.append(point)

  spline.append("EndSpline")

  splines.append(points_str)

  return spline


def area(p):
  return 0.5 * abs(
    sum(x0 * y1 - x1 * y0 for ((x0, y0), (x1, y1)) in segments(p)))


def segments(p):
  return zip(p, p[1:] + [p[0]])


import tqdm

colors = np.unique(im.reshape(-1, im.shape[-1]), axis=0)
for color in tqdm.tqdm(colors):
  mask = np.where(np.all(im == color, 2), 255, 0)
  mask = mask.astype(np.uint8)
  #mask = cv2.erode(mask, kernel)
  #mask = cv2.dilate(mask, kernel)
  ret, labels = cv2.connectedComponents(mask)
  (h, w) = mask.shape[:2]
  color = (color) / 255 * 1.6

  for label in range(1, ret):
    component = np.zeros(labels.shape[:2], dtype=np.uint8)
    component[labels == label] = 255
    a, _ = cv2.findContours(component, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

    mask3 = cv2.dilate(component, kernel)

    mask2 = component
    #mask2 = cv2.erode(component, kernel)
    #mask2 = cv2.dilate(mask2, kernel)
    #mask2 = cv2.dilate(component, kernel)

    # outline
    points = a[0]
    points = points.reshape(points.shape[0], 2)

    points = points.tolist()
    points = [(p[0], p[1]) for p in points]
    points.append(points[0])

    old_len = len(points)
    points = optimize.optimize(points, 4)
    points = [(round(p[0]), round(p[1])) for p in points]

    #if len(points) < 2:
    #  continue

    #if area(points) < 4:
    #  continue

    # scanlines
    (ys, _) = np.where(mask2 == 255)
    if len(ys) == 0: continue
    lines = []
    for y in range(ys[0] + 2, ys[-1], 2):
      (xs, ) = np.where(mask2[y, :] == 255)
      if len(xs) == 0: continue
      x = xs[0]
      while x < xs[-1]:
        while x <= xs[-1] and mask2[y, x] == 0:
          x += 1

        if x > xs[-1]: continue

        if mask2[y, x] == 0: continue

        start = x
        while x < xs[-1] and mask2[y, x + 1] == 255:
          x += 1

        if x <= start:
          x = start + 1
          continue

        line = [(start + 1, y), (x - 1, y)]
        lines.append(line)

        x += 1

    def extend(line, line2):
      line.append(line[-1])
      line[-2] = add(line[-2], norm(sub(line[-3], line[-2])))
      line2.insert(0, line2[0])
      line2[1] = add(line2[1], norm(sub(line2[2], line2[1])))
      line.extend(line2)

    if lines:
      i = 0
      while i < len(lines):
        line = lines[i]
        while True:
          changed = False
          for line2 in lines[i:]:
            if line2[0][1] == line[-1][1]: continue
            if line2[-1][1] == line[-1][1]: continue
            if abs(line2[0][1] - line[-1][1]) > 2: continue
            if abs(line2[-1][1] - line[-1][1]) > 2: continue
            dista = abs(line[-1][0] - line2[0][0])
            distb = abs(line[-1][0] - line2[-1][0])
            mida = (line[-1][0] + line2[0][0]) // 2
            midb = (line[-1][0] + line2[-1][0]) // 2
            if line[-2][0] > line[-1][0] and dista <= 8 and mask3[
                round(line[-1][1]), mida] != 0:
              extend(lines[i], line2.copy())
              lines.remove(line2)
              changed = True
              break
            if line[-2][0] < line[-1][0] and distb <= 8 and mask3[
                round(line[-1][1]), midb] != 0:
              extend(lines[i], line2[::-1])
              lines.remove(line2)
              changed = True
              break
          if not changed: break
        i += 1

    for line in lines:
      line = optimize.optimize(line, 2)
      layer1.extend(create_spline(splines1, line, color, 0.6))

    # add outline

    if area(points) < 8:
      continue

    for i in range(len(points) - 1):
      if sim(points[i], points[i + 1]):
        direction = random.choice(gradients)
        points[i +
               1] = [points[i][0] + direction[0], points[i][1] + direction[1]]

    layer2.extend(create_spline(splines2, points, [-2] * 3, 0.5))
    layer3.extend(create_spline(splines3, points, color, 0.5))

    # outline inner components
    component_a = component.copy()
    for contour in a:
      cv2.drawContours(component_a, [contour], 0, 255, -1)

    inner_components = np.bitwise_and(255 - component, component_a)
    ret2, labels2 = cv2.connectedComponents(inner_components)

    for label2 in range(1, ret2):
      inner = np.zeros(labels2.shape[:2], dtype=np.uint8)
      inner[labels2 == label2] = 255
      if inner[0, 0] and inner[-1, -1] and inner[0, -1] and inner[-1, 0]:
        continue
      a, _ = cv2.findContours(inner, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

      points = a[0]
      points = points.reshape(points.shape[0], 2)

      points = points.tolist()
      points = [(p[0], p[1]) for p in points]
      points.append(points[0])

      old_len = len(points)
      points = optimize.optimize(points, 2)
      points = [(round(p[0]), round(p[1])) for p in points]
      if len(points) < 2:
        continue

      for i in range(len(points) - 1):
        if sim(points[i], points[i + 1]):
          direction = random.choice(gradients)
          points[i + 1] = [
            points[i][0] + direction[0], points[i][1] + direction[1]
          ]

      layer2.extend(create_spline(splines2, points, [-2] * 3, 0.5))
      layer3.extend(create_spline(splines3, points, color, 0.5))

layer1 = [s for s in layer1 if s]
layer2 = [s for s in layer2 if s]
layer2 = [s for s in layer2 if s]

layer1.append("- Set 0")
layer2.append("- Set 255")
layer3.append("- Set 0")

open("test.txt", "w+").write("\n".join(layer1 + layer2 + layer3))
