#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

struct Point {
  double x, y;
  double dx, dy;
};

struct IT {
  int i;
  double t;
};

struct ITT {
  int i;
  double t;
  double t2;
};

struct ITPoint {
  int i;
  double t;
  Point p;
};

double dist(Point a, Point b) {
  return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

bool sim(Point a, Point b) {
  return round(a.x) == round(b.x) && round(a.y) == round(b.y);
}

double length(Point p) { return sqrt(pow(p.x, 2) + pow(p.y, 2)); }

double cubic(double a, double b, double c, double d, double x) {
  return a * pow(x, 3) + b * pow(x, 2) + c * x + d;
}

Point catmull_rom(Point p0, Point p1, Point p2, Point p3, double t) {
  double alpha = 0.5;

  double a1, b1, c1, d1, a2, b2, c2, d2;

  a1 = -p0.x + 3 * p1.x - 3 * p2.x + p3.x;
  a2 = -p0.y + 3 * p1.y - 3 * p2.y + p3.y;

  b1 = 2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x;
  b2 = 2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y;

  c1 = -p0.x + p2.x;
  c2 = -p0.y + p2.y;

  d1 = 2 * p1.x;
  d2 = 2 * p1.y;

  double t3 = pow(t, 3);
  double t2 = pow(t, 2);

  Point p = {(a1 * t3 + b1 * t2 + c1 * t + d1) * alpha,
             (a2 * t3 + b2 * t2 + c2 * t + d2) * alpha,
             (3 * t2 * a1 + 2 * t * b1 + c1) * alpha,
             (3 * t2 * a2 + 2 * t * b2 + c2) * alpha};

  return p;
}

Point catmull_rom(const std::vector<Point> points, int start, double t) {
  Point p0, p1, p2, p3;
  p0 = start ? points[start - 1] : points[0];
  p1 = points[start + 0];
  p2 = points[start + 1];
  p3 = start < points.size() - 2 ? points[start + 2] : p2;

  return catmull_rom(p0, p1, p2, p3, t);
}

double cubic_solve(double a, double b, double c, double d, Point p0, Point p1,
                   Point p2, Point p3, Point target) {
  if (a == 0) {
    return 0;
  }

  if (d == 0) {
    return 0;
  }

  b /= a;
  c /= a;
  d /= a;

  double disc, q, r, dum1, s, t, term1, r13;
  q = (3.0 * c - (b * b)) / 9.0;
  r = -(27.0 * d) + b * (9.0 * c - 2.0 * (b * b));
  r /= 54.0;
  disc = q * q * q + r * r;
  term1 = b / 3.0;

  // One real root, two are complex.
  if (disc > 0) {
    s = r + sqrt(disc);
    s = (s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, (1.0 / 3.0));
    t = r - sqrt(disc);
    t = (t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, (1.0 / 3.0));
    double root1 = -term1 + s + t;
    term1 += (s + t) / 2.0;
    double real = term1;

    double y1 = dist(catmull_rom(p0, p1, p2, p3, root1), target);
    double y2 = dist(catmull_rom(p0, p1, p2, p3, real), target);
    if (y1 < y2) {
      return root1;
    } else {
      return real;
    }
  }

  // All roots real, at least two are equal.
  if (disc == 0) {
    r13 = (r < 0) ? -pow(-r, (1.0 / 3.0)) : pow(r, (1.0 / 3.0));
    double root1 = -term1 + 2.0 * r13;
    double root2 = -(r13 + term1);
    double y1 = dist(catmull_rom(p0, p1, p2, p3, root1), target);
    double y2 = dist(catmull_rom(p0, p1, p2, p3, root2), target);
    if (y1 < y2) {
      return root1;
    } else {
      return root2;
    }
  }

  // All roots real
  q = -q;
  dum1 = q * q * q;
  dum1 = acos(r / sqrt(dum1));
  r13 = 2.0 * sqrt(q);

  double root1 = -term1 + r13 * cos(dum1 / 3.0);
  double root2 = -term1 + r13 * cos((dum1 + 2.0 * M_PI) / 3.0);
  double root3 = -term1 + r13 * cos((dum1 + 4.0 * M_PI) / 3.0);

  double y1 = dist(catmull_rom(p0, p1, p2, p3, root1), target);
  double y2 = dist(catmull_rom(p0, p1, p2, p3, root2), target);
  double y3 = dist(catmull_rom(p0, p1, p2, p3, root3), target);

  if (y1 < y2 && y1 < y3) {
    return root1;
  } else if (y2 < y3) {
    return root2;
  } else {
    return root3;
  }
}

ITPoint catmull_rom_solve(const std::vector<Point> &points, int start,
                          Point target) {
  double alpha = 0.5;
  Point p0, p1, p2, p3;
  double a1, b1, c1, d1, a2, b2, c2, d2;

  p0 = start > 0 ? points[start - 1] : points[0];
  p1 = points[start + 0];
  p2 = points[start + 1];
  p3 = start < points.size() - 2 ? points[start + 2] : p2;

  a1 = -p0.x + 3 * p1.x - 3 * p2.x + p3.x;
  b1 = 2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x;
  c1 = -p0.x + p2.x;
  d1 = 2 * p1.x;

  a2 = -p0.y + 3 * p1.y - 3 * p2.y + p3.y;
  b2 = 2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y;
  c2 = -p0.y + p2.y;
  d2 = 2 * p1.y;

  double root_x =
      cubic_solve(a1, b1, c1, d1 - target.x / alpha, p0, p1, p2, p3, target);
  double root_y =
      cubic_solve(a2, b2, c2, d2 - target.y / alpha, p0, p1, p2, p3, target);
  Point p_x = catmull_rom(points, start, root_x);
  Point p_y = catmull_rom(points, start, root_y);

  if (dist(p_x, target) < dist(p_y, target)) {
    return {start, root_x, p_x};
  } else {
    return {start, root_y, p_y};
  }
}

// find t
Point find_on_spline(const std::vector<Point> &points, Point point, int start) {
  std::vector<ITPoint> candidate_t;

  int first = std::max(0, start - 4);
  int last = std::min((int)points.size() - 1, start + 4);

  for (int i = first; i < last; i++) {
    candidate_t.push_back(catmull_rom_solve(points, i, point));
  }

  ITPoint closest = candidate_t[0];
  double closest_distance = dist(point, closest.p);
  int closest_index = 0;

  for (int i = 1; i < candidate_t.size(); i++) {
    double candidate_distance = dist(point, candidate_t[i].p);
    if (candidate_distance < closest_distance) {
      closest_distance = candidate_distance;
      closest_index = i;
      closest = candidate_t[i];
    }
  }

  return closest.p;
}

std::vector<Point> copy_points(const std::vector<Point> &points) {
  std::vector<Point> new_points(points.size());
  memcpy(new_points.data(), points.data(), sizeof(Point) * points.size());
  return new_points;
}

double polygonArea(double X[], double Y[], int n) {
  // Initialize area
  double area = 0.0;

  // Calculate value of shoelace formula
  int j = n - 1;
  for (int i = 0; i < n; i++) {
    area += (X[j] + X[i]) * (Y[j] - Y[i]);
    j = i; // j is previous vertex to i
  }

  // Return absolute value
  return abs(area / 2.0);
}

double calculate_error(const std::vector<Point> samples,
                       const std::vector<Point> &points, int start) {

  std::vector<Point> samples_points;
  for (Point point : samples) {
    samples_points.push_back(find_on_spline(points, point, start));
  }

  double area = 0;
  for (int i = 0; i < samples.size() - 1; i++) {
    Point p0 = samples[i];
    Point p1 = samples[i + 1];

    Point p2 = samples_points[i];
    Point p3 = samples_points[i + 1];

    double x[] = {p0.x, p1.x, p3.x, p2.x};
    double y[] = {p0.y, p1.y, p3.y, p2.y};

    area += polygonArea(x, y, 4);
  }

  return area;
}

double arc_length(const std::vector<Point> &points, int start, double res = 1) {
  double l = 0;
  Point last_point = catmull_rom(points, start, 0);
  for (double t = res / length({last_point.dx, last_point.dy}); t < 1;) {
    Point p = catmull_rom(points, start, t);
    l += dist(last_point, p);
    last_point = p;
    t += res / length(p);
  }
  l += dist(catmull_rom(points, start, 1), last_point);
  return l;
}

std::vector<double> equidistant(const std::vector<Point> &points, int start,
                                double res) {
  std::vector<double> sample_t;
  Point last_point = catmull_rom(points, start, 0);
  for (double t = res / length({last_point.dx, last_point.dy}); t < 1;) {
    Point p = catmull_rom(points, start, t);
    last_point = p;
    sample_t.push_back(t);
    t += res / length(p);
  }
  return sample_t;
}

// A - B - (C) - D - E
std::vector<Point> get_samples(std::vector<Point> points, int start,
                               double res = 1) {
  std::vector<Point> samples;

  // [A - B)
  if (start > 1) {
    samples.push_back(points[start - 2]);
    for (double t : equidistant(points, start - 2, res)) {
      samples.push_back(catmull_rom(points, start - 2, t));
    }
  }
  // [B - C)
  samples.push_back(points[start - 1]);
  for (double t : equidistant(points, start - 1, res)) {
    samples.push_back(catmull_rom(points, start - 1, t));
  }
  // [C - D]
  samples.push_back(points[start]);
  for (double t : equidistant(points, start, res)) {
    samples.push_back(catmull_rom(points, start, t));
  }
  samples.push_back(points[start + 1]);
  // (D - E]
  if (start < points.size() - 2) {
    for (double t : equidistant(points, start + 1, res)) {
      samples.push_back(catmull_rom(points, start + 1, t));
    }
    samples.push_back(points[start + 2]);
  }

  return samples;
}

IT optimize(std::vector<IT> &errors, std::vector<Point> &points, int start = 1,
            int end = -1, double move = 0, int recalculate = 0) {
  for (int i = start; i < points.size() + end; i++) {
    if (recalculate == -1 || recalculate != -2 && abs(i - recalculate) > 5) {
      continue;
    }

    double d = arc_length(points, i - 1) + arc_length(points, i);
    if (i > 1) {
      d += arc_length(points, i - 2);
    }
    if (i < points.size() - 2) {
      d += arc_length(points, i + 1);
    }

    std::vector<Point> new_points = copy_points(points);
    new_points.erase(new_points.begin() + i);

    if (i < points.size() + end - 1) {
      new_points[i] = catmull_rom(points, i, move);
    }

    std::vector<Point> samples = get_samples(points, i);
    double error = calculate_error(samples, new_points, i);
    errors[i] = {i, error};
    // if (error < d * 0.1) {
    //   points = new_points;
    //   i++;
    // }
  }

  IT smallest = errors[start];
  int in = start;

  for (int i = start + 1; i < points.size() + end; i++) {
    if (errors[i].t < smallest.t) {
      in = i;
      smallest = errors[i];
    }
  }

  return smallest;

  return {};
}

void step(std::vector<Point> &points, int direction) {
  if (direction > 0) {
    Point p1 = catmull_rom(points, 0, 0);
    Point p2 = {p1.x + p1.dx, p1.y + p1.dy};
    p2 = find_on_spline(points, p2, 0);

    double d = dist(p1, p2);

    if (d < 2) {
      points.insert(points.begin() + 1, p2);
    }
  } else {
    Point p1 = catmull_rom(points, points.size() - 2, 1);
    Point p2 = {p1.x - p1.dx, p1.y - p1.dy};
    p2 = find_on_spline(points, p2, points.size() - 2);

    double d = dist(p1, p2);

    if (d < 2) {
      points.insert(points.begin() + points.size() - 1, p2);
    }
  }
}

int main() {
  std::ifstream in("original_spline");

  std::vector<Point> points;

  double x, y;
  while (in >> x >> y) {
    points.push_back({x, y});
  }

  size_t starting_size = points.size();
  size_t current_size = starting_size;
  printf("%zu points\n", starting_size);

  step(points, 1);
  step(points, -1);

  // Greedy
  /*
  std::vector<IT> errors(points.size());
  while (true) {
    optimize(errors, points, 2, -2, 0.5);
    if (current_size != points.size()) {
      current_size = points.size();
    } else {
      break;
    }
  }
  printf("done first step\n");

  for (int i = 0; i <= 8; i++) {
    printf("step %d\n", i);
    optimize(errors, points, 2, -2, i / 8.0);
  }
  */


  // Globally optimal
  std::vector<std::vector<IT>> errors(9);
  for (int i = 0; i < errors.size(); i++) {
    errors[i].resize(points.size());
    optimize(errors[i], points, 2, -2, i / (errors.size() - 1.0), -2);
  }

  printf("done first step\n");

  int last_delete = -1;
  std::vector<IT> best_error(errors.size());

  while (points.size() > 4) {
    for (int j = 0; j < errors.size(); j++) {
      best_error[j] = optimize(errors[j], points, 2, -2,
                               j / (errors.size() - 1.0), last_delete);
    }

    IT smallest = best_error[0];
    int in = 0;
    for (int j = 1; j < best_error.size(); j++) {
      if (best_error[j].t < smallest.t) {
        smallest = best_error[j];
        in = j;
      }
    }

    printf("%f\n", smallest.t);
    if (smallest.t < 1) {
      Point new_point =
          catmull_rom(points, smallest.i, in / (errors.size() - 1.0));
      points.erase(points.begin() + smallest.i);
      for (int j = 0; j < errors.size(); j++) {
        errors[j].erase(errors[j].begin() + smallest.i);
        for (int k = smallest.i; k < errors[j].size(); k++) {
          errors[j][k].i -= 1;
        }
      }

      last_delete = smallest.i;
      points[smallest.i] = new_point;
    } else {
      break;
    }
  }

  if (points.size() != starting_size) {
    printf("%d\n", (int)(points.size() - starting_size));
  } else {
    printf("0\n");
  }

  std::ofstream out("optimized_spline");
  for (Point point : points) {
    out << point.x << " " << point.y << "\n";
  }

  return 0;
}
