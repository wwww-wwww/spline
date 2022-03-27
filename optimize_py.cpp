#include "Python.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string.h>
#include <vector>

// sampling resolution
#define RES 10

struct Point {
  double x;
  double y;

  Point operator+(const Point &a) { return {x + a.x, y + a.y}; }
  Point operator-(const Point &a) { return {x - a.x, y - a.y}; }
  Point operator+(const double &a) { return {x + a, y + a}; }
  Point operator-(const double &a) { return {x - a, y - a}; }
  Point operator*(const double &a) { return {x * a, y * a}; }
};

double dist(Point p0, Point p1) {
  return sqrt(pow(p1.x - p0.x, 2) + pow(p1.y - p1.x, 2));
}

double shoelace(Point p0, Point p1, Point p2, Point p3) {
  return abs(((p0.x * p1.y - p0.y * p1.x) + (p1.x * p2.y - p1.y * p2.x) +
              (p2.x * p3.y - p2.y * p3.x) + (p3.x * p0.y - p3.y * p0.x)) /
             2);
}

Point catmull_rom(Point points[4], double t) {
  double alpha = 0.5;

  Point p0 = points[0];
  Point p1 = points[1];
  Point p2 = points[2];
  Point p3 = points[3];

  Point a = p1 * 3.0 - p2 * 3.0 + p3 - p0;
  Point b = p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3;
  Point c = p2 - p0;
  Point d = p1 * 2.0;

  return (a * pow(t, 3) + b * pow(t, 2) + c * t + d) * alpha;
}

template <typename T>
double error_remove(std::vector<Point> &points, int remove, T &remove_op) {
  int start = 0;
  int center = 0;
  int end = 0;

  std::vector<Point> line1;
  std::vector<Point> line2;

  for (int i = remove - 3; i <= remove + 3; i++) {
    if (i >= 0 && i < (int)points.size()) {
      if (i == remove - 3 && line1.size() == 0) {
        start = 1;
      }

      if (i == remove) {
        center = line1.size();
      }

      line1.push_back(points[i]);
      line2.push_back(points[i]);

      if (i == remove + 2) {
        end = line1.size() - 1;
      }
    }
  }

  remove_op(line1, line2, start, center, end);

  if (end == 0) {
    end = line1.size() - 1;
  }

  std::vector<Point> line1_spls;
  std::vector<Point> line2_spls;

  auto add_samples = [](std::vector<Point> &out_points,
                        std::vector<Point> points, int i, int count) {
    out_points.push_back(points[i]);

    Point line[4] = {
        i ? points[i - 1] : points[0],
        points[i],
        points[i + 1],
        i + 2 < (int)points.size() ? points[i + 2] : points[i + 1],
    };

    for (int j = 1; j < count; j++) {
      Point p = catmull_rom(line, (double)j / count);
      out_points.push_back(p);
    }
  };

  for (int i = start; i < end; i++) {
    add_samples(line1_spls, line1, i, RES);
  }
  line1_spls.push_back(line1[end]);

  for (int i = start; i < end - 1; i++) {
    add_samples(line2_spls, line2, i, i == center - 1 ? RES * 2 : RES);
  }
  line2_spls.push_back(line2[end - 1]);

  Point last_p1 = line1_spls[0];
  Point last_p2 = line2_spls[0];

  double total_area = 0;
  for (int i = 1; i < (int)line1_spls.size(); i++) {
    total_area += shoelace(last_p1, last_p2, line2_spls[i], line1_spls[i]);
    last_p1 = line1_spls[i];
    last_p2 = line2_spls[i];
  }

  return total_area;
}

static PyObject *optimize_py(PyObject *self, PyObject *args) {
  PyObject *py_points;
  float error;

  PyArg_ParseTuple(args, "Of", &py_points, &error);

  size_t n_points = PyObject_Length(py_points);
  std::vector<Point> points;

  for (size_t i = 0; i < n_points; i++) {
    PyObject *point = PyList_GetItem(py_points, i);
    PyObject *py_x = PyTuple_GetItem(point, 0);
    PyObject *py_y = PyTuple_GetItem(point, 1);

    double x = PyFloat_AsDouble(py_x);
    double y = PyFloat_AsDouble(py_y);

    points.push_back({x, y});
  }

  // delete middle
  int margin = 1;
  if (points.size() - margin * 2 > 0) {
    auto remove_op = [](std::vector<Point> &line1, std::vector<Point> &line2,
                        int start, int center,
                        int end) { line2.erase(line2.begin() + center); };

    std::vector<double> errors;
    errors.resize(points.size() - margin * 2);
    for (int i = margin; i < (int)points.size() - margin; i++) {
      errors[i - margin] = error_remove(points, i, remove_op);
    }

    int min;
    float err = 0;
    while (err < error) {
      min = std::min_element(errors.begin(), errors.end()) - errors.begin();
      err = errors[min];

      if (err < error) {
        errors.erase(errors.begin() + min);
        points.erase(points.begin() + min + margin);
        for (int i = min - 2; i < min + 2; i++) {
          if (i >= 0 && i < (int)errors.size()) {
            errors[i] = error_remove(points, i + margin, remove_op);
          }
        }
      }
    }
  }

  // delete middle and shift surrounding points to the 1/3 and 2/3 positions
  margin = 2;
  if (points.size() - margin * 2 > 0) {
    auto remove_op = [](std::vector<Point> &line1, std::vector<Point> &line2,
                        int start, int center, int end) {
      line2.erase(line2.begin() + center);
      line2[center - 1] = catmull_rom(line1.data() + center - 1, 1.0 / 3);
      line2[center] = catmull_rom(line1.data() + center, 2.0 / 3);
    };

    std::vector<double> errors;
    errors.resize(points.size() - margin * 2);
    for (int i = margin; i < (int)points.size() - margin; i++) {
      errors[i - margin] = error_remove(points, i, remove_op);
    }

    int min;
    float err = 0;
    while (err < error / 2) {
      min = std::min_element(errors.begin(), errors.end()) - errors.begin();
      err = errors[min];

      if (err < error / 2) {
        errors.erase(errors.begin() + min);
        points.erase(points.begin() + min + margin);
        for (int i = min - 2; i < min + 2; i++) {
          if (i >= 0 && i < (int)errors.size()) {
            errors[i] = error_remove(points, i + margin, remove_op);
          }
        }
      }
    }
  }

  PyObject *result = PyList_New(points.size());

  for (size_t i = 0; i < points.size(); i++) {
    PyObject *point = Py_BuildValue("(ff)", points[i].x, points[i].y);
    PyList_SetItem(result, i, point);
  }

  return result;
}

static PyMethodDef methods[] = {
    {"optimize", optimize_py, METH_VARARGS, NULL},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef optimize = {PyModuleDef_HEAD_INIT, "optimize",
                                      "optimize", -1, methods};

PyMODINIT_FUNC PyInit_optimize(void) { return PyModule_Create(&optimize); }
