from setuptools import setup, find_packages, Extension

with open("requirements.txt") as fh:
  install_requires = fh.read()

setup(
  name='optimize',
  version='1.0',
  ext_modules=[Extension(
    'optimize',
    ['optimize_py.cpp'],
  )],
  install_requires=install_requires,
  packages=find_packages(),
  entry_points={
    "console_scripts": ["svg2jxl=svg2jxl.svg2jxl:main"],
  },
)
