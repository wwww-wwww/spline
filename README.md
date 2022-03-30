# svg2jxl

## Installation

```
python3 setup.py build
python3 setup.py install --user
```


## Usage

```
python3 svg2jxl.py [OPTIONS] INPUT OUTPUT
```

## Options

```
-h, --help                                    show this help message and exit
-s SCALE, --scale SCALE                       Scale for working resolution
-u UPSAMPLE, --upsample UPSAMPLE              JXL tree Upsample
-e ERROR, --error ERROR                       Allowed error, multiplied by scale^2
-b BITDEPTH, --bitdepth BITDEPTH              Bitdepth for colors and JXL
-t THICKNESS, --thickness THICKNESS           Thickness for black lines
-tc THICKNESSC, --thicknessc THICKNESSC       Thickness for colored lines
-sat SATURATION, --saturation SATURATION      Saturation for black lines
-satc SATURATIONC, --saturationc SATURATIONC  Saturation for colored lines
```

## Example
```
svg2jxl image.svg tree.txt -s 1.5 -e 1.3 -b 4 -t 0.5
jxl_from_tree tree.txt tree.jxl
```
