# topogtools

Tools related to changing ocean model topography and regenerating dependent model inputs.

Below if a list of included tools and short documentation for each.

## bulldozer

bulldozer.py is a simple tool to modify the MOM bathymetry/topography file, adding or removing land points.

It is used as follows:

```bash
echo -e "112, 246, 0.0, 50.0" | ./bulldozer.py topog.nc --new_topog new_topog.nc
```

The input is taken from stdin. The format of the input is a list of
comma-separated line with the following values:

'i index', 'j index', 'original depth', 'new depth'

For example to run on a single point:

```bash
echo -e "112, 246, 0.0, 50.0" | ./bulldozer.py test/topog.nc
```

To run a whole file:

```bash
cat file.csv | ./bulldozer.py test/topog.nc
```

Where the contents of file could look something like this:

112, 246, 0.0, 50.0
113, 246, 0.0, 50.0

## unmask

unmask.py removes the masked regions of a variable. For example it can fill in land with value from the nearest ocean point.

This is useful when bulldozer.py is used to modify the land-sea mask, potentially creating undefined values in masked initial conditions. e.g. an OASIS restart has 0's on land by default so when we change the mask some of these may become ocean points and cause the model to crash. In this case unmask.py would be used to remove all land points.

Example use:

```bash
./unmask.py test/test_data/i2o.nc test/test_data/kmt.nc kmt --output_file test/test_data/new_i2o.nc --flip_mask
```

## topog2mask

topog2mask.py takes a topog file and outputs a mask.

Example use:

```bash
./topog2mask.py test/test_data/topog.nc test/test_data/new_kmt.nc
```

