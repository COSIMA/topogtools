# topogtools

Tools related to changing ocean model topography and regenerating dependent model inputs.

Below if a list of included tools and short documentation for each.

## bulldozer

Bulldozer is a simple tool to modify the MOM bathymetry/topography file, adding or removing land points.

It is used as follows:

```bash
echo "112, 246, 0.0, 50.0" | ./bulldozer.py topog.nc --new_topog new_topog.nc
```

The input is taken from stdin. The format of the input is a list of
comma-separated line with the following values:

'i index', 'j index', 'original depth', 'new depth'

For example to run on a single point:

echo "112, 246, 0.0, 50.0" | ./bulldozer.py test/topog.nc

To run a whole file:

cat file.csv | ./bulldozer.py test/topog.nc

Where the contents of file could look something like this:

112, 246, 0.0, 50.0
113, 246, 0.0, 50.0

