# topogtools

Tools related to changing ocean model topography and regenerating dependent model inputs.

Below if a list of included tools and short documentation for each.

## bulldozer

Bulldozer is a simple tool to modify the MOM bathymetry/topography file, adding or removing land points.

It is used as follows:

```bash
./bulldozer.py orig_topog.nc new_topog.nc changes.csv
```

It takes `orig_topog.nc` as input and applies `changes.csv` to create `new_topog.nc`.

`changes.csv` is a comma separated list of changes, one per line. It must contain a header like the following:

```
i index, j index, original depth, new depth
```
