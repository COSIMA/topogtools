# topogtools

Tools related to changing ocean model topography and regenerating dependent model inputs.

Below if a list of included tools and short documentation for each.

## bulldozer

Bulldozer is a simple tool to modify the MOM bathymetry/topography file, adding or removing land points.

It is used in a two step process:
1. __mark__: mark out the changes to be made. Each change is written as a line into an output file.
2. __drive__: apply the changes created above to a topography file. The result will be a new file.


