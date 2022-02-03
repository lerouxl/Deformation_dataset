# Deformation_dataset

## Description
---
This script is used generate a dataset of random cubes with a physique deformation label. 
As this dataset is intended for an AI needing a fix number of edges, all generated mesh will have the same number of edges.
The deformation is computed using Fenics and is made to look like additive manufacturing deformation.

## Use:
This code is using the fenics library, if you are on windows, I recommand you to use the fenics docker images.
1. Install fenics (Or download the docker image)
1. Clone this repository
1. Install the required library `python3 -m pip install -r requirements.txt`
1. Edit `dataset_generator.py` to have the dataset as you want
1. Launch `dataset_generator.py`

## Examples
Here is an example of a mesh and his deformation. The vtk files can be opened using Paraview.
![Glyph](https://github.com/hy-son/Deformation_dataset/blob/main/imgs/magnitude.JPG)
![Glyph](https://github.com/hy-son/Deformation_dataset/blob/main/imgs/glyph.JPG)

## Thanks
Many thanks to Dr Pierre Kerfriden for coding the inherent strain code deforming the mesh.
Many thanks to Pierre L for his [stl to msh code](https://github.com/Gost65/STL_to_MSH)

