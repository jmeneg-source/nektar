# `nektar-workbook` image

This is an image based on the `scipy-workbook` image which contains all of the
Nektar++ Python bindings, common pre- and post-processing tools, and some of the
core solvers (e.g. incompressible Navier-Stokes and the advection-diffusion
reaction solver). The code is compiled with virtually all options enabled,
including HDF5, MPI and OpenCascade for mesh generation. It can either be used
standalone, or as part of a Jupyterhub installation for multi-user access.

## Building

The image requires as build context the path to the Nektar++ source tree.  Build
the image using a command similar to:

```sh
docker build -t nektarpp/nektar-notebook -f Dockerfile ~/nektar++
```

## About Nektar++

Nektar++ is an open-source framework, distributed under the MIT license, for the
spectral/hp element method. These images contain useful environments for users
and developers, and are built automatically with new tags, and commits to the
main branch.

For more information on Nektar++, [see our website](https://www.nektar.info) and
[source code released on our GitLab
instance](https://gitlab.nektar.info/nektar/nektar).


