# `nekmesh` image

This image provides a full installation of the `NekMesh` executable with all
mesh generation capabilities enabled via the OpenCASCADE CAD framework. It is
designed to be a smaller image than the `nektarpp/nektar` image if you are
interested in solely running `NekMesh` as a standalone executable.

## Building

The image is built using `nektarpp/nektar-env` and requires as build context the
path to the Nektar++ source tree. It supports several additional build
arguments:

- `ENV_IMAGE` is used to select the environment to build against. This is used
  by the CI to e.g. consistently build against the correct commits. By default
  this is set to `nektarpp/nektar-env:default`.
- `INSTALL_PREFIX` can be set to adjust the install prefix, which is
  `/usr/local` by default.
  
Then build the image using a command similar to:

```sh
docker build -t nektarpp/nekmesh -f Dockerfile ~/nektar++
```

Note that this is a multi-stage build: after the initial build phase is
completed, the build tree is erased along with development headers, and the
libraries/executables installed into a fresh Debian 12 container (with
appropriate packages for runtime libraries). If you want to keep the build
files, you should instead only build up to the end of the build stage using
`--target build`:

```sh
docker build --target build -t nekmesh-build -f Dockerfile ~/nektar++
```

## About Nektar++

Nektar++ is an open-source framework, distributed under the MIT license, for the
spectral/hp element method. These images contain useful environments for users
and developers, and are built automatically with new tags, and commits to the
main branch.

For more information on Nektar++, [see our website](https://www.nektar.info) and
[source code released on our GitLab
instance](https://gitlab.nektar.info/nektar/nektar).
