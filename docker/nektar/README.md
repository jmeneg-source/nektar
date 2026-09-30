# `nektar` image

This is an image containing all of the Nektar++ solvers, pre- and
post-processing tools, and Python bindings for Nektar++. The code is compiled
with virtually all options enabled, including HDF5, MPI and OpenCascade for mesh
generation. Note that to limit this image's size, it does _not_ contain
developer headers or associated packages: for this you can instead use the
`nektarpp/nektar-dev` image.

This image provides a full installation of Nektar++ with the following options
enabled:

- `NEKTAR_USE_MPI`
- `NEKTAR_USE_FFT`
- `NEKTAR_BUILD_PYTHON`
- `NEKTAR_USE_HDF5`
- `NEKTAR_USE_MESHGEN`
- `NEKTAR_USE_PETSC`
- `NEKTAR_USE_CCM`
- `NEKTAR_USE_VTK`

This image is based on either:

- Debian 12 for `latest`;
- Debian 11 for `v5.5.0` and all prior tags.

## Building

The image is built using `nektarpp/nektar-env` or similar and requires as build
context the path to the Nektar++ source tree. It supports several additional
build arguments:

- `ENV_IMAGE` is used to select the environment to build against. This is used
  by the CI to e.g. consistently build against the correct commits. By default
  this is set to `nektarpp/nektar-env:default`.
- `BUILD_DEMOS` can be set to `ON` to build demos, which are disabled by
  default;
- `BUILD_SOLVERS` can be set to `OFF` to disable build of solvers, which are
  enabled by default;
- `BUILD_ARGS` can be set to `ON` to build doxygen documentation, which are
  enabled by default;
- `INSTALL_PREFIX` can be set to adjust the install prefix, which is
  `/usr/local` by default.
  
Then build the image using a command similar to:

```sh
docker build -t nektarpp/nektar -f Dockerfile ~/nektar++
```

Note that this is a multi-stage build: after the initial build phase is
completed, the build tree is erased along with development headers, and the
libraries/executables installed into a fresh Debian 12 container (with
appropriate packages for runtime libraries). If you want to keep the build
files, you should instead only build up to the end of the build stage using
`--target build`:

```sh
docker build --target build -t nektar-build -f Dockerfile ~/nektar++
```

## About Nektar++

Nektar++ is an open-source framework, distributed under the MIT license, for the
spectral/hp element method. These images contain useful environments for users
and developers, and are built automatically with new tags, and commits to the
main branch.

For more information on Nektar++, [see our website](https://www.nektar.info) and
[source code released on our GitLab
instance](https://gitlab.nektar.info/nektar/nektar).
