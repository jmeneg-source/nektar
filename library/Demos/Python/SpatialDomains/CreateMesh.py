###############################################################################
##
## File: CreateMeshGraph_UnitTest.py
##
## For more information, please see: http://www.nektar.info
##
## The MIT License
##
## Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
## Department of Aeronautics, Imperial College London (UK), and Scientific
## Computing and Imaging Institute, University of Utah (USA).
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included
## in all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
## OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Description: Unit tests for the Equation class.
##
###############################################################################

import itertools
import os.path
import re
import shutil
import tempfile

import NekPy.SpatialDomains as SD

NX = 4
NY = 5

FLOAT = r"-?\d\.\d+e[+-]\d\d"
VERTICES = re.compile(
    r'<\s*V\s+ID="\d+"\s*>\s*'
    + FLOAT
    + r"\s+"
    + FLOAT
    + r"\s+"
    + FLOAT
    + r"\s*</\s*V\s*>",
    re.I,
)
SEGMENTS = re.compile(r'<\s*E\s+ID="\d+"\s*>\s*\d+\s+\d+\s*</\s*E\s*>', re.I)
QUADS = re.compile(r'<\s*Q\s+ID="\d+"\s*>\s*\d+\s+\d+\s+\d+\s+\d+\s*</\s*Q\s*>', re.I)
COMPOSITES = re.compile(r'<\s*C\s+ID="\d+"\s*>.*</\s*C\s*>', re.I)
DOMAINS = re.compile(r'<\s*D\s+ID="\d+"\s*>.*</\s*D\s*>', re.I)


def check_count(item_type, actual_count, expected_count):
    if actual_count != expected_count:
        print(
            "ERROR: Found {} {}; expected {}".format(
                actual_count, item_type, expected_count
            )
        )
        return False
    else:
        print("Found {} {}, as expected".format(actual_count, item_type))
        return True


def main():
    print("Creating {} by {} mesh".format(NX - 1, NY - 1))
    mesh = SD.MeshGraph()
    mesh.SetSpaceDimension(2)
    mesh.SetMeshDimension(2)
    composites = mesh.GetComposites()
    domains = mesh.GetDomain()
    x_bound_lower = []
    x_bound_upper = []
    y_bound_lower = []
    y_bound_upper = []

    print("Generating points...")
    for i, (x, y) in enumerate(itertools.product(range(NX), range(NY))):
        mesh.AddGeom(SD.PointGeom(2, i, x, y, 0.0))
    print("Generating vertical segments...")
    for i, (x, y) in enumerate(itertools.product(range(NX), range(NY - 1))):
        mesh.AddGeom(SD.SegGeom(i, 2, [mesh.points[x * NY + y], mesh.points[x * NY + y + 1]]))
        if x == 0:
            x_bound_lower.append(mesh.segments[i])
        if x == NX - 1:
            x_bound_upper.append(mesh.segments[i])
    print("Generating horizontal segments...")
    horiz_start = NX * (NY - 1)
    for i, (x, y) in enumerate(
        itertools.product(range(NX - 1), range(NY)), start=horiz_start
    ):
        mesh.AddGeom(SD.SegGeom(i, 2, [mesh.points[x * NY + y], mesh.points[(x + 1) * NY + y]]))
        if y == 0:
            y_bound_lower.append(mesh.segments[i])
        if y == NY - 1:
            y_bound_upper.append(mesh.segments[i])
    print("Generating quads...")
    for i, (x, y) in enumerate(
        itertools.product(range(NX - 1), range(NY - 1))
    ):
        mesh.AddGeom(SD.QuadGeom(
            i,
            [
                mesh.segments[x * (NY - 1) + y],
                mesh.segments[horiz_start + x * NY + y],
                mesh.segments[(x + 1) * (NY - 1) + y],
                mesh.segments[horiz_start + x * NY + y + 1],
            ],
        ))
        assert mesh.quads[i].IsValid()
    print("Generating domain...")
    composites[0] = SD.Composite(list(mesh.quads.values()))
    comp_map = SD.CompositeMap()
    comp_map[0] = composites[0]
    domains[0] = comp_map
    print("Generating boundaries...")
    composites[1] = SD.Composite(x_bound_lower)
    composites[2] = SD.Composite(x_bound_upper)
    composites[3] = SD.Composite(y_bound_lower)
    composites[4] = SD.Composite(y_bound_upper)
    # Write out the mesh then read it in as a string
    tmpdir = tempfile.mkdtemp()
    try:
        writer = SD.MeshGraphIO.Create("Xml")
        writer.SetMeshGraph(mesh)
        outfile = os.path.join(tmpdir, "output.xml")
        print("Writing mesh to temporary file " + outfile)
        writer.Write(outfile, True, SD.FieldMetaDataMap())
        with open(outfile) as f:
            xml = f.read()
    finally:
        shutil.rmtree(tmpdir)
    # Check the mesh is of hte expected size
    passing = check_count("vertices", len(VERTICES.findall(xml)), NX * NY)
    passing = check_count("segments", len(SEGMENTS.findall(xml)),
        (NX - 1) * NY + NX * (NY - 1)
    ) and passing
    passing = check_count("quads", len(QUADS.findall(xml)), (NX - 1) * (NY - 1)) and passing
    passing = check_count("composites", len(COMPOSITES.findall(xml)), 5) and passing
    passing = check_count("domains", len(DOMAINS.findall(xml)), 1) and passing
    if passing:
        print("Test successful!")
    else:
        print("Test unsuccessful")


if __name__ == "__main__":
    main()
