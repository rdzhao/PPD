# ProteinPocketDetection

## Introduction
This is the source code for paper "Protein Pocket Detection via Convex Hull Surface Evolution and Associated Reeb Graph".
This only supports Linux OS.

## Dependencies
Make sure all the following libraries installed in your system.

[Boost](https://www.boost.org/)
[OpenEXR](http://www.openexr.com/)
[TBB](https://www.threadingbuildingblocks.org/)
[CGAL](https://www.cgal.org/)
[OpenVDB](http://www.openvdb.org/)

## Building
1. Create a new folder under root directory
2. Command "cmake .."
3. Command "make"

## Usage
./CavityDetection obj_file xyzr_file minimum_depth minimum_area minimum_volume
./CavityDetection: the executable
obj_file: molecular surface model obj file.
xyzr_file: atoms information file, xyz coordinate and radius per line for each atom.
minimum_depth: minimum required pocket depth
minimum_area: minimum required pocket cross-section area
minimum_volume: minimum required pocket volume.
