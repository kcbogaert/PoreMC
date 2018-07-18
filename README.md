# PoreMC
Monte Carlo simulation of nanopore growth in 2D transition metal dichalcogenides

## How to cite
Please cite the following work if you want to use PoreMC.

<insert citation here>

## Usage
### Inputs
**NumIterations** is a guess of the number of iterations used for initializing the area, perimeter, and heatmap outputs and determining x-axis scale in video figure.

**MetalProb** is a number between 0 and 1 which corresponds to the probability that any edge metal atom will be ejected.

**ChalcProb** is same as MetalProb, but for chalcogen atoms.

### Example Inputs
```
[Area, Perimeter, Lattice, HeatMap, AFit,PFit] = PoreMC(1000, 0.03, 0.03);
[Area, Perimeter, Lattice, HeatMap, AFit,PFit] = PoreMC(100, 0.3, 0.3);
```

### Outputs
**Area** is an array. First column is iteration number, second column is number of vacant atomic columns.

**Perimeter** is same as Area, but for the number of edge atoms.

**Lattice** is an array map of where the remaining atoms are located.

**HeatMap** is an array map that labels which iteration caused each ejection to occur.

**AFit** is a quadratic fit of the Area.

**PFit** is a linear fit of the Perimeter.

**MP4 Video** (left) visualizing the evolution of Lattice and (right) plotting Area and Perimeter vs. iteration number.

## Authors
This software was primarily written by [Kevin Bogaert](https://scholar.google.com/citations?user=6m28RhcAAAAJ&hl=en) who is advised by [Prof. Silvija Gradecak](https://dmse.mit.edu/faculty/profile/gradecak).

## License
PoreMC is released under the MIT License.
