# Gravo MG C++ library
[[Paper]](https://graphics.tudelft.nl/~klaus/papers/Gravo_MG.pdf) [[Project page]](https://rubenwiersma.nl/gravomg)

![](https://rubenwiersma.nl/assets/img/publications/gravomg/teaser_gravomg.png)

C++ library implementing Gravo MG. Gravo MG is a Geometric Multigrid Method for solving linear systems on curved surfaces. For more information, check out our [project page](https://rubenwiersma.nl/gravomg).

## Linking in CMake
You can add this library to your project by simply including the files in your folder for dependencies (e.g., by using git submodules) and then adding the subdirectory to your CMake file:
```cmake
add_subdirectory(deps/gravomg)
target_link_libraries(project_name PRIVATE gravomg)
```

## Usage
Create the solver and construct a hierarchy:
```cpp
GravoMG::MultigridSolver(positions, neighbors, mass)
solver->buildHierarchy();
```

Solve a linear system
```cpp
Eigen::SparseMatrix<double> lhs;
Eigen::MatrixXd rhs, x;
solver->solve(lhs, rhs, x);
```

## Citations
Please cite our paper if this code contributes to an academic publication:

```bib
@Article{WiersmaNasikun2023GravoMG,
author = {Ruben Wiersma, Ahmad Nasikun, Elmar Eisemann, Klaus Hildebrandt},
journal = {SIGGRAPH 2023},
title = {A Fast Geometric Multigrid Method for Curved Surfaces},
year = {2023},
month = jul,
number = {4},
volume = {41},
doi = {10.1145/3588432.3591502},
publisher = {ACM}
}
```