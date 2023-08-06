# Gravo MG C++ library
Library that can be included from any C++ project to use Gravo MG. Gravo MG is a Geometric Multigrid Method for solving linear systems on curved surfaces. For more information, check out our [project page](https://rubenwiersma.nl/gravomg).

## Linking in CMake
You can add this library to your project by simply including the files in your folder for dependencies (e.g., by using git submodules) and then adding the subdirectory to your CMake file:
```cmake
add_subdirectory(deps/gravomg)
target_link_libraries(project_name PRIVATE gravomg)
```

## Usage
Create the solver and construct a hierarchy:
```cpp
MGBS::MultigridSolver(positions, neighbors, mass)
solver->buildHierarchy();
```

Solve a linear system
```cpp
Eigen::SparseMatrix<double> lhs;
Eigen::MatrixXd rhs, x;
solver->solve(lhs, rhs, x);
```