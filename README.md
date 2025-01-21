# README

## Project Overview

This project implements an **analytic placement algorithm** to solve a standard cell placement problem inspired by the FastPlace method described in the provided IEEE paper. The algorithm includes cell placement, quadratic optimization, cell spreading, and visualization using Python. The project leverages several programming components, including:

- **C++ modules** for parsing benchmark files, matrix computations, and placement algorithms.
- **Python scripts** for generating graphical representations of cell and pad placements.
- **Generative AI assistance**, primarily ChatGPT, for code generation and debugging.

### Key Features
1. **Quadratic Placement Algorithm:**
   - Converts wirelength minimization into a quadratic programming problem.
   - Uses a conjugate gradient solver for efficient optimization.
   
2. **Cell Spreading:**
   - Employs a bin-based approach to distribute cells uniformly over the chip area while maintaining relative order.

3. **Visualization:**
   - Python GUI to plot cell and pad placements before and after spreading.
   - Zoom-in and zoom-out capabilities for detailed inspection.

4. **File Handling:**
   - Reads benchmark files (.net, .are, .kiaPad).
   - Outputs wirelength values, cell coordinates, and graphical plots.

5. **Hybrid Net Model:**
   - Incorporates clique and star models to represent nets, reducing non-zero entries in the connectivity matrix for faster computation.

---

## Code Breakdown

### 1. **C++ Code**

#### File: `main.cpp`
This is the main entry point for the program. It:
- Parses input benchmark files using `suraj_parser.cpp`.
- Initializes matrices and data structures for placement.
- Outputs intermediate results, including parsed cell and pad information.

#### File: `suraj_parser.cpp` & `suraj_parser.h`
These files handle the parsing of benchmark files:
- `.net` files: Describe the netlist and connectivity.
- `.are` files: Specify cell areas (ignored in this project).
- `.kiaPad` files: Define pad locations.
- The parser creates the connectivity matrix (`matQ`) and initializes vectors for optimization.

#### Key Functions:
- **`conjugateGradient`**: Solves the quadratic minimization problem.
- **`cellSpreading`**: Implements the grid-based cell spreading to minimize overlaps.
- **`writeCoordinatesToFile`**: Outputs the final cell coordinates.

#### Implementation Bugs Introduced by Generative AI:
1. Incorrect indexing in cell spreading logic.
   - Fixed by carefully analyzing array bounds and ensuring all bins are utilized.
2. Overlooked edge cases in the conjugate gradient solver.
   - Addressed by adding convergence checks and boundary conditions.

### 2. **Python Code**

#### File: `plot_placement.py`
This script visualizes the placement results:
- Reads `.kiaPad` files to extract cell and pad coordinates.
- Plots the coordinates using `matplotlib`.
- Saves images of the placement before and after cell spreading.

#### Usage:
```bash
python3 plot_placement.py <input_file> <output_image>
```

---

## How the Code Uses FastPlace Concepts

### Quadratic Placement
- The connectivity matrix (`matQ`) is constructed using the clique and star hybrid model.
- A conjugate gradient solver minimizes the quadratic objective, producing initial cell placements.
- The quadratic objective function represents the total wirelength as a sum of weighted distances between connected cells. Minimizing this function involves solving a linear system of equations derived from the connectivity matrix (`matQ`).
- The conjugate gradient method is particularly efficient for this task because the `matQ` matrix is symmetric, positive-definite, and sparse, aligning well with the method's strengths. The algorithm iteratively refines the solution, converging quickly without requiring explicit matrix inversion, which would be computationally expensive.
- Each iteration of the conjugate gradient method evaluates the residual error and adjusts the cell placements to reduce the quadratic cost, ensuring both accuracy and computational efficiency.
- **Hybrid Net Model:** This project uses a mixed approach by combining the clique and star net models. For smaller nets (2-3 pins), the clique model is used to minimize complexity, while the star model is employed for larger nets (4+ pins) to reduce non-zero entries in the connectivity matrix. This hybrid approach ensures a balance between accuracy and computational efficiency, as described in the FastPlace methodology.

![image](https://github.com/user-attachments/assets/2cdd32db-837a-4cfc-8540-c1961fad9c7a)

- The connectivity matrix (`matQ`) is constructed using the clique and star hybrid model.
- A conjugate gradient solver minimizes the quadratic objective, producing initial cell placements.
- The quadratic objective function represents the total wirelength as a sum of weighted distances between connected cells. Minimizing this function involves solving a linear system of equations derived from the connectivity matrix (`matQ`).
- The conjugate gradient method is particularly efficient for this task because the `matQ` matrix is symmetric, positive-definite, and sparse, aligning well with the method's strengths. The algorithm iteratively refines the solution, converging quickly without requiring explicit matrix inversion, which would be computationally expensive.
- Each iteration of the conjugate gradient method evaluates the residual error and adjusts the cell placements to reduce the quadratic cost, ensuring both accuracy and computational efficiency.
- The connectivity matrix (`matQ`) is constructed using the clique and star hybrid model.
- A conjugate gradient solver minimizes the quadratic objective, producing initial cell placements.

### Cell Spreading
- A 5x5 grid divides the chip area.
- The `cellSpreading` function translates the concept of bin utilization as described in the FastPlace paper. Initially, the chip area is divided into uniform bins, each representing a portion of the placement region. The function calculates the utilization of each bin based on the number of cells within it. Cells are shifted from over-utilized bins to neighboring less-utilized bins while maintaining relative ordering.
- To prevent cells from reverting to their previous positions during subsequent iterations, spreading forces are calculated and applied. These forces are modeled as pseudo-nets connecting cells to virtual pins along the chip boundary, as described in the FastPlace methodology.

  ![image](https://github.com/user-attachments/assets/e4fac33b-3ceb-47ab-bf1c-6cd84c10cc4c)

- The function also accounts for movement control parameters, ensuring that cells do not move excessively between iterations, thereby stabilizing the placement process.
- Cells are redistributed based on bin utilization, ensuring an even distribution.

  ![image](https://github.com/user-attachments/assets/fb99ef33-6f42-4d84-8516-1f09f552cb65)

### Wirelength Optimization
- The wirelength is computed using:
  - Pre-optimization coordinates.
  - Post-optimization coordinates, after applying cell spreading.

---

## Workflow

### Input
1. Benchmark files (.net, .are, .kiaPad) placed in the working directory.
2. Command-line input specifying the benchmark name (without extension).

### Execution
1. Compile the program using the provided `Makefile`:
   ```bash
   make
   ```
2. Run the program:
   ```bash
   ./parser <benchmark_name>
   ```
3. Generate plots using:
   ```bash
   python3 plot_placement.py coordinates_before_spreading.kiaPad placement_before.png
   python3 plot_placement.py coordinates_after_spreading.kiaPad placement_after.png
   ```

### Output
- **Wirelength** before and after spreading (console).
- **Coordinate files**:
  - `coordinates_before_spreading.kiaPad`
  - `coordinates_after_spreading.kiaPad`
- **Placement Plots**:
  - `placement_before.png`
  - `placement_after.png`

---

## Challenges and AI Contributions

### Assistance from ChatGPT
- **Matrix Operations**: ChatGPT provided the initial implementation for matrix multiplication and conjugate gradient methods, ensuring compatibility with sparse, symmetric matrices as required by the FastPlace methodology.
- **Hybrid Net Model**: ChatGPT helped automate the construction of the connectivity matrix (`matQ`) using both the clique and star models, reducing non-zero entries to improve computational efficiency.
- **Iterative Local Refinement**: Suggestions from ChatGPT inspired the design of localized wirelength minimization using heuristic adjustments, as described in the FastPlace paper. This significantly reduced final wirelength after placement iterations.
- **Visualization**: Python script for plotting placements before and after cell spreading was AI-generated, offering seamless integration with the placement pipeline and clear graphical outputs for analysis.
- **Matrix Operations**: ChatGPT provided the initial implementation for matrix multiplication and conjugate gradient methods.
- **Visualization**: Python script for plotting was AI-generated.

### Bugs and Resolutions
1. **Matrix Initialization Issues:**
   - AI code failed to handle large matrices efficiently.
   - Solution: Optimized memory allocation and utilized sparse matrix techniques.
2. **File Parsing Errors:**
   - Misinterpretation of `.kiaPad` indexing caused incorrect pad placements.
   - Solution: Manual debugging and adding unit tests.

---

## Future Enhancements
1. Integrate a dynamic GUI for real-time placement adjustments.
2. Extend the algorithm to handle non-uniform cell sizes.
3. Automate benchmarking and comparison with other placement tools.

Example Output for toy01
Upon compiling and running the program with the toy01 benchmark, the following results are produced:

Compilation Commands:
g++ -g -c suraj_parser.cpp
g++ -g -c main.cpp
g++ -g -o placer suraj_parser.o main.o

Execution Command:
./placer toy01  # Replace 'toy01' with your input circuit file prefix
Execution Output:

Reading circuit file toy01
numCellPins, numhyper, numCellsAndPads, numCells_noPads = 15, 6, 6, 3
matQ:
4.25    -1      0       -1.25
-1      4.25    0       -1.25
0       0       1.25    -1.25
-1.25   -1.25   -1.25   6.25

dx Vector:
Node 0: -2
Node 1: -7
Node 2: 0
Node 3: -6.25

dy Vector:
Node 0: -5
Node 1: -5
Node 2: 0
Node 3: -12.5

Converged in 4 iterations.
x is: -1.83333
x is: -2.78571
x is: -2.40476
x is: -2.40476
Converged in 3 iterations.
y is: -3.09524
y is: -3.09524
y is: -4.04762
y is: -4.04762

Converged in 4 iterations.
Converged in 3 iterations.
Total wirelength is: 1.94015

Coordinates written to coordinates_before_spreading.kiaPad
Cell coordinates before spreading:
Cell 0: (-1.83333, -3.09524)
Cell 1: (-2.78571, -3.09524)
Cell 2: (-2.40476, -4.04762)
Cell 3: (-2.40476, -4.04762)

Coordinates written to coordinates_after_spreading.kiaPad
Cell coordinates after spreading:
Cell 0: (-1.88095, -3.14286)
Cell 1: (-2.7381, -3.14286)
Cell 2: (-2.39524, -4)
Cell 3: (-2.39524, -4)

Wirelength after spreading: 3.2606
Visualization Commands: The generated coordinate files can be visualized using the following commands:

python3 plot_placement.py coordinates_before_spreading.kiaPad placement_before.png
python3 plot_placement.py coordinates_after_spreading.kiaPad placement_after.png

![image](https://github.com/user-attachments/assets/c9b2a9a8-2212-456d-8eee-71bf99e83d46)

---

## References
1. FastPlace paper: "Efficient Analytical Placement using Cell Shifting, Iterative Local Refinement, and a Hybrid Net Model" by Viswanathan & Chu.
2. EE5301 Project Assignment Details.

---

### Acknowledgments
Special thanks to Prof. Bazargan for project guidance and the insightful FastPlace methodology. Generative AI tools like ChatGPT played a crucial role in coding and debugging.


