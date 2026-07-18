# RAMSES Eigenanalysis Tool

A comprehensive MATLAB-based tool for performing eigenanalysis on power system dynamic data extracted from the RAMSES simulator, part of the [STEPSS](https://stepss.sps-lab.org/) power system simulation platform. This tool provides multiple numerical methods for computing eigenvalues and eigenvectors of power system models in descriptor form.

## Overview

The RAMSES Eigenanalysis tool is designed to analyze the small-signal stability of power systems by computing eigenvalues and eigenvectors of the system Jacobian matrix. It supports both differential-algebraic equation (DAE) systems and provides various computational methods optimized for different system sizes and requirements.

## Features

- **Multiple Analysis Methods**:
  - **QZ Method**: Standard eigenvalue decomposition using MATLAB's `eig()` function with algebraic variable elimination
  - **ARPACK Method**: Sparse descriptor approach using Krylov-Schur algorithm via MATLAB's `eigs()` function
  - **Arnoldi Method**: Iterative shift-invert eigenvalue computation for large systems (commented in code)
  - **JDQR Method**: Jacobi-Davidson QR method for targeted eigenvalue computation

- **Interactive Analysis**: Post-processing tools for examining results including:
  - Dominant eigenvalue identification and filtering
  - Mode shape analysis
  - Participation factor calculations
  - Damping factor computation and plotting
  - Interactive result exploration

- **Data Processing**: 
  - Automatic Jacobian matrix extraction from RAMSES simulation data
  - Support for descriptor form matrices (A, E)
  - Flexible input format handling

## Installation

**Requirements:** MATLAB R2016a or later (base MATLAB only — sparse matrices are supported natively and no additional toolboxes are required), and [PyRAMSES](https://stepss.sps-lab.org/pyramses/overview/) — the Python interface to the RAMSES simulator — for extracting Jacobian matrices from simulation data.

1. Clone or download this repository
2. Add the repository folder to your MATLAB path:
   ```matlab
   addpath('path/to/stepss-eigenanalysis')
   addpath('path/to/stepss-eigenanalysis/scripts')
   ```

## Usage

### Basic Usage

The main function `ssa()` performs the eigenanalysis:

```matlab
% Basic usage with default settings
ssa('jac_val.dat', 'jac_eqs.dat', 'jac_var.dat', 'jac_struc.dat')

% With custom parameters
ssa('jac_val.dat', 'jac_eqs.dat', 'jac_var.dat', 'jac_struc.dat', real_limit, damp_ratio, method)
```

### Function Parameters

- `jac_val.dat`: Matrix values in coordinate format
- `jac_eqs.dat`: Equation descriptions (differential/algebraic)
- `jac_var.dat`: Variable descriptions (differential/algebraic)  
- `jac_struc.dat`: Decomposed power system structure (optional)
- `real_limit`: Real part threshold for dominant eigenvalues (default: -∞)
- `damp_ratio`: Damping factor limit drawn on the eigenvalue plot (default: 1.0)
- `method`: Analysis method - 'QZ' or 'ARP' (default: 'QZ')

### Example Workflow

1. **Extract Jacobian Data**: Use PyRAMSES to generate the required data files
   ```python
   import pyramses
   ram = pyramses.sim()
   case = pyramses.cfg('cmd.txt')
   ram.execSim(case)
   ```

2. **Run Eigenanalysis**: Execute in MATLAB
   ```matlab
   ssa('jac_val.dat', 'jac_eqs.dat', 'jac_var.dat', 'jac_struc.dat')
   ```

3. **Interactive Analysis**: Use the interactive menu to explore results:
   - View dominant eigenvalues
   - Analyze mode shapes
   - Calculate participation factors
   - Examine system structure

### Analysis Methods

#### QZ Method (Default)
- Suitable for small to medium systems (< 50,000 states)
- Provides complete eigenvalue spectrum
- Supports mode shape and participation factor analysis
- Uses algebraic variable elimination

#### ARPACK Method
- Optimized for large sparse systems
- Uses descriptor form directly (no elimination)
- Limited post-processing capabilities
- Memory efficient for large systems

For very large systems (> 50,000 states), use the ARPACK method or uncomment the Arnoldi iterations in `ssa.m`.

## Output Files

The tool generates several output files:

- `modal_reduction.mat`: MATLAB workspace with all computed results
- `eigs.fig`: Figure showing eigenvalue plot
- Console output with timing and convergence information

The dynamic (state) Jacobian and the computed eigenvalues/eigenvectors are also assigned to the base MATLAB workspace.

## Example Files

The `example/` folder contains sample data files, a Jupyter notebook (`simply_load_and_run.ipynb`) for extracting the Jacobian data with PyRAMSES, and a `Readme.md` describing the complete workflow through to the MATLAB eigenanalysis.

## Important Notes

- The Jacobian extraction requires synchronous reference frame in RAMSES settings (`$OMEGA_REF SYN`)
- The dynamic Jacobian is automatically created in the workspace
- Mode shapes and participation factors are only available with the QZ method
- For very large systems, consider using the commented Arnoldi method in `ssa.m`; its sparse solver backend (`scripts/eigs_solver.m`) can optionally use external solvers (KLU, PARDISO) that require separately installed MEX interfaces

## File Structure

```
stepss-eigenanalysis/
├── ssa.m                      # Main eigenanalysis function
├── scripts/                   # Analysis functions
│   ├── init.m                 # Data initialization
│   ├── calc_Jdyn.m            # Dynamic (state) Jacobian construction
│   ├── eigenvals_eig.m        # QZ method implementation
│   ├── eigenvals_eig_descr.m  # ARPACK descriptor method implementation
│   ├── eigenvals_eigs.m       # Arnoldi (shift-invert) method implementation
│   ├── eigenvals_jdqr.m       # JDQR method implementation
│   ├── eigs_solver.m          # Sparse linear solver backend for the Arnoldi method
│   ├── analyze_results.m      # Result analysis and filtering
│   ├── loop_analysis.m        # Interactive result exploration
│   └── ...
├── example/                   # Sample data and workflows
├── LICENSE                    # Apache License 2.0
├── NOTICE                     # Notes on proprietary STEPSS components
└── README.md                  # This documentation
```

## Documentation

- [Eigenanalysis user guide](https://stepss.sps-lab.org/user-guide/eigenanalysis/) on the STEPSS documentation site
- In-source help: run `help ssa` in MATLAB, or see the function headers in `scripts/`
- [`example/Readme.md`](example/Readme.md): step-by-step example workflow

## License

The RAMSES Eigenanalysis Tool is distributed under the **Apache License 2.0** — see [LICENSE](LICENSE). Copyright © Petros Aristidou. The [NOTICE](NOTICE) file describes the licensing of proprietary components of the wider STEPSS suite, which are not included in this repository.

## Authors

Developed and maintained by the [Sustainable Power Systems Laboratory (SPS-L)](https://sps-lab.org/) at the Cyprus University of Technology, under the direction of Dr. Petros Aristidou.

The bundled Jacobi-Davidson QR implementation (`scripts/jdqr.m`) is by Gerard Sleijpen (Copyright 1998).

For questions or support, please contact info@sps-lab.org.

## Citation

If you use this tool in your research, please cite the relevant papers from the Sustainable Power Systems Lab and acknowledge the use of this eigenanalysis tool.
