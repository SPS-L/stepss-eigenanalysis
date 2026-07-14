# RAMSES Eigenanalysis Tool

A comprehensive MATLAB-based tool for performing eigenanalysis on power system dynamic data extracted from the RAMSES simulator. This tool provides multiple numerical methods for computing eigenvalues and eigenvectors of power system models in descriptor form.

## Overview

The RAMSES Eigenanalysis tool is designed to analyze the small-signal stability of power systems by computing eigenvalues and eigenvectors of the system Jacobian matrix. It supports both differential-algebraic equation (DAE) systems and provides various computational methods optimized for different system sizes and requirements.

## Features

- **Multiple Analysis Methods**:
  - **QZ Method**: Standard eigenvalue decomposition using MATLAB's `eig()` function with algebraic variable elimination
  - **ARPACK Method**: Sparse descriptor approach using Krylov-Schur algorithm via MATLAB's `eigs()` function
  - **Arnoldi Method**: Iterative eigenvalue computation for large systems (commented in code)
  - **JDQR Method**: Jacobi-Davidson QR method for targeted eigenvalue computation

- **Interactive Analysis**: Post-processing tools for examining results including:
  - Dominant eigenvalue identification and filtering
  - Mode shape analysis
  - Participation factor calculations
  - Damping ratio analysis
  - Interactive result exploration

- **Data Processing**: 
  - Automatic Jacobian matrix extraction from RAMSES simulation data
  - Support for descriptor form matrices (A, E)
  - Flexible input format handling

## Installation and Setup

### Prerequisites

1. **MATLAB**: R2016a or later with the following toolboxes:
   - MATLAB Base
   - MATLAB Sparse Matrices support

2. **PyRAMSES**: The Python interface to RAMSES simulator
   - Install from: [stepss.sps-lab.org/pyramses](https://stepss.sps-lab.org/pyramses/)
   - Required for extracting Jacobian matrices from simulation data

### Installation

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
- `damp_ratio`: Damping ratio threshold (default: 1.0)
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

## Output Files

The tool generates several output files:

- `modal_reduction.mat`: MATLAB workspace with all computed results
- `eigs.fig`: Figure showing eigenvalue plot
- Console output with timing and convergence information

## Example Files

The `example/` folder contains sample data files and a Jupyter notebook (`simply_load_and_run.ipynb`) demonstrating the complete workflow from PyRAMSES data extraction to MATLAB eigenanalysis.

## System Requirements

- **Small Systems** (< 50,000 states): QZ method recommended
- **Large Systems** (> 50,000 states): ARPACK method or uncomment Arnoldi iterations in `ssa.m`
- **Memory**: Depends on system size; ARPACK method is more memory-efficient

## Important Notes

- The Jacobian extraction requires synchronous reference frame in RAMSES settings (`$OMEGA_REF SYN`)
- The dynamic Jacobian is automatically created in the workspace
- For very large systems, consider using the commented Arnoldi method in `ssa.m`

## File Structure

```
stepss-eigenanalysis/
├── ssa.m                    # Main eigenanalysis function
├── scripts/                 # Analysis functions
│   ├── init.m              # Data initialization
│   ├── eigenvals_eig.m     # QZ method implementation
│   ├── eigenvals_eigs.m    # Arnoldi method implementation
│   ├── eigenvals_eig_descr.m # ARPACK method implementation
│   ├── analyze_results.m   # Result analysis and filtering
│   ├── loop_analysis.m     # Interactive result exploration
│   └── ...
├── example/                 # Sample data and workflows
└── README.md               # This documentation
```

## Acknowledgments

We would like to express our gratitude to the [Sustainable Power Systems Lab (SPSL)](https://sps-lab.org) for their invaluable support and contributions to the development of this project. For questions or support, please contact info@sps-lab.org.

## License

See the [LICENSE](LICENSE) file for details.

## Citation

If you use this tool in your research, please cite the relevant papers from the Sustainable Power Systems Lab and acknowledge the use of this eigenanalysis tool.
