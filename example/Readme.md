# Procedure

This procedure assumes you have stepss (stepss.sps-lab.org) installed.

1. Load the python file simply_load_and_run.ipynb and execute the cell. This will load the data and extract the Jacobian in the form of 4 files ('jac_val.dat','jac_eqs.dat','jac_var.dat','jac_struc.dat').
2. Load MATLAB in the folder and run ssa('jac_val.dat','jac_eqs.dat','jac_var.dat','jac_struc.dat'). You will get the eigenvalues and their plot.
3. Use the interactive commands to get more information.

## Comments

- The command to extract the Jacobian works only if you use synchronous reference frame in the settings ($OMEGA_REF SYN ;).
- The dynamic Jacobian is created automatically in the workspace.
- You can limit the eigenvalues you see by giving more arguments to the ssa script. See ssa.m for documentation.
- This only works for small systems. In large ones (e.e., more than 50k states), you would need to run the ARNOLDI method in the ssa.m which is currently commented out.
