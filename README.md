# recon_code
Code to perform density field reconstruction given an object catalog
and two random files.  The goal of density field reconstruction is to
sharpen the baryon acoustic oscillation feature for measuring distances.
This code implements the standard, lowest-order algorithm presented in
Eisenstein et al. (2007) with periodic boundary conditions using a
multigrid relaxation technique with a full multigrid V-cycle based on
damped Jacobi iteration.
Inputs are a catalog of objects and two random catalogs which serve to
specify the selection function and act as the "shifted" fields. Output
are `shifted' versions of the object and second random catalogs.
