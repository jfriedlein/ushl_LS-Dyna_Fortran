# ushl_LS-Dyna_Fortran
Basics to implement user-defined shell elements (ushl, uel) in LS-Dyna with Fortran 

Focus here is on 2D elements for plane strain and axisymmetry, not on actual shell elements

## Generalised interface for use of separate element subroutines for 2D plane strain and axisymmetry
If you use a separate subroutine for the element formulation or e.g. AceGen to generate the element routine, I can recommend the general interface stated in "ushl_e101_generalInterface.f".

