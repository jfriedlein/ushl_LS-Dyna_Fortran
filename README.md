# ushl_LS-Dyna_Fortran
Basics to implement user-defined shell elements (ushl, uel) in LS-Dyna with Fortran 

Focus here is on 2D elements for plane strain and axisymmetry, not on actual shell elements

## Generalised interface for use of separate element subroutines for 2D plane strain and axisymmetry
If you use a separate subroutine (e.g. a function that computes the force vector and stiffness matrix based on the displacement vector) for the element formulation or e.g. AceGen to generate the element routine, I can recommend the general interface stated in "ushl_e101_generalInterface.f".

## Numerical examples
The folder "numericalExamples_LS-Dyna" contains examples to test user-shell elements for plane strain and axisymmetry.
The folder "userLoading_LS-Dyna" contains the modified subroutine loadud to apply constant pressure on an edge of a shell element to enable plane strain or axisymmetric pressure loading of user-shell elements.


## todo
- "UEL_helper_Fortran_LS-Dyna" is included as submodule in git, but is therefore not downloaded when using "Code"->"Download ZIP". This is a common (and annoying) limitation/bug in git (status 2024). Please manually download the submodule or clone the entire repo (including submodules)
- function "isNan(*)" is only available for ifort compiler, e.g. not for pgi. Maybe use ( a /= 1 ) or ( a == NaN )?
