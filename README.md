Dicontinuous Galerkin solver for a smooth periodic convection diffusion
type PDE with possibly time-dependent uniform velocity field.
The basis is a tensor product basis of normalized Legendre polynomials;
due to the orthogonality of Legendre polynomials and choice of normalization,
the mass matrix is the identity matrix.

This is intended to be run in a web browser and is compiled down
to WebAssembly from C++. 

Build dependencies: emscripten, eigen

