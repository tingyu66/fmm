#ifndef types_h
#define types_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>

//! Namespace of exafmm
namespace exafmm {
  // Basic type definitions
  typedef double real_t;                                        //!< Floating point type
  typedef std::complex<real_t> complex_t;                       //!< Complex type

  //! Structure of bodies
  struct Body {
    real_t X[3];                                                //!< Position
    real_t q;                                                   //!< Charge
    real_t p;                                                   //!< Potential
    real_t F[3];                                                //!< Force
  };
  typedef std::vector<Body> Bodies;                             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int NCHILD;                                                 //!< Number of child cells
    int NBODY;                                                  //!< Number of descendant bodies
    Cell * CHILD;                                               //!< Pointer of first child cell
    Body * BODY;                                                //!< Pointer of first body
    real_t X[3];                                                //!< Cell center
    real_t R;                                                   //!< Cell radius
    std::vector<complex_t> M;                                   //!< Multipole expansion coefs
    std::vector<complex_t> L;                                   //!< Local expansion coefs
  };
  typedef std::vector<Cell> Cells;                              //!< Vector of cells
}
#endif
