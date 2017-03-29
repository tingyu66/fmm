#include "build_tree.h"
#include "kernel.h"
#include "timer.h"
#include "traversal.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  const int numBodies = 1000;                                   // Number of bodies
  P = 10;                                                       // Order of expansions
  ncrit = 64;                                                   // Number of bodies per leaf cell
  theta = 0.4;                                                  // Multipole acceptance criterion

  printf("--- %-16s ------------\n", "FMM Profiling");          // Start profiling
  //! Initialize bodies
  start("Initialize bodies");                                   // Start timer
  Bodies bodies(numBodies);                                     // Initialize bodies
  real_t average = 0;                                           // Average charge
  srand48(0);                                                   // Set seed for random number generator
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies
    for (int d=0; d<3; d++) {                                   //  Loop over dimension
      bodies[b].X[d] = drand48() * 2 * M_PI - M_PI;             //   Initialize positions
    }                                                           //  End loop over dimension
    bodies[b].q = drand48() - .5;                               //  Initialize charge
    average += bodies[b].q;                                     //  Accumulate charge
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  average /= bodies.size();                                     // Average charge
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies
    bodies[b].q -= average;                                     // Charge neutral
  }                                                             // End loop over bodies
  stop("Initialize bodies");                                    // Stop timer

  //! Build tree
  start("Build tree");                                          // Start timer
  Cell * cells = buildTree(bodies);                             // Build tree
  stop("Build tree");                                           // Stop timer

  //! FMM evaluation
  start("Upward pass");                                         // Start timer
  initKernel();                                                 // Initialize kernel
  upwardPass(cells);                                            // Upward pass for P2M, M2M
  stop("Upward pass");                                          // Stop timer
  start("Traversal");                                           // Start timer
  traversal(cells, cells);                                      // Traversal for M2L, P2P
  stop("Traversal");                                            // Stop timer
  start("Downward pass");                                       // Start timer
  downwardPass(cells);                                          // Downward pass for L2L, L2P
  stop("Downward pass");                                        // Stop timer

  //! Direct N-Body
  start("Direct N-Body");                                       // Start timer
  const int numTargets = 10;                                    // Number of targets for checking answer
  Bodies jbodies = bodies;                                      // Save bodies in jbodies
  int stride = bodies.size() / numTargets;                      // Stride of sampling
  for (int b=0; b<numTargets; b++) {                            // Loop over target samples
    bodies[b] = bodies[b*stride];                               //  Sample targets
  }                                                             // End loop over target samples
  bodies.resize(numTargets);                                    // Resize bodies
  Bodies bodies2 = bodies;                                      // Backup bodies
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies
    bodies[b].p = 0;                                            //  Clear potential
    for (int d=0; d<3; d++) bodies[b].F[d] = 0;                 //  Clear force
  }                                                             // End loop over bodies
  direct(bodies, jbodies);                                      // Direct N-Body
  stop("Direct N-Body");                                        // Stop timer

  //! Verify result
  real_t pSum = 0, pSum2 = 0, FDif = 0, FNrm = 0;
  for (int b=0; b<int(bodies.size()); b++) {                    // Loop over bodies & bodies2
    pSum += bodies[b].p * bodies[b].q;                          // Sum of potential for bodies
    pSum2 += bodies2[b].p * bodies2[b].q;                       // Sum of potential for bodies2
    FDif += (bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +// Difference of force
      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]) +// Difference of force
      (bodies[b].F[2] - bodies2[b].F[2]) * (bodies[b].F[2] - bodies2[b].F[2]);// Difference of force
    FNrm += bodies[b].F[0] * bodies[b].F[0] + bodies[b].F[1] * bodies[b].F[1] +// Value of force
      bodies[b].F[2] * bodies[b].F[2];
  }                                                             // End loop over bodies & bodies2
  real_t pDif = (pSum - pSum2) * (pSum - pSum2);                // Difference in sum
  real_t pNrm = pSum * pSum;                                    // Norm of the sum
  printf("--- %-16s ------------\n", "FMM vs. direct");         // Print message
  printf("%-20s : %8.5e s\n","Rel. L2 Error (p)", sqrt(pDif/pNrm));// Print potential error
  printf("%-20s : %8.5e s\n","Rel. L2 Error (F)", sqrt(FDif/FNrm));// Print force error
  return 0;
}
