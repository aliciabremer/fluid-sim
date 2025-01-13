#ifndef THREADED_FLUID_SOLVER_H
#define THREADED_FLUID_SOLVER_H
#include <vector>
#include <array>
#include <utility>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <boost/asio.hpp>
#include <boost/asio/thread_pool.hpp>
#include "fluidSolver.h"
using namespace std;

namespace fluid {

class ThreadedFluidSolver: public FluidSolver {

    boost::asio::thread_pool pool;
    size_t num_threads;

  public:
    ThreadedFluidSolver(int width, int height, int dx, float viscosity, size_t num_threads);
    ~ThreadedFluidSolver();
    void step(double dt_ms) override;
    void linearSolverJacobi(std::vector<std::array<float, DIM>> &arr,
                                   std::vector<std::array<float, DIM>> &prevArr,
                                   std::vector<std::array<float, DIM>> &tmp,
                                   float a, float denom, float scale, size_t dims) override;
    
};

}

#endif

