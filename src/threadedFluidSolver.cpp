#include <threadedFluidSolver.h>
#include <cmath>
#include <iostream>
using namespace std;

fluid::ThreadedFluidSolver::ThreadedFluidSolver(int width, int height, int dx, float viscosity, size_t num_threads):
    FluidSolver{width, height, dx, viscosity}, num_threads{num_threads}, pool{num_threads} {
        
}

fluid::ThreadedFluidSolver::~ThreadedFluidSolver() {
}

void fluid::ThreadedFluidSolver::step(double dt_ms) {
    std::swap(velocity, prevVelocity);

    size_t div = numX / num_threads + 1;
    for (size_t n = 0; n < num_threads; n++) {
        boost::asio::post(pool,
        [this, dt_ms, n, div]() mutable
        {
            addAllForce(div*n, min(numX, div*(n+1)), dt_ms);
        });
    }

    pool.wait();

    for (size_t n = 0; n < num_threads; n++) {
        boost::asio::post(pool,
        [this, dt_ms, n, div]() mutable
        {
            transportAllFluid(div*n, min(numX, div*(n+1)), dt_ms);
        });
    }

    pool.wait(); 
    
    diffuse(dt_ms);
    project();

}

void fluid::ThreadedFluidSolver::linearSolverJacobi(std::vector<std::array<float, DIM>> &arr,
                                   std::vector<std::array<float, DIM>> &prevArr,
                                   std::vector<std::array<float, DIM>> &tmp,
                                   float a, float denom, float scale, size_t dims)
{
    for (int k = 0; k < ITER; k++) {
        size_t div = (numX-2) / num_threads + 1;
        for (size_t n = 0; n < num_threads; n++) {
            boost::asio::post(pool,
            [this, n, &arr, div, &prevArr, &tmp, a, denom, scale, dims]() mutable
            {
                iterationsLinearSolverJacobi(1+div*n, min(numX-1, div*(n+1)), arr,
                                            prevArr, tmp, a, denom, scale, dims);
            });
        }

        pool.wait();
        
        std::swap(tmp, arr);
        updateBoundary(arr, scale);
    }
}