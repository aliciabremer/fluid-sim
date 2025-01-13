#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H
#include <vector>
#include <array>
#include <utility>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
using namespace std;

namespace fluid {

class FluidSolver {
    int width;
    int height;
    int dx;

    GLuint vbo, vao, offsetVBO, colorVBO, vertexShader, fragmentShader, shaderProgram, posAttrib, colAttrib;
    std::vector<glm::vec2> translations;
    std::vector<glm::vec3> colors;


    std::pair<float, float> traceParticle(size_t x, size_t y, float timestep);
    float linearInterp(std::pair<float, float> pos, size_t dim);

  public:
    static const int DIM = 2;
    static const int ITER = 10;

    FluidSolver(int width, int height, int dx, float viscosity);
    virtual ~FluidSolver();
    virtual void step(double dt_ms);
    virtual void initGL();
    virtual void render();
    virtual void linearSolverJacobi(std::vector<std::array<float, DIM>> &arr,
                                   std::vector<std::array<float, DIM>> &prevArr,
                                   std::vector<std::array<float, DIM>> &tmp,
                                   float a, float denom, float scale, size_t dims);
    void setExternalForce(float x_pos, float y_pos, std::vector<float> f);
  
  protected:
    size_t numX;
    size_t numY;
    float viscosity;
    std::vector<std::array<float, DIM>> tmp;
    std::vector<std::array<float, DIM>> pressure;
    std::vector<std::array<float, DIM>> divergence;
    std::vector<std::array<float, DIM>> velocity;
    std::vector<std::array<float, DIM>> prevVelocity;
    std::vector<std::array<float, DIM>> externalForce;

    size_t getCoord(size_t x, size_t y);
    void addAllForce(size_t x_begin, size_t x_end, float dt);
    void addForce(size_t x, size_t y, float dt);
    void transportAllFluid(size_t x_begin, size_t x_end, float dt);
    void transportFluid(size_t x, size_t y, float dt);
    void diffuse(float dt);
    void project();
    void updateBoundary(std::vector<std::array<float, DIM>> &arr, float s);

    void iterationsLinearSolverJacobi(size_t x_begin, size_t x_end,
                                   std::vector<std::array<float, DIM>> &arr,
                                   std::vector<std::array<float, DIM>> &prevArr,
                                   std::vector<std::array<float, DIM>> &tmp,
                                   float a, float denom, float scale, size_t dims);
};

}

#endif
