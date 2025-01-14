#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <thread>
#include <iostream>
#include "fluidSolver.h"
#include "threadedFluidSolver.h"
using namespace std;

int main(int argc, char* argv[])
{

    int width = 1000;
    int height = 1000;
    int dx = 10;
    float viscosity = 10.f;
    size_t num_threads = 8;

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    GLFWwindow *window = glfwCreateWindow(width, height, "OpenGL", nullptr, nullptr); // Windowed
    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    glewInit();

    glEnable(GL_DEPTH_TEST);

    // GLuint vertexBuffer;
    // glGenBuffers(1, &vertexBuffer);

    // printf("%u\n", vertexBuffer);

    fluid::FluidSolver *solver = new fluid::FluidSolver(width, height, dx, viscosity);
    // fluid::FluidSolver *solver = new fluid::ThreadedFluidSolver(width, height, dx, viscosity, num_threads);
    double dt_ms = 0;
    bool setExternalForce = false;
    double xposition, yposition;
    double prev_xposition, prev_yposition;


    while(!glfwWindowShouldClose(window))
    {
        // get start time
        auto t1 = std::chrono::high_resolution_clock::now();

        if (setExternalForce) {
            solver->setExternalForce(xposition, yposition, {0.1f*float(xposition-prev_xposition), -0.1f*float(yposition-prev_yposition)});
        }

        // update for previous
        solver->step(dt_ms);
        solver->render();

        glfwSwapBuffers(window);
        glfwPollEvents();

        if (setExternalForce) {
            solver->setExternalForce(xposition, yposition, {0.f, 0.f});
        }

        

        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GL_TRUE);

        if(glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            prev_xposition = xposition;
            prev_yposition = yposition;
            glfwGetCursorPos(window, &xposition, &yposition);
            if (!setExternalForce) {
                prev_xposition = xposition-1;
                prev_yposition = yposition-1;
            }
            setExternalForce = true;
        }
        else {
            setExternalForce = false;
        }
        
        // get end time
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> diff = t2 - t1;
        dt_ms = diff.count();
    }

    delete(solver);
    glfwTerminate();
}