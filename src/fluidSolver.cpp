#include <fluidSolver.h>
#include <iostream>
#include <boost/asio.hpp>
using namespace std;

fluid::FluidSolver::FluidSolver(int width, int height, int dx, float viscosity):
    width{width}, height{height}, dx{dx},
    numX{size_t(width/dx)}, numY{size_t(height/dx)}, viscosity{viscosity},
    tmp((width/dx)*(height/dx), {0}),
    pressure((width/dx)*(height/dx), {0}),
    divergence((width/dx)*(height/dx), {0}),
    velocity((width/dx)*(height/dx), {0}),
    prevVelocity((width/dx )* (height/dx), {0}),
    externalForce((width/dx )* (height/dx), {0}) {
    
    initGL();
}

void fluid::FluidSolver::step(double dt_ms) {
    std::swap(velocity, prevVelocity);
    
    addAllForce(0, numX, dt_ms);

    transportAllFluid(0, numX, dt_ms);
    
    diffuse(dt_ms);
    project();
    
}

void fluid::FluidSolver::initGL() {

    // following:
    // https://learnopengl.com/Advanced-OpenGL/Instancing

    // shaders
    const char* vertexSource = R"glsl(
        #version 330 core

        layout (location = 0) in vec2 aPos;
        layout (location = 1) in vec3 aColor;
        layout (location = 2) in vec2 aOffset;

        out vec3 fcolor;

        void main()
        {
            fcolor = aColor;
            gl_Position = vec4(aPos + aOffset, 0.0, 1.0);
        }
    )glsl";
    
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);

    GLint status;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
    if (status != GL_TRUE) {
        cerr << "Error compiling vertex shader" << endl;
    }

    const char* fragmentSource = R"glsl(
        #version 330 core

        in vec3 fcolor;

        out vec4 outColor;

        void main()
        {
            outColor = vec4(fcolor, 1.0);
        }
    )glsl";
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &status);
    if (status != GL_TRUE) {
        cerr << "Error compiling fragment shader" << endl;
    }

    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);

    glLinkProgram(shaderProgram);

    // instance arrays
    float offsetX = 2.f/numX;
    float offsetY = 2.f/numY;

    for (size_t i = 0; i < numX; i++) {
        for (size_t j = 0; j < numY; j++) {
            glm::vec2 t;
            t.x = i*offsetX - 1.f;
            t.y = j*offsetY - 1.f;
            translations.push_back(t);

            glm::vec3 c;
            c.r = 0.f;
            c.b = 0.f;
            c.g = 0.f;
            colors.push_back(c);
        }
    }
    
    // offset instance arrays
    glGenBuffers(1, &offsetVBO);
    glBindBuffer(GL_ARRAY_BUFFER, offsetVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2) * translations.size(), &translations[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // color instance arrays
    glGenBuffers(1, &colorVBO);
    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3) * colors.size(), &colors[0], GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // make vertices
    float vertices[] = {
        // positions    
        0.f,  0.,
        0.f, offsetY,
        offsetX, 0.f,

        offsetX,  0.f,
        offsetX,  offsetY,
        0.f, offsetY
    };

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
      
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    // position attribute
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);

    // also set instance data
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // also set other instance data
    glEnableVertexAttribArray(2);
    glBindBuffer(GL_ARRAY_BUFFER, offsetVBO);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // instanced attribute
    glVertexAttribDivisor(1, 1);
    glVertexAttribDivisor(2, 1);    

}

fluid::FluidSolver::~FluidSolver() {
    glDeleteProgram(shaderProgram);
    glDeleteShader(fragmentShader);
    glDeleteShader(vertexShader);

    glDeleteBuffers(1, &vbo);

    glDeleteVertexArrays(1, &vao);
}

void fluid::FluidSolver::render() {
    // clear stuff on the screen
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // #pragma omp parallel for
    for (size_t i = 0; i < numX; i++) {
        for (size_t j = 0; j < numY; j++) {
            size_t ind = getCoord(i,j);

            glm::vec3 c;
            c.r = 0.f;
            c.b = pressure[ind][0]*10000;
            c.g = 0.f;
            colors[ind] = c;
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec3) * colors.size(), &colors[0]);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glUseProgram(shaderProgram);
    glBindVertexArray(vao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, numX*numY);
    glBindVertexArray(0);

}

void fluid::FluidSolver::setExternalForce(float x_pos, float y_pos, std::vector<float> f) {

    size_t ind = getCoord(size_t(x_pos/dx), numY-1-size_t(y_pos/dx));
    for (int d = 0; d < DIM; d++) {
        externalForce[ind][d] = f[d];
    }
}

inline size_t fluid::FluidSolver::getCoord(size_t x, size_t y) {
    return x*numX + y;
}

void fluid::FluidSolver::addAllForce(size_t x_begin, size_t x_end, float dt) {
    for (size_t i = x_begin; i < x_end; i++) {
        for (size_t j = 0; j < this->numY; j++) {
            this->addForce(i, j, dt);
        }
    }
}

void fluid::FluidSolver::addForce(size_t x, size_t y, float dt) {
    size_t ind = this->getCoord(x, y);
    for (size_t d = 0; d < this->DIM; d++) {
        this->prevVelocity[ind][d] += dt*this->externalForce[ind][d];
    }
}

void fluid::FluidSolver::transportAllFluid(size_t x_begin, size_t x_end, float dt) {
    for (size_t i = x_begin; i < x_end; i++) {
        for (size_t j = 0; j < numY; j++) {
            transportFluid(i, j, dt);
        }
    }
}

void fluid::FluidSolver::transportFluid(size_t x, size_t y, float dt) {
    size_t ind = getCoord(x, y);

    std::pair<float, float> new_x0 = traceParticle(x, y, -dt);
    for (size_t d = 0; d < DIM; d++) {
        velocity[ind][d] = linearInterp(new_x0, d);
    }
}

void fluid::FluidSolver::diffuse(float dt) {
    float a = dt*viscosity*numX*numY;
    linearSolverJacobi(velocity, prevVelocity, tmp, a, 1.f+4*a, -1, DIM);
}

void fluid::FluidSolver::project() {
    for (int i = 1; i < numX-1; i++) {
        for (int j = 1; j < numY-1; j++) {
            pressure[getCoord(i,j)][0] = 0;
            divergence[getCoord(i,j)][0] = -0.5*(velocity[getCoord(i+1,j)][0] - velocity[getCoord(i-1,j)][0] 
                                              + velocity[getCoord(i,j+1)][1] - velocity[getCoord(i, j-1)][1]) / numX;
        }
    }

    updateBoundary(divergence, 1);
    updateBoundary(pressure, 1);
    linearSolverJacobi(pressure, divergence, tmp, 1, 4, 1, 1);

    for (int i = 1; i < numX-1; i++) {
        for (int j = 1; j < numY-1; j++) {
            size_t ind = getCoord(i,j);
            velocity[ind][0] = velocity[ind][0] - 0.5*(pressure[getCoord(i+1, j)][0] - pressure[getCoord(i-1, j)][0])*numX;
            velocity[ind][1] = velocity[ind][1] - 0.5*(pressure[getCoord(i, j+1)][0] - pressure[getCoord(i, j-1)][0])*numY;
        }
    }
    

    updateBoundary(velocity, -1);

}

void fluid::FluidSolver::updateBoundary(std::vector<std::array<float, DIM>> &arr, float s) {
    for (int i = 1; i < numX-1; i++) {
        for (int d= 0; d < DIM; d++) {
            arr[getCoord(i, 0)][d] =  s*arr[getCoord(i, 1)][d];
            arr[getCoord(i, numY-1)][d] = s*arr[getCoord(i, numY-2)][d];
        }
    }

    for (int j = 1; j < numY-1; j++) {
        for (int d= 0; d < DIM; d++) {
            arr[getCoord(0, j)][d] =  s*arr[getCoord(1, j)][d];
            arr[getCoord(numX-1, j)][d] = s*arr[getCoord(numX-2, j)][d];
        }
    }

    for (int d = 0; d < DIM; d++) {
        arr[getCoord(0,0)][d] = (arr[getCoord(0, 1)][d] + arr[getCoord(1, 1)][d]) / 2.f;
        arr[getCoord(numX-1,0)][d] = (arr[getCoord(numX-1, 1)][d] + arr[getCoord(numX-2, 0)][d]) / 2.f;
        arr[getCoord(0, numY-1)][d] = (arr[getCoord(1, numY-1)][d] + arr[getCoord(0, numY-2)][d]) / 2.f;
        arr[getCoord(numX-1, numY-1)][d] = (arr[getCoord(numX-2, numY-1)][d] + arr[getCoord(numX-1, numY-2)][d]) / 2.f;
    }
}

std::pair<float, float> fluid::FluidSolver::traceParticle(size_t x, size_t y, float timestep) {
    size_t ind = getCoord(x, y);

    // use second order Runge-Kutta

    // x and y but in positions on grid
    // these are yn of equation
    float x_dist = (float(x) + 0.5) * dx;
    float y_dist = (float(y) + 0.5) * dx;

    // compute new position on grid
    float new_x = x_dist + 0.5*timestep*prevVelocity[ind][0];
    float new_y = y_dist + 0.5*timestep*prevVelocity[ind][1];

    // get positions on grid but in coordinates
    size_t new_x_ind = size_t(new_x/dx);
    size_t new_y_ind = size_t(new_y/dx);

    // check coords
    if (new_x_ind < 0) {
        new_x_ind = 0;
    }
    if (new_y_ind < 0) {
        new_y_ind = 0;
    }
    if (new_x_ind > numX-1) {
        new_x_ind = numX-1;
    }
    if (new_y_ind > numY-1) {
        new_y_ind = numY-1;
    }

    size_t new_ind = getCoord(new_x_ind, new_y_ind);

    // k2
    return std::pair<float, float>{
        x_dist + timestep*prevVelocity[new_ind][0],
        y_dist + timestep*prevVelocity[new_ind][1],
    };
}

float fluid::FluidSolver::linearInterp(std::pair<float, float> pos, size_t dim) {
    // get new indices
    std::pair<size_t, size_t> ind1 = make_pair(
        size_t(pos.first/(float)dx - 0.5),
        size_t(pos.second/(float)dx - 0.5)
    );
    if (ind1.first < 0) {
        ind1.first = 0;
    }
    if (ind1.second < 0) {
        ind1.second = 0;
    }
    if (ind1.first > numX-2) {
        ind1.first = numX-2;
    }
    if (ind1.second > numY-2) {
        ind1.second = numY-2;
    }
    std::pair<size_t, size_t> ind2 = make_pair(ind1.first, ind1.second+1);
    std::pair<size_t, size_t> ind3 = make_pair(ind1.first+1, ind1.second);
    std::pair<size_t, size_t> ind4 = make_pair(ind1.first+1, ind1.second+1);
    
    // linearly interpolate in the x-dir
    float v_1 = (
            ((float(ind3.first)+0.5)*dx - pos.first)*prevVelocity[getCoord(ind1.first, ind1.second)][dim] +
            (pos.first - (float(ind1.first)+0.5)*dx)*prevVelocity[getCoord(ind3.first, ind3.second)][dim]
        ) / dx;
    
    float v_2 = (
            ((float(ind4.first)+0.5)*dx - pos.first)*prevVelocity[getCoord(ind2.first, ind2.second)][dim] +
            (pos.first - (float(ind2.first)+0.5)*dx)*prevVelocity[getCoord(ind4.first, ind4.second)][dim]
        ) / dx;
    
    return (((float(ind2.second)+0.5)*dx - pos.second)*v_1 + (pos.second - (float(ind1.second)+0.5)*dx)*v_2) / dx;
}

void fluid::FluidSolver::linearSolverJacobi(std::vector<std::array<float, DIM>> &arr,
                                   std::vector<std::array<float, DIM>> &prevArr,
                                   std::vector<std::array<float, DIM>> &tmp,
                                   float a, float denom, float scale, size_t dims)
{
    for (int k = 0; k < ITER; k++) {
        iterationsLinearSolverJacobi(1, numX-1, arr, prevArr, tmp,
                                    a, denom, scale, dims);
        
        std::swap(tmp, arr);
        updateBoundary(arr, scale);
    }
}

void fluid::FluidSolver::iterationsLinearSolverJacobi(size_t x_begin, size_t x_end,
                                   std::vector<std::array<float, DIM>> &arr,
                                   std::vector<std::array<float, DIM>> &prevArr,
                                   std::vector<std::array<float, DIM>> &tmp,
                                   float a, float denom, float scale, size_t dims) {
    for (int i = x_begin; i < x_end; i++) {
        for (int j = 1; j < numY-1; j++) {
            size_t ind = getCoord(i,j);
            for (int d = 0; d < dims; d++) {
                tmp[ind][d] = (prevArr[ind][d] + a*(arr[getCoord(i+1,j)][d]
                                                +  arr[getCoord(i-1, j)][d]
                                                +  arr[getCoord(i, j+1)][d] 
                                                +  arr[getCoord(i, j-1)][d]))/denom;
            }
        }
    }
}
