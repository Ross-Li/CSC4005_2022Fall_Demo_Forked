#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

// Declare variables for number of bodies and number of iteration
int n_body;
int n_iteration;

///
/// @brief Function to generate data for every body
/// @param[in,out] m,x,y,vx,vy Pointers pointing to the start of arrays
/// @param[in] n_body The number of bodies
/// 
void generate_data(double *m, double *x,double *y, double *vx, double *vy, int n_body) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n_body; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}


/// 
/// @brief Function to update the positions of every body for 1 iteration
/// @param[in,out] x,y,vx,vy Pointers to the start of the arrays 
/// @param[in,out] n_body The number of bodies
///
void update_position(double *x, double *y, double *vx, double *vy, int n_body) {
    //TODO: update position 
    for (int i = 0; i < n_body; i++) {
        x[i] = x[i] + vx[i] * dt;
        y[i] = y[i] + vy[i] * dt;
    }
    //* If the body hit the wall...
    //* Assume that they just go out of the picture
}

///
/// @brief Function to update the velocities of every body for 1 iteration
/// @param[in] m,x,y Pointers to the start of the arrays of  
/// @param[in,out] vx,vy Pointers to the head of arrays of 
/// @param[in] n_body The number of bodies
/// 
void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n_body) {
    //TODO: calculate force and acceleration, update velocity
    // For each body i...
    for (int i = 0; i < n_body; i++) {
        // For every body j other than the current body i...
        for (int j = 0; j < n_body; j++) {
            if (i == j)
                continue;

            double deltaX = x[j] - x[i];
            double deltaY = y[j] - y[i];
            double distance = sqrt(deltaX * deltaX + deltaY * deltaY); 

            //* If 2 bodies collide, we just assmue that they go on as their original directions
            //* Hense we ignore the following computations because the  
            if (distance == 0) 
                continue;

            // Calculate the acceleration along the X axis and Y axis of body i under the force of body j
            double accX = gravity_const * m[j] * deltaX / pow(distance, 3);
            double accY = gravity_const * m[j] * deltaY / pow(distance, 3);

            // Update the velocity of body i with the influence of body j
            vx[i] = vx[i] + accX * dt;
            vy[i] = vy[i] + accY * dt;            
        }
    }
}


void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("seq", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        update_velocity(m, x, y, vx, vy, n_body);
        update_position(x, y, vx, vy, n_body);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        
        //! Comment this if you don't want your terminal to be filled
        // printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(x, y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    
}


int main(int argc, char *argv[]){
    
    // Get input from terminal
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

    printf("Student ID: 119020026\n"); // replace it with your student id
    printf("Name: Li Peilin\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation Sequential Implementation\n");
    
    return 0;

}


