#include "asg2.h"
#include <stdio.h>
#include <mpi.h>


int rank;
int world_size;


void master() {
	//TODO: procedure run in master process
	// MPI_Scatter...
	// MPI_Gather...
	// the following code is not a necessary, please replace it with MPI implementation.
	
	Point* p = data;
	for (int index = 0; index < total_size; index++){
		compute(p);
		p++;
	}

	//TODO END

}


void slave() {
	//TODO: procedure run in slave process
	// MPI_Scatter...
	// MPI_Gather...

	//TODO END
}


int main(int argc, char *argv[]) {
	if ( argc == 4 ) {
		X_RESN = atoi(argv[1]);
		Y_RESN = atoi(argv[2]);
		max_iteration = atoi(argv[3]);
	} else {
		X_RESN = 1000;
		Y_RESN = 1000;
		max_iteration = 300;
	}

	if (rank == 0) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("MPI");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, X_RESN, 0, Y_RESN);
		glutDisplayFunc(plot);
	}

	/* computation part begin */
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (rank == 0) {
		initData();
		master();
	} else {
		slave();
	}

	MPI_Finalize();
	/* computation part end */

	if (rank == 0){
		glutMainLoop();
	}

	return 0;
}
