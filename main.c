#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "file_operations.h"
#include "graphics.h"
#include "quad.h"

void updateParticles(double delta_t, particle_t *particles, int N) {
   //Set constants
   double *forcex=(double*)malloc(N*sizeof(double));
   double *forcey=(double*)malloc(N*sizeof(double));
   const double G = 100.0/N;
   const double eps = 0.001;
   double abs_r;
   double r_x, r_y;
   double x;
   double y;
   double m_j;
   double m_i;
   double k;
   
   for(int i=0; i<N; i++){
      x = particles[i].x_pos;
      y = particles[i].y_pos;
      m_i = particles[i].mass;
      
      // For each particle i, calculate the sum of the forces acting on it
      // two for loops!!!
      for(int j=i+1; j<N; j++){
         
            m_j = particles[j].mass;
            
            // Calculate the distance betweem particles i and j.
            // USE SQRTF????
            abs_r = sqrt((x-particles[j].x_pos)*(x-particles[j].x_pos)+(y-particles[j].y_pos)*(y-particles[j].y_pos));
            r_x = x-particles[j].x_pos;
            r_y = y-particles[j].y_pos;
            // Plumber spheres
            // use dummy variable???
            k = -G*m_i*m_j/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
            forcex[i] += k*r_x;
            forcey[i] += k*r_y;
            forcex[j] += -k*r_x;
            forcey[j] += -k*r_y;
         
      }   
        
   }
   // Using the force, update the velocity and position.
   for(int i=0;i<N;i++){
      m_i = 1/particles[i].mass;
      particles[i].vel_x += delta_t*forcex[i]*m_i;
      particles[i].vel_y += delta_t*forcey[i]*m_i;
      particles[i].x_pos += delta_t*particles[i].vel_x;
      particles[i].y_pos += delta_t*particles[i].vel_y;
   }
   free(forcex);
   free(forcey);
}

 
int main(int argc, const char* argv[]) { 
 // read in N filename nsteps delta_t graphics
 // N number of stars/particles to simulate 
 // filename is the filename of the file to read the initial configuration from 
 // nsteps is the number of timesteps
 // theta_max is the threshold value to be used in the Barnes-Hut algorithm
 // graphics is 1 or 0 meaning graphics on/off

// check if the parameters in the command line are correct, otherwise error message with instructions.	
  	if(argc != 7) {
      printf("Please give in: N filename nsteps delta_t theta_max graphics.\n");
      return -1;
    }
 
// read in N, check if N is 1 or larger otherwise error message.
  	int N = atoi(argv[1]);
    printf("N = %d\n", N);
    if(N < 1) {
      printf("Error: (N < 1).\n");
      return -1;
    }
 // read in theta_max
   double theta_max = atof(argv[5]);
      printf("theta_max = %f\n", theta_max);
   
 // read in filename and open filename. 	
   FILE *ptr_file;
 

  	ptr_file = fopen(argv[2], "r");
 
  	if(!ptr_file){
  		printf("File does not exist." );
  		return 1;}
 	  // store data of opened file 
   //
 		
  	fclose(ptr_file);
   
 	int nsteps = atoi(argv[3]);
 	double delta_t = atof(argv[4]);
 	
   int graphics = atoi(argv[6]);
  
 double *values =(double*)malloc(5*N*sizeof(double));
 read_doubles_from_file(atoi(argv[1])*5, values, argv[2]);
 
 //Allocate memory for particles  
 particle_t *particles = (particle_t*)malloc(N*sizeof(particle_t));
 
 //Set the particle data  
 int i = 0;
 int j = 0;  
 while(j<N){
    particles[j].x_pos = values[i];
    particles[j].y_pos = values[i+1];
    particles[j].mass = values[i+2];
    particles[j].vel_x = values[i+3];
    particles[j].vel_y = values[i+4];
    j++;
    i=j*5;
 }
   
 // assignment 4
   

 double epsilon=0.001;
 const double G=100.0/N;
  if(graphics ==0){ 
 for(int t=0;t<nsteps;t++) {
    p_qtree * head=(p_qtree *) malloc(sizeof(p_qtree));
    (*head).nw = NULL;
    (*head).ne = NULL;
    (*head).sw = NULL; 
    (*head).se = NULL; 
    (*head).width = 1.0;
    (*head).centerX = 0.5;
    (*head).centerY = 0.5;
    (*head).mass = 0;
    (*head).massCenterX = 0.5;
    (*head).massCenterY = 0.5;
    force_t * force = (force_t*)calloc(1,sizeof(force_t));
   
   // insert(&head, particles[0]);
   for(int k=0;k<N;k++)
   {
       insert(&head, particles[k]);
   }
   
  massification(&head);
         
   for(int i=0;i<N;i++){
     getForce(&head, particles[i],*force,theta_max,G,epsilon);
      
      double m_i = 1/particles[i].mass;
      particles[i].vel_x += delta_t*(*force).x*m_i;
      particles[i].vel_y += delta_t*(*force).y*m_i;
      particles[i].x_pos += delta_t*particles[i].vel_x;
      particles[i].y_pos += delta_t*particles[i].vel_y;  
   }
   
   delete(&head);
   free(force);
   }
  }
   
   
 
   if(graphics==5) {
      for(int t=0;t<nsteps;t++) {
         // dont use function?
         updateParticles(delta_t, particles, N);
      }
   }
   else if(graphics ==1) {
      int L = 1;
      int W = 1;
      int windowWidth = 600;
      int windowHeight = 600;
      SetCAxes(0,1);
      InitializeGraphics("",windowWidth,windowHeight);
      double x, y, circleRadius;
         
        for(int t=0;t<nsteps;t++) {
           
                          p_qtree * head=(p_qtree *) malloc(sizeof(p_qtree));
    (*head).nw = NULL;
    (*head).ne = NULL;
    (*head).sw = NULL; 
    (*head).se = NULL; 
    (*head).width = 1;
    (*head).centerX = 0.5;
    (*head).centerY = 0.5;
    (*head).mass = 0;
    (*head).massCenterX = 0.5;
    (*head).massCenterY = 0.5;
    
   
    insert(&head, particles[0]);
          // printf("x0: %lf\n",particles[0].x_pos);
          // printf("y0: %lf\n",particles[0].y_pos);
          // printf("x1: %lf\n",particles[1].x_pos);
           //printf("y1: %lf\n",particles[1].y_pos);
           
           
            for(int k=1;k<N;k++)
   {
       insert(&head, particles[k]);
   }
           
  massification(&head);
  printTree(&head);
           
           
           ClearScreen();           
           for(int i=0;i<N;i++) {
              x = particles[i].x_pos;
              //printf("%lf\n", x);
              y = particles[i].y_pos;
              //printf("%lf\n", y);
              circleRadius = 0.005;
              DrawCircle(x, y, L, W, circleRadius, 0.1);          
           }
           Refresh();
           usleep(8000);
           

           
   for(int i=0;i<N;i++){
      force_t * force = (force_t*)calloc(1,sizeof(force_t));
     getForce(&head, particles[i], *force,theta_max,G,epsilon);
      
      double m_i = 1/particles[i].mass;
      particles[i].vel_x += delta_t*(*force).x*m_i;
      printf("%lf\n\n", delta_t*(*force).x*m_i);
      particles[i].vel_y += delta_t*(*force).y*m_i;
      particles[i].x_pos += delta_t*particles[i].vel_x;
      particles[i].y_pos += delta_t*particles[i].vel_y;  
      free(force);
   }
           
   
   delete(&head);
   
   }
   
         
    
     FlushDisplay();
     CloseDisplay();
   } 

 i = 0;
 j = 0;  
 while(j<N){
    values[i] = particles[j].x_pos;
    values[i+1] = particles[j].y_pos;
    values[i+2] = particles[j].mass;
    values[i+3] = particles[j].vel_x;
    values[i+4] = particles[j].vel_y;
    j++;
    i=j*5;
 }
   write_doubles_to_file(5*N,values,"result.gal");
   
  return 0;
 
}
