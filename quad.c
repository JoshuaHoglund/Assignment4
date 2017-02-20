
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "quad.h"
#include <math.h>

double dist(double x1, double x2, double y1, double y2) {
	double d = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	return d;
}

force_t* getForce(p_qtree ** node, particle_t p, double thetamax, double G, double eps) {
	double theta = (**node).width/dist(p.x_pos, (**node).centerX, p.y_pos, (**node).centerY);
	//If (**node).nw==NULL all other children are also NULL
	
	
	if ((**node).nw==NULL && p.x_pos!=(**node).massCenterX) {
		force_t* force = calloc(1,sizeof(force_t));
		double r_x = p.x_pos - (**node).massCenterX;
		double r_y = p.y_pos - (**node).massCenterY;
		
		double abs_r = sqrt(r_x*r_x + r_y*r_y);
		double r3=1/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
		(*force).x = -G*p.mass*(**node).mass*r_x*r3;
		(*force).y = -G*p.mass*(**node).mass*r_y*r3;
		/*if (sqrt((*force).x*(*force).x)<0.001) {
			(*force).x = 0;
			(*force).y = 0;
		}*/
		return force;
	}
	
	
	
	else if (theta>thetamax && p.x_pos!=(**node).massCenterX) {
		force_t * force = calloc(1,sizeof(force_t));
		force_t * nwforce = getForce((&(**node).nw),p, thetamax, G, eps);
		force_t * neforce = getForce((&(**node).ne),p, thetamax, G, eps);
		force_t * swforce = getForce((&(**node).sw),p, thetamax, G, eps);
		force_t * seforce = getForce((&(**node).se),p, thetamax, G, eps);
		
		(*force).x = (*nwforce).x + (*neforce).x + (*swforce).x + (*seforce).x;
		(*force).y = (*nwforce).y + (*neforce).y + (*swforce).y + (*seforce).y;
		free(nwforce);
		free(neforce);
		free(swforce);
		free(seforce);
		
		return force;
	}
	
	else if(p.x_pos!=(**node).massCenterX) {
		force_t * force = calloc(1,sizeof(force));
		double r_x = p.x_pos - (**node).massCenterX;
		double r_y = p.y_pos - (**node).massCenterY;
		double abs_r = sqrt(r_x*r_x + r_y*r_y);
		double r3=1/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
		(*force).x += -G*p.mass*(**node).mass*r_x*r3;
		(*force).y += -G*p.mass*(**node).mass*r_y*r3;
		return force;
	}
	
	
	
	else {
		force_t * force = calloc(1,sizeof(force_t));
		return force;
	}
}

void printTree(p_qtree ** node) {
	if ((**node).nw==NULL) {
		if ((**node).mass==0) {
			printf("Empty \n");
		} else {
		printf("x: %lf \n y: %lf \n mass: %lf \n",(**node).massCenterX,(**node).massCenterY,(**node).mass);
		}
	}
	else {
		printf("new head \n");
		printf("nw \n");
		printTree(&(**node).nw);
		printf("ne \n");
		printTree(&(**node).ne);
		printf("sw \n");
		printTree(&(**node).sw);
		printf("se \n");
		printTree(&(**node).se);
		printf("end of head \n");
	}
}

void delete(p_qtree ** node) {
	if ((**node).nw==NULL) {
		free(*node);
		(*node)=NULL;
	}
	else {
		delete(&(**node).nw);
		delete(&(**node).ne);
		delete(&(**node).sw);
		delete(&(**node).se);
		free(*node);
		(*node)=NULL;
	}	
}

int compass(double px, double py, double centerX, double centerY) {
	int res; // nw=1 ne=2 sw=3 se=4
	// west
			if (px<centerX) {
				// south
				if(py<centerY) {
					res = 3;
				}
				// north
				else {
					res = 1;
				}
			}
			// east south
			else if(py<centerY){
				res = 4;
			}
			// east north
			else {
				res = 2;
			}
			return res;
}

void assignHome(int home, particle_t p , p_qtree * nw, p_qtree * ne, 
				p_qtree* sw, p_qtree* se) {	
	switch(home) {
			case 1:
			    (*nw).p=p;
				(*nw).mass = p.mass;
				break;
			case 2:
			    (*ne).p=p;
				(*ne).mass = p.mass;
				break;
			case 3:
				(*sw).p=p;
				(*sw).mass = p.mass;
				break;
			case 4:
				(*se).p=p;
				(*se).mass = p.mass;
				break;
		}
}

void massification(p_qtree ** node) {
	double mass = (**node).mass;
	
	
	if ((**node).nw==NULL && mass==0) {
		return;
	}
	else if((**node).nw==NULL){
		(**node).massCenterX = (**node).p.x_pos;
		(**node).massCenterY = (**node).p.y_pos;
	}
	else {
		massification(&(**node).nw);
		massification(&(**node).sw);
		massification(&(**node).ne);
		massification(&(**node).se);
		(**node).mass = (*(**node).nw).mass + (*(**node).se).mass 
		 + (*(**node).sw).mass +(*(**node).ne).mass; 
		(**node).massCenterX = ( (*(**node).nw).mass * (*(**node).nw).massCenterX +
				(*(**node).ne).mass * (*(**node).ne).massCenterX +
				(*(**node).sw).mass * (*(**node).sw).massCenterX +
				(*(**node).se).mass * (*(**node).se).massCenterX)/(**node).mass;
		(**node).massCenterY = ( (*(**node).nw).mass * (*(**node).nw).massCenterY +
				(*(**node).ne).mass * (*(**node).ne).massCenterY +
				(*(**node).sw).mass * (*(**node).sw).massCenterY +
				(*(**node).se).mass * (*(**node).se).massCenterY)/(**node).mass;
				
								
		}
}

void nullifyChildren(p_qtree ** node) {
	(**node).nw = NULL;
	(**node).ne = NULL;
	(**node).sw = NULL;
	(**node).se = NULL;
}

void insert(p_qtree ** node, particle_t p) {
	double width = (**node).width;
	double centerX = (**node).centerX;
	double centerY = (**node).centerY;
	double mass = (**node).mass;
	// external node, if (**node).nw==NULL all other children are also NULL
	if ((**node).nw==NULL) {
		if(mass!=0) {
		(**node).nw = (p_qtree *) calloc(1,sizeof(p_qtree));
		(**node).ne = (p_qtree *) calloc(1,sizeof(p_qtree));
		(**node).sw = (p_qtree *) calloc(1,sizeof(p_qtree));
		(**node).se = (p_qtree *) calloc(1,sizeof(p_qtree));
		nullifyChildren(&(**node).nw);
		nullifyChildren(&(**node).ne);
		nullifyChildren(&(**node).sw);
		nullifyChildren(&(**node).se);
		(*(**node).nw).width = 0.5*width; 
		(*(**node).ne).width = 0.5*width;
		(*(**node).sw).width = 0.5*width;
		(*(**node).se).width = 0.5*width;
		(*(**node).nw).mass = 0; 
		(*(**node).ne).mass = 0;
		(*(**node).sw).mass = 0;
		(*(**node).se).mass = 0;
		(*(**node).nw).centerX = centerX-0.25*width;
		(*(**node).ne).centerX = centerX+0.25*width;
		(*(**node).sw).centerX = centerX-0.25*width;
		(*(**node).se).centerX = centerX+0.25*width;
		
		(*(**node).nw).centerY = centerY+0.25*width;
		(*(**node).ne).centerY = centerY+0.25*width;
		(*(**node).sw).centerY = centerY-0.25*width;
		(*(**node).se).centerY = centerY-0.25*width;		
		
		int home1 = compass(p.x_pos, p.y_pos, centerX, centerY);
		int home2 = compass((**node).p.x_pos, (**node).p.y_pos, centerX, centerY);
		if (home1==home2) {
				switch(home1) {
						case 1:
							insert(&(**node).nw,p);
							insert(&(**node).nw,(**node).p);
							break;
						case 2:
							insert(&(**node).ne,p);
							insert(&(**node).ne,(**node).p);
							break;
						case 3:
							insert(&(**node).sw,p);
							insert(&(**node).sw,(**node).p);
							break;
						case 4:
							insert(&(**node).se,p);
							insert(&(**node).se,(**node).p);
							break;
				}
				
			}
			else {
				assignHome(home1, p, (**node).nw, (**node).ne, (**node).sw, (**node).se);
				assignHome(home2, (**node).p, (**node).nw, (**node).ne, (**node).sw, (**node).se);
			}			
		}		
		else {
			(**node).p=p;
			(**node).mass = p.mass;
		}		
	}
	else {
		int home = compass(p.x_pos, p.y_pos, centerX, centerY);
		switch(home) {
			case 1:
				insert(&(**node).nw,p);
				break;
			case 2:
				insert(&(**node).ne,p);
				break;
			case 3:
				insert(&(**node).sw,p);
				break;
			case 4:
				insert(&(**node).se,p);
				break;
				}
	}
}



