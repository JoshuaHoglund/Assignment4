
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "quad.h"
#include <math.h>

double dist(double x1, double x2, double y1, double y2) {
	double d = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	return d;
}

force_t getForce(p_qtree ** node, particle_t p, force_t force, double thetamax, double G, double eps) {
	double theta = (**node).width/dist(p.x_pos, (**node).centerX, p.x_pos, (**node).centerY);
	if ((**node).nw==NULL) {
		double r_x = p.x_pos - (**node).centerX;
		double r_y = p.y_pos - (**node).centerY;
		double abs_r = sqrt(r_x*r_x + r_y*r_y);
		force.x = -G*p.mass*(**node).mass*r_x/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
		force.y = -G*p.mass*(**node).mass*r_y/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
		return force;
	}
	if (theta>thetamax) {
		force.x = getForce((&(**node).nw).x_pos, p, force, thetamax, G, eps) + getForce((&(**node).ne).x_pos, p, force, thetamax, G, eps) + getForce((&(**node).sw).x_pos, p, force, thetamax, G, eps) + getForce((&(**node).se).x_pos, p, force, thetamax, G, eps);
		force.y = getForce(&(**node).nw).y + getForce(&(**node).ne).y + getForce(&(**node).sw).y + getForce(&(**node).se).y;
		return force;
	}
	else {
		double r_x = p.x - (**node).centerX;
		double r_y = p.y - (**node).centerY;
		double abs_r = sqrt(r_x*r_x + r_y*r_y);
		force.x = -G*p.mass*(**node).mass*r_x/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
		force.y = -G*p.mass*(**node).mass*r_y/((abs_r+eps)*(abs_r+eps)*(abs_r+eps));
		return force;
	}
}

void delete(p_qtree ** node) {
	if ((**node).nw==NULL) {
		free(*p_qtree);
		(*p_qtree)=NULL;
	}
	else {
		delete(&(**node).nw);
		delete(&(**node).ne);
		delete(&(**node).sw);
		delete(&(**node).se);
		free(*p_qtree);
		(*p_qtree)=NULL;
	}	
}

int compass(double px, double py, centerX, centerY) {
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
				p_qtree* sw, p_qtree se) {	
	switch(home) {
			case 1:
			    (*nw).particle=p;
				(*nw).mass = p.mass;
				break;
			case 2:
			    (*ne).particle=p;
				(*ne).mass = p.mass;
				break;
			case 3:
				(*sw).particle=p;
				(*sw).mass = p.mass;
				break;
			case 4:
				(*se).particle=p;
				(*se).mass = p.mass;
				break;
		}
}

void massification(p_qtree ** node) {
	double mass = (**node).mass;
	
	
	if (mass==0) {
		return;
	}
	else if((**node).nw==NULL){
		(**node).massCenterX = (**node).particle.x;
		(**node).massCenterY = (**node).particle.y;
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


void insert(p_qtree ** node, particle_t p) {
	p_qtree * nw = (**node).nw;
	p_qtree * sw = (**node).sw;
	p_qtree * ne = (**node).ne;
	p_qtree * se = (**node).se;
	double width = (**node).width;
	double centerX = (**node).centerX;
	double centerY = (**node).centerY;
	double mass = (**node).mass;
	
	// external node
	if (nw==NULL) {
		nw = (p_qtree *) malloc(sizeof(p_qtree));
		ne = (p_qtree *) malloc(sizeof(p_qtree));
		sw = (p_qtree *) malloc(sizeof(p_qtree));
		se = (p_qtree *) malloc(sizeof(p_qtree));
		(*nw).width = 0.5*width; 
		(*ne).width = 0.5*width;
		(*sw).width = 0.5*width;
		(*se).width = 0.5*width;
		(*nw).mass = 0; 
		(*ne).mass = 0;
		(*sw).mass = 0;
		(*se).mass = 0;
		(*nw).centerX = centerX-0.25*width;
		(*ne).centerX = centerX+0.25*width;
		(*sw).centerX = centerX-0.25*width;
		(*se).centerX = centerX+0.25*width;
		
		(*nw).centerY = centerY+0.25*width;
		(*ne).centerY = centerY+0.25*width;
		(*sw).centerY = centerY-0.25*width;
		(*se).centerY = centerY-0.25*width;
		
		
		if (mass==0) {
			assignHome(compass(p.x,p.y,centerX,centerY), p, nw, ne, sw, se);
		}
		else {
			int home1 = compass(p.x, p.y, centerX, centerY);
			int home2 = compass((**node).particle.x, (**node).particle.y, centerX, centerY);
			if (home1==home2) {
				switch(home1) {
						case 1:
							insert(&nw,p);
							insert(&nw,(**node).particle);
							break;
						case 2:
							insert(&ne,p);
							insert(&ne,(**node).particle);
							break;
						case 3:
							insert(&sw,p);
							insert(&sw,(**node).particle);
							break;
						case 4:
							insert(&se,p);
							insert(&se,(**node).particle);
							break;
				}
				
			}
			else {
				assignHome(home1, p, nw, ne, sw, se);
				assignHome(home2, (**node).particle, nw, ne, sw, se);
			}
		}
	}
	else {
		int home = compass(p.x, p.y, centerX, centerY);
		switch(home) {
			case 1:
				insert(&nw,p);
				break;
			case 2:
				insert(&ne,p);
				break;
			case 3:
				insert(&sw,p);
				break;
			case 4:
				insert(&se,p);
				break;
				}
	}
}



