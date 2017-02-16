//header file for quad



typedef struct force{
 double x;
 double y;
} force_t;

typedef struct particle
{
   double             x_pos;
   double             y_pos;
   double 	     mass;
   double      vel_x;
   double       vel_y;
} particle_t;

struct Quadtree {
	p_qtree * nw;
	p_qtree * ne;
	p_qtree * se;
	p_qtree * sw;
	particle p;
	double width;
	double centerX;
	double centerY;
	double massCenterX;
	double massCenterY;
	double mass=0;	
}p_qtree;

void insert(p_qtree ** node, particle p);

force_t getForce(p_qtree ** node, particle p, force_t force, double thetamax, double G, double epsilon);

void delete(p_qtree ** node);

void massification(p_qtree ** node);
