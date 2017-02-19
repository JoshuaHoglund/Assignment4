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

typedef struct Quadtree {
	struct Quadtree * nw;
	struct Quadtree * ne;
	struct Quadtree * se;
	struct Quadtree * sw;
	particle_t p;
	double width;
	double centerX;
	double centerY;
	double massCenterX;
	double massCenterY;
	double mass;	
}p_qtree;

void insert(p_qtree ** node, particle_t p);

force_t* getForce(p_qtree ** node, particle_t p, double thetamax, double G, double epsilon);

void delete(p_qtree ** node);

void nullifyChildren(p_qtree **node);

void massification(p_qtree ** node);

void printTree(p_qtree ** node);
