#include <iostream>
#include <queue>
#include <set>
#include <math.h>
#include <GL/glut.h>
#define PI 3.14159265

using namespace std;

//Leaf Venation
void growth(bool sqhex);


//Voronoi function
void voronoi(bool sqhex);

// Notation for working with points
typedef pair<double, double> point;
#define x first
#define y second

// Arc, event, and segment datatypes
struct arcbst;
struct seg;

struct event {
   double x;
   point p;
   arcbst *a;
   bool valid;

   event(double xx, point pp, arcbst *aa)
      : x(xx), p(pp), a(aa), valid(true) {}
};

struct arcbst {
   point p1,p2;
   arcbst *left, *right, *parent;
   event *e;
   bool tp;//node type -> tells you if it's internal or leaf node ... could be done checking whether left or right children are null, however, code is cleaner this way (not really worried about an extra bool) ----- tp = true means the node is an internal node ... tp = false means its an leaf node
   bool gt;//in case its an internal node this bools tells you whether you are looking for the greater or the lesser intersection (with regard to the y coordinate)
   seg *s;

   arcbst(point pp, bool nodetype, arcbst *b,  point pp2, bool gtt=0)
    : p1(pp), left(0), right(0), parent(b), e(0), s(0), tp(nodetype), p2(pp2), gt(gtt) {}
};

vector<seg*> output;  // Array of output segments.

struct seg {
   point start, end;//the extreme points of the segment (AKA edge)
   bool st, ed;//these just tells the function whether the start or the end points were defined already

   seg(point pp) : start(pp), end(0,0), st(true), ed(false) { output.push_back(this); }//one can start the edge with one of the extremes defined
   seg() : start(0,0), end(0,0), st(false), ed(false) { output.push_back(this); }//or one can start the edge without any of the extremes defined
   
   // Set the end point and mark as "done."
   void finishst(point p) { if (st) return; start = p; st = true; }
   void finishend(point p) { if (ed) return; end = p; ed = true; }
};

arcbst *root = 0; // First node of the binary search tree that represents the beachline

// Function declarations
point getmin(arcbst *i);
arcbst* findarc(point p);
arcbst* getnext(arcbst *i, bool full);
arcbst* getprev(arcbst *i, bool full);
void arcinsert(point p);
void arcdelete(arcbst *i);
void finishtraversal(arcbst *ptr, double l);
void printtree(arcbst* ptr, double l);
void openscadoutput();
void checkintersection(arcbst *i, arcbst *temp, double x0);
void circleevent();
bool circle(point a, point b, point c, double *x, point *o);
void check_circle_event(arcbst *i, double x0);
bool intersect(point p, arcbst *i, point *result = 0);
point intersection(arcbst *i, double l);
void finish_edges(bool sqhex);
void deleteroot(arcbst *root);

// "Greater than" comparison, for reverse sorting in priority queue.
struct gt {
   bool operator()(point a, point b) {return a.x==b.x ? a.y>b.y : a.x>b.x;}
   bool operator()(event *a, event *b) {return a->x>b->x;}
};

// Bounding box coordinates.
double X0 = 0, X1 = 0, Y0 = 0, Y1 = 0;
double uX0 = 0, uX1 = 0, uY0 = 0, uY1 = 0;

//DCEL
//edges

struct veinnode;
struct vert;
struct face;

struct hedge{
	hedge *next,*prev,*twin;
	vert *origin;
	face *f;
	hedge() : next(0), prev(0), twin(0), origin(0), f(0) {}
};

vector<vert*> dcelist;
vector<face*> dcefacelist;
vector<vert*> veinvtcs;
vector<face*> veinfaces;
vector<vert*> auxinvtcs;
vector<face*> auxinfaces;

struct vert{
	point p;
	hedge *edge;
	vert(point pp)
	: p(pp) { dcelist.push_back(this); 
	dcelist.back()->edge=0; }
	vert() {}
};

struct face{
	hedge *edge;
	point site; //site of the corresponding face
	veinnode *vein;
	bool del;
	vector<point> auxin;
	face(hedge *halfedge) : edge(halfedge), site(0,0),del(0), vein(0) {}
};

struct auxinnode;

struct veinnode{
	point p;
	vector<veinnode*> next;
	veinnode *prev;
	int order;
	vector<auxinnode*> auxin;
	double width;//according to Murray's law what's the diameter of this vein segment???
	veinnode(point pp) : p(pp) , width(0), prev(0), order(10) {}
}*veinroot(0),*veinroot2(0),*veinroot3(0);

struct auxinnode{
	point p;
	auxinnode *next, *prev;
	vector<veinnode*> vein;
	bool del;
	auxinnode(point pp) : p(pp), next(0), prev(0), del(0) {}
}*auxinroot(0);

void updateveinnodes(auxinnode *i);
void auxinadd(point p);
auxinnode* auxindel(auxinnode *temp);
veinnode* findvein(point p, veinnode *ptr);
bool veinintersection(point p0, point p1, veinnode *ptr);
void leafvenationsim(double denst, int numdiv, bool sqhex);
bool voronoistructure(bool veinorauxin, bool sqhex);
void cleanup();
void updateveingrowth(veinnode *ptr, double unitymod);
void printvein(veinnode *ptr);
void printveinout(veinnode *ptr);
void defineveinorder(veinnode *ptr, int order);
void veinlengtharea(veinnode *ptr, int order);
void printveinout2(veinnode *ptr);
double voronoicellarea(hedge *orgn);
void defineveinwidth(veinnode *ptr);
void deltreebeyondcutoff(veinnode *ptr);
double computeAverageDistTree();
vector<double> computeAverageDistSite(point p);

hedge *cche=0,*che=0; //counterclockwise half-edge and clockwise half-edge
int facecounter=0;
bool flag;

void addface(hedge *temp);
bool findinst(point p);
int findit(point p);
int findfaceit(hedge *temp);
void starthedge(int i);
void getcounterclockhe(hedge *temp, int i,bool change);
int createdcel();
void addvert(point p);
void deletefaceedges(hedge *temp0, hedge *temp1);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Border treatment functions~~~~~~~~~~~~~~~~
point checkborder(point p0, point p1,bool sqhex);
void increaseborder(point *p);
int outofbordercheck(int i, bool sqhex);
//bool outofbordercheck2(point p);
bool orderx(seg* i, seg* j);
bool ordery(seg* i, seg* j);
bool orderxinv(seg* i, seg* j);
bool orderyinv(seg* i, seg* j);
int signof(double a);
double crossproduct2d(point p1, point p2, point p3);//the result is in the magnitude in the z direction (point p1 is defined as the origin of the two vectors being croosproducted)
bool linesegintersection(point p1, point p2, point p3, point p4);//do the line segments p1-p2 and p3-p4 intersect?
bool checkpointhex(point p);
bool checkpointsq(point p);
bool pointsegintersection(point p1, point p2, point p3);//does point p1 intersect with segment p2-p3?
double rounding(double x, int n);//round number up to the nth decimal digit
bool roundingpoint(point p1, point p2, int n);//round number up to the nth decimal digit


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Binary search tree data structure implementation~~~~~~~~~~~~~~~

struct bstnode{
	bstnode *left;
	bstnode *right;
	double ycoord;//no points aside from the six border points have the exact same y coordinate, and they should NOT be deleted so no worries
	vector<int> index;//indexes of output vector in which the point appear
	vector<bool> strorend;//is the point the start point or the end point? 1 - start// 2- end
	bstnode(double yy,int i,bool stoe) : left(0) , right(0), ycoord(yy) { index.push_back(i); strorend.push_back(stoe); };
}*bstroot(0);

void bstinsert(double ycoord, int ind, bool stoe);
void deletebst(bstnode* ptr);

struct bstdcelnode{
	bstdcelnode *left;
	bstdcelnode *right;
	double ycoord;//y coordinates, again will be used to generate the bst
	vector<int> index;//index of the vertice in the dcelist vector 
	bstdcelnode(double yy,int i) : left(0) , right(0), ycoord(yy) { index.push_back(i); };
	bstdcelnode() {};
}*bstdcelroot(0);

void bstdcelinsert(double ycoord, int ind, point p0);
void inOrder2(bstdcelnode *ptr);
void deletebstdcel(bstdcelnode* ptr);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Adding and Deleting Sites of Voronoi Tesselation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void deletesite(point p);
void deletesite2(point p, bool veinorauxin);
void getnearestneighbors(hedge *temp0,bool clear);
void getnearestneighbors2(hedge *temp0,bool clear,bool nnonly);
void createsite2(point p, bool veinorauxin);
void createsite(point p);
void delaunayneighbors(point p);

vector<seg*> tempoutput;





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Point location query functions and data structure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Slab Method~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//bool findinst2(point p);

struct slabbsty;

struct slabbstx{
	slabbstx *right;
	slabbstx *left;
	slabbstx *next,*prev;
	double xcoord,xcoordf;
	vector<int> index;
	slabbsty *slaby;
	slabbstx(double xx, int i) : xcoord(xx) , right(0) , left(0), slaby(0), next(0), prev(0) {index.push_back(i);};
	slabbstx() : xcoord(0) , right(0), left(0), slaby(0), next(0), prev(0) {};
}*slabbstxroot(0), *slabbstxroot2(0);

void slabbstxinsert(double xcoord,int ind, point p0,bool root1vs2);
void slabbstxaddxf(bool root1vs2);
void deleteslabbstx(slabbstx* ptr);
void deleteslabbsty(slabbsty* ptr);

struct slabbsty{
	slabbsty *up;
	slabbsty *down;
	face *fup, *fdown;
	double a,b,ycoord;
	slabbsty(double aa, double bb, double yy) : a(aa), b(bb), ycoord(yy), up(0), down(0), fup(0), fdown(0) {};
}*slabbstyroot(0);

void createslabstruct(bool root1vs2);
void slabbstyinsert(double ycoord, slabbsty** rootptr, double a, double b, hedge* temphedge);

struct slabybuild{
	double a,b,xcoord;
	hedge *edge;
};

void inOrder(slabbsty *ptr,double xcoord, double xcoordf);
face* plquery(point p0, bool root1vs2);

void sitefacelinkage(bool p2ortsites);

///////////////Grid data structure/////////////////////
struct griddiv{
	vector<point> vetcs;
	bool ok;
	int quota;
	griddiv() : ok(0), quota(0) {}
};

vector<griddiv*> grid;

void initiatesqgrid(int numdiv);
void checkgrid();
void genrddist(double denst, int numdiv, bool fullcheck);//generates random distribution of points in each square partition of the grid
int gridindget(point p, int numdiv);
void checkveingrid(veinnode *ptr, int numdiv);
bool orderxpt(point i, point j);
int getoccupiedboxes(veinnode *ptr,int numdiv,bool ven1or2);
bool checknneighbors(int xind, int yind, int numdiv2, int knn);
