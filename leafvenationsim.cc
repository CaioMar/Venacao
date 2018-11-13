///////////////////////////////////// leafvenationsim3.cc ///////////////////////////////////////////
// Program developed to generate open venation patterns and interdigitated patterns //

#include "leafvenationsim2.hh"
#include "mtrand.h"
#include "gl2ps.h"
#include <fstream>
#include <algorithm>
#include <GL/glut.h>

GLdouble pointRadius(0.01); 
priority_queue<point,  vector<point>,  gt> points; // site events
priority_queue<event*, vector<event*>, gt> events; // circle events
vector<point> points2, tempsites;
vector<point> veinpoints, veinpoints2, veinpoints3, auxinpoints;
double scalefactor(1);
int mosaic;//do you want to assign colors to each face?
double col;//variable used to assign the color to each face in the diagram
int tentativa(0);
point hexvert[8],hexvert2[8]; //border points of the hexagon
point bintersection, squarevert[6];
point pproblem(0,0);////////////////////////Voronoi site for the border face (could be any point outside the border)//////////////////////
double moveh(0),movev(0);
face *queryface;
vector<seg*> testout;
double winsy,winsx;
int iterdcel(0),checker(0),auxincutoff(0),ordercutoffprim(0),ordercutoffcomp(0);
point center;
int ven2multiplier(3), growthstepmax(20),maxtries(300), printorder(0);
double maxsteps2, bdmultiplier, kdmultiplier, denstmultiplier, veindistmultiplier,venationlength(0), bifurcationcount(0);
double finallength(0),birthdist(4),killdist(4),veindist(1.84),seed(1),lstepgrowth(2), minwidth(0.5), chanheight(0.250), murrayexp(0), initialdiam(0),murrayexp2(0);
ofstream output1,output2,output3;
bool squarehex, v2, background, diagbool, removeendchannels, autooutput;

/*
void growth(bool sqhex){
	double unitymod(rounding(sqrt(pow(center.x,2)+pow(center.y,2)),12));
	point r1;
	if(sqhex){
		r1.x=0;
		r1.y=0;
	}else{
		r1.x=1;
		r1.y=0;
	}
	//growth of voronoi sites
	for(int i=0;i<dcefacelist.size();i++){
		dcefacelist[i]->site.x=rounding(dcefacelist[i]->site.x+lstepgrowth*((dcefacelist[i]->site.x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		dcefacelist[i]->site.y=rounding(dcefacelist[i]->site.y+lstepgrowth*((dcefacelist[i]->site.y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	for(int i=0;i<veinfaces.size();i++){
		veinfaces[i]->site.x=rounding(veinfaces[i]->site.x+lstepgrowth*((veinfaces[i]->site.x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		veinfaces[i]->site.y=rounding(veinfaces[i]->site.y+lstepgrowth*((veinfaces[i]->site.y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	for(int i=0;i<dcelist.size();i++){
		dcelist[i]->p.x=rounding(dcelist[i]->p.x+lstepgrowth*((dcelist[i]->p.x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		dcelist[i]->p.y=rounding(dcelist[i]->p.y+lstepgrowth*((dcelist[i]->p.y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	for(int i=0;i<veinvtcs.size();i++){
		veinvtcs[i]->p.x=rounding(veinvtcs[i]->p.x+lstepgrowth*((veinvtcs[i]->p.x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		veinvtcs[i]->p.y=rounding(veinvtcs[i]->p.y+lstepgrowth*((veinvtcs[i]->p.y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	if(sqhex){
	for(int i=0;i<8;i++){
		hexvert[i].x=rounding(hexvert[i].x+lstepgrowth*((hexvert[i].x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		hexvert[i].y=rounding(hexvert[i].y+lstepgrowth*((hexvert[i].y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	}else{
	for(int i=0;i<6;i++){
		squarevert[i].x=rounding(squarevert[i].x+lstepgrowth*((squarevert[i].x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		squarevert[i].y=rounding(squarevert[i].y+lstepgrowth*((squarevert[i].y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	}
	for(int i=0;i<veinpoints.size();i++){
		veinpoints[i].x=rounding(veinpoints[i].x+lstepgrowth*((veinpoints[i].x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		veinpoints[i].y=rounding(veinpoints[i].y+lstepgrowth*((veinpoints[i].y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	for(int i=0;i<grid.size();i++){
		for(int j=0;j<4;j++){
			grid[i]->vetcs[j].x=rounding(grid[i]->vetcs[j].x+lstepgrowth*((grid[i]->vetcs[j].x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
			grid[i]->vetcs[j].y=rounding(grid[i]->vetcs[j].y+lstepgrowth*((grid[i]->vetcs[j].y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
		}
	}
	for(int i=0;i<auxinpoints.size();i++){
		auxinpoints[i].x=rounding(auxinpoints[i].x+lstepgrowth*((auxinpoints[i].x-center.x)/unitymod+(center.x-r1.x)/unitymod),12);
		auxinpoints[i].y=rounding(auxinpoints[i].y+lstepgrowth*((auxinpoints[i].y-center.y)/unitymod+(center.y-r1.y)/unitymod),12);
	}
	//take the center of the hexagon to its next coordinates
	updateveingrowth(veinroot,unitymod);
	updateveingrowth(veinroot2,unitymod);
//	cout<<"Ok"<<endl;
	if(slabbstxroot){
	deleteslabbstx(slabbstxroot);
	slabbstxroot=0;   
	for(int i=0;i<dcelist.size();i++){
		if(dcelist[i]) slabbstxinsert(dcelist.at(i)->p.x,i,dcelist.at(i)->p,1);
   	}
   	slabbstxaddxf(1);/////////adds the x coordinate where each slab ends to each slab in the of the data structure
   	createslabstruct(1);//////creates all of the binary search trees in the y direction and adds them to each slab of the data structure
	}
//	cout<<"Ok2"<<endl;
	if(slabbstxroot2){
	deleteslabbstx(slabbstxroot2);
	slabbstxroot2=0;   
	for(int i=0;i<veinvtcs.size();i++){
		if(veinvtcs[i]) slabbstxinsert(veinvtcs.at(i)->p.x,i,veinvtcs.at(i)->p,0);
   	}
   	slabbstxaddxf(0);/////////adds the x coordinate where each slab ends to each slab in the of the data structure
   	createslabstruct(0);//////creates all of the binary search trees in the y direction and adds them to each slab of the data structure
	}

	center.x=rounding(center.x+lstepgrowth*((center.x-r1.x)/unitymod),12);
	center.y=rounding(center.y+lstepgrowth*((center.y-r1.y)/unitymod),12);

//	cout<<"New center: "<<center.x<<" "<<center.y<<endl;
}
*/

void growth(bool sqhex){
	point r1, lprime;
	if(sqhex){
		r1.x=0;r1.y=0;
	}else{
		r1.x=1;r1.y=0;
	}
	double unitymod(rounding(sqrt(pow(r1.x-center.x,2)+pow(r1.y-center.y,2)),12)),modulus;
//	modulus = sqrt(pow(r1.x-center.x,2)+pow(r1.y-center.y,2));
	lprime.x = lstepgrowth*(r1.x-center.x)/(unitymod);
	lprime.y = lstepgrowth*(r1.y-center.y)/(unitymod);
	//growth of voronoi sites
	for(int i=0;i<dcefacelist.size();i++){
		dcefacelist[i]->site.x=rounding(dcefacelist[i]->site.x+lstepgrowth*((dcefacelist[i]->site.x-center.x)/(unitymod)) - lprime.x,12);
		dcefacelist[i]->site.y=rounding(dcefacelist[i]->site.y+lstepgrowth*((dcefacelist[i]->site.y-center.y)/(unitymod)) - lprime.y,12);
	}
	for(int i=0;i<veinfaces.size();i++){
		veinfaces[i]->site.x=rounding(veinfaces[i]->site.x+lstepgrowth*((veinfaces[i]->site.x-center.x)/(unitymod)) - lprime.x,12);
		veinfaces[i]->site.y=rounding(veinfaces[i]->site.y+lstepgrowth*((veinfaces[i]->site.y-center.y)/(unitymod)) - lprime.y,12);
	}
	for(int i=0;i<dcelist.size();i++){
		dcelist[i]->p.x=rounding(dcelist[i]->p.x+lstepgrowth*((dcelist[i]->p.x-center.x)/(unitymod)) - lprime.x,12);
		dcelist[i]->p.y=rounding(dcelist[i]->p.y+lstepgrowth*((dcelist[i]->p.y-center.y)/(unitymod)) - lprime.y,12);
	}
	for(int i=0;i<veinvtcs.size();i++){
		veinvtcs[i]->p.x=rounding(veinvtcs[i]->p.x+lstepgrowth*((veinvtcs[i]->p.x-center.x)/(unitymod)) - lprime.x,12);
		veinvtcs[i]->p.y=rounding(veinvtcs[i]->p.y+lstepgrowth*((veinvtcs[i]->p.y-center.y)/(unitymod)) - lprime.y,12);
	}
	if(sqhex){
	for(int i=0;i<8;i++){
		hexvert[i].x=rounding(hexvert[i].x+lstepgrowth*((hexvert[i].x-center.x)/(unitymod)) - lprime.x,12);
		hexvert[i].y=rounding(hexvert[i].y+lstepgrowth*((hexvert[i].y-center.y)/(unitymod)) - lprime.y,12);
	}
	}else{
	for(int i=0;i<6;i++){
		squarevert[i].x=rounding(squarevert[i].x+lstepgrowth*((squarevert[i].x-center.x)/(unitymod)) - lprime.x,12);
		squarevert[i].y=rounding(squarevert[i].y+lstepgrowth*((squarevert[i].y-center.y)/(unitymod)) - lprime.y,12);
//		cout<<"vertex "<<i<<":"<<squarevert[i].x<<" "<<squarevert[i].y<<endl;
	}
	}
	for(int i=0;i<veinpoints.size();i++){
		veinpoints[i].x=rounding(veinpoints[i].x+lstepgrowth*((veinpoints[i].x-center.x)/(unitymod)) - lprime.x,12);
		veinpoints[i].y=rounding(veinpoints[i].y+lstepgrowth*((veinpoints[i].y-center.y)/(unitymod)) - lprime.y,12);
	}
	if(veinpoints3.size() != 0){
		for(int i=0;i<veinpoints3.size();i++){
			veinpoints3[i].x=rounding(veinpoints3[i].x+lstepgrowth*((veinpoints3[i].x-center.x)/(unitymod)) - lprime.x,12);
			veinpoints3[i].y=rounding(veinpoints3[i].y+lstepgrowth*((veinpoints3[i].y-center.y)/(unitymod)) - lprime.y,12);
		}
	}
	for(int i=0;i<grid.size();i++){
		for(int j=0;j<4;j++){
			grid[i]->vetcs[j].x=rounding(grid[i]->vetcs[j].x+lstepgrowth*((grid[i]->vetcs[j].x-center.x)/(unitymod)) - lprime.x,12);
			grid[i]->vetcs[j].y=rounding(grid[i]->vetcs[j].y+lstepgrowth*((grid[i]->vetcs[j].y-center.y)/(unitymod)) - lprime.y,12);
		}
	}
	for(int i=0;i<auxinpoints.size();i++){
		auxinpoints[i].x=rounding(auxinpoints[i].x+lstepgrowth*((auxinpoints[i].x-center.x)/(unitymod)) - lprime.x,12);
		auxinpoints[i].y=rounding(auxinpoints[i].y+lstepgrowth*((auxinpoints[i].y-center.y)/(unitymod)) - lprime.y,12);
	}
	//take the center of the hexagon to its next coordinates
	updateveingrowth(veinroot,unitymod);
	updateveingrowth(veinroot2,unitymod);
//	cout<<"Ok"<<endl;
	if(slabbstxroot){
	deleteslabbstx(slabbstxroot);
	slabbstxroot=0;   
	for(int i=0;i<dcelist.size();i++){
		if(dcelist[i]) slabbstxinsert(dcelist.at(i)->p.x,i,dcelist.at(i)->p,1);
   	}
   	slabbstxaddxf(1);/////////adds the x coordinate where each slab ends to each slab in the of the data structure
   	createslabstruct(1);//////creates all of the binary search trees in the y direction and adds them to each slab of the data structure
	}
//	cout<<"Ok2"<<endl;
	if(slabbstxroot2){
	deleteslabbstx(slabbstxroot2);
	slabbstxroot2=0;   
	for(int i=0;i<veinvtcs.size();i++){
		if(veinvtcs[i]) slabbstxinsert(veinvtcs.at(i)->p.x,i,veinvtcs.at(i)->p,0);
   	}
   	slabbstxaddxf(0);/////////adds the x coordinate where each slab ends to each slab in the of the data structure
   	createslabstruct(0);//////creates all of the binary search trees in the y direction and adds them to each slab of the data structure
	}
	
	center.x=rounding(center.x - lprime.x,12);
	center.y=rounding(center.y - lprime.y,12);

//	cout<<"New center: "<<center.x<<" "<<center.y<<endl;
}

double computeAverageDistTree(){
	double avdist(0);//stores the average distance between trees in this variable
	vector<double> tempdist;

	output1.open("disthistogram");
	//Computes the average distance between complementary tree with sites stored in the veinpoints2 vector and the
	//primary tree with sites stored in the veinpoints vector
	for(int i=0;i<veinpoints2.size();i++){
		tempdist = computeAverageDistSite(veinpoints2[i]);
		avdist += tempdist[0];//computeAverageDistSite(veinpoints2[i]);//computes harmonic average distance between this site and its face neighbors in the primary tree
	//	cout<<"av tot: "<<avdist<<" "<<veinpoints2[i].x<<" "<<veinpoints2[i].y<<endl;
		output1 << tempdist[0] <<'\t'<< tempdist[1]<< '\t'<< tempdist[2]<<'\t'<<tempdist[3]<<'\t'<<tempdist[4]<<'\t'<<veinpoints2[i].x<<'\t'<<veinpoints2[i].y<<'\t'<<tempdist[5]<<'\t'<<tempdist[6]<<'\t'<<tempdist[7]<<'\t'<<tempdist[8]<<endl;
	}
	avdist = avdist/veinpoints2.size();
	output1.close();

	output1.open("venconnections");
	output1 <<"startx"<<'\t'<<"starty"<<'\t'<<"endx"<<'\t'<<"endy"<<'\t'<<"disttoroot"<<
endl;
	saveVenationStructure(veinroot,0);
	output1.close();

	output1.open("ven2connections");
	output1 <<"startx"<<'\t'<<"starty"<<'\t'<<"endx"<<'\t'<<"endy"<<'\t'<<"disttoroot"<<endl;
	saveVenationStructure(veinroot2,0);
	output1.close();
	return avdist;//returns average distance between trees
}

vector<double> computeAverageDistSite(point p){
	double harmavdist(0),arthavdist(0),geomavdist(1), tempdist(0);//stores the harmonic mean in this variable
	bool flagpoint(0);
	double NN(0);
	int totalNeighbors, NNindx(0);
	vector<double> averages;
	averages.push_back(harmavdist);
	averages.push_back(arthavdist);
	averages.push_back(geomavdist);
	averages.push_back(double(findvein(p, veinroot2)->order));
	averages.push_back(0);
	averages.push_back(0);
	averages.push_back(0);
	averages.push_back(0);
	averages.push_back(0);

	//Adds the point of interest to the veinpoints vector
	veinpoints.push_back(p);

	//creates the voronoi diagram with the current veinpoints vector
	cleanup();
	voronoistructure(1,squarehex);
	
	//gets the face pointer of the point of interest
	queryface=plquery(p,1);

	//Adds all the face neighbors of the point of interest to the tempsites vector	
	hedge *temp1(queryface->edge->next), *temp2;
	tempsites.clear();
	temp2=queryface->edge->twin;
	if(temp2->f->site!=pproblem) tempsites.push_back(temp2->f->site);
	while(temp1!=queryface->edge){
		temp2=temp1->twin;
		if(temp2->f->site!=pproblem){
			flagpoint=true;
			for(int i=0;i<tempsites.size();i++){
				if(temp2->f->site==tempsites.at(i)||temp2->f->site==queryface->site) flagpoint=false;
			}
			if(flagpoint) tempsites.push_back(temp2->f->site);
		}
		temp1=temp1->next;	
	}	

	
	if(tempsites.size() == 0) return averages;	

	//Adds the point of interest to the tempsites vector
	tempsites.push_back(p);
	
	//Removes the point of interest from the veinpoints vector
	veinpoints.erase(veinpoints.begin()+veinpoints.size()-1);

	totalNeighbors = tempsites.size() - 1;
	//Determines the harmonic mean of the distances between all face neighbor sites and the point of interest 
	NN = sqrt(pow(tempsites[0].x - tempsites[tempsites.size()-1].x,2)+pow(tempsites[0].y - tempsites[tempsites.size()-1].y,2));
	for(int i = 0; i < tempsites.size() - 1; i++){
		tempdist = sqrt(pow(tempsites[i].x - tempsites[tempsites.size()-1].x,2)+pow(tempsites[i].y - tempsites[tempsites.size()-1].y,2));
		harmavdist += 1.0/tempdist;
		arthavdist += tempdist;
		geomavdist *= tempdist;
		if(tempdist < NN){
			NN = tempdist;
			NNindx = i;
		}
	}
	//	cout<<"av inside: "<<avdist<<" "<<tempsites[i].x<<" "<<tempsites[i].y<<" "<<(sqrt(pow(tempsites[i].x - tempsites[tempsites.size()-1].x,2)+pow(tempsites[i].y - tempsites[tempsites.size()-1].y,2)))<<endl;
		
	harmavdist = totalNeighbors / float(harmavdist);
	arthavdist = arthavdist / float(totalNeighbors);
	geomavdist = pow(geomavdist, 1.0/float(totalNeighbors));
	
	averages[0] = harmavdist;
	averages[1] = arthavdist;
	averages[2] = geomavdist;
	averages[4] = double(findvein(tempsites[NNindx], veinroot)->order);
	averages[5] = double(tempsites[NNindx].x);
	averages[6] = double(tempsites[NNindx].y);
	averages[7] = double(findvein(p, veinroot2)->width);
	averages[8] = double(findvein(tempsites[NNindx], veinroot)->width);

	return averages;//returns harmonic mean
}



void updateveingrowth(veinnode *ptr, double unitymod){
	if(!ptr) return;
	point r1,lprime;
	if(squarehex){
		r1.x=0;
		r1.y=0;
	}else{
		r1.x=1;
		r1.y=0;
	}
	double modulus;
	lprime.x = lstepgrowth*(r1.x-center.x)/(unitymod);
	lprime.y = lstepgrowth*(r1.y-center.y)/(unitymod);
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			updateveingrowth(ptr->next[i],unitymod);
		}
	}
	modulus = sqrt(pow(ptr->p.x-center.x,2)+pow(ptr->p.y-center.y,2));
	ptr->p.x=rounding(ptr->p.x+lstepgrowth*((ptr->p.x-center.x)/(unitymod)) - lprime.x,12);
	ptr->p.y=rounding(ptr->p.y+lstepgrowth*((ptr->p.y-center.y)/(unitymod)) - lprime.y,12);
}


void printvein(veinnode *ptr){
	MTRand mt(seed); 
   	double a,b,c,d,col;
   	bool removemin(0),flagrem(0);
//	cout<<"DEBUG printvein 1"<<endl;
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			removemin=0;flagrem=0;
//			if(!ptr->next[i]->next.empty()){ if(ptr->next[i]->next[0]->next.empty() && !ptr->next[i]->next[1]) flagrem=1;}
//			if(v2 && (ptr->next[i]->next.empty() || flagrem) && removeendchannels) removemin=1;
			//print line
			if(!removemin){
			if(v2){
			if(!veinintersection(ptr->p, ptr->next[i]->p, veinroot)){
			printvein(ptr->next[i]);	
			col=1.0*ptr->next[i]->order;
   			glLoadIdentity();             // Reset the projection matrix
			glTranslatef(moveh,movev,0.0);
		   	a = (2*ptr->p.x/(center.y*2)-1)*scalefactor;
		   	b = (2*ptr->p.y/(center.y*2)-1)*scalefactor;
		   	c = (2*ptr->next[i]->p.x/(center.y*2)-1)*scalefactor;
		   	d = (2*ptr->next[i]->p.y/(center.y*2)-1)*scalefactor;
//			cout<<"vein width: "<<ptr->next[i]->width<<endl;
   			glLineWidth(ptr->next[i]->width);
		   	gl2psLineWidth(ptr->next[i]->width);
			glBegin(GL_LINES);
			if(printorder)	glColor3d(pow(cos(50*col),5),pow(sin(200*col),10), pow(cos(55*col),8));
	//	      	glColor3f(0.0f, 0.2f, 0.0f);  // Red
		      	glVertex2d(a, b);       // Center of circle
		      	glVertex2d(c, d);       // Center of circle
		   	glEnd();
			}
			}else{
			printvein(ptr->next[i]);
			//if(ptr->next[i]->order>4) ptr->next[i]->order=4;
			col=1.0*ptr->next[i]->order;
			glLoadIdentity();             // Reset the projection matrix
			glTranslatef(moveh,movev,0.0);
		   	a = (2*ptr->p.x/(center.y*2)-1)*scalefactor;
		   	b = (2*ptr->p.y/(center.y*2)-1)*scalefactor;
		   	c = (2*ptr->next[i]->p.x/(center.y*2)-1)*scalefactor;
		   	d = (2*ptr->next[i]->p.y/(center.y*2)-1)*scalefactor;
//			cout<<"vein width: "<<ptr->next[i]->width<<endl;
   			glLineWidth(ptr->next[i]->width);
		   	gl2psLineWidth(ptr->next[i]->width);
			glBegin(GL_LINES);
	//	      	glColor3f(0.0f, 0.2f, 0.0f);  // Red
			if(printorder)	glColor3d(pow(cos(50*col),5),pow(sin(200*col),10), pow(cos(55*col),8));
		      	glVertex2d(a, b);       // Center of circle
		      	glVertex2d(c, d);       // Center of circle
		   	glEnd();
			}
			}
		}
	}
}
/*
void printveinout2(veinnode *ptr){
	double length, angle;
	point center2;
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			printveinout(ptr->next[i]);
			//print line
			if(!ptr->width != minwidth && (ptr->next[i]->width != minwidth)){
			angle = atan2(ptr->next[i]->p.y - ptr->p.y,ptr->next[i]->p.x - ptr->p.x);
			length = sqrt(pow(ptr->next[i]->p.x - ptr->p.x,2)+pow(ptr->next[i]->p.y - ptr->p.y,2));
			center2.x = ptr->p.x + (ptr->next[i]->p.x - ptr->p.x)/2;
			center2.y = ptr->p.y + (ptr->next[i]->p.y - ptr->p.y)/2;
			if(ptr->width - ptr->next[i]->width < 4 && ptr->next[i]->width > 0.1) {
			output1<<"hull(){"<<endl;
			output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0) << ", r1 = " <<  ptr->width/2  << ", r2 = "<<ptr->next[i]->width/2<<", center = true, $fn = 18);" <<endl;
			output1<< " translate([" <<ptr->next[i]->p.x <<","<<ptr->next[i]->p.y<<",0]) sphere(r = "<< ptr->next[i]->width/2<< ", $fn=25);"<<endl<<"}"<<endl;
			}else{
			output1<<"hull(){"<<endl;
			output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0) << ", r = " << ptr->next[i]->width/2 <<", center = true,  $fn=18);" <<endl;
			output1<< " translate([" <<ptr->next[i]->p.x<< ","<<ptr->next[i]->p.y<<",0]) sphere(r = "<< ptr->next[i]->width/2<< ", $fn=25);"<<endl<<"}"<<endl;
			}
			}
		}
	}
}*/
void printveinout(veinnode *ptr){
	double length, length2, angle, angle2;
	point center2, center3;
	bool removemin(0),flagrem;
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			printveinout(ptr->next[i]);
			angle = atan2(ptr->next[i]->p.y - ptr->p.y,ptr->next[i]->p.x - ptr->p.x);
			length = sqrt(pow(ptr->next[i]->p.x - ptr->p.x,2)+pow(ptr->next[i]->p.y - ptr->p.y,2));
			center2.x = ptr->p.x + (ptr->next[i]->p.x - ptr->p.x)/2; center2.y = ptr->p.y + (ptr->next[i]->p.y - ptr->p.y)/2;
			removemin=0; flagrem=0;
			if(!ptr->next[i]->next.empty()){ if(ptr->next[i]->next[0]->next.empty() && !ptr->next[i]->next[1]) flagrem=1;}
			if(v2 && (ptr->next[i]->next.empty() || flagrem) && removeendchannels) removemin=1;
			if(!removemin){
				if(ptr->width - ptr->next[i]->width < 4 && ptr->next[i]->width > 0.1) {
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0)+0.002 << ", r1 = " <<  ptr->width/2  << ", r2 = "<<ptr->next[i]->width/2+0.001<<", center = true, $fn = 15);" <<endl;
				}else{
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0)+0.002 << ", r = " << ptr->next[i]->width/2  +0.001<<", center = true,  $fn = 15);" <<endl;
				}
			for(int j=0;j<ptr->next[i]->next.size();j++){
				if(ptr->next[i]->next[j]){
				removemin=0;
				if(v2 && (ptr->next[i]->next[j]->next.empty()) && removeendchannels) removemin=1;
				if(!removemin){
				if(&ptr->next[i]->next[j]->p.y){
				//print line
				angle2 = atan2(ptr->next[i]->next[j]->p.y - ptr->next[i]->p.y,ptr->next[i]->next[j]->p.x - ptr->next[i]->p.x);
				length2 = sqrt(pow(ptr->next[i]->next[j]->p.x - ptr->next[i]->p.x,2)+pow(ptr->next[i]->next[j]->p.y - ptr->next[i]->p.y,2));
				center3.x = ptr->next[i]->p.x + (ptr->next[i]->next[j]->p.x - ptr->next[i]->p.x)/2; center3.y = ptr->next[i]->p.y + (ptr->next[i]->next[j]->p.y - ptr->next[i]->p.y)/2;
				if(ptr->width - ptr->next[i]->width < 4 && ptr->next[i]->width > 0.1) {
					output1<<"difference(){"<<endl;
					output1<<"intersection(){"<<endl;
				//	output1<<"union(){"<<endl;
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0)+10 << ", r = " <<  ptr->next[i]->width/2  <<", center = true, $fn = 15);" <<endl;
				//	output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0)+1 << ", r = " <<  ptr->next[i]->width/2  <<", center = true, $fn = 15);" <<endl;
				//	output1<<"}"<<endl;
					output1<<"hull(){"<<endl;
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0) << ", r = " <<  ptr->next[i]->width/2  <<", center = true, $fn = 15);" <<endl;
					output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0) << ", r = " <<  ptr->next[i]->width/2  <<", center = true, $fn = 15);" <<endl;
					output1<<"}"<<endl;
					output1<<"}"<<endl;
					output1<<"union(){"<<endl;
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0)+0.001  << ", r = " << ptr->next[i]->width/2+0.001<<", center = true, $fn = 15);" <<endl;
					output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0)+0.001 << ", r = " <<  ptr->next[i]->width/2 +0.001 << ", center = true, $fn = 15);" <<endl;
					output1<<"}"<<endl;
					output1<<"}"<<endl;
	//				output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0)+0.002 << ", r1 = " <<  ptr->next[i]->width/2 +0.001 << ", r2 = "<<ptr->next[i]->next[j]->width/2<<", center = true, $fn = 15);" <<endl;
				}else{
					output1<<"difference(){"<<endl;
					output1<<"intersection(){"<<endl;
				//	output1<<"union(){"<<endl;
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0)+10 << ", r = " <<  ptr->next[i]->width/2  <<", center = true, $fn = 15);" <<endl;
				//	output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0)+1 << ", r = " <<  ptr->next[i]->next[j]->width/2  <<", center = true, $fn = 15);" <<endl;
				//	output1<<"}"<<endl;
					output1<<"hull(){"<<endl;
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0) << ", r = " <<  ptr->next[i]->width/2  <<", center = true, $fn = 15);" <<endl;
					output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0) << ", r = " <<  ptr->next[i]->next[j]->width/2 <<", center = true, $fn = 15);" <<endl;
					output1<<"}"<<endl;
					output1<<"}"<<endl;
					output1<<"union(){"<<endl;
					output1 << "translate([" << center2.x << "," << center2.y << ",0]) rotate([0,90," << (angle * 180)/PI << "]) cylinder (h = " << length*(1+0)+0.001  << ", r = " << ptr->next[i]->width/2+0.001<<", center = true, $fn = 15);" <<endl;
					output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0)+0.001 << ", r = " <<  ptr->next[i]->next[j]->width/2 +0.001 << ", center = true, $fn = 15);" <<endl;
					output1<<"}"<<endl;
					output1<<"}"<<endl;
	//				output1 << "translate([" << center3.x << "," << center3.y << ",0]) rotate([0,90," << (angle2 * 180)/PI << "]) cylinder (h = " << length2*(1+0)+0.002 << ", r = " << ptr->next[i]->next[j]->width/2  +0.001<<", center = true,  $fn = 15);" <<endl;
	//				output1<< " translate([" <<ptr->next[i]->p.x<< ","<<ptr->next[i]->p.y<<",0]) sphere(r = "<< ptr->next[i]->width/2<< ", $fn=25);"<<endl<<"}"<<endl;
				}
				}
				}
			}
			}
			}
		}
	}
}

void redefineveinwidth(veinnode *ptr){
	//if(!ptr->next.empty()){
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->width!=minwidth){
			redefineveinwidth(ptr->next[i]);
			//ptr->width=ptr->width+pow(ptr->next[i]->width,1);
		}else{
			ptr=0;
		}
	}
//	cout<<"width defined before "<<ptr->width<<endl;
//	ptr->width=pow(ptr->width,1);
//	cout<<"width defined not empty "<<ptr->width<<endl;
	//}else{
//		ptr->width=minwidth;
//	cout<<"width defined empty "<<ptr->width<<endl;
//	}
}

void veinlengtharea(veinnode *ptr,int ordem){
	bool removemin(0);
	if(!ptr->next.empty()){
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			removemin=0;
			if(v2 && ptr->next[i]->next.empty() && removeendchannels) removemin=1;
			if(!removemin){
			veinlengtharea(ptr->next[i],ordem);
			if(ptr->order==ordem||ordem==0) venationlength = venationlength + sqrt(pow(ptr->next[i]->p.x - ptr->p.x,2)+pow(ptr->next[i]->p.y - ptr->p.y,2));
			}
		}
	}
	}
}


//deletes tree branches with order beyond cutoff
void deltreebeyondcutoff(veinnode *ptr, int order){
	if(ptr->order >= order){
		if(!ptr->next.empty()){
			for(int i=0;i<ptr->next.size();i++){
				if(ptr->next[i]){
					deltreebeyondcutoff(ptr->next[i], order);
				}
			}
			
			ptr->next.clear();
			vector<veinnode*> ().swap(ptr->next);
			delete ptr;
		}else{
			delete ptr;
		}
	}else{
		if(!ptr->next.empty()){
			for(int i=0;i<ptr->next.size();i++){
				if(ptr->next[i]){
					if(ptr->next[i]->order >= order){ deltreebeyondcutoff(ptr->next[i], order); ptr->next[i]=0; }else{ deltreebeyondcutoff(ptr->next[i], order); }
				}
			}
		}
	}
}


void bifurcationcountmethod(veinnode *ptr){
	if(!ptr->next.empty()){
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			bifurcationcountmethod(ptr->next[i]);
			if(i!=0) bifurcationcount++;
		}
	}
	}
}

void defineveinwidth(veinnode *ptr){
	if(!ptr->next.empty()){
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			defineveinwidth(ptr->next[i]);
			ptr->width=ptr->width + pow(ptr->next[i]->width,murrayexp);
		}
	}
	ptr->width=pow(ptr->width,(1.0/murrayexp));
	}else{
		ptr->width=minwidth;
	}
}

void saveVenationStructure(veinnode *ptr, double disttoroot){
	double tempdist(0);
	if(!ptr->next.empty()){
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			tempdist = sqrt(pow(ptr->p.x - ptr->next[i]->p.y,2)+pow(ptr->p.y - ptr->next[i]->p.y,2));
			tempdist += disttoroot;
			saveVenationStructure(ptr->next[i], tempdist);
			output1 << ptr->p.x<<'\t'<<ptr->p.y<<'\t'<<ptr->next[i]->p.x<<'\t'<<ptr->next[i]->p.y<<'\t'<<tempdist<<endl;
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void defineveinorder(veinnode *ptr,int ordem){
	bool veininter(0),flag(0),flag2(0);
//	int veinordergt(0);
	int contando(0);
	int greater;
	if(v2&&ordem<3) {
//		cout<<"width comparison"<<veinroot2->width<<" "<<ptr->width<<endl;
		if(veinroot2->width-ptr->width < veinroot2->width*0.35){ ordem=1;}
//		if((ptr->prev->width-ptr->width)<ptr->prev->width*0.05) ordem--;
	//	if((ptr->prev->width-ptr->width)<ptr->prev->width*0.95) ordem++;
		//if((ptr->prev->width-ptr->width)<ptr->prev->width*0.99) ordem++;
	}
	if(!ptr->next.empty()){
//	if(ptr->next.size()>1){
	while(ptr->next.size()!=0){
	ptr->order=ordem;
	flag=0;
	contando++;
//	cout<<"contando "<<contando<<endl;
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
	//		defineveinorder(ptr->next[i]);
			if(!flag){ greater = i;} flag=1;
			if(ptr->next[greater]->width < ptr->next[i]->width){ defineveinorder(ptr->next[greater],ordem+1); greater = i; }else{if(greater!=i)defineveinorder( ptr->next[i],ordem+1);}
			//if(veinordergt < ptr->next[i]->order) veinordergt = ptr->next[i]->order;
	//		cout<<"width: "<< ptr->next[i]->width<<endl<<" order "<<ptr->next[i]->order<<endl;
		}
	}
	ptr=ptr->next[greater];
	}
	}
	ptr->order=ordem;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void murrayfrombase(veinnode *ptr,double diameter){
	if(!ptr->next.empty()){
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
				murrayfrombase(ptr->next[i],(diameter*(ptr->next[i]->width/ptr->width)));
		}
	}
	}
	if(diameter>minwidth)
	ptr->width=diameter;
	else
	ptr->width=minwidth;
}

void display();


void openscadoutput(){
   double dim(center.y*2*(1.1)),length, angle;
      cout<<dim<<endl;
      cout<<dim<<endl;
	    output1.open("main.scad");
	    printf("Writing 'main.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"difference() {"<< endl<<"translate([" << center.x <<","<< center.y<<",0]) cube(["<<finallength<<","<< finallength <<","<< initialdiam + 1 <<"], center = true);"<<endl<<"union() {"<< endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x));
		length = 2*sqrt(pow(-veinroot->p.y/angle,2)+pow(-veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-veinroot->p.y/angle<< ","<<0<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	    output1<<"//second venation starts here"<<endl;
	    printveinout(veinroot2);
		if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y+((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x));
		length = 2*sqrt(pow((finallength-veinroot2->p.y)/angle,2)+pow((finallength-veinroot2->p.y),2));
		output1<< " translate([" <<veinroot2->p.x+(finallength-veinroot2->p.y)/angle<< ","<<finallength<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r=" <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	   
	    //output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y<<",0]) sphere(r = "<< veinroot2->width/2<< ", $fn=15);"<<endl;
	    output1 <<"}"<<endl<<"}"<<endl;
	    output1.close();
		printf("Done!\n");
		output1.open("venation1.scad");
	    printf("Writing 'venation1.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"union() {"<< endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x));
		length = 2*sqrt(pow(-veinroot->p.y/angle,2)+pow(-veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-veinroot->p.y/angle<< ","<<0<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	    output1<<"}"<<endl;
	    output1.close();
		printf("Done!\n");
		output1.open("venation2.scad");
	    printf("Writing 'venation2.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"union() {"<< endl;
	    printveinout(veinroot2);
		if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y+((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x));
		length = 2*sqrt(pow((finallength-veinroot2->p.y)/angle,2)+pow((finallength-veinroot2->p.y),2));
		output1<< " translate([" <<veinroot2->p.x+(finallength-veinroot2->p.y)/angle<< ","<<finallength<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r=" <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    output1 <<"}"<<endl;
	    output1.close();
		printf("Done!\n");
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~OPEN GL FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void keyboard(unsigned char key, int x, int y) {
   FILE *fp;
   int state=GL2PS_OVERFLOW, buffsize=0;
   double dim(center.y*2*(1.1)),length, angle;
	
   switch (key) {
      case 27:     // ESC key
         exit(0);
         break;
      case 'b':
      cout<<dim<<endl;
      cout<<dim<<endl;
	    output1.open("main.scad");
	    printf("Writing 'main.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"difference() {"<< endl<<"translate([" << center.x <<","<< center.y<<",0]) cube(["<<finallength<<","<< finallength <<","<< initialdiam + 1 <<"], center = true);"<<endl<<"union() {"<< endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x));
		length = 2*sqrt(pow(-veinroot->p.y/angle,2)+pow(-veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-veinroot->p.y/angle<< ","<<0<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	    output1<<"//second venation starts here"<<endl;
	    printveinout(veinroot2);
		if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y+((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x));
		length = 2*sqrt(pow((finallength-veinroot2->p.y)/angle,2)+pow((finallength-veinroot2->p.y),2));
		output1<< " translate([" <<veinroot2->p.x+(finallength-veinroot2->p.y)/angle<< ","<<finallength<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r=" <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	   
	    //output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y<<",0]) sphere(r = "<< veinroot2->width/2<< ", $fn=15);"<<endl;
	    output1 <<"}"<<endl<<"}"<<endl;
	    output1.close();
		printf("Done!\n");
		output1.open("venation1.scad");
	    printf("Writing 'venation1.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"union() {"<< endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x));
		length = 2*sqrt(pow(-veinroot->p.y/angle,2)+pow(-veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-veinroot->p.y/angle<< ","<<0<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	    output1<<"}"<<endl;
	    output1.close();
		printf("Done!\n");
		output1.open("venation2.scad");
	    printf("Writing 'venation2.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"union() {"<< endl;
	    printveinout(veinroot2);
		if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y+((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = tan(atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x));
		length = 2*sqrt(pow((finallength-veinroot2->p.y)/angle,2)+pow((finallength-veinroot2->p.y),2));
		output1<< " translate([" <<veinroot2->p.x+(finallength-veinroot2->p.y)/angle<< ","<<finallength<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << length << ", r=" <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    output1 <<"}"<<endl;
	    output1.close();
		printf("Done!\n");
	    break;
	  case 'd':
      cout<<dim<<endl;
      cout<<dim<<endl;
	    output1.open("main.scad");
	    printf("Writing 'main.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"difference() {"<< endl<<"translate([" << center.x <<","<< center.y<<",0]) cube(["<<finallength<<","<< finallength <<","<< initialdiam + 1 <<"], center = true);"<<endl<<"union() {"<< endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x);
		length = 2*sqrt(pow(veinroot->next[0]->p.x - veinroot->p.x,2)+pow(veinroot->next[0]->p.y - veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-length*sin(angle)<< ","<<veinroot->p.y-length*cos(angle)<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << 3*length << ", r=" <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	    output1<<"//second venation starts here"<<endl;
	    printveinout(veinroot2);
	    if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y+((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x);
		length = 2*sqrt(pow(veinroot2->next[0]->p.x - veinroot2->p.x,2)+pow(veinroot2->next[0]->p.y - veinroot2->p.y,2));
		output1<< " translate([" <<veinroot2->p.x-length*sin(angle)<< ","<<veinroot2->p.y-length*cos(angle)<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << 3*length << ", r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    //output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y<<",0]) sphere(r = "<< veinroot2->width/2<< ", $fn=15);"<<endl;
	    output1 <<"}"<<endl<<"translate([" << center.x <<","<< center.y<<","<<(initialdiam+1)/2 <<"]) cube(["<<finallength<<","<< finallength <<","<< initialdiam + 1 <<"], center = true);"<<endl<<"}"<<endl;
	    output1.close();
		printf("Done!\n");
		output1.open("venation1.scad");
	    printf("Writing 'venation1.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"union() {"<< endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x);
		length = 2*sqrt(pow(veinroot->next[0]->p.x - veinroot->p.x,2)+pow(veinroot->next[0]->p.y - veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-length*sin(angle)<< ","<<veinroot->p.y-length*cos(angle)<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << 3*length << ", r= " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	    output1<<"}"<<endl;
	    output1.close();
		printf("Done!\n");
		output1.open("venation2.scad");
	    printf("Writing 'venation2.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"union() {"<< endl;
	    printveinout(veinroot2);
	    if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y+((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r= " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x);
		length = 2*sqrt(pow(veinroot2->next[0]->p.x - veinroot2->p.x,2)+pow(veinroot2->next[0]->p.y - veinroot2->p.y,2));
		output1<< " translate([" <<veinroot2->p.x-length*sin(angle)<< ","<<veinroot2->p.y-length*cos(angle)<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << 3*length << ", r=" <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    output1 <<"}"<<endl;
	    output1.close();
		printf("Done!\n");
	    break;
	  case 'c':
	    output1.open("output2.scad");
	    printf("Writing 'output.scad'... ");
	    output1 <<"translate([" << -center.x <<","<< -center.y<<",0])"<<endl;
	    output1 <<"difference() {"<< endl <<"union() {"<< endl << "translate([" << center.x <<","<< center.y<<",-0.6]) cube(["<<finallength+10<<","<<finallength+10<<",0.9], center = true);"<<endl;
		if(!diagbool){
		output1<< " translate([" <<veinroot->p.x<< ","<<veinroot->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ", r1=0.375, r2 = " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = atan2(veinroot->next[0]->p.y - veinroot->p.y,veinroot->next[0]->p.x - veinroot->p.x);
		length = 2*sqrt(pow(veinroot->next[0]->p.x - veinroot->p.x,2)+pow(veinroot->next[0]->p.y - veinroot->p.y,2));
		output1<< " translate([" <<veinroot->p.x-length*sin(angle)<< ","<<veinroot->p.y-length*cos(angle)<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << 3*length << ", r1=0.375, r2 = " <<  veinroot->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    printveinout(veinroot);
	   // output1<<"//second venation starts here"<<endl;
	   // output1<< " translate([" <<veinroot3->p.x<< ","<<veinroot3->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r1=0.375, r2 = " <<  veinroot3->width/2 <<", center = true, $fn = 15);" <<endl;
	   // printveinout(veinroot3);
	    output1<<"//complementary venation starts here"<<endl;
	    printveinout(veinroot2);
	    if(!diagbool){
		output1<< " translate([" <<veinroot2->p.x<< ","<<veinroot2->p.y-((dim/2-center.y)/2)<<",0]) rotate([0,90,90]) cylinder (h = " << (dim/2-center.y)+1 << ",r2=0.375, r1 = " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}else{
		angle = atan2(veinroot2->next[0]->p.y - veinroot2->p.y,veinroot2->next[0]->p.x - veinroot2->p.x)+PI;
		length = 2*sqrt(pow(veinroot2->next[0]->p.x - veinroot2->p.x,2)+pow(veinroot2->next[0]->p.y - veinroot2->p.y,2));
		output1<< " translate([" <<veinroot2->p.x-length*sin(angle)<< ","<<veinroot2->p.y-length*cos(angle)<<",0]) rotate([0,90,"<<(angle * 180)/PI<<"]) cylinder (h = " << 3*length << ", r1=0.375, r2 = " <<  veinroot2->width/2 <<", center = true, $fn = 15);" <<endl;
		}
	    output1 <<"}"<<endl<< "translate([" << center.x <<","<< center.y<<",-3]) cube(["<<dim*1.0001<<","<<dim*1.0001<<",3], center = true);"<<endl<<"}"<<endl;
	    output1.close();
	    printf("Done!\n");
	    break;
      case 's':
       fp = fopen("out.eps", "wb");
       printf("Writing 'out.eps'... ");
       while(state == GL2PS_OVERFLOW){
         buffsize += 1024*1024;
         gl2psBeginPage("leafvenation", "Caio Martins", NULL, GL2PS_EPS, GL2PS_SIMPLE_SORT,
                     GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT,
                     GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, "out.eps");
         display();
         state = gl2psEndPage();
       }
       fclose(fp);
       printf("Done!\n");
       break;
    }
}


double voronoicellarea(hedge *orgn){
	double area(0);
	for(hedge *i=orgn->next;i!=orgn->prev;i=i->next){
		area=area+crossproduct2d(orgn->origin->p,i->origin->p,i->next->origin->p);
	}
	return area;
}


void initGL() {
   // Set "clearing" or background color
   glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Black and opaque
}


void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
   // Compute aspect ratio of the new window
   winsy=height; winsx=width;
   if (height == 0) height = 1;                // To prevent divide by 0
   GLfloat aspect = (GLfloat)width / (GLfloat)height;
 
   // Set the viewport to cover the new window
   glViewport(0, 0, width, height);
 
   // Set the aspect ratio of the clipping area to match the viewport
   glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
   glLoadIdentity();             // Reset the projection matrix
   if (width >= height) {
     // aspect >= 1, set the height from -1 to 1, with larger width
      gluOrtho2D(-1.0 * aspect, 1.0 * aspect, -1.0, 1.0);
   } else {
      // aspect < 1, set the width to -1 to 1, with larger height
     gluOrtho2D(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect);
   }
}

void idle() {
   glutPostRedisplay();   // Post a re-paint request to activate display()
}
 


void display(){
   glClear(GL_COLOR_BUFFER_BIT);
   glMatrixMode(GL_MODELVIEW);    // To operate on the model-view matrix
   glLoadIdentity();              // Reset model-view matrix
   point p,p0,p1;
   int numSegments = 100,tot,i,j;
   GLdouble angle;
   double a,b,c,d;
     hedge *temp1,*temp2;
   //-----------------------------------------COLORING FACES (ONLY IF REQUIRED)--------------------------//
/*   hedge *temp1,*temp2;
   for(j=0;j<dcefacelist.size()-1;j++){
//   cout<<"step1"<<endl;
   temp1=dcefacelist.at(j)->edge;
   temp2=temp1->next;
   p=temp1->origin->p;
   a = (2*p.x/(center.y*2)-1)*scalefactor;
   b = (2*p.y/(center.y*2)-1)*scalefactor;
  // cout << "Pontos: " << a << " " << b << endl;
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
//   glTranslatef(a,b,0.0);
//   cout<<"step2"<<endl;
   glBegin(GL_TRIANGLE_FAN);
   col=1.0*j/dcefacelist.size();
  // cout<<"step21"<<endl;
  // cout<<"COLOR: "<< j <<" " <<col<<endl;
   glColor3d(pow(cos(50*col),5),pow(sin(100*col),2), pow(cos(55*col),8));  // Red
  // myfile<<output.at(j)->start.x<<" "<<output.at(j)->start.y<<" "<<output.at(j)->end.x<<" "<<output.at(j)->end.y<<endl;
      glVertex2d(a, b);       // Center of circle
 //  cout<<"step22"<<endl;
      while(temp2!=temp1) { // Last vertex same as first vertex
  // cout<<"step23"<<endl;
//	cout<<"infinite loop"<<endl;   
	 p=temp2->origin->p;
//   cout<<"step231"<<endl;
   a = (2*p.x/(center.y*2)-1)*scalefactor;
   b = (2*p.y/(center.y*2)-1)*scalefactor;
   hedge *temp1,*temp2;ndl;
         glVertex2d(a, b);
 //  cout<<"step233"<<endl;
	 temp2=temp2->next;
 //  cout<<"step24"<<endl;
      }
  // cout<<"step3"<<endl;
   glEnd();
   }
*/   
   if(squarehex){ ////////////////////////////////////////HEXAGON/////////////////////////////////////////////
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
   glBegin(GL_TRIANGLE_FAN);
   glColor3d(0.0,0.9, 0.0);  // dark gray
   p0=hexvert[0];
   a = (2*p0.x/(center.y*2)-1)*scalefactor;
   b = (2*p0.y/(center.y*2)-1)*scalefactor;
   glVertex2d(a, b);       // Center of circle
   for(int j=1;j<6;j++){
      p1=hexvert[j];
      c = (2*p1.x/(center.y*2)-1)*scalefactor;
      d = (2*p1.y/(center.y*2)-1)*scalefactor;
      glVertex2d(c, d);       // Center of circle
   }
   glEnd();
   ///////////////////////////
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
   glLineWidth(3);
   gl2psLineWidth(3);
   glBegin(GL_LINE_LOOP);
   p0=hexvert[0];
   a = (2*p0.x/(center.y*2)-1)*scalefactor;
   b = (2*p0.y/(center.y*2)-1)*scalefactor;
      glColor3f(0.0f, 0.2f, 0.0f);  // Red
   glVertex2d(a, b);       // Center of circle
   for(int j=1;j<6;j++){
      p1=hexvert[j];
      c = (2*p1.x/(center.y*2)-1)*scalefactor;
      d = (2*p1.y/(center.y*2)-1)*scalefactor;
      glColor3f(0.0f, 0.2f, 0.0f);  // Red
      glVertex2d(c, d);       // Center of circle
   }
   glEnd();
   }else{//////////////////////////////SQUARE/////////////////////////////
   if(!printorder){
   if(!background){
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
   glBegin(GL_TRIANGLE_FAN);
   glColor3d(0.0,0.9, 0.0);  // dark gray
   p0=squarevert[0];
   a = (2*p0.x/(center.y*2)-1)*scalefactor;
   b = (2*p0.y/(center.y*2)-1)*scalefactor;
   glVertex2d(a, b);       // Center of circle
   for(int j=1;j<6;j++){
      p1=squarevert[j];
      c = (2*p1.x/(center.y*2)-1)*scalefactor;
      d = (2*p1.y/(center.y*2)-1)*scalefactor;
      glVertex2d(c, d);       // Center of circle
   }
   glEnd();
   }
   }
   ///////////////////////////
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
   glLineWidth(3);
   gl2psLineWidth(3);
 
  glBegin(GL_LINE_LOOP);
   p0=squarevert[0];
   a = (2*p0.x/(center.y*2)-1)*scalefactor;
   b = (2*p0.y/(center.y*2)-1)*scalefactor;
      glColor3f(0.0f, 0.1f, 0.0f);  // Red
   glVertex2d(a, b);       // Center of circle
   for(int j=1;j<6;j++){
      p1=squarevert[j];
      c = (2*p1.x/(center.y*2)-1)*scalefactor;
      d = (2*p1.y/(center.y*2)-1)*scalefactor;
      glColor3f(0.0f, 0.1f, 0.0f);  // Red
      glVertex2d(c, d);       // Center of circle
   }
   glEnd();   

   }
   v2=0;
   printvein(veinroot); 
   //v2=1;
         glColor3f(1.0f, 0.0f, 0.0f);  // Red
   printvein(veinroot2);

  /*for(j=0;j<dcefacelist.size();j++)
  // cout<<"step4"<<endl;
//  if(j==9){
   temp1=dcefacelist.at(j)->edge;
   temp2=temp1->next;
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
   glBegin(GL_LINE_LOOP);
   p0=temp1->origin->p;
   a = (2*p0.x/20-1)*scalefactor;
   b = (2*p0.y/20-1)*scalefactor;
   glVertex2d(a, b);       // Center of circle
   while(temp2!=temp1){
      p1=temp2->origin->p;
      c = (2*p1.x/20-1)*scalefactor;
      d = (2*p1.y/20-1)*scalefactor;
      glColor3f(0.0f, 0.0f, 0.0f);  // Red
      glVertex2d(c, d);       // Center of circle
      temp2=temp2->next;
   }
   glEnd();
   }
*/
  //------------------------------------------DRAWING POINTS---------------------------------------------//
 /*     for(j=0;j<veinpoints.size();j++){
  // cout<<"step4"<<endl;
//  if(j==9){

   temp1=dcefacelist.at(j)->edge;
   temp2=temp1->next;
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(moveh,movev,0.0);
   glBegin(GL_LINE_LOOP);
   p0=temp1->origin->p;
   a = (2*p0.x/(center.x*2)-1)*scalefactor+0.2;
   b = (2*p0.y/(center.y*2)-1)*scalefactor;
      glColor3f(0.0f, 0.0f, 0.0f);  // Red
   glVertex2d(a, b);       // Center of circle
   while(temp2!=temp1){
      p1=temp2->origin->p;
      c = (2*p1.x/(center.x*2)-1)*scalefactor+0.2;
      d = (2*p1.y/(center.y*2)-1)*scalefactor;
      glColor3f(0.0f, 0.0f, 0.0f);  // Red
      glVertex2d(c, d);       // Center of circle
      temp2=temp2->next;
   }
   glEnd();*/
/* 
 p=veinpoints[j]; 
   a = (2*p.x/(center.y*2)-1)*scalefactor;
   b = (2*p.y/(center.y*2)-1)*scalefactor;
  // cout << "Pontos: " << a << " " << b << endl;
   glLoadIdentity();             // Reset the projection matrix
   glTranslatef(a+moveh,b+movev,0.0);
   glBegin(GL_TRIANGLE_FAN);
//   if(j!=points2.size()-1){
	 glColor3f(1.0f, 0.0f, 0.0f);  // Red
//}else{
//	 glColor3f(1.0f, 0.0f, 0.0f);  // Red
//}      
glVertex2d(0,0);       // Center of circle
      for (int i = 0; i <= numSegments; i++) { // Last vertex same as first vertex
         angle = i * PI * 2.0 / numSegments;  // 360 deg for all segments
         glVertex2d(cos(angle) * pointRadius * 1.0, sin(angle) * pointRadius * 1.0);
      }
   glEnd();
*/
 // }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
   glutSwapBuffers();   // Double buffered - swap the front and back buffers
}

void specialKeys(int key, int x, int y) {
   switch (key) {
      case GLUT_KEY_RIGHT:    // Right: increase x speed
         moveh -= 1.5f; break;
      case GLUT_KEY_LEFT:     // Left: decrease x speed
         moveh += 1.5f; break;
      case GLUT_KEY_UP:       // Up: increase y speed
         movev -= 1.5f; break;
      case GLUT_KEY_DOWN:     // Down: decrease y speed
         movev += 1.5f; break;
      case GLUT_KEY_PAGE_UP:  // Page-Up: increase ball's radius
         scalefactor *= 1.1f;
         break;
      case GLUT_KEY_PAGE_DOWN: // Page-Down: decrease ball's radius
         scalefactor *= 0.9f;
         break;
      case GLUT_KEY_HOME:
	 scalefactor *=1.01f; break;
      case GLUT_KEY_END:
	 scalefactor *=0.99f; break;
   }
}

// Callback handler for mouse event
void mouse(int button, int state, int xcoord, int ycoord) {
   if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) { // Pause/resume
	point query;
	double xcoord2,ycoord2;
   if (winsy == 0) winsy = 1;                // To prevent divide by 0
   double aspect(winsx/winsy);
	xcoord2=2.0*xcoord/winsx-1;  ycoord2=(2.0*ycoord/winsy-1)*(-1);
   if (winsx >= winsy) {
        query.x=center.y*(aspect*xcoord2/scalefactor+1); query.y=center.y*(ycoord2/scalefactor+1);
   } else {
        query.x=center.y*(xcoord2/scalefactor+1); query.y=center.y*(ycoord2/(aspect*scalefactor)+1);
   }
	cout<<"The query point is: "<<query.x<<" "<<query.y<<endl;
//	queryface=plquery(query);
//	cout<<"The Voronoi site of this is: "<< queryface->site.x<< " "<<queryface->site.y<<endl;
//	deletesite2(query,true);  
  // 	cout<<"size of the face array: "<<dcefacelist.size()<<endl; 
   }
   if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) { // Pause/resume
//	point query;
//	double xcoord2,ycoord2;
 // if (winsy == 0) winsy = 1;                // To prevent divide by 0
//   double aspect(winsx/winsy);
//	xcoord2=2.0*xcoord/winsx-1;  ycoord2=(2.0*ycoord/winsy-1)*(-1);
 //  if (winsx >= winsy) {
 //       query.x=10*(aspect*xcoord2/scalefactor+1); query.y=10*(ycoord2/scalefactor+1);
  // } else {
 //       query.x=10*(xcoord2/scalefactor+1); query.y=10*(ycoord2/(aspect*scalefactor)+1);
  // }
//	cout<<"The query point is: "<<query.x<<" "<<query.y<<endl;
//	queryface=plquery(query);
//	cout<<"The Voronoi site of this is: "<< queryface->site.x<< " "<<queryface->site.y<<endl;
//	createsite2(query);  
  //	cout<<"size of the face array: "<<dcefacelist.size()<<endl; 
   growth(squarehex);
   }
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END OF OPENGL FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

/////////////////////////////Stores all neighbor points of the query into the tempsites array///////////////////////////////
void getnearestneighbors(hedge *temp0, bool clear){
	bool flagpoint;
	hedge *temp1(temp0->next), *temp2, *temp3;
	if(clear) tempsites.clear();
	while(temp1!=temp0){
		temp2=temp1->twin; temp3=temp2->next;
		if(temp2->f->site!=pproblem){
		flagpoint=true;
		for(int i=0;i<tempsites.size();i++){
			if(temp2->f->site==tempsites.at(i)||temp2->f->site==temp0->f->site) flagpoint=false;
		}
		if(flagpoint) tempsites.push_back(temp2->f->site);
		while(temp3!=temp2){
			if(temp3->twin->f->site!=pproblem){
			flagpoint=true;
			for(int i=0;i<tempsites.size();i++){
			if(temp3->twin->f->site==tempsites.at(i)||temp3->twin->f->site==temp0->f->site) flagpoint=false;
			}
			if(flagpoint) tempsites.push_back(temp3->twin->f->site);
			}
			temp3=temp3->next;
		}
		}
		temp1=temp1->next;
	}
	temp2=temp1->twin; temp3=temp2->next;
	if(temp2->f->site!=pproblem){
	flagpoint=true;
	for(int i=0;i<tempsites.size();i++){
		if(temp2->f->site==tempsites.at(i)||temp2->f->site==temp0->f->site) flagpoint=false;
	}
	if(flagpoint) tempsites.push_back(temp2->f->site);
	while(temp3!=temp2){
		if(temp3->twin->f->site!=pproblem){
		flagpoint=true;
		for(int i=0;i<tempsites.size();i++){
			if(temp3->twin->f->site==tempsites.at(i)||temp3->twin->f->site==temp0->f->site) flagpoint=false;
		}
		if(flagpoint) tempsites.push_back(temp3->twin->f->site);
		}
		temp3=temp3->next;
	}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////Delete site function//////////////////////////////////////////////////////////////////////////////////////////////////////
void deletesite2(point p, bool veinorauxin)
{
	queryface=plquery(p,1);
	if(veinorauxin){
	for(int i=0;i<veinpoints.size();i++){
		if(queryface->site==veinpoints[i]){
			veinpoints.erase(veinpoints.begin()+i);
			break;
		}
	}
	}else{
	for(int i=0;i<auxinpoints.size();i++){
		if(queryface->site==auxinpoints[i]){
			auxinpoints.erase(auxinpoints.begin()+i);
			break;
		}
	}
	}
	if(points.size()>0) while(points.size()!=0) points.pop();
	//////////////////////////////inserts all sites from tempsites array into the priority queue for subsquent voronoi diagram construction/////////////////////////////
	if(veinorauxin){
	for(int i=0;i<veinpoints.size();i++){
		points.push(veinpoints[i]);
	}
	}else{
	for(int i=0;i<auxinpoints.size();i++){
		points.push(auxinpoints[i]);
	}
	}
	cleanup();
	if(voronoistructure(veinorauxin, squarehex)) return;
}

bool roundingpoint(point p1, point p2, int n){
	if(rounding(p1.x,n)==rounding(p2.x,n)&&rounding(p1.y,n)==rounding(p2.y,n)){ 
		return true; 
	}else{ 
		return false;
	}
}

void createsite2(point p,bool veinorauxin)
{
	if(veinorauxin){
	for(int i=0;i<veinpoints.size();i++){
		if(roundingpoint(veinpoints[i], p, 6)) { cout<<"Site already exists. Please, choose different coordinates."<<endl; return; }
	}
	veinpoints.push_back(p);
	}else{
	for(int i=0;i<auxinpoints.size();i++){
		if(roundingpoint(auxinpoints[i], p, 6)) { cout<<"Site already exists. Please, choose different coordinates."<<endl; return; }
	}
	auxinpoints.push_back(p);
	}
	if(points.size()>0) while(points.size()!=0) points.pop();
	//////////////////////////////inserts all sites from tempsites array into the priority queue for subsquent voronoi diagram construction/////////////////////////////
	if(veinorauxin){
	for(int i=0;i<veinpoints.size();i++){
		points.push(veinpoints[i]);
	}
	}else{
	for(int i=0;i<auxinpoints.size();i++){
		points.push(auxinpoints[i]);
	}
	}
	cleanup();
	if(voronoistructure(veinorauxin, squarehex)) return;
//	ptr->index.clear();
//	vector<bool> ().swap(ptr->strorend);
}
//////////////////////////////////////////////////////////


void initiatesqgrid(int nb){
	double lgth((squarevert[1].x-squarevert[0].x)/nb);
	point tmp;
	int k,j;
	griddiv *temp;
	for(int i=0;i<nb*nb;i++){
		temp = new griddiv();
		grid.push_back(temp);
		k=(int) i/nb;
		j=i%nb;
		tmp.x=squarevert[0].x + j*lgth;
		tmp.y=squarevert[0].y + k*lgth;
		grid[i]->vetcs.push_back(tmp);
		tmp.x+=lgth;
		grid[i]->vetcs.push_back(tmp);
		tmp.x-=lgth;
		tmp.y+=lgth;
		grid[i]->vetcs.push_back(tmp);
		tmp.x+=lgth;
		grid[i]->vetcs.push_back(tmp);
		grid[i]->ok=1;
	}
}


void cleanup(){
	if(bstdcelroot)	deletebstdcel(bstdcelroot);
	bstdcelroot=0;
	vector<hedge*> deltemp;
	if(!dcefacelist.empty()){
	for(int i=0;i<dcefacelist.size();i++){
		deltemp.clear();
		for(hedge *j=dcefacelist[i]->edge;j!=dcefacelist[i]->edge->prev;j=j->next){
			deltemp.push_back(j);
		}
		deltemp.push_back(dcefacelist[i]->edge->prev);
		for(int j=0;j<deltemp.size();j++){
			delete deltemp[j];
			deltemp[j]=0;
		}
		vector<point> ().swap(dcefacelist[i]->auxin);
		delete dcefacelist[i];
		dcefacelist[i]=0;
	}
	dcefacelist.clear();
	}
	if(!dcelist.empty()){
	for(int i=0;i<dcelist.size();i++){
		delete dcelist[i];
		dcelist[i]=0;
	}
	dcelist.clear();
	}
	if(slabbstxroot) deleteslabbstx(slabbstxroot);
	slabbstxroot=0;   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void deletefaceedges(hedge *temp0, hedge *temp1){
	if(temp1!=temp0){deletefaceedges(temp0,temp1->next); }
	delete temp1;
	temp1=0;
}

void createslabstruct(bool root1vs2)
{
	slabbstx *current;
	hedge *temphedge, *tempref;
	slabybuild slabytemp;
	double a,b,xhalf,yhalf;
	int size;
	point pf,pi;
//	current=slabbstxroot;
	if(root1vs2){ current = slabbstxroot; }else{ current = slabbstxroot2; dcelist.swap(veinvtcs);}
	vector<slabybuild> slabhistory;
	while(current->left){
		current=current->left;
	}
	while(current){
		xhalf=(current->xcoordf-current->xcoord)/2+current->xcoord;
		size=slabhistory.size();
		for(int i=0;i<slabhistory.size();i++){
			if(slabhistory.at(i).xcoord==current->xcoord){ slabhistory.erase(slabhistory.begin()+i); i--;} //<==== try to modify this line ...  erase method is too slow
		}			
		for(int i=0;i<slabhistory.size();i++){
			yhalf=slabhistory.at(i).a*xhalf+slabhistory.at(i).b;
			slabbstyinsert(yhalf,&current->slaby,slabhistory.at(i).a,slabhistory.at(i).b,slabhistory.at(i).edge);	
		}
		for(int i=0;i<current->index.size();i++){	
			temphedge=tempref=dcelist.at(current->index.at(i))->edge->twin;	
			if(temphedge->origin->p.x>current->xcoord){
				pf=temphedge->origin->p;
				pi=dcelist.at(current->index.at(i))->p;
				a=(pf.y-pi.y)/(pf.x-pi.x);
				b=pf.y-a*pf.x;	
				yhalf=a*xhalf+b;
				slabbstyinsert(yhalf,&current->slaby,a,b,temphedge);		
				slabytemp.a=a; slabytemp.b=b; slabytemp.xcoord=pf.x; slabytemp.edge=temphedge;
				slabhistory.push_back(slabytemp);
			}
			while(temphedge->next->twin!=tempref){
				temphedge=temphedge->next->twin;
				if(temphedge->origin->p.x>current->xcoord){
					pf=temphedge->origin->p;
					pi=dcelist.at(current->index.at(i))->p;
					a=(pf.y-pi.y)/(pf.x-pi.x);
					b=pf.y-a*pf.x;	
					yhalf=a*xhalf+b;
					slabbstyinsert(yhalf,&current->slaby,a,b,temphedge);		
					slabytemp.a=a; slabytemp.b=b; slabytemp.xcoord=pf.x; slabytemp.edge=temphedge;
					slabhistory.push_back(slabytemp);
				}
			}
		}
		current=current->next;
	}
	if(!root1vs2){ dcelist.swap(veinvtcs); }
}


void slabbstyinsert(double ycoord,slabbsty** rootptr, double a, double b, hedge* temphedge)
{
	if(!*rootptr){
		*rootptr = new slabbsty(a,b,ycoord);
		(*rootptr)->fup=temphedge->twin->f;
		(*rootptr)->fdown=temphedge->f;
		return;
	}
	slabbsty *temproot;
	temproot = *rootptr;
	while (temproot){
		if(ycoord<temproot->ycoord){
			if(temproot->down){
				temproot=temproot->down;
			}else{
				temproot->down = new slabbsty(a, b, ycoord);
				temproot->down->fup=temphedge->twin->f;
				temproot->down->fdown=temphedge->f;
				return;
			}
		}
		if(ycoord>temproot->ycoord){
			if(temproot->up){
				temproot=temproot->up;
			}else{
				temproot->up = new slabbsty(a, b, ycoord);
				temproot->up->fup=temphedge->twin->f;
				temproot->up->fdown=temphedge->f;
				return;
			}
		}
		if(ycoord==temproot->ycoord){ 
			return;
		}
	}
}

void inOrder2(bstdcelnode *ptr){
	if(ptr->left) inOrder2(ptr->left);
	cout<<"bstdcel tree point: "<< ptr->ycoord<<endl;
	if(ptr->right) inOrder2(ptr->right);
}

void inOrder(slabbsty *ptr, double xi, double xf)
{
	if(ptr->down) inOrder(ptr->down, xi, xf);
	seg *tempseg;
	tempseg=new seg();
	tempseg->start.x=xi;
	tempseg->end.x=xf;
	tempseg->start.y=ptr->a*xi+ptr->b;
	tempseg->end.y=ptr->a*xf+ptr->b;
	testout.push_back(tempseg);
	if(ptr->up) inOrder(ptr->up,xi,xf);
}

face* plquery(point p, bool root1vs2)
{
	slabbstx *temproot;
	if(root1vs2){ temproot = slabbstxroot; }else{ temproot = slabbstxroot2; }
	slabbsty *tempry;
//	cout<<"checking... step 1"<<endl;
	while(temproot){
		if(p.x<temproot->xcoord) if(temproot->left) {temproot=temproot->left;}else{temproot=temproot->prev; break;}
		if(p.x>temproot->xcoord) if(temproot->right){ temproot=temproot->right;}else{break;}
		if(p.x==temproot->xcoord) break;
	}	
	tempry = temproot->slaby;
//	inOrder(tempry, temproot->xcoord, temproot->xcoordf);
//	cout<<"checking... step 2"<<endl;
	while(tempry){
		if(p.y>tempry->a*p.x+tempry->b){ 
			if(tempry->up) {
				tempry=tempry->up; 
			}else{
				return tempry->fup;
			}
		}
		if(p.y<tempry->a*p.x+tempry->b) {
			if(tempry->down) {
				tempry=tempry->down;
			}else{
				return tempry->fdown;
			}
		}
	}
}

void sitefacelinkage(bool veinorauxin)
{	
//	if(veinorauxin){
	for(int i=0;i<veinpoints.size();i++){
		queryface=plquery(veinpoints[i], veinorauxin);
		queryface->site=veinpoints[i];
	}
//	}else{
//	for(int i=0;i<auxinpoints.size();i++){
//		queryface=plquery(auxinpoints[i]);
//		queryface->site=auxinpoints[i];
//	}
//	}
}

bool voronoistructure(bool veinorauxin, bool sqhex){
	if(points.size()>0) while(points.size()!=0) points.pop();
	if(veinorauxin){
	for(int i=0;i<veinpoints.size();i++){
		points.push(veinpoints[i]);
	}
	}else{
	for(int i=0;i<auxinpoints.size();i++){
		points.push(auxinpoints[i]);
	}
	}
	voronoi(sqhex);
	if(createdcel()==2){
	for(int i=0;i<dcelist.size();i++){
		if(dcelist[i]) 	slabbstxinsert(dcelist.at(i)->p.x,i,dcelist.at(i)->p,1);
   	}
	slabbstxaddxf(1);/////////adds the x coordinate where each slab ends to each slab in the of the data structure
	createslabstruct(1);//////creates all of the binary search trees in the y direction and adds them to each slab of the data structure
	sitefacelinkage(veinorauxin);///////adds each voronoi site point coordinates to its corresponding voronoi face in the data structure
	return true; //everything went well
	}else{
	return false; //a face may not have been detected or maybe its an issue with the voronoi code itself <- issues still must be fixed
	}
}

veinnode* findvein(point p, veinnode *ptr){
	veinnode *best(ptr), *temp;
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			temp=findvein(p, ptr->next[i]);
			if(pow(best->p.x-p.x,2)+pow(best->p.y-p.y,2)>pow(temp->p.x-p.x,2)+pow(temp->p.y-p.y,2)){
				best=temp;	
			}
		}
	}
	return best;
}

bool veinintersection(point p0, point p1, veinnode *ptr){
	bool check(0), temp(0);
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			temp = veinintersection(p0, p1, ptr->next[i]);
			if(temp){ check = temp; break; }
			temp = linesegintersection(p0, p1, ptr->p, ptr->next[i]->p);
			if(temp){ check = temp; break; }
		}
	}
	return check;
}

void auxinadd(point p){
	auxinnode *temp;
	temp = new auxinnode(p);
	if(auxinroot){
		temp->next = auxinroot;
		auxinroot->prev=temp;
	}
	auxinroot = temp;
}

auxinnode* auxindel(auxinnode *temp){
	auxinnode *temp0;
	bool ok(1);
	if(temp->next){ if(temp->prev){ temp->prev->next=temp->next; temp->next->prev=temp->prev; temp0=temp->prev; }else{ auxinroot=temp->next; auxinroot->prev=0; temp0=auxinroot; } }else{ if(temp->prev){ temp->prev->next=0; temp0=temp->prev; }else{ ok=0; auxinroot=0;} }
	delete temp;
	temp=0;
	if(ok){
		return temp0;
	}else{
		return temp;
	}
}

/////////ok=0 every square grid cell that is too close to a vein node and sets quota to 0
void checkgrid(){
	point tmp;
	double dist, distf(0);
	for(int i = 0; i < grid.size(); i++){
		for(int k = 0; k < grid[i]->vetcs.size() ; k++){
			grid[i]->quota = 0;
	//		cout<<"x "<<grid[i]->vetcs[k].x<<" y "<<grid[i]->vetcs[k].y<<endl;
			tmp = grid[i]->vetcs[k];
			if(signof(tmp.x - squarevert[0].x)*(tmp.x - squarevert[0].x) < 0.0000001){ tmp.x+=0.000001; }
			if(signof(tmp.y - squarevert[0].y)*(tmp.y - squarevert[1].y) < 0.0000001){ tmp.y+=0.000001; }
			if(signof(tmp.x - squarevert[1].x)*(tmp.x - squarevert[1].x) < 0.0000001){ tmp.x-=0.000001; }
			if(signof(tmp.y - squarevert[2].y)*(tmp.y - squarevert[2].y) < 0.0000001){ tmp.y-=0.000001; }
	//		cout<<"x "<<tmp.x<<" y "<<tmp.y<<endl;
			tmp = plquery(tmp,1)->site;
			dist = sqrt(pow(grid[i]->vetcs[k].x-tmp.x,2)+pow(grid[i]->vetcs[k].y-tmp.y,2));
			if(distf == 0 || dist > distf){
				distf = dist;
			}
		}
		if( distf < birthdist ) { grid[i]->ok = 0; }else{ grid[i]->ok = 1;}
	}
}

int gridindget(point p, int numdiv){
	int xind,yind,gridind;
	xind = floor((p.x - squarevert[0].x)/(grid[0]->vetcs[1].x - grid[0]->vetcs[0].x));
	yind = floor((p.y - squarevert[0].y)/(grid[0]->vetcs[2].y - grid[0]->vetcs[0].y));
	gridind = xind + numdiv*yind;
	return gridind;
}

bool orderxpt(point i,point j) {return (i.x<j.x);}


//If a vein passes through a square just ok=0 it
void checkveingrid(veinnode *ptr, int numdiv){
	int gridind,gridind2,gridind3;
	gridind = gridindget(ptr->p,numdiv); grid[gridind]->ok=0;
	if(!ptr->next.empty()){
	point p1,p2,temppoint;
	vector< vector<point> > segpoints;
	vector<point> intersections;
	double diffx, diffy, divlengthx(grid[0]->vetcs[1].x - grid[0]->vetcs[0].x),divlengthy(grid[0]->vetcs[2].y - grid[0]->vetcs[0].y),maxy,miny,maxx,minx,gtx,gty;
	for(int i=0;i<ptr->next.size();i++){
		if(ptr->next[i]){
			if(ptr->width != minwidth && (ptr->next[i]->width != minwidth)){
			checkveingrid(ptr->next[i], numdiv);
			gridind2 = gridindget(ptr->next[i]->p,numdiv); grid[gridind2]->ok=0;
			if(gridind!=gridind2){
				/////////////////Chunk of code bellow stores grid line segments that intersect with p2-p1 segment//////////////
				p2=ptr->next[i]->p; p1=ptr->p; //p2 - grindind2 ... p1 - gridind
	//			cout<<"p2 "<<p2.x<<" "<<p2.y<<" p1 "<<p1.x<<" "<<p1.y<<endl;
				
				maxy=p2.y; 	if(maxy < p1.y){ maxy=p1.y; miny=p2.y; }else{ miny=p1.y; }
				maxx=p2.x; 	if(maxx < p1.x){ maxx=p1.x; minx=p2.x; }else{ minx=p1.x; }
				if(signof(grid[gridind2]->vetcs[0].x - grid[gridind]->vetcs[0].x)*(grid[gridind2]->vetcs[0].x - grid[gridind]->vetcs[0].x) > 0.0000001){
	//				cout<<"g2 "<<grid[gridind2]->vetcs[0].x<<" g1 "<<grid[gridind]->vetcs[1].x<<endl;
				if(p2.x > p1.x){
					diffx = grid[gridind2]->vetcs[0].x - grid[gridind]->vetcs[1].x;
					gtx=grid[gridind2]->vetcs[0].x;
				}else{
					diffx = grid[gridind]->vetcs[0].x - grid[gridind2]->vetcs[1].x;
					gtx=grid[gridind]->vetcs[0].x;
				}
				for(int j = 0; j < round((signof(diffx)*diffx)/divlengthx)+1; j++){
					vector<point> tempinter;
					temppoint.x = gtx - j*divlengthx; temppoint.y = maxy;
	//				cout<<"tmp "<<temppoint.x<<" "<<temppoint.y<<endl;
					tempinter.push_back(temppoint);
					temppoint.x = gtx - j*divlengthx; temppoint.y = miny;
	//				cout<<"tmp "<<temppoint.x<<" "<<temppoint.y<<endl;
					tempinter.push_back(temppoint);
					segpoints.push_back(tempinter);
	//				cout<<"segp size "<<segpoints.size()<<endl;
				}
				}
				if(signof(grid[gridind2]->vetcs[0].y - grid[gridind]->vetcs[0].y)*(grid[gridind2]->vetcs[0].y - grid[gridind]->vetcs[0].y) > 0.0000001){
	//				cout<<"g2 "<<grid[gridind2]->vetcs[0].y<<" g1 "<<grid[gridind]->vetcs[2].y<<endl;
				if(p2.y > p1.y){
					diffy = grid[gridind2]->vetcs[0].y - grid[gridind]->vetcs[2].y;
					gty=grid[gridind2]->vetcs[0].y;
				}else{
					diffy = grid[gridind2]->vetcs[0].y - grid[gridind]->vetcs[2].y;
					gty=grid[gridind]->vetcs[0].y;
				}
				for(int j = 0; j < round((signof(diffy)*diffy)/divlengthy)+1; j++){
					vector<point> tempinter;
					temppoint.x = maxx; temppoint.y = gty - j*divlengthy;
					tempinter.push_back(temppoint);
					temppoint.x = minx; temppoint.y = gty - j*divlengthy;
					tempinter.push_back(temppoint);
					segpoints.push_back(tempinter);
				}	
				}
				//////////////////////Now we store the intersections//////////////////////
	//				cout<<"segp size2 "<<segpoints.size()<<endl;
				for(int j = 0; j < segpoints.size() ; j++){
					if(linesegintersection(p1,p2,segpoints[j][0],segpoints[j][1])){
	//					cout<<" segpoints 1 "<< segpoints[j][0].x<<" "<<segpoints[j][0].y<<" segpoints2 "<<segpoints[j][1].x<<" "<<segpoints[j][1].y<<endl;
						intersections.push_back(bintersection);
	//					cout<<"intersection.x "<<intersections[j].x<<" intersections.y "<<intersections[j].y<<endl;
					}
				}
				/////////////finally by making small displacements in x and y we can find the indices of both boxes which form the line segment the intersection belongs to///////////////
				for(int j = 0; j < intersections.size(); j++){
					temppoint.x = intersections[j].x + 0.0000001; temppoint.y = intersections[j].y + 0.0000001;
					gridind3 = gridindget(temppoint,numdiv); grid[gridind3]->ok=0;
					temppoint.x = intersections[j].x - 0.0000001; temppoint.y = intersections[j].y - 0.0000001;
					gridind3 = gridindget(temppoint,numdiv); grid[gridind3]->ok=0;
				}
				//////////////////We are not taking into account the case where the intersection is one of the boxes vertices, since the chance that happens is too slim to even worry about it///////////
			}
			segpoints.clear(); intersections.clear();/////clears both vectors, since they may be used for in another iteration of the ptr->next[i] "for"
		}
	}
	}
	}
}

bool bboxpic(0),bboxflag(0);

int getoccupiedboxes(veinnode *ptr, int numdiv, bool ven1or2){
	int j,k,N(0);
	vector< vector<int> > bbox;
	//intiates bbox
	for(int i=0;i<numdiv;i++){
		vector<int>temp;
		for(int l=0;l<numdiv;l++){
				temp.push_back(0);
		}
		bbox.push_back(temp);
	}
	//fills bbox and increases N
	for(int i=0;i<grid.size();i++){
		k=(int) i/numdiv;
		j=i%numdiv;
		bbox[k][j] = 1*grid[i]->ok;
		if(!grid[i]->ok) N++;
	}	
	//writes on bbox file
	FILE *fl;
	if(ven1or2) fl=fopen("bbox.dat","w"); else fl=fopen("bbox2.dat","w");
	for(int i = numdiv - 1;i>-1;i--){
		for(int l = 0;l<numdiv;l++){
			fprintf(fl,"%d\t",bbox[i][l]);
		}
		fprintf(fl,"\n");
	}
	fclose(fl);
	if(bboxpic){
	FILE *fl;
	fl=fopen("bbox3.dat","w"); 
	for(int i = numdiv - 1;i>-1;i--){
		for(int l = 0;l<numdiv;l++){
			fprintf(fl,"%d\t",bbox[i][l]);
		}
		fprintf(fl,"\n");
	}
	fclose(fl);
	}
	return N;
}


int getoccupiedboxes2(veinnode *ptr, int numdiv){
	int j,k,N(0);
	vector< vector<int> > bbox;
	//intiates bbox
	for(int i=0;i<numdiv;i++){
		vector<int>temp;
		for(int l=0;l<numdiv;l++){
				temp.push_back(0);
		}
		bbox.push_back(temp);
	}
	//fills bbox and increases N
	for(int i=0;i<grid.size();i++){
		k=(int) i/numdiv;
		j=i%numdiv;
		bbox[k][j] = 1*grid[i]->ok;
		if(!grid[i]->ok) N++;
	}	
	//writes on bbox file
	FILE *fl;
	fl=fopen("bbox2.dat","w");
	for(int i = numdiv - 1;i>-1;i--){
		for(int l = 0;l<numdiv;l++){
			fprintf(fl,"%d\t",bbox[i][l]);
		}
		fprintf(fl,"\n");
	}
	fclose(fl);
	return N;
}

///////////adds auxins to the square grid cells in case they are not occupied by a vein edge or node
void genrddist(double denst, int numdiv, bool fullcheck){
	MTRand mt(seed); ////Starts the pseudo-random number generator
	double area((grid[0]->vetcs[1].x - grid[0]->vetcs[0].x)*(grid[0]->vetcs[2].y - grid[0]->vetcs[0].y));
	int fails(0), tries(0), xind, yind, gridint;
	point tmp, p;
	bool sourceok;
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
			xind = floor((i->p.x - squarevert[0].x)/(grid[0]->vetcs[1].x - grid[0]->vetcs[0].x));
			yind = floor((i->p.y - squarevert[0].y)/(grid[0]->vetcs[2].y - grid[0]->vetcs[0].y));
			gridint = xind + numdiv*yind;
			grid[gridint]->quota++;
	}
	for(int i = 0; i < grid.size() ; i++){
		if(grid[i]->ok){
			for(int j = grid[i]->quota; j < round(area*denst) + fails; j++){
				sourceok=1;
				p.x = mt()*(grid[0]->vetcs[1].x - grid[0]->vetcs[0].x) + grid[i]->vetcs[0].x;
				p.y = mt()*(grid[0]->vetcs[2].y - grid[0]->vetcs[0].y) + grid[i]->vetcs[0].y;
				tmp = plquery(p,1)->site;
				if(sqrt(pow(p.x-tmp.x,2)+pow(p.y-tmp.y,2))<birthdist) sourceok=0;
				if(fullcheck){
					 tmp = plquery(p,0)->site;
					if(sqrt(pow(p.x-tmp.x,2)+pow(p.y-tmp.y,2))<birthdist) sourceok=0;
				}
				if(sourceok){
					for( auxinnode *l = auxinroot ; l ; l = l->next ){
					if(sqrt(pow(p.x - l->p.x,2)+pow(p.y - l->p.y,2))<birthdist){ sourceok=0; break;}
					}
				}
				if(sourceok){ auxinadd(p); }else{ fails++; tries++;  }
				if(tries==maxtries){break;}
			}
			fails=0; tries=0;
		}
	}
}

////////////////////////checks all k-neighbors and see of they are populated or not... if one is populated return 0 else return 1//////////////////////////////////////////////////
bool checknneighbors(int xind, int yind, int numdiv2, int knn){
	bool check(1),check2(1);
	for(int i = 0; i < (2*knn + 1) ; i++){
		for(int l = 0; l < (2*knn + 1) ; l++){
				if((yind + l - knn) < 0 || (yind + l - knn) >= numdiv2 || (xind + i - knn) < 0 || (xind + i - knn) >= numdiv2 ){check2 = 0; break;}
		}
	}
	if(check2){
	for(int i = 0; i < (2*knn + 1) ; i++){
		for(int l = 0; l < (2*knn + 1) ; l++){
			if(i!=l){
				if(!grid[(xind + i - knn)+numdiv2*(yind + l - knn)]->ok){
					check = 0; break;
				}
			}
		}
	}
	}else{
		int sel;
		if(round(1.0*xind/numdiv2) == 0){ sel = 0; } if(round(1.0*xind/numdiv2) == 1){ sel = 1; } if(round(1.0*yind/numdiv2) == 0){ sel = 2; } if(round(1.0*yind/numdiv2) == 1){ sel = 3; } 
	//	sel = round(1.0*xind/numdiv2) + 2*round((1.0*yind/numdiv2));
		switch (sel){
			case 0:
				for(int i = 0; i < knn + 1 ; i++){
					for(int l = 0; l < (2*knn + 1) ; l++){
						if(i!=l){
						if((yind + l - knn) >= 0 && (yind + l - knn) < numdiv2){
						if(!grid[(xind + i)+numdiv2*(yind + l - knn)]->ok){
							check = 0; break;
						}
						}
						}
					}
				}
				break;
			case 1:
				for(int i = 0; i < knn + 1 ; i++){
					for(int l = 0; l < (2*knn + 1) ; l++){
						if(i!=l){
						if((yind + l - knn) >= 0 && (yind + l - knn) < numdiv2){
						if(!grid[(xind + i - knn)+numdiv2*(yind + l - knn)]->ok){
							check = 0; break;
						}
						}
						}
					}
				}
				break;
			case 2:
				for(int i = 0; i < (2*knn + 1) ; i++){
					for(int l = 0; l < (knn + 1) ; l++){
						if(i!=l){
						if((xind + i - knn) >= 0 && (xind + i - knn) < numdiv2){
						if(!grid[(xind + i - knn)+numdiv2*(yind + l)]->ok){
							check = 0; break;
						}
						}
						}
					}
				}
				break;
			case 3:
				for(int i = 0; i < (2*knn + 1) ; i++){
					for(int l = 0; l < (knn + 1) ; l++){
						if(i!=l){
						if((xind + i - knn) >= 0 && (xind + i - knn) < numdiv2){
						if(!grid[(xind + i - knn)+numdiv2*(yind + l - knn)]->ok){
							check = 0; break;
						}
						}
						}
					}
				}
				break;
		}
	}
	return check;
}


//////////////////////////////////Main Leaf Venation Simulation Function///////////////////////////////////////
void leafvenationsim(double denst, int numdiv, bool sqhex){
	
   //	MTRand mt(seed); ////Starts the pseudo-random number generator
   	
   	///////////Initiates first vein seed//////////////
	point p; veinnode *tvein,*tvein2,*tvein3;
	if(diagbool){ p.x=squarevert[0].x+(center.x-squarevert[0].x)*0.1; p.y=0.0001; }else{ p.x=center.x; p.y=0.0001;}
	//p.x=0.5*center.x+1; p.y=0.0001;
	points.push(p); veinpoints.push_back(p); veinroot=new veinnode(p);

		///////////Initiates more variables and arrays used throughout the simulation////////
	vector<point> newnodes; 
	int newauxins(0),growthstep(0), tries,xind,yind;
	bool sourceok,nodeok,pointok;
	double modulus;
	
	
	
	/////////////////Leaf venation simulation stars here////////////////////
	if(voronoistructure(1, squarehex)){//////////////generates first Voronoi structure//////////


	
	///////////Starts step counting////////////////
	while(growthstep<growthstepmax){
	//birthdist=birthdist*(1.003);/////////
	cout<<"step "<<growthstep<<endl;
	newauxins=0;
	tries=0;
	
	///////Checks grid divisions and eliminates partitions inside vein nodes range
	checkgrid();
	
	///////excludes all grid cells near the border with 'extreme' x coordinates
	for(int i = 0; i < grid.size(); i++){
		yind=(int) i/numdiv;
		xind=i%numdiv;
		if(xind==0 || xind==1 || xind==numdiv-2 || xind==numdiv-1 || yind == numdiv - 1 || yind == numdiv - 2 || yind==0 || yind ==1){
			grid[i]->ok=0;
		}
	}
	
	///////Adds new auxin sources randomly throughout the system and checks whether they are valid sources or not. In case they aren't they are not included in the auxin vector 
	genrddist(denst, numdiv, 0);
	
	///////Adds source sites to facelist data structure
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
		plquery(i->p,1)->auxin.push_back(i->p);
	}
	///////Clears new nodes vector before this step's new nodes are computed
	newnodes.clear();
	
	///////Computes new vein nodes and add them to the vein node data structure under veinroot pointer
	for(int i=0;i<dcefacelist.size();i++){
		p.x=0;p.y=0;
		if(!dcefacelist[i]->auxin.empty()){
		for(int j=0;j<dcefacelist[i]->auxin.size();j++){
			modulus=sqrt(pow(dcefacelist[i]->auxin[j].x-dcefacelist[i]->site.x,2)+pow(dcefacelist[i]->auxin[j].y-dcefacelist[i]->site.y,2));
			p.x=p.x+(dcefacelist[i]->auxin[j].x-dcefacelist[i]->site.x)/modulus;
			p.y=p.y+(dcefacelist[i]->auxin[j].y-dcefacelist[i]->site.y)/modulus;
		}
		modulus=sqrt(pow(p.x,2)+pow(p.y,2));
		///growthmod
		p.x=dcefacelist[i]->site.x+veindist*p.x/modulus;
		p.y=dcefacelist[i]->site.y+veindist*p.y/modulus;
		if(sqhex){pointok=checkpointhex(p);}else{pointok=checkpointsq(p);}
		if(pointok){
			nodeok=1;
			if(grid[gridindget(p,numdiv)]->ok){
			for(int j=0;j<newnodes.size();j++){
				if(sqrt(pow(newnodes[j].x-p.x,2)+pow(newnodes[j].y-p.y,2))<sqrt(pow(veindist,2))){ nodeok=0; break; }
			}
			if(nodeok){
				newnodes.push_back(p);
				tvein=findvein(dcefacelist[i]->site,veinroot);
				tvein2=new veinnode(p);
				tvein2->prev=tvein;
				tvein->next.push_back(tvein2);
			}
			}
		}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	veinpoints.insert(veinpoints.end(),newnodes.begin(),newnodes.end());
	cleanup();
	if(voronoistructure(1, squarehex)){
	nodeok=1;
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
		if(!nodeok){ i=i->prev; nodeok=1;}
		queryface=plquery( i->p, 1 );
		if(sqrt(pow(i->p.x - queryface->site.x,2)+pow(i->p.y - queryface->site.y,2))<killdist){i=auxindel(i); if(!i){ break;} if(!i->prev){nodeok=0;}}
	}
	growth(sqhex);
	growthstep++;
	}else{//////////////////Problem Occured - tries to circumvent Voronoi problem//////////////////////
		break;
		for(int i=0;i<newnodes.size();i++){
			tvein=findvein(newnodes[i],veinroot);
			for(int j=0;j<tvein->prev->next.size();j++){
				if(tvein==tvein->prev->next[j]){ tvein->prev->next[j]=0; tvein->prev->next.erase(tvein->prev->next.begin()+j); break;}
			}
			delete tvein;
			tvein=0;
			veinpoints.erase(veinpoints.begin()+(veinpoints.size()-1));
		}
		for( auxinnode *i = auxinroot ; i ; i = i->next ){
			if( i->prev ) { delete i->prev; i->prev = 0; }
			if( !i->next ) { auxinroot=i; }
		}
		if(auxinroot) delete auxinroot;
		auxinroot=0;
		growthstep++;
		cleanup();
		lstepgrowth=lstepgrowth*(-1);
		while(!voronoistructure(1, squarehex)){growth(sqhex); cleanup();birthdist=birthdist*(0.998);}
		lstepgrowth=lstepgrowth*(-1);
	}//////////////////////////////////////////////////////////////////////////////////////////////////
	}
	

	///////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////SIMULATION INTERMISSION///////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////Necessary in order clear dcelist/dcefacelist vectors for the new venation generation////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	//   for(int i=0;i<10;i++){growth();}
	dcefacelist.swap(veinfaces);
	dcelist.swap(veinvtcs);
        cleanup();
	slabbstxroot2 = 0;
	for(int i=0;i<veinvtcs.size();i++){
		if(veinvtcs[i]) slabbstxinsert(veinvtcs.at(i)->p.x,i,veinvtcs.at(i)->p,0);
   	}
   	slabbstxaddxf(0);/////////adds the x coordinate where each slab ends to each slab in the of the data structure
   	createslabstruct(0);//////creates all of the binary search trees in the y direction and adds them to each slab of the data structure
	sitefacelinkage(0);///////adds each voronoi site point coordinates to its corresponding voronoi face in the data structure
	cleanup();
	veinpoints3.swap(veinpoints);
 	veinpoints.clear();
	
	///////Generates new seed for the second part of the simulation//////////////
	if(diagbool){  p.x=squarevert[0].x+(center.x-squarevert[0].x)*1.9; p.y=2*center.y-0.0001; }else{ p.x=center.x; p.y=2*center.y-0.0001;}

	//p.x=squarevert[2].x - 0.0001; p.y=2*center.y-0.0001;
	//p.x=center.x; p.y=2*center.y-0.0001;
	//p.x=1.5*center.x; p.y=2*center.y-0.0001;
	points.push(p); veinpoints.push_back(p); veinroot2=new veinnode(p);
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
			if( i->prev ) { delete i->prev; i->prev = 0; }
			if( !i->next ) { auxinroot=i; }
	}
	if(auxinroot) delete auxinroot;
	auxinroot=0;
	
	veindist*=100;
	birthdist*= bdmultiplier;
	killdist*= kdmultiplier;
	denst*= denstmultiplier;
	for(int i=0;i<10;i++){growth(0);}

	//defines widths and vein order	


	////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////
	///////////////////////////SIMULATION PART 2///////////////////
	//////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////
	if(voronoistructure(1, squarehex)){//////////////generates first Voronoi structure of the second venation//////////
	growthstep=0;
	veindist=sqrt(pow(veinroot->next[0]->p.x-veinroot->p.x,2)+pow(veinroot->next[0]->p.y-veinroot->p.y,2))*veindistmultiplier;
	
	
	
	///////////Starts step counting////////////////
	while(growthstep < maxsteps2){
	//birthdist=birthdist*(1.003);
	cout<<"step "<<growthstep<<endl;
	newauxins=0;
	tries=0;
	
	///////Checks grid divisions and eliminates partitions inside vein nodes range
	checkgrid();
	
	///////Adds new auxin sources randomly throughout the system and checks whether they are valid sources or not. In case they aren't they are not included in the auxin vector 
	genrddist(denst, numdiv, 1);
	
	///////Adds source sites to facelist data structure
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
		queryface=plquery(i->p,1);
		//if(pow(i->p.x-queryface->site.x,2)+pow(i->p.y-queryface->site.y,2)<pow(5*veindist,2)){	queryface->auxin.push_back(i->p); }
		queryface->auxin.push_back(i->p);
	}
	
	///////Clears new nodes vector before this step's new nodes are computed
	newnodes.clear();
	
	///////Computes new vein nodes and add them to the vein node data structure under veinroot pointer
	for(int i=0;i<dcefacelist.size();i++){
		p.x=0;p.y=0;
		if(!dcefacelist[i]->auxin.empty()){
		for(int j=0;j<dcefacelist[i]->auxin.size();j++){
			modulus=sqrt(pow(dcefacelist[i]->auxin[j].x-dcefacelist[i]->site.x,2)+pow(dcefacelist[i]->auxin[j].y-dcefacelist[i]->site.y,2));
			p.x=p.x+(dcefacelist[i]->auxin[j].x-dcefacelist[i]->site.x)/modulus;
			p.y=p.y+(dcefacelist[i]->auxin[j].y-dcefacelist[i]->site.y)/modulus;
		}
		modulus=sqrt(pow(p.x,2)+pow(p.y,2));
		///growthmod
		p.x=dcefacelist[i]->site.x+veindist*p.x/modulus;
		p.y=dcefacelist[i]->site.y+veindist*p.y/modulus;
		if(sqhex){pointok=checkpointhex(p);}else{pointok=checkpointsq(p);}
		if(pointok){
			nodeok=1;
			if(!veinintersection(dcefacelist[i]->site, p, veinroot)){
			if(nodeok){
				newnodes.push_back(p);
				tvein=findvein(dcefacelist[i]->site,veinroot2);
				tvein2=new veinnode(p);
				tvein2->prev=tvein;
				tvein->next.push_back(tvein2);
			}
			}
		}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	veinpoints.insert(veinpoints.end(),newnodes.begin(),newnodes.end());
	cleanup();
	if(voronoistructure(1,squarehex)){
	nodeok=1;
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
		if(!nodeok){ i=i->prev; nodeok=1;}
		queryface=plquery( i->p, 1 );
		if(sqrt(pow(i->p.x - queryface->site.x,2)+pow(i->p.y - queryface->site.y,2))<killdist){i=auxindel(i); if(!i){ break;} if(!i->prev){nodeok=0;}}
	}
	growth(sqhex);
	growthstep++;


	}else{//////////////////Problem Occured - tries to circumvent Voronoi problem//////////////////////
		break;
		for(int i=0;i<newnodes.size();i++){
			tvein=findvein(newnodes[i],veinroot2);
			for(int j=0;j<tvein->prev->next.size();j++){
				if(tvein==tvein->prev->next[j]){ tvein->prev->next[j]=0; tvein->prev->next.erase(tvein->prev->next.begin()+j); break;}
			}
			delete tvein;
			tvein=0;
			veinpoints.erase(veinpoints.begin()+(veinpoints.size()-1));
		}
		for( auxinnode *i = auxinroot ; i ; i = i->next ){
			if( i->prev ) { delete i->prev; i->prev = 0; }
			if( !i->next ) { auxinroot=i; }
		}
		if(auxinroot) delete auxinroot;
		auxinroot=0;
		growthstep++;
		cleanup();
		lstepgrowth=lstepgrowth*(-1);
		while(!voronoistructure(1,squarehex)){growth(sqhex); cleanup();birthdist=birthdist*(0.998);}
		lstepgrowth=lstepgrowth*(-1);
	}//////////////////////////////////////////////////////////////////////////////////////////////////////
	}
	
	/////////////////////////////////////////////////////
	//////////////////////////////////////////////////////
	////////////////////////////////////////////////
	////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
			if( i->prev ) { delete i->prev; i->prev = 0; }
			if( !i->next ) { auxinroot=i; }
	}
	if(auxinroot) delete auxinroot;
	auxinroot=0;
	int numdiv2;
	vector<int> array;
	array.reserve(10000);
	growthstep=0;xind=0;
	bool threshold(0);
	
	////2nd main loop
	while(growthstep<ven2multiplier*growthstepmax){
	//birthdist=birthdist*(1.003);
	cout<<"step "<<growthstep<<endl;
	
	//did the number of auxins cross the 100 source threshold? If so, check if the number of auxins is still greater than 100. If not don't do anything.
	if(xind>auxincutoff||!threshold){
		
	newauxins=0;
	tries=0;
	
	///////Checks grid divisions and eliminates partitions inside vein nodes range
	//checkgrid();
	numdiv2 = 100;
	grid.clear();
	array.clear(); 	array.reserve(10000);
	initiatesqgrid(numdiv2);
	checkveingrid(veinroot,numdiv2);
//	denst*= denstmultiplier;
	
	//cout<<"ok"<<getoccupiedboxes(veinroot,numdiv2)<<endl;
	for(int i = 0; i < grid.size(); i++){
		for(int k = 0; k < grid[i]->vetcs.size() ; k++){
			grid[i]->quota = 0;
		}
//		cout<<"gok"<<endl;
		yind=(int) i/numdiv2;
		xind=i%numdiv2;
//		cout<<"gok2"<<endl;
		
		if(grid[i]->ok){
//		if(){
//				if(grid[xind+numdiv2*(yind+1)]->ok || grid[xind+numdiv2*(yind-1)]->ok || grid[(xind+1)+numdiv2*(yind+1)]->ok || grid[(xind+1)+numdiv2*(yind)]->ok || grid[(xind+1)+numdiv2*(yind-1)]->ok){
//					array.push_back(i);
//				}
//			}
//			cout<<"ok=0 it"<<xind<<" "<<yind<<endl;
			//checks all neighbors and see if they are not populated by and vein edge or node: in case they aren't stores index of the cell in the array to ok=0 (treating it as if it were populated
	
		if(checknneighbors(xind,yind,numdiv2,9)){
				array.push_back(i);
		}
		if(!checknneighbors(xind,yind,numdiv2,3)){
				array.push_back(i);
		}
		}
	}
//	cout<<"array size "<<array.size()<<endl;
	for(int i = 0; i < array.size(); i++){
		grid[array[i]]->ok=0;
	}
	for(int i = 0; i < grid.size(); i++){
		yind=(int) i/numdiv2;
		xind=i%numdiv2;
		if(xind==0 || xind ==1|| yind==0 || xind == numdiv2-2 || xind==numdiv2-1){
			grid[i]->ok=1;
		}
	}
	if(!bboxflag){bboxpic=1; getoccupiedboxes(veinroot,numdiv2,1); bboxflag=1; bboxpic=0;}
//	cout<<"ok"<<getoccupiedboxes2(veinroot,numdiv2)<<endl;
	
	///////Adds new auxin sources randomly throughout the system and checks whether they are valid sources or not. In case they aren't they are not included in the auxin vector 
	genrddist(denst, numdiv2, 1);
	xind=0;
	///////Adds source sites to facelist data structure
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
		queryface=plquery(i->p,1);
		if(sqrt(pow(i->p.x-queryface->site.x,2)+pow(i->p.y-queryface->site.y,2))<sqrt(pow(6*birthdist,2))){	queryface->auxin.push_back(i->p); }
	//	queryface->auxin.push_back(i->p);
		xind++;
	}
	cout<<"auxins added: "<<xind<<endl; if(xind>auxincutoff) threshold=1;
	///////Clears new nodes vector before this step's new nodes are computed
	newnodes.clear();
	
	///////Computes new vein nodes and add them to the vein node data structure under veinroot pointer
	for(int i=0;i<dcefacelist.size();i++){
		p.x=0;p.y=0;
		if(!dcefacelist[i]->auxin.empty()){
		for(int j=0;j<dcefacelist[i]->auxin.size();j++){
			modulus=sqrt(pow(dcefacelist[i]->auxin[j].x-dcefacelist[i]->site.x,2)+pow(dcefacelist[i]->auxin[j].y-dcefacelist[i]->site.y,2));
			p.x=p.x+(dcefacelist[i]->auxin[j].x-dcefacelist[i]->site.x)/modulus;
			p.y=p.y+(dcefacelist[i]->auxin[j].y-dcefacelist[i]->site.y)/modulus;
		}
		modulus=sqrt(pow(p.x,2)+pow(p.y,2));
		///growthmod
		p.x=dcefacelist[i]->site.x+veindist*p.x/modulus;
		p.y=dcefacelist[i]->site.y+veindist*p.y/modulus;
		if(sqhex){pointok=checkpointhex(p);}else{pointok=checkpointsq(p);}
		if(pointok){
			nodeok=1;
			if(grid[gridindget(p,numdiv2)]->ok){
			if(!veinintersection(dcefacelist[i]->site, p, veinroot)){
			if(nodeok){
				newnodes.push_back(p);
				tvein=findvein(dcefacelist[i]->site,veinroot2);
				tvein2=new veinnode(p);
				tvein2->prev=tvein;
				tvein->next.push_back(tvein2);
			}
			}
			}
		}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////
	
	veinpoints.insert(veinpoints.end(),newnodes.begin(),newnodes.end());
	cleanup();
	if(voronoistructure(1,squarehex)){
	nodeok=1;
	for( auxinnode *i = auxinroot ; i ; i = i->next ){
		if(!nodeok){ i=i->prev; nodeok=1;}
		queryface=plquery( i->p, 1 );
		if(sqrt(pow(i->p.x - queryface->site.x,2)+pow(i->p.y - queryface->site.y,2))<killdist){i=auxindel(i); if(!i){ break;} if(!i->prev){nodeok=0;}}
	}


	}else{//////////////////Problem Occured - tries to circumvent Voronoi problem//////////////////////
		break;
		for(int i=0;i<newnodes.size();i++){
			tvein=findvein(newnodes[i],veinroot2);
			for(int j=0;j<tvein->prev->next.size();j++){
				if(tvein==tvein->prev->next[j]){ tvein->prev->next[j]=0; tvein->prev->next.erase(tvein->prev->next.begin()+j); break;}
			}
			delete tvein;
			tvein=0;
			veinpoints.erase(veinpoints.begin()+(veinpoints.size()-1));
		}
		for( auxinnode *i = auxinroot ; i ; i = i->next ){
			if( i->prev ) { delete i->prev; i->prev = 0; }
			if( !i->next ) { auxinroot=i; }
		}
		if(auxinroot) delete auxinroot;
		auxinroot=0;
		growthstep++;
		cleanup();
		lstepgrowth=lstepgrowth*(-1);
		while(!voronoistructure(1,squarehex)){growth(sqhex); cleanup();birthdist=birthdist*(0.998);}
		lstepgrowth=lstepgrowth*(-1);
	}//////////////////////////////////////////////////////////////////////////////////////////////////////
	}
		growth(sqhex);
		growthstep++;
		if(squarevert[1].x-squarevert[0].x > finallength*0.9) break;
	}/////////////////////////////////////////////////////////////////////////////////////////////
	
	
	}else{cout<<"Problem occurred. Try again."<<endl;}

	}else{cout<<"Problem occurred. Try again."<<endl;}
	
	////////////End of the simulation/////////////////////
	cout<<veinpoints.size()<<" vein nodes to be computed"<<endl;
	cleanup();
	voronoistructure(1,squarehex);
	////////////Grows leaf a little further to produce a smoother result/////////////////
//	for(int i=0;i<10;i++){growth(sqhex);}
	////////////Grows leaf a little further to produce a smoother result/////////////////
	while(squarevert[1].x-squarevert[0].x < finallength*0.9){		
		growth(sqhex);
		growthstep++;
		cout<<"step "<<growthstep<<endl;
//		 break;
	}

	/////////////////////////Saves all vein nodes from primary and complementary pattern into the veinpoints and veinpoints2 vectors respectively/////////////////	
	cleanup();	
	veinpoints2.swap(veinpoints);
	cout << "SIZE 1: "<<veinpoints2.size()<<endl;
	veinpoints.clear();
	veinpoints.swap(veinpoints3);
	cout << "SIZE 2: "<<veinpoints.size()<<endl;
	voronoistructure(1,squarehex);


	//////////////Defines the width of each vein segment according to murray's law//////////////
	v2=0;
	defineveinwidth(veinroot);
	defineveinorder(veinroot,1);

	int ord(ordercutoffprim);
	cout << "DEBUG 1"<<endl;	
//	deltreebeyondcutoff(veinroot, ord);
	cout << "DEBUG 2"<<endl;	
	minwidth=0.00001;
	murrayfrombase(veinroot,initialdiam);
//	minwidth=0.25;// v2=1;
        murrayexp=murrayexp2;
	v2=1;
	defineveinwidth(veinroot2);
	defineveinorder(veinroot2,1);
	ord = ordercutoffcomp;
//	deltreebeyondcutoff(veinroot2, ord);
	v2=0;
//	minwidth=0.3; 
	murrayfrombase(veinroot2,initialdiam);
//	redefineveinwidth(veinroot2);
//	defineveinwidth(veinroot2);
	//////////////////////////////////////////////////////////////////////////////////////////////

	cout << "DEBUG 3"<<endl;	
	/////////////////////////////////Calculates Vein density of both venations (VLA) per order////////////////////////////////
	double templength(0), firstordervla(0), totalvla(0), secondordervla(0), thirdordervla(0);
	output1.open("Venationdensities");
	output1<<"Primary Venation VLAs"<<endl<<endl;
	for(int l=1;l<4;l++){
	veinlengtharea(veinroot,l);	
	output1<<"Length(mm) of "<<l<<" order veins: "<<venationlength<<" corresponding vein density VLA(mm-1): "<<venationlength/(finallength*finallength)<<endl;
	templength = templength + venationlength;
	if(l==1) firstordervla += venationlength;
	if(l==2) secondordervla += venationlength;
	if(l==3) thirdordervla += venationlength;
	venationlength=0;
	}
	veinlengtharea(veinroot,0);	
	totalvla = totalvla + venationlength;
	output1<<"Total venation Length(mm) (including minor order veins) "<<venationlength<<" corresponding vein density VLA(mm-1): "<<venationlength/(finallength*finallength)<<endl;
	venationlength = venationlength - templength;
	output1<<"Minor order venation Length(mm)"<<venationlength<<" corresponding vein density VLA(mm-1): "<<venationlength/(finallength*finallength)<<endl<<endl;
	output1<<"Complementary Venation VLAs"<<endl<<endl;
	venationlength=0; templength=0;
	for(int l=1;l<4;l++){
	veinlengtharea(veinroot2,l);	
	output1<<"Length(mm) of "<<l<<" order veins: "<<venationlength<<" corresponding vein density VLA(mm-1): "<<venationlength/(finallength*finallength)<<endl;
	templength = templength + venationlength;
	if(l==1) firstordervla += venationlength;
	if(l==2) secondordervla += venationlength;
	if(l==3) thirdordervla += venationlength;
	venationlength=0;
	}
	veinlengtharea(veinroot2,0);
	totalvla = totalvla + venationlength;	
	output1<<"Total venation Length(mm) (including minor order veins) "<<venationlength<<" corresponding vein density VLA(mm-1): "<<venationlength/(finallength*finallength)<<endl;
	venationlength = venationlength - templength;
	output1<<"Minor order venation Length(mm) "<<venationlength<<" corresponding vein density VLA(mm-1): "<<venationlength/(finallength*finallength)<<endl<<endl;
	output1<<"Total Venation VLAs"<<endl<<endl;
	for(int l=1;l<6;l++){
	if(l==1) output1<<"Length(mm) of "<<l<<" order veins: "<<firstordervla<<" corresponding vein density VLA(mm-1): "<<firstordervla/(finallength*finallength)<<endl;
	if(l==2) output1<<"Length(mm) of "<<l<<" order veins: "<<secondordervla<<" corresponding vein density VLA(mm-1): "<<secondordervla/(finallength*finallength)<<endl;
	if(l==3) output1<<"Length(mm) of "<<l<<" order veins: "<<thirdordervla<<" corresponding vein density VLA(mm-1): "<<thirdordervla/(finallength*finallength)<<endl;
    if(l==4) output1<<"Total venation Length(mm) (including minor order veins) "<<totalvla<<" corresponding vein density VLA(mm-1): "<<totalvla/(finallength*finallength)<<endl;
    if(l==5) output1<<"Minor order venation Length(mm) "<<totalvla-firstordervla-secondordervla-thirdordervla<<" corresponding vein density VLA(mm-1): "<<(totalvla-firstordervla-secondordervla-thirdordervla)/(finallength*finallength)<<endl;
	}
	bifurcationcount=0; bifurcationcountmethod(veinroot);
	output1<<"Number of birfurcation in the primary venation: "<<bifurcationcount<<endl;
	cout<<"Number of birfurcation in the primary venation: "<<bifurcationcount<<endl;
	bifurcationcount=0; bifurcationcountmethod(veinroot2);
	output1<<"Number of birfurcation in the complementary venation: "<<bifurcationcount<<endl;
	cout<<"Number of birfurcation in the complementary venation: "<<bifurcationcount<<endl;
	output1.close();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "DEBUG 4"<<endl;	
	/////////////////////////////////////fractal dimension of first vein pattern determination/////////////////////////////////////////
/*	FILE *fb;
	fb=fopen("fractaldim.dat","w");
	int ndiv[500],N[500],saturation(1);
	double meanN(0), meandelta(0),covar(0),var(0),fractaldim,zerocoeff;
	bool flag01(0),flag02(0),flag03(0);
	for(int i=1;i<30;i++){
		grid.clear();
		ndiv[i]=i+15;
		initiatesqgrid(ndiv[i]);
		checkveingrid(veinroot,ndiv[i]);
		N[i]=getoccupiedboxes(veinroot,ndiv[i],1);
//		cout<<"total boxes "<<ndiv[i]*ndiv[i]<<" Occupied (N) "<<N[i]<<" box length "<<(squarevert[1].x-squarevert[0].x)/ndiv[i]<<endl;
		fprintf(fb,"%f\t%f\n",log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])),log(N[i]));
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	meanN=log(N[1]); meandelta = log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]));
	for(int i=2;i<30;i++){
//		if(flag02&&signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag03=true;}else{flag02=false;}
//		if(flag01&&signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag02=true;}else{flag01=false;}
//		if(signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag01=true;}
//		if(flag03){flag01=false; flag02=false; flag03=false; break;}
		saturation++;
		meanN += log(N[i]);
		meandelta += log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i]));
	}
	meanN = meanN/(1.0*saturation);
	meandelta = meandelta/(1.0*saturation);
	cout<<"saturation "<<saturation<<endl;
	for(int i=30;i<saturation;i++){
		covar+=(log(N[i]) - meanN)*(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - meandelta);
		var+=pow((log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - meandelta),2);
	}
	fractaldim=covar/var;
	zerocoeff=meanN - fractaldim*meandelta;
	output1.open("fractaldimensions");
	cout<<"Fractal dimension "<<fractaldim<<endl<<"Linear coefficient "<<zerocoeff<<endl;
    output1<<"Box counting fractal dimension of primary venation: "<<fractaldim<<endl<<"Linear coefficient: "<<zerocoeff<<endl<<"Min and max x: "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]))<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[saturation]))<<endl;
    output2.open("fractaldimensions2");
    output2<<fractaldim<<" "<<zerocoeff<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]))<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[saturation]))<<" ";
	fclose(fb);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << "DEBUG 5"<<endl;	
	
	
	
	
		/////////////////////////////////////fractal dimension of second vein pattern determination/////////////////////////////////////////
	fb=fopen("fractaldimcomp.dat","w");
	saturation=1; meanN=0; meandelta=0; covar=0; var=0; fractaldim=0;
	for(int i=1;i<30;i++){
		grid.clear();
		ndiv[i]=i+15;
		initiatesqgrid(ndiv[i]);
		checkveingrid(veinroot2,ndiv[i]);
		N[i]=getoccupiedboxes(veinroot2,ndiv[i],0);
//		cout<<"total boxes "<<ndiv[i]*ndiv[i]<<" Occupied (N) "<<N[i]<<" box length "<<(squarevert[1].x-squarevert[0].x)/ndiv[i]<<endl;
		fprintf(fb,"%f\t%f\n",log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])),log(N[i]));
	}
	//////we do this to get a nice resolution of the vein pattern when transforming it in a matrix, which will be used in future monte carlo simulations////
//	grid.clear();
//	initiatesqgrid(numdiv);
//	checkveingrid(veinroot2,numdiv);
//	getoccupiedboxes(veinroot2,numdiv,0);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	meanN=log(N[1]); meandelta = log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]));
	for(int i=2;i<30;i++){
//		if(flag02&&signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag03=true;}else{flag02=false;}
//		if(flag01&&signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag02=true;}else{flag01=false;}
//		if(signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag01=true;}
//		if(flag03){flag01=false; flag02=false; flag03=false; break;}
		saturation++;
		meanN += log(N[i]);
		meandelta += log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i]));
		
	}
	meanN = meanN/(1.0*saturation);
	meandelta = meandelta/(1.0*saturation);
	cout<<"saturation "<<saturation<<endl;
	for(int i=30;i<saturation;i++){
		covar+=(log(N[i]) - meanN)*(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - meandelta);
		var+=pow((log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - meandelta),2);
	}
	fractaldim=covar/var;
	zerocoeff=meanN - fractaldim*meandelta;
	cout<<"Fractal dimension "<<fractaldim<<endl<<"Linear coefficient "<<zerocoeff<<endl;
	output1<<"Box counting fractal dimension of primary venation: "<<fractaldim<<endl<<"Linear coefficient: "<<zerocoeff<<endl<<"Min and max x: "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]))<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[saturation]))<<endl;
	output2<<fractaldim<<" "<<zerocoeff<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]))<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[saturation]))<<" ";
	fclose(fb);




			/////////////////////////////////////fractal dimension of complete pattern determination/////////////////////////////////////////
	fb=fopen("fractaldimtotal.dat","w");
	saturation=1; meanN=0; meandelta=0; covar=0; var=0; fractaldim=0;
	for(int i=1;i<30;i++){
		grid.clear();
		ndiv[i]=i+15;
		initiatesqgrid(ndiv[i]);
		checkveingrid(veinroot2,ndiv[i]);
		checkveingrid(veinroot,ndiv[i]);
		N[i]=getoccupiedboxes(veinroot2,ndiv[i],0);
//		cout<<"total boxes "<<ndiv[i]*ndiv[i]<<" Occupied (N) "<<N[i]<<" box length "<<(squarevert[1].x-squarevert[0].x)/ndiv[i]<<endl;
		fprintf(fb,"%f\t%f\n",log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])),log(N[i]));
		cout<<log(1/((squarevert[1].x-squarevert[0].x)/(ndiv[i]-15)))<<" "<<(squarevert[1].x-squarevert[0].x)/(ndiv[i]-15)<<endl;
	}
	//////we do this to get a nice resolution of the vein pattern when transforming it in a matrix, which will be used in future monte carlo simulations////
//	grid.clear();
//	initiatesqgrid(numdiv);
//	checkveingrid(veinroot2,numdiv);
//	getoccupiedboxes(veinroot2,numdiv,0);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	meanN=log(N[1]); meandelta = log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]));
	for(int i=16;i<30;i++){
//		if(flag02&&signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag03=true;}else{flag02=false;}
//		if(flag01&&signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag02=true;}else{flag01=false;}
//		if(signof(log(N[i])-log(N[i-1]))*(log(N[i])-log(N[i-1]))/(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i-1]))) < 1.1){flag01=true;}
//		if(flag03){flag01=false; flag02=false; flag03=false; break;}
		saturation++;
		meanN += log(N[i]);
		meandelta += log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i]));
		
	}
	meanN = meanN/(1.0*saturation);
	meandelta = meandelta/(1.0*saturation);
	cout<<"saturation "<<saturation<<endl;
	for(int i=30;i<saturation;i++){
		covar+=(log(N[i]) - meanN)*(log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - meandelta);
		var+=pow((log(1/((squarevert[1].x-squarevert[0].x)/ndiv[i])) - meandelta),2);
	}
	fractaldim=covar/var;
	zerocoeff=meanN - fractaldim*meandelta;
	cout<<"Fractal dimension "<<fractaldim<<endl<<"Linear coefficient "<<zerocoeff<<endl;
	output1<<"Box counting fractal dimension of primary venation: "<<fractaldim<<endl<<"Linear coefficient: "<<zerocoeff<<endl<<"Min and max x: "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]))<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[saturation]))<<endl;
	output2<<fractaldim<<" "<<zerocoeff<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[1]))<<" "<<log(1/((squarevert[1].x-squarevert[0].x)/ndiv[saturation]))<<endl;
	output1.close();
	output2.close();
	fclose(fb);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout << "DEBUG 6"<<endl;	
	system("gnuplot script.p");
	system("gnuplot script2.p");
    system("gnuplot script5.p");

*/	cout << "DEBUG 7"<<endl;	
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
int main(int argc, char* argv[])
{
   // Creates vertices of the hexagonal border storing them in this array
   double hexsidelength(5);
   double hwidth(hexsidelength*sqrt(3)/2); //half-width of the hexagon
   double fheight(hexsidelength/2);
   hexvert[0].x=hwidth;hexvert[0].y=0;//A point at the base of the hexagon ("leaf")
   hexvert[1].x=2*hwidth;hexvert[1].y=fheight;//B next point of the hexagon vertex set in the counterclockwise direction
   hexvert[2].x=2*hwidth;hexvert[2].y=3*fheight;//C next point
   hexvert[3].x=hwidth;hexvert[3].y=4*fheight;//D next point
   hexvert[4].x=0;hexvert[4].y=3*fheight;//E next point
   hexvert[5].x=0;hexvert[5].y=fheight;//F next point
   hexvert[6].x=hwidth;hexvert[6].y=0;//A again - completing the loop -------------- all edges of the border are defined with this vertex 
   hexvert[7].x=2*hwidth;hexvert[7].y=fheight;//B again <- this guy is used simplify if a point is inside an hexagon or not

   double squareside;//  square side length


   
   double denst;
   int numdiv;

   if(argc!=30){
   cout << "Using default parameters. Usage:" << endl;
   }else{
   seed=atof(argv[1]);lstepgrowth=atof(argv[2]);ven2multiplier=atoi(argv[3]);growthstepmax=atoi(argv[4]);
   veindist=atof(argv[5]);birthdist=atof(argv[6]);killdist=atof(argv[7]);maxtries=atoi(argv[8]);
   squarehex=atoi(argv[9]);denst=atof(argv[10]);numdiv=atoi(argv[11]);squareside=atof(argv[12]);
   finallength=atof(argv[13]);maxsteps2=atoi(argv[14]);bdmultiplier=atof(argv[15]);kdmultiplier=atof(argv[16]);
   denstmultiplier=atof(argv[17]);veindistmultiplier=atof(argv[18]);auxincutoff=atoi(argv[19]);murrayexp=atof(argv[20]);printorder=atoi(argv[21]);
   background=atoi(argv[22]);diagbool=atoi(argv[23]);initialdiam=atof(argv[24]);murrayexp2=atof(argv[25]);removeendchannels=atoi(argv[26]);autooutput=atoi(argv[27]);
   ordercutoffprim=atoi(argv[28]);ordercutoffcomp=atoi(argv[29]);
   }
   
	output1.open("parameters");
	output1<<seed<<" "<<lstepgrowth<<" "<<ven2multiplier<<" "<<growthstepmax<<" "<<veindist<<" "<<birthdist<<" "<< killdist<<" "<<maxtries<<" "<<squarehex<<" "<<denst<<" "<<numdiv<<" "<<squareside<<" ";
	output1<<finallength<<" "<<maxsteps2<<" "<<bdmultiplier<<" "<<kdmultiplier<<" "<<denstmultiplier<<" "<<veindistmultiplier<<" "<<auxincutoff<<" "<<murrayexp<<" "<<printorder<<" "<<background<<" ";
	output1<<diagbool<<" "<<initialdiam<<" "<<murrayexp2<<" "<<removeendchannels<<" "<<autooutput<<" "<<ordercutoffprim<<" "<<ordercutoffcomp;
	output1.close();
	
   double minx(1),miny(0),maxx(squareside+minx),maxy(squareside+miny);
   squarevert[0].x=minx; squarevert[0].y=miny;
   squarevert[1].x=maxx; squarevert[1].y=miny;
   squarevert[2].x=maxx; squarevert[2].y=maxy;
   squarevert[3].x=minx; squarevert[3].y=maxy;
   squarevert[4].x=minx; squarevert[4].y=miny;
   squarevert[5].x=maxx; squarevert[5].y=miny;
   
    
   
   lstepgrowth *= sqrt(2)/2;
   
   if(squarehex){
      	center.x=hwidth;center.y=2*fheight;
   }else{
	center.x=squareside/2+minx;center.y=squareside/2+miny;
   }

   initiatesqgrid(numdiv);

   leafvenationsim(denst, numdiv, squarehex);

   cout << "Average distance between sites: " << computeAverageDistTree() << endl;

   if(autooutput){ 
      openscadoutput();
   }else{
   // The Chunk of code below prints the output onto the screen (Leaf) - OpenGL/GLUT/Glew
   glutInit(&argc, argv);          // Initialize GLUT
   glutInitDisplayMode(GLUT_DOUBLE);  // Enable double buffered mode
   glutInitWindowSize(800, 600);   // Set the window's initial width & height
   winsy=600; winsx=800;
   glutInitWindowPosition(500, 50); // Position the window's initial top-left corner
   glutCreateWindow("Voronoi Diagram");  // Create window with the given title
   glutDisplayFunc(display);       // Register callback handler for window re-paint event
   glutReshapeFunc(reshape);       // Register callback handler for window re-size event
   glutSpecialFunc(specialKeys);   // Register callback handler for special-key event
   glutKeyboardFunc(keyboard);   // Register callback handler for special-key event
   glutMouseFunc(mouse);	// Register callback handler for mouse event
   glutIdleFunc(idle);             // Register callback handler if no other event   
   initGL();                       // Our own OpenGL initialization
   glutMainLoop();                 // Enter the event-processing loop
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   };
   return 0;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END OF MAIN FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//




//Is the point inside the square border or not?
bool checkpointsq(point p)
{
	double cp;
	int signcpp1,signcpreference;//sign of the crossproduct between edge vector and vector to the query point ... AND ... sign of the crossproduct between edge vector and another edge vector (REFERENCE)
	for(int j=0;j<4;j++){
		cp=signof(crossproduct2d(squarevert[j],p,squarevert[j+1]));
		if(cp!=signof(crossproduct2d(squarevert[j],squarevert[j+2],squarevert[j+1]))) return false;
		if(cp==0) return false;
	}
	return true;
}

//Is the point inside the hexagonal border or not?
bool checkpointhex(point p)
{
	double cp;
	int signcpp1,signcpreference;//sign of the crossproduct between edge vector and vector to the query point ... AND ... sign of the crossproduct between edge vector and another edge vector (REFERENCE)
	for(int j=0;j<6;j++){
		cp=signof(crossproduct2d(hexvert[j],p,hexvert[j+1]));
		if(cp!=signof(crossproduct2d(hexvert[j],hexvert[j+2],hexvert[j+1]))) return false;
		if(cp==0) return false;
	}
	return true;
}

//Displace vertex p1 that is outside the hexagonal or square border to a new point that intersects the border
point checkborder(point p0, point p1,bool sqhex)
{
//p1  --- > the point that will be changed (end)
//p0 the one that won't change (start)
	point tempint; //temporary intersection
	int countintersection(0);
	if(sqhex){
	for(int j=0;j<6;j++){
		if(linesegintersection(p0,p1,hexvert[j],hexvert[j+1])){ 
			if(countintersection==0) {
				tempint=bintersection;//bintersection - border intersection ia a point that is assigned a new value whenever we call the function linesegintersection
				countintersection++;//this is counting how many intersections the segment has with the border - value ranges from 0 to 2
			}else{ countintersection++; break;}
		 }
	}  
	}else{
	for(int j=0;j<4;j++){
		if(linesegintersection(p0,p1,squarevert[j],squarevert[j+1])){ 
			if(countintersection==0) {
				tempint=bintersection;//bintersection - border intersection ia a point that is assigned a new value whenever we call the function linesegintersection
				countintersection++;//this is counting how many intersections the segment has with the border - value ranges from 0 to 2
			}else{ countintersection++; break;}
		}
	}
	}
	if(countintersection==1){ 
		return tempint;//if there is just one intersection return the value stored previosly in the tempint variable
	}else{//in case there are 2 intersections calculate the one which will replace the p1 vertex (in other words calculate the distance from p1 to each intersection and move p1 to the closest one)
		double length1(pow(p1.x-tempint.x,2)+pow(p1.y-tempint.y,2)), length2(pow(p1.x-bintersection.x,2)+pow(p1.y-bintersection.y,2)); //check distance from p1 to each one of the edge intersections
		if(length1<length2){
			return tempint;
		}else
			return bintersection;
	}
}	

//Sign function
int signof(double a){ return (a==0) ? 0 : ((a<0) ? -1 : 1);}

//Cross product in "2D" --- z=0
double crossproduct2d(point p1, point p2, point p3){ return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);}//point p1 is the origin of both vectors

//Do the line segments p1-p2 and p3-p4 intersect?
bool linesegintersection(point p1, point p2, point p3, point p4)
{
	bintersection.x=0;
	bintersection.y=0;
	double a1((p1.y-p2.y)/(p1.x-p2.x)),b1(p1.y-a1*p1.x),xi;
	if(p3.x!=p4.x){
	double a2((p3.y-p4.y)/(p3.x-p4.x)),b2(p3.y-a2*p3.x);
	point temppoint;
	if(a1==a2&&b1!=b2) return false;//they are parallel segments ....... finish the situation where b1 and b2 are equal in the future!
	xi=(b2-b1)/(a1-a2); // value of x at the intersection of the two lines
	bintersection.x=xi;
	bintersection.y=a1*xi+b1;
	if(p1.x<p2.x){temppoint=p1; p1=p2; p2=temppoint;} //swap point so p1.x will always be greater
	if(p3.x<p4.x){temppoint=p3; p3=p4; p4=temppoint;} //swap point so p3.x will always be greater
	if(xi<p1.x && xi>p2.x && xi<p3.x && xi>p4.x){ return true;}else return false;
	}else{
	point temppoint;
	xi=p3.x; // value of x at the intersection of the two lines
	double yi(a1*xi+b1);
	bintersection.x=xi;
	bintersection.y=yi;
	if(p1.x<p2.x){temppoint=p1; p1=p2; p2=temppoint;} //swap point so p1.x will always be greater
	if(p3.y<p4.y){temppoint=p3; p3=p4; p4=temppoint;} //swap point so p3.y will always be greater
	if(xi<p1.x && xi>p2.x && yi<p3.y && yi>p4.y){ return true;}else return false;
	}
}

//Does point p1 intersect with segment p2-p3?
bool pointsegintersection(point p1, point p2, point p3)
{
	
	//check if they are colinear---- if they are colinear the crossproduct must be 0
	//rounding
	if(rounding(crossproduct2d(p2,p1,p3),5)!=0) return false;

	//check if point p1 is between points p2 and p3	
	point temppoint;
	if(p2.x!=p3.x){
	if(p2.x<p3.x){temppoint=p2; p2=p3; p3=temppoint;} //swap point so p2.x will always be greater
	if(p1.x<=p2.x && p1.x>=p3.x) return true;
	}else{
	if(p2.y<p3.y){temppoint=p2; p2=p3; p3=temppoint;} //swap point so p2.x will always be greater
	if(p1.y<=p2.y && p1.y>=p3.y) return true;
	}
	
	//point p1 is colinear but it's out of the segment frontier
	return false;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~Rounds number up to the nth decimal digit~~~~~~~~~~~~~~~~~~~~
double rounding(double x, int n){ return sqrt(pow(round(x*pow(10,n))/pow(10,n),2));}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//What kind of edge is the index i edge? //Edge types are described in the comments below
int outofbordercheck(int i,bool sqhex)
{
	point a1,a2;
	a1=output.at(i)->start;// the two vetices of the edge
	a2=output.at(i)->end;// the two vetices of the edge
	bool b1,b2;
	if(sqhex){
		b1=checkpointhex(a1);
        	b2=checkpointhex(a2);
	}else{
		b1=checkpointsq(a1);
        	b2=checkpointsq(a2);
	}
	if(b1&&b2) return 1;
	if(!b1&&b2) return 2;
	if(b1&&!b2) return 3;
	if(sqhex){
	for(int j = 0;j<6;j++){
		if(linesegintersection(a1,a2,hexvert[j],hexvert[j+1])){ return 4; }//does the edge intersect with any of the boundary segments????
	}
	}else{
	for(int j = 0;j<4;j++){
		if(linesegintersection(a1,a2,squarevert[j],squarevert[j+1])){ return 4; }//does the edge intersect with any of the boundary segments????
	}
	}
        return 0;
}//do the colinear case later
//If it returns 0, the edge is completely out of the boundary so it can be deleted.
//If it returns 1, the edge is completely inside the box
//If it returns 2, the start vertex of the edge is outside the box and the end vertex is inside
//If it returns 3, the end vertex of the edge is outside the box and the start vertex is inside
//If it returns 4, both vertices are outside but part of the edge is inside

////////////////////////////////////FORTUNE'S ALGORITHM FUNCTIONS///////////////////////////////////////////
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Voronoi functions chunk~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Behold! Voronoi main function approaching! Let there be Voronoi Tesselation!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
void voronoi(bool sqhex)
{
   // Process the queues; select the top element with smaller x coordinate.
   while (!points.empty())
      if (!events.empty() && events.top()->x <= points.top().x){
      	 //circle event processing
	 circleevent();
      }else{
	 //site event processing
         arcinsert(points.top());
	 points.pop();
      }

   // After all points are processed, do the remaining circle events.
   while (!events.empty())
      circleevent();

   // Creates border vertices and edges (the border is hexagonal but it can be easily modified) and deletes unnecessary vertices
   finish_edges(sqhex);
   deleteroot(root);
   root=0;
}

//////////////Deletes all remaining arc pointers and reclaims their memory back///////////////////////////
void deleteroot(arcbst *root){
	if(root->left) deleteroot(root->left);
	if(root->right) deleteroot(root->right);
	delete root;
}

void circleevent()
{   
   event *e(events.top());// Get the next event from the queue.
   events.pop();
   if (e->valid) {
      arcdelete(e->a);// Remove the associated arc from the arcbst
   }
   delete e;
}

void arcdelete(arcbst *i){
	bool left;//boolean that says whether j (the other node that share the same parent with i) is the right or the left child of k (i's parent node)
	bool pointatbase(0);
	double nextx(-1);
	arcbst *j(0), *k(i->parent), *a(0), *temp(0);//create null pointers that will be defined later ... k is defined as i's parent node
	point p1, p2;//declare 2 point variables that will be used later -> they will be the focuses that define temp's breakpoint later
      	seg *s = new seg(i->e->p);//here I create the new edge that will be leaving this vertex (the vertex was previously calculated in a circleevent call and stored in i->e->p)
	if(k->left == i){ j = k->right; a = getprev(i, 0); left = 0; }else{ j = k->left; a = getnext(i, 0); left = 1;}//
	p1 = getprev(i, 1)->p1;
	p2 = getnext(i, 1)->p1;
	if(!a->s->st) a->s->finishst(i->e->p); else{ if(!a->s->ed) a->s->finishend(i->e->p); }
	if(!k->s->st) k->s->finishst(i->e->p); else{ if(!k->s->ed) k->s->finishend(i->e->p); }		
	temp = new arcbst(p1,1,a->parent,p2,0);
	double dist1, dist2;
	dist1= pow(intersection(temp, i->e->x).x - i->e->p.x, 2) + pow(intersection(temp, i->e->x).y - i->e->p.y, 2);
	temp->gt = 1;
	dist2= pow(intersection(temp, i->e->x).x - i->e->p.x, 2) + pow(intersection(temp, i->e->x).y - i->e->p.y, 2);
	if(rounding(sqrt(dist1),7) != rounding(sqrt(dist2),7)) {	
		if(rounding(sqrt(dist1),6) < rounding(sqrt(dist2),6)) { temp->gt = 0; }else{ temp->gt = 1; }	
	}else{
//		cout<<"PPPPROBLEM"<<endl;
		pointatbase = 1;
		temp->gt = 0;
		if(!points.empty()){
      		if (!events.empty() && events.top()->x <= points.top().x)
			nextx = events.top()->x;
		else
			nextx = points.top().x;
		}else
			nextx = events.top()->x;
		if(nextx == -1) nextx = i->e->x + 0.5;
	}
	if(!left){
		if(k == a->right){
		temp->right = j; j->parent = temp;
		}else{
		temp->right = a->right; a->right->parent = temp;
		k->parent->left = j; j->parent = k->parent;		
		}
		if(a->parent){ if(a->parent->left == a) a->parent->left = temp; else a->parent->right = temp; }else{ root = temp; }
		temp->left = a->left; a->left->parent = temp;
		if(pointatbase){
		if(temp->left->tp){
			checkintersection(temp->left, temp, nextx);
			if(getprev(temp, 0)){ if(intersection(temp, nextx).y < intersection(getprev(temp, 0), nextx).y) temp->gt = 1; }
		}else{
			if(getprev(temp, 0)){ if(intersection(temp, nextx).y < intersection(getprev(temp, 0), nextx).y) temp->gt = 1; }
		}
		}
		while(j->tp){//while i is an internal node
			j = j->left;//just move down the tree going left as far as possible
		}
		check_circle_event(j, i->e->x); check_circle_event(getprev(j,1), i->e->x);
	}else{
		if(k == a->left){
		temp->left = j; j->parent = temp;
		}else{
		temp->left = a->left; a->left->parent = temp;
		k->parent->right = j; j->parent = k->parent;		
		}
		if(a->parent){ if(a->parent->left == a) a->parent->left = temp; else a->parent->right = temp; }else{ root = temp; }
		temp->right = a->right; a->right->parent = temp;
		if(pointatbase){
		if(temp->left->tp){
			checkintersection(temp->left, temp, nextx);
			if(getprev(temp, 0)){ if(intersection(temp, nextx).y < intersection(getprev(temp, 0), nextx).y) temp->gt = 1; }
		}else{
			if(getprev(temp, 0)){ if(intersection(temp, nextx).y < intersection(getprev(temp, 0), nextx).y) temp->gt = 1; }
		}
		}
		while(j->tp){//while i is an internal node
			j = j->right;//just move down the tree going left as far as possible
		}
		check_circle_event(j, i->e->x); check_circle_event(getnext(j,1), i->e->x);
	}		
	temp->s = s;
	delete a; delete k; delete i;
}

void checkintersection(arcbst *i, arcbst *temp, double x0){
	while(i->right->tp){
		i = i->right;
	}
	if(intersection(temp, x0).y < intersection(i, x0).y){ temp->gt = 1; return;}
}

//returns the leftmost node point data of an arcbst tree(or subtree)!
point getmin(arcbst* i){
	while(i->tp){//while i is an internal node
		i = i->left;//just move down the tree going left as far as possible
	}
	return i->p1;//returns this arcnode's point p1
}

arcbst *nextparent(0), *prevparent(0), *parent(0);//store next and prev neighbor internal nodes of an arc: important to deal with the edges that are being generated!

//This function inserts the new site encountered by the directrix into the arc bst
void arcinsert(point p){

	//in case root is null
	if(!root){
		root = new arcbst(p, 0, 0, pproblem);//if the root is null all we need to do is add the new site to the tree... our first arc!
		return;
	}

	//gets the correct leaf node in which the new site is going to be inserted;
	arcbst *i(findarc(p)), *j(0);
	
	//in case the arc i is not an internal node of the tree do the following (YES... although unlikely the new site CAN have the same y coordinate as one of the breakpoints!)
	if(!i->tp){
		// the subtree that will replace the arc node depends on whether the new node focus y coordinate is greater than the old node focus y coordinate
   			// Invalidate any old event which was stored previously in i: remember that if there
			if (i->e){//<<== check if this is really true!
					i->e->valid = false;
			}
			i->e = 0;// will be deleted later when this event is queried in the priority queue
		if( p.y > i->p1.y ){
			i->tp = 1; i->p2 = p;					   // before	p1	    old	bst	p1			after	p1	    new bst     b1
			i->gt = 0;						   //	______________									     '/'  '\'	
			i->left = new arcbst(i->p1, 0, i, pproblem);		   //		    p2								p2           p1    b2
			i->right = new arcbst(i->p1, 1, i, p, 1);		   //  label:							__________________		 '/''\'
			i->right->left = new arcbst(p, 0, i->right, pproblem);	   //  b1 - "smaller" breakpoint (y coodinate is smaller)					 p2  p1
			i->right->right = new arcbst(i->p1, 0, i->right, pproblem);//  b2 - "greater" breakpoint ... p1 - old arc focus ... p2 - new arc focus	   
			i->s = new seg(); i->right->s = i->s; //creates a segment between both breakpoints (b1 and b2)
			check_circle_event(i->left, p.x); check_circle_event(i->right->right, p.x);//tests both extreme leaf nodes of that subtree that has been created for new circle events;
		}else{ // in truth the code above would work for the case bellow too without any modifications --- it would work even if p1 and p2 had the same y coordinate, but I'll keep the code bellow for now
			i->tp = 1;  i->p2 = p;					// before	p1	    old	bst	p1			after	p1	    new bst     b2
			i->gt = 1;						//	______________									     '/'  '\'	
			i->right= new arcbst(i->p1, 0, i, pproblem);		//	   p2    						  p2     		     b1    p1
			i->left = new arcbst(i->p1, 1, i, p, 0);		//  label:							__________________	   '/''\'
			i->left->right = new arcbst(p, 0, i->left, pproblem);	//  b1 - "smaller" breakpoint (y coodinate is smaller)					   p1  p2
			i->left->left = new arcbst(i->p1, 0, i->left, pproblem);//  b2 - "greater" breakpoint ... p1 - old arc focus ... p2 - new arc focus	   
			i->s = new seg(); i->left->s = i->s; //creates a segment between both breakpoints (b1 and b2)
			check_circle_event(i->right, p.x); check_circle_event(i->left->left, p.x);//tests both extreme leaf nodes of that subtree that has been created for new circle events;
		}
	}else{
		//THIS PART OF THE ARCINSERT FUNCTION IS ADMITTEDLY NOT WORKING!!! HOWEVER THE CHANCES THIS PART IS EVER USED ARE SO SLIM THAT IT'S JUST NOT WORTH THE EFFORT SOLVING THIS
		//I TESTED THE VORONOI FUNCTION WITH 1,000,000 SITES AND THIS PART HAS NOT BEEN REQUESTED EVEN ONCE! AND I AM NOT GOING ANYWHERE NEAR THAT RANGE OF SITES IN MY APPLICATIONS
		//IF YOU, WHO ARE READING THIS, WANT TO FIX THIS GLITCH, BY ALL MEANS, SUIT YOURSELF!
		cout<<"ARC INSERT FUNCTION PROBLEM"<<endl;
		seg* newseg(0);
		//what the fuck do I do if the new point has the same y coordinate as a breakpoint? --- although very unlikely it is possible!	
		//first I need to know which point (p1 or p2) is related to which side (right or left)
		if(i->p1 == getmin(i->right)){//the function getmin gets the leftmost node of a tree(or subtree) so if the point p1 is equal to the leftmost node of the right subtree p1 corresponds to 'right'
			j = new arcbst(p, 1, i, i->p2, 0);//p2 here corresponds to the left of the i breakpoint so I include it in the new breakpoint being created as well as p;
			i->p2 = p;//then I change the point p2 in the old breakpoint to p;
		}else{//in case p1 does not correspond to 'right', p2 does.
			j = new arcbst(p, 1, i, i->p1, 0);//p1 here corresponds to the left of the i breakpoint so I include it in the new breakpoint being created as well as p;
			i->p1 = p;//then I change the point p1 in the old breakpoint to p;
			
		}
		//i the 'old' breakpoint corresponds to the 'greater' parabola intersection while the j we've just created corresponds to the 'lesser' one	
		i->gt = 1;
		j->right = new arcbst(p, 0, j, pproblem); //here I add the right child of the 'lesser' breakpoint j and assign the focus p that's being inserted into the tree
		j->left = i->left; //finally I include both nodes we've created into the 'old' tree ... j will now be where i->left was previously 
		i->left = j;
//		i->s = new seg(); i->left->s = i->s; //creates a segment between both breakpoints (b1 and b2)
	}
}

//fuction that prints the arc bst onto the screen (IN ASCENDING ORDER): mainly used for debugging
void printtree(arcbst *ptr, double l){
	if(ptr->left) printtree(ptr->left, l);
	if(ptr){ if(ptr->tp) cout<<"breakpoint y: "<<intersection(ptr, l).y<<endl; else cout<<"leaf y:"<<ptr->p1.y<<endl;}
	if(ptr->right) printtree(ptr->right, l);
}

//finds the arc in which a leaf node is in the arc bst
arcbst* findarc(point p){
	arcbst *i(root); 
	while(i->tp){////while the node is an internal node of the bst do the following
		if(p.y == intersection(i, p.x).y){//in case the y coordinate of the new site p is exactly the same as a break point, get out of the loop and return current breakpoint node
			break;
		}
		if(p.y > intersection(i, p.x).y){///I can handle both points of the arc node (p1 and p2 have values on them) since I know i is not a leaf node
			///the above line stores the nearest previous internal neighbor node of the point p
			i = i->right;	
			
		}
		if(!i->tp) break;
		if(p.y < intersection(i, p.x).y){///I can handle both points of the arc node (p1 and p2 have values on them) since I know i is not a leaf node
			i = i->left;	
		}
	}
	return i;//return the arc(internal or leaf) as well as previous and next internal nodes through the global variables prevparent and nextparent
}
///////////////////////////////////////////////////////////

//Given an arc in the bst what's the next leafnode (full=1) or what's the next breakpoint (full=0)
arcbst* getnext(arcbst *i, bool full){
	if(i->parent){
	while(i == i->parent->right){
		if(i->parent->parent){ i = i->parent; }else{ return 0; }
	}
	}
	if(!full) return i->parent;
	arcbst *j(i->parent->right);
	while(j->tp){//while i is an internal node
		j = j->left;//just move down the tree going left as far as possible
	}
	return j;//returns this arcnode
}

//Given an arc in the bst what's the previous leafnode (full=1) or what's the previous breakpoint (full=0)
arcbst* getprev(arcbst *i, bool full){
	if(i->parent){
	while(i == i->parent->left){
		if(i->parent->parent){ i = i->parent; }else{ return 0; }
	}
	}
	if(!full) return i->parent;
	arcbst *j(i->parent->left);
	while(j->tp){//while i is an internal node
		j = j->right;//just move down the tree going left as far as possible
	}
	return j;//returns this arcnode
}

// Look for a new circle event for arc i.
void check_circle_event(arcbst *i, double x0)
{
   	arcbst *nx(getnext(i, 1)), *pr(getprev(i, 1));//done so prev and next are only computed once
   	if (!nx || !pr)//if, by any reason, the next arc or prev arc are null, stop the function: there are no circle events! 
      		return;//stoping function-> no circle events to be computed
   
   	// Invalidate any old event which was stored previously in i: only true when the function was called processing a circle event
   	if (i->e){
        	i->e->valid = false;
  	}
	i->e = 0;// will be deleted later when this event is queried in the priority queue

   	double xx;
   	point o;

   	if (circle(pr->p1, i->p1, nx->p1, &xx,&o)){
		if( xx > x0) {
	    	// Create new event.
     	    	i->e = new event(xx, o, i);
      	    	events.push(i->e);
		}
   	}
}

// Find the rightmost point on the circle through a,b,c. This function was implemented by .... link:
bool circle(point a, point b, point c, double *xx, point *o)
{
   // Check that bc is a "right turn" from ab.
   if ((b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y) > 0)
      return false;

   // Algorithm from O'Rourke 2ed p. 189.
   double A = b.x - a.x,  B = b.y - a.y,
          C = c.x - a.x,  D = c.y - a.y,
          E = A*(a.x+b.x) + B*(a.y+b.y),
          F = C*(a.x+c.x) + D*(a.y+c.y),
          G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

   if (rounding(G,14) == 0) return false;  // Points are co-linear.

   // Point o is the center of the circle.
   o->x = (D*E-B*F)/G;
   o->y = (A*F-C*E)/G;

   // o.x plus radius equals max x coordinate.
   *xx = o->x + sqrt( pow(a.x - o->x, 2) + pow(a.y - o->y, 2) );
   return true;
}

// Where do two parabolas or lines intersect given the location of the directrix?
point intersection(arcbst *i, double l)
{
   point p0(i->p1), p1(i->p2);
   bool gt(i->gt); //defines which intersection we are looking for, the one with greater y coordinate or the one with the lesser y value
   point res, p = p0;
   double y1, y2;
   double z0(2*(p0.x - l));
   double z1(2*(p1.x - l));

///////////Finds the intersection y coordinate first
   if (pow(p0.x - p1.x,2) < pow(10, -20)){
      res.y = (p0.y + p1.y) / 2;
   }else if (pow(p1.x - l,2) < pow(10,-20)){
      res.y = p1.y;
   }else if (pow(p0.x - l,2) < pow(10,-20)) {
      res.y = p0.y;
      p = p1;
   } else {
      // Use the quadratic formula.
      double a(1/z0 - 1/z1);
      double b(-2*(p0.y/z0 - p1.y/z1));
      double c((p0.y*p0.y + p0.x*p0.x - l*l)/z0 - (p1.y*p1.y + p1.x*p1.x - l*l)/z1);

      //calculate the y coordinate of both intersections (in case they exist)
      y1 = ( -b - sqrt(b*b - 4*a*c) ) / (2*a);
      y2 = ( -b + sqrt(b*b - 4*a*c) ) / (2*a);

      if(gt){
        if(y1 < y2) y1 = y2;
        res.y = y1;
      }else{
        if(y1 > y2) y1 = y2;
        res.y = y1;
      }
	
   }

   // Plug back into one of the parabola equations and gets the intersections x coordinate
   res.x = (p.x*p.x + (p.y-res.y)*(p.y-res.y) - l*l)/(2*p.x-2*l);
   return res;
}

/////end of voronoi functions///////////////////

//~~~~~~~~~~~~~~~~~~~~~~SORTING FUNCTIONS (used to sort output vector instances)~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool orderx(seg* i,seg* j) {return (i->start.x<j->start.x);}
bool ordery(seg* i,seg* j) {return (i->start.y<j->start.y);}
bool orderxinv(seg* i,seg* j) {return (i->start.x>j->start.x);}
bool orderyinv(seg* i,seg* j) {return (i->start.y>j->start.y);}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void finishtraversal(arcbst *ptr, double l){
	if(ptr){
	if(ptr->left) finishtraversal(ptr->left, l);
	if(ptr->right) finishtraversal(ptr->right, l);
	if(ptr->tp){ if(ptr->s){ if(!ptr->s->st){ ptr->s->finishst(intersection(ptr, 2*l)); }else{ if(!ptr->s->ed) ptr->s->finishend(intersection(ptr, 2*l)); }} }
	}
}

//This is where border treatment is done
void finish_edges(bool sqhex)
{
   point temporary;
   seg *tempseg;
   int size=output.size();

   // Advance the sweep line so no parabolas can cross the bounding box.
   double l = 6*center.y;

     // printtree(root, 2*l);
   // Extend each remaining segment to the new parabola intersections.
   finishtraversal(root, l);

//Erases all edges and part of edges that are outside the border
   int erasedquant(0),minusmax(0);//How many instances were erased from output vector? -> erasedquant --- and how many were included? -> minusmax
   for(int i=0;i<output.size()-minusmax;i++){
	//each of the possible cases is explained thoroughly on the comments after the outofbordercheck function in this file
	switch(outofbordercheck(i,sqhex)){
	     case 0 :
		delete output.at(i);
		output.erase(output.begin()+i);//<<erase function sux gotta change this
		i--;
		erasedquant++;
		break;
	     case 2 :
		output.at(i)->start=checkborder(output.at(i)->end,output.at(i)->start,sqhex);
		tempseg=new seg(output.at(i)->start);
		minusmax++;
		break;
             case 3 :
		output.at(i)->end=checkborder(output.at(i)->start,output.at(i)->end,sqhex);
		tempseg=new seg(output.at(i)->end);
		minusmax++;
		break;
	     case 4:
		output.at(i)->start=checkborder(output.at(i)->end,output.at(i)->start,sqhex);
		tempseg=new seg(output.at(i)->start);
		output.at(i)->end=checkborder(output.at(i)->start,output.at(i)->end,sqhex);
		tempseg=new seg(output.at(i)->end);
		minusmax+=2;
		break;
	     default :
		break;
	}//nothing is done if it's equal to 1
   }

//Includes border hexagon or square vertices in the output vector
if(sqhex){
   for(int i=0;i<6;i++){ tempseg=new seg(hexvert[i]); }
}else{
   for(int i=0;i<4;i++){ tempseg=new seg(squarevert[i]); }
}

// Sorting border points in counter-clockwise order starting from A
vector<seg*> tempborder;
int newsize(0);
if(sqhex){
for(int i=0;i<6;i++){
   for (int j=size-erasedquant;j<output.size();j++){
	if(pointsegintersection(output.at(j)->start,hexvert[i],hexvert[i+1])){//does the start point intersect with edge hexvert(i->i+1)???
		tempborder.push_back(output.at(j));//in case it does, store it in this temp vector for future sorting
		output.erase(output.begin()+j);//erase point from output <=====bad funcion erase is too slow gotta change this
		j--;//with point erased, the output vector size will change to size-1, so it's necessary to subtract 1 from the counter so no instance will be skipped
	}
   }
   if(i<3){
	   sort(tempborder.begin()+newsize,tempborder.end(),ordery);//y coordinate of points must grow if i is under 3 for an hexagon (remember we are tracking the border in the counterclockwise direction starting from point A)
   }else
	   sort(tempborder.begin()+newsize,tempborder.end(),orderyinv);//y coordinate of points must now decrease if i is over 3 for the same hexagon
   newsize=tempborder.size();
}
}else{
for(int i=0;i<4;i++){
   for (int j=size-erasedquant;j<output.size();j++){
	if(pointsegintersection(output.at(j)->start,squarevert[i],squarevert[i+1])){//does the start point intersect with edge hexvert(i->i+1)???
		tempborder.push_back(output.at(j));//in case it does, store it in this temp vector for future sorting
		output.erase(output.begin()+j);//erase point from output <=====bad funcion erase is too slow gotta change this
		j--;//with point erased, the output vector size will change to size-1, so it's necessary to subtract 1 from the counter so no instance will be skipped
	}
   }
   if(i==0)   sort(tempborder.begin()+newsize,tempborder.end(),orderx);//x coordinate of points must grow for first border edge located at the bottom of the square
   if(i==1)   sort(tempborder.begin()+newsize,tempborder.end(),ordery);//y coordinate of points must grow for the second border edge (the next one in counterclockwise order)
   if(i==2)   sort(tempborder.begin()+newsize,tempborder.end(),orderxinv);//x coordinate of points must now decrease for third border edge
   if(i==3)   sort(tempborder.begin()+newsize,tempborder.end(),orderyinv);//y coordinate of points must now decrease for last border edge
   newsize=tempborder.size();
}
}

//Connecting the dots! (border)
   for (int j=0;j<tempborder.size();j++){
	if(j!=tempborder.size()-1){
	tempborder.at(j)->end=tempborder.at(j+1)->start;
	}else{
	tempborder.at(j)->end=tempborder.at(0)->start;
	}
   } 


//Inserting temporary vector used to sort points back into the output vector   
   output.insert(output.end(),tempborder.begin(),tempborder.end());

//Delete tempborder vector;
   tempborder.clear();

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BORDER FUNCTIONS END HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void deletebst(bstnode *ptr){
	if(ptr->left) {deletebst(ptr->left);}
	if(ptr->right) {deletebst(ptr->right);}
	ptr->strorend.clear();
	ptr->index.clear();
	vector<bool> ().swap(ptr->strorend);
	vector<int> ().swap(ptr->index);
	if(ptr)	delete ptr;
	ptr=0;
}

void deletebstdcel(bstdcelnode *ptr){
	if(ptr->left) {deletebstdcel(ptr->left);}
	if(ptr->right) {deletebstdcel(ptr->right);}
	ptr->index.clear();
	vector<int> ().swap(ptr->index);
	if(ptr)	delete ptr;
	ptr=0;
}

void deleteslabbstx(slabbstx *ptr){
	if(ptr->left) {deleteslabbstx(ptr->left);}
	if(ptr->right) {deleteslabbstx(ptr->right);}
	if(ptr->slaby) {deleteslabbsty(ptr->slaby);}
	ptr->index.clear();
	vector<int> ().swap(ptr->index);
	if(ptr)	delete ptr;
	ptr=0;
}

void deleteslabbsty(slabbsty *ptr){
	if(ptr->down) {deleteslabbsty(ptr->down);}
	if(ptr->up) { deleteslabbsty(ptr->up); }
//	delete (*ptr)->fup;
//	(*ptr)->fup=0;
//	delete (*ptr)->down;
//	(*ptr)->down=0;
	if(ptr)	delete ptr;
	ptr=0;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~BINARY SEARCH TREES FUNTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void bstinsert(double ycoord,int ind,bool stoe)
{
	if(!bstroot){
		bstroot = new bstnode(ycoord,ind,stoe);
		return;
	}
	bstnode *temproot;
	temproot = bstroot;
	while (temproot){
		if(ycoord<temproot->ycoord){
			if(temproot->left){
				temproot=temproot->left;
			}else{
				temproot->left = new bstnode(ycoord, ind, stoe);
				return;
			}
		}
		if(ycoord>temproot->ycoord){
			if(temproot->right){
				temproot=temproot->right;
			}else{
				temproot->right = new bstnode(ycoord, ind, stoe);
				return;
			}
		}
		if(ycoord==temproot->ycoord){ 
			bool indcheck(true);
			for(int i=0;i<temproot->index.size();i++){
				if(temproot->index.at(i)==ind) indcheck=false;
			}
			if(indcheck){
				temproot->index.push_back(ind);
				temproot->strorend.push_back(stoe);
			}
			return;
		}	
	}
}

//////////the overwrite boolean determines whether the function will overwrite a node in case someone tries to input coordinates that match the coordinates of an existing node
void bstdcelinsert(double ycoord,int ind, point p0){
	if(!bstdcelroot){
		bstdcelroot = new bstdcelnode(ycoord,ind);
		return;
	}
	bstdcelnode *temproot;
	temproot = bstdcelroot;
	while (temproot){
		if(ycoord<temproot->ycoord){
			if(temproot->left){
				temproot=temproot->left;
			}else{
				temproot->left = new bstdcelnode(ycoord, ind);
				return;
			}
		}
		if(ycoord>temproot->ycoord){
			if(temproot->right){
				temproot=temproot->right;
			}else{
				temproot->right = new bstdcelnode(ycoord, ind);
				return;
			}
		}
		////there might be a case in which one tries to insert the exact same coordinates into the tree/////
		////the code bellow will handle that situation: the old pointer in the dcelist with the given coordinates will be set as NULL before inserting a new node with the same coordinates into the binary search tree//////////////////////////////////////////////////////
		if(ycoord==temproot->ycoord){ 
	//	cout<<"bla"<<endl;
//			bool indcheck(true);
			for(int i=0;i<temproot->index.size();i++){
				if(dcelist[temproot->index[i]]){ if(dcelist[temproot->index[i]]->p==p0){ delete dcelist[temproot->index[i]]; dcelist[temproot->index[i]]=0; break;} }
			}
//			if(indcheck){
			temproot->index.push_back(ind);
//			}
			return;
		}
	}
}


void slabbstxinsert(double xcoord,int ind, point p0,bool root1vs2){
	slabbstx *temproot;
	if(root1vs2){ temproot = slabbstxroot; }else{ temproot = slabbstxroot2; }//There are 2 roots because there are 2 voronoi diagrams: one linked to the source channel and the other linked to the collecting channel
	if(!temproot){
		if(root1vs2){
			slabbstxroot = new slabbstx(xcoord,ind);
			return;
		}else{ 
			slabbstxroot2 = new slabbstx(xcoord,ind);
			return;
		}
	}
	while (temproot){
		if(xcoord<temproot->xcoord){
			if(temproot->left){
				temproot=temproot->left;
			}else{
				temproot->left = new slabbstx(xcoord, ind);
				return;
			}
		}
		if(xcoord>temproot->xcoord){
			if(temproot->right){
				temproot=temproot->right;
			}else{
				temproot->right = new slabbstx(xcoord, ind);
				return;
			}
		}
		if(xcoord==temproot->xcoord){ 
			bool indcheck(true);
			for(int i=0;i<temproot->index.size();i++){
				if(root1vs2){ if(dcelist.at(temproot->index.at(i))->p==p0) indcheck=false; }else{  if(veinvtcs.at(temproot->index.at(i))->p==p0) indcheck=false; }
			}
			if(indcheck){
				temproot->index.push_back(ind);
			}
			return;
		}
	}
}

void slabbstxaddxf(bool root1vs2)
{
	slabbstx *tempi;
	slabbstx *tempf, *tempmin;
	double xi,xf;//the left (initial) side (x) of the slab
	int size;
	if(root1vs2) { size = dcelist.size(); }else{ size = veinvtcs.size(); }
	for(int i=0;i<size;i++){
		if(root1vs2){
		if(dcelist[i]){
		tempi = slabbstxroot;
	//	cout<<"TEST"<<endl;
		xi = dcelist.at(i)->p.x;
	//	cout<<"Problem"<<endl;
		tempf=0;	
		while(tempi){
			if(xi<tempi->xcoord){
				tempf=tempi;
				tempi=tempi->left;
			}
			if(xi>tempi->xcoord)
				tempi=tempi->right;
			if(xi==tempi->xcoord){
				//if(!tempi->next){
				if(tempi->right){
					tempmin=tempi->right;
					while(tempmin->left){
						tempmin=tempmin->left;
					}
					tempi->next=tempmin;
					tempi->next->prev=tempi;
					tempi->xcoordf=tempmin->xcoord;
					break;
				}else{
					if(tempf){
						tempi->next=tempf;
						tempi->next->prev=tempi;
						tempi->xcoordf=tempf->xcoord;
						break;
					}else{
						tempi->next=0;	
						break;
					}
				}
				//}
			}
		}
		}
		}else{
		if(veinvtcs[i]){
		tempi = slabbstxroot2;
	//	cout<<"TEST"<<endl;
		xi = veinvtcs[i]->p.x;
	//	cout<<"Problem"<<endl;
		tempf=0;	
		while(tempi){
			if(xi<tempi->xcoord){
				tempf=tempi;
				tempi=tempi->left;
			}
			if(xi>tempi->xcoord)
				tempi=tempi->right;
			if(xi==tempi->xcoord){
				//if(!tempi->next){
				if(tempi->right){
					tempmin=tempi->right;
					while(tempmin->left){
						tempmin=tempmin->left;
					}
					tempi->next=tempmin;
					tempi->next->prev=tempi;
					tempi->xcoordf=tempmin->xcoord;
					break;
				}else{
					if(tempf){
						tempi->next=tempf;
						tempi->next->prev=tempi;
						tempi->xcoordf=tempf->xcoord;
						break;
					}else{
						tempi->next=0;	
						break;
					}
				}
				//}
			}
		}
		}
		}
	}
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//The functions below are poorly commented --- I need to improve this
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DOUBLY CONNECTED EDGE LIST FUNCTION SET START HERE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void addvert(point p){
	vert *temp= new vert(p);
	bstdcelinsert(p.y,dcelist.size()-1,p);
};

/////////bool del parameter tells you whether the node of tree will be deleted or not -> del=true means the node will be deleted and del=false means it will not
/////////must take into account the fact that now there might be null elements inside the dcelist vector/////////////////
bool findinst(point p){//find if the point is NOT in the dcelist vector /// true->not in the list .... false->in the list
	if(dcelist.empty()) return true;
	bstdcelnode *temproot;
	temproot=bstdcelroot;
	while(temproot){
		if(p.y>temproot->ycoord){
			if(temproot->right){
				temproot=temproot->right;
			}else
				return true;
		}
		if(p.y<temproot->ycoord){
			if(temproot->left){
				temproot=temproot->left;
			}else
				return true;
		}
		if(p.y==temproot->ycoord) { break; }
	}
	for(int i=0;i<temproot->index.size();i++){
		if(dcelist[temproot->index[i]]){ //fix this!!!!!!!!!!!!! THINK!!!! IT MIGHT WORK FOR THE FIRST TIME YOU CLICK IT, BUT IT SURE WON'T WORK AFTERWARDS!!!
		if(p==dcelist[temproot->index[i]]->p){
//		if(rounding(p.x,6)==rounding(dcelist.at(temproot->index.at(i))->p.x,6) && rounding(p.y,6)==rounding(dcelist.at(temproot->index.at(i))->p.y,6) ){
			iterdcel=temproot->index.at(i);//iterdcel returns the index of the vertex in the dcel vector
			return false;
		}
		}else{
	//		cout<<"blablabls"<<endl;	
			temproot->index.erase(temproot->index.begin()+i);
			return true;
		}
	}
	return true;
};

//I NEED TO OPTIMIZE THIS FUNCTION instead of linear search I need to implement a bst
int findfaceit(hedge *temp){//find the index location of the point in the dcelist
	if(temp->f){
	for(int i=0;i<dcefacelist.size();i++){
	        if(dcefacelist.at(i)==temp->f){
			return i;
		};
	};
	return -1;
	}else{
	hedge *temp2;
	temp2=temp->next;
	while(!temp2->f){
		temp2=temp2->next;
	};
	for(int i=0;i<dcefacelist.size();i++){
		 if(dcefacelist.at(i)==temp2->f){
			return i;
		 };
	};
	return -1;
	};
};

void addface(hedge *temp){
//	cout<<
	if(!flag){
	temp->f=new face(temp);
	dcefacelist.push_back(temp->f);
//	++facecounter;
	}else{
	int b=findfaceit(temp);
//	cout<<"Índice da face no vetor dcefacelist: "<<b<<endl;
//	cout<<"Tentativa de deletar face."<<endl;
	if(b!=-1){
//		if(dcefacelist[b]==dcefacelist[b]->edge->f){
		hedge *temp1,*temp2;
		temp1=temp;
		temp2=temp1->next;
//		cout<<"Origin: "<<temp1->origin->p.x<<" "<<temp1->origin->p.y<<endl;
//		delete temp1->f;
		temp1->f=0;
		while(temp2!=temp1){
//			cout<<"Origin: "<<temp2->origin->p.x<<" "<<temp2->origin->p.y<<endl;
//			if(!temp2->next) break;
	//		delete temp2->f;
			temp2->f=0;
			temp2=temp2->next;	
		}
		temp1=temp->next->twin;
		temp2=temp1->next;
//		delete temp1->f;
		temp1->f=0;
		while(temp2!=temp1){
//			cout<<"Origin: "<<temp2->origin->p.x<<" "<<temp2->origin->p.y<<endl;
			if(!temp2->next) break;
	//		delete temp2->f;
			temp2->f=0;
			temp2=temp2->next;	
		}
		temp2=temp1->prev;
//		delete temp1->f;
		temp1->f=0;
		while(temp2!=temp1){
//			cout<<"Origin: "<<temp2->origin->p.x<<" "<<temp2->origin->p.y<<endl;
			if(!temp2->prev) break;
	//		delete temp2->f;
			temp2->f=0;
			temp2=temp2->prev;	
		}
//		delete dcefacelist[b];
		dcefacelist[b]=0;
//		cout<<"face deletada"<<endl;
		dcefacelist.erase(dcefacelist.begin()+b);
//		delete temp->f;
//		temp->f=0;
//		}
	}else{
//		cout<<"face nao deletada"<<endl;
	}
	temp->f=new face(temp);
	dcefacelist.push_back(temp->f);
	};
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~Fundamental function of the dcelist construction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Given a half-edge and a point, what is the most 'counterclockwise' halfedge and what is the most 'clockwise' halfedge?
void getcounterclockhe(hedge *temp, int i,bool change){//CONCLUIDA
	double angle,angle2,angle3;
	double x,y,a,b,c,d;
//	cout<<"get1"<<endl;
	hedge *temp2,*temp3,*temp4;
	temp2=temp;
	temp3=temp;
	temp4=temp;
//	cout<<"get2"<<endl;
	if(change){
	a=output.at(i)->end.x;
	b=output.at(i)->end.y;
	c=output.at(i)->start.x;
	d=output.at(i)->start.y;
	}else{
	a=output.at(i)->start.x;
	b=output.at(i)->start.y;
	c=output.at(i)->end.x;
	d=output.at(i)->end.y;
	};
//	cout<<"get3"<<endl;
	x=a-c;
	y=b-d;
//	cout<<"get31"<<endl;
	angle=atan2(y,x);

	x=((temp->origin->p.x-c)*(-1)*cos(angle))+((-1)*(temp->origin->p.y-d)*sin(angle));

	y=((temp->origin->p.y-d)*(-1)*cos(angle))+((temp->origin->p.x-c)*sin(angle));
	angle2=atan2(y,x);
	angle3=angle2;
//	cout<<"get4"<<endl;
	while(temp2->next!=NULL){
		if(temp2->next->twin!=temp){
		temp2=temp2->next->twin;
		x=((temp2->origin->p.x-c)*(-1)*cos(angle))+((-1)*(temp2->origin->p.y-d)*sin(angle));
		y=((temp2->origin->p.y-d)*(-1)*cos(angle))+((temp2->origin->p.x-c)*sin(angle));
		if(atan2(y,x)<angle2){
			angle2=atan2(y,x);
			temp3=temp2;			
		};
		}else{
		break;
		};
	};
//	cout<<"get5"<<endl;
	temp2=temp;
	cche=temp3;
	while(temp2->twin->prev!=NULL){
		if(temp2->twin->prev!=temp){
		temp2=temp2->twin->prev;
		x=((temp2->origin->p.x-c)*(-1)*cos(angle))+((-1)*(temp2->origin->p.y-d)*sin(angle));
		y=((temp2->origin->p.y-d)*(-1)*cos(angle))+((temp2->origin->p.x-c)*sin(angle));
		if(atan2(y,x)>angle3){
			angle3=atan2(y,x);
			temp4=temp2;
		};	
		}else{
		break;
		};
	};
//	cout<<"get6"<<endl;
	che=temp4->twin;
};



int createdcel(){
	int it1,it2,eulercheck;
	bool b1,b2,valid,lala;
	hedge *temp1,*temp2, *temp3, *cchestart, *chestart, *ccheend, *cheend;
	for(int i=0;i<output.size();i++){
		b1=b2=false;
//		cout<<"flag1"<<endl;
		if(findinst(output.at(i)->start)){
			addvert(output.at(i)->start);	
 			it1=dcelist.size()-1;
			dcelist.at(it1)->edge=new hedge();
			dcelist.at(it1)->edge->origin=dcelist.at(it1);	
			b1=true;
		}else{
			it1=iterdcel;
		};
//		cout<<"flag2"<<endl;
		if(findinst(output.at(i)->end)){
			addvert(output.at(i)->end);	
			it2=dcelist.size()-1;
			dcelist.at(it2)->edge=new hedge();
			dcelist.at(it2)->edge->origin=dcelist.at(it2);
			b2=true;
		}else{
			it2=iterdcel;
		};
//			cout<<"output[i]: "<<output[i]->start.x<<" "<<output[i]->start.y<<" "<<output[i]->end.x<<" "<<output[i]->end.y<<endl;	
//			if(output[i]->start==output[i]->end){cout<<"OMG!!!WTF?!!!"<<endl;}
//			cout<<"it1 and it2 "<<it1<<" "<<it2<<" "<<dcelist[it1]->p.x<<" "<<dcelist[it1]->p.y<<" "<<dcelist[it2]->p.x<<" "<<dcelist[it2]->p.y<<endl;
//		cout<<"flag3"<<endl;
/////////////////////////////////////////////TWINS//////////////////////////////////////////////////////////////////////
		if(b1&&b2){
			dcelist.at(it1)->edge->twin=dcelist.at(it2)->edge;
			dcelist.at(it2)->edge->twin=dcelist.at(it1)->edge;
		};
		if(!b1&&b2){
			dcelist.at(it2)->edge->twin=new hedge();
			temp1=dcelist.at(it1)->edge->twin;
			getcounterclockhe(temp1, i,true);
			dcelist.at(it2)->edge->twin->prev=cche;
			dcelist.at(it2)->edge->twin->prev->next=dcelist.at(it2)->edge->twin;
			dcelist.at(it2)->edge->twin->twin=dcelist.at(it2)->edge;
			dcelist.at(it2)->edge->twin->origin=dcelist.at(it1);
			dcelist.at(it2)->edge->next=che;
			dcelist.at(it2)->edge->next->prev=dcelist.at(it2)->edge;
//			dcelist[it2]->edge->prev=dcelist[it2]->edge->twin;
//			dcelist[it2]->edge->prev->next=dcelist[it2]->edge;
		};
		if(b1&&!b2){
			dcelist.at(it1)->edge->twin=new hedge();
			temp1=dcelist.at(it2)->edge->twin;
			getcounterclockhe(temp1, i,false);
			dcelist.at(it1)->edge->twin->prev=cche;
			dcelist.at(it1)->edge->twin->prev->next=dcelist.at(it1)->edge->twin;
			dcelist.at(it1)->edge->twin->twin=dcelist.at(it1)->edge;
			dcelist.at(it1)->edge->twin->origin=dcelist.at(it2);
			dcelist.at(it1)->edge->next=che;
			dcelist.at(it1)->edge->next->prev=dcelist.at(it1)->edge;
//			dcelist[it1]->edge->prev=dcelist[it1]->edge->twin;
//			dcelist[it1]->edge->prev->next=dcelist[it1]->edge;
		};
//		cout<<"flag5"<<endl;
		if(!b1&&!b2){
			temp1=dcelist.at(it1)->edge->twin;		//need to get cche	
			getcounterclockhe(temp1,i,true);//peguei o cche e o che do starting point
			cchestart=cche;
			chestart=che;
			temp2=dcelist.at(it2)->edge->twin;			//need to get cche
			getcounterclockhe(temp2,i,false);//peguei o cche e o che do end point
			ccheend=cche;
			cheend=che;
			cchestart->next= new hedge();
			cchestart->next->prev=cchestart;
			cchestart->next->next=cheend;//need to change to cche
			cchestart->next->next->prev=cchestart->next;
			cchestart->next->origin=dcelist.at(it1);
			cchestart->next->twin= new hedge();
			cchestart->next->twin->twin=cchestart->next;
			cchestart->next->twin->origin=dcelist.at(it2);
			cchestart->next->twin->prev=ccheend;
			cchestart->next->twin->prev->next=cchestart->next->twin;
			cchestart->next->twin->next=chestart;
			cchestart->next->twin->next->prev=cchestart->next->twin;	
		};
	};
//	cout<<"faces to be added"<<endl;
	bool stopthisshit(0);
	for(int i=0;i<dcelist.size();i++){
//		cout<<"step1"<<endl;
		if(!dcelist[i]->edge->f){
			dcelist[i]->edge->f=new face(dcelist[i]->edge);
			dcefacelist.push_back(dcelist[i]->edge->f);
//		cout<<"step11"<<endl;
			for(hedge *j=dcelist[i]->edge->next;j!=dcelist[i]->edge;j=j->next){
				j->f=j->prev->f;
			}
//		cout<<"step12"<<endl;
		}
//		cout<<"step2"<<endl;
		for(hedge *k=dcelist[i]->edge->prev->twin;k!=dcelist[i]->edge;k=k->prev->twin){
			if(!k->f){
				k->f=new face(k);
				dcefacelist.push_back(k->f);
//		cout<<"step21"<<endl;
				if(k->next){
				for(hedge *j=k->next;j!=k;j=j->next){
					if(j->prev&&j&&j->next){
					j->f=j->prev->f;
				}else{
				cout<<i<<endl;		stopthisshit=1; break;
				}
				}
				}
//		cout<<"step22"<<endl;
			}
		if(stopthisshit) break;
		}
		if(stopthisshit) break;
//		cout<<"step3"<<endl;
	}
//	cout<<"faces to be added ok"<<endl;
//	cout<<"Dcel completed. Now displaying each detected face vertices: "<<endl;
//	for(int i=0;i<dcefacelist.size();i++){
//		cout<<"Face "<<i<<endl;
//		temp1=dcefacelist[i]->edge;
//		cout<<"Origin: "<<temp1->origin->p.x<<" "<<temp1->origin->p.y<<endl;
///		temp2=temp1->next;
//		while(temp2!=temp1){
//			cout<<"Origin: "<<temp2->origin->p.x<<" "<<temp2->origin->p.y<<endl;
//			temp2=temp2->next;
//		}
//	}
//	cout<<dcefacelist.size()<<" faces were created. "<<endl;
//	cout<<"There are "<<dcelist.size()<<" vertices."<<endl;
	eulercheck=dcelist.size()-output.size()+dcefacelist.size();
//	cout<<"Euler formula check: "<< dcelist.size()-output.size()+dcefacelist.size()<<endl;
	if(eulercheck==2){
	for(int i=0;i<output.size();i++){
		delete output[i];
		output[i]=0;
	}
	output.clear();
	}
	return eulercheck;
};
