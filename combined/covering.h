#include <iostream>

using namespace std;

struct Point
{	
	double x, y;

	void print();
	Point();
	Point(double xx, double yy);
	Point (const Point & p);
};

bool good_rectangle_cover(vector<double> & radii,Point lowleft,Point upleft,Point upright,Point lowright,vector<pair<double,Point> > & ans);

bool bad_rectangle_cover_bounded_radii(vector<double> & diskList, Point bl, Point tl, Point tr, Point br
	, vector<pair<double, Point> > & output);

bool eq(vector<double> & diskList, vector<pair<double, Point> > & output, Point point0, double scale);

bool coverEqTriangles(vector<double> & diskList, vector<pair<double, Point> > & output);

double greedySplit(vector<double> & input, vector<double> & out1, vector<double> & out2);

bool checkRectangle(vector<pair<double, Point> > & solution, Point bl, Point tl, Point tr, Point br);

//bool checkEq(vector<pait<double, Point> > & solution, Point point0, double scale);
