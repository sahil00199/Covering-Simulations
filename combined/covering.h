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

bool eq(vector<double> & diskList, vector<pair<double, Point> > & output, Point point0, double scale);

bool coverEqTriangles(vector<double> & diskList, vector<pair<double, Point> > & output);

double greedySplit(vector<double> & input, vector<double> & out1, vector<double> & out2);

