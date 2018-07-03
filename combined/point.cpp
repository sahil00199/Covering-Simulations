#include <iostream>
#include "covering.h"

using namespace std;

void Point::print()
{
	cout<<"("<<x<<","<<y<<") ";

	return;
}

Point::Point()
{
	x = -1;
	y = -1;
}

Point::Point(double xx, double yy)
{
	x = xx;
	y = yy;
}

Point::Point (const Point & p)
{
	x = p.x;
	y = p.y;
}

