#include <iostream>
#include <vector>
#include <cmath>
#include "covering.h"

using namespace std;


bool right(vector<double> & diskList, vector<pair<double, Point> > & output, Point point0, double scale, int configuration)
{
	//base cases
	//no disk
	if (diskList.size() == 0) return false;
	
	double r1 = diskList[0];
	
	//1 disk
	if (diskList.size() == 1)
	{
		if (r1 < scale/sqrt(2.0)) return false;
		else
		{
			output.resize(0);
			output.push_back(make_pair(r1, Point(scale/2.0, scale/2.0)));
		}
	}

	//2 disks
	else if (diskList.size() == 2)
	{
		double r2 = diskList[1];
		if (r1*r1 + r2*r2 < scale*scale*0.5) return false;
		else
		{
			output.resize(0);
			output.push_back(make_pair(r1, Point(0.5*scale, sqrt(r1*r1 - 0.25*scale*scale))));
			output.push_back(make_pair(r2, Point(scale, 0.5*scale + sqrt(r1*r1 - 0.25*scale*scale))));
		}
	}

	//other cases
	else
	{
		//case 1
		if (r1 >= scale/sqrt(2.0))
		{
			output.resize(0);
			output.push_back(make_pair(r1, Point(scale/2.0, scale/2.0)));
		}

		//case 2
		else if (r1 >= scale/2.0)
		{

		}
	}

	//adjust for non-standard configuration
	if (configuration == 1)
	{
		for (int i = 0 ; i < output.size() ; i ++)
		{
			output[i].second.x -= scale;
			double newX = (output[i].second.y - output[i].second.x)/sqrt(2.0);
			double newY = (-output[i].second.x - output[i].second.y)/sqrt(2.0);
			output[i].second.x = newX;
			output[i].second.y - newY;
		}
	}

	//adjust for point0
	for (int i = 0 ; i < output.size() ; i ++)
	{
		output[i].second.x += point0.x;
		output[i].second.y += point0.y;
	}


	return true;
}
