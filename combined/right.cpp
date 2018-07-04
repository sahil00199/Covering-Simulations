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
			//recurse on a right triangle in configuration 1
			Point newOrigin (scale/2.0 + sqrt(r1*r1 - 0.25*scale*scale), scale/2.0 + sqrt(r1*r1 - 0.25*scale*scale));
			double newScale = (scale - 2.0*sqrt(r1*r1 - 0.25*scale*scale))/sqrt(2.0);
			vector<double> newDiskList(diskList.size() - 1);
			for (int i = 1 ; i < diskList.size() ; i ++)
			{
				newDiskList[i-1] = diskList[i];
			}
			vector<pair<double, Point> > answer;
			bool possible = right(newDiskList, answer, newOrigin, newScale, 1);

			if (!possible) return false;

			//place the first disk and then copy the result of the recursion
			output.resize(0);
			output.push_back(make_pair(r1, Point(scale/2.0, sqrt(r1*r1 - 0.25*scale*scale))));
			for (int i = 0 ; i < answer.size() ; i ++)
			{
				output.push_back(answer[i]);
			}
		}

		//case 3&4
		else if (r1 >= scale*0.5/sqrt(2.0))
		{
			double r2 = diskList[1];
			//case 3
			if (r1 + r2 < scale/sqrt(2.0)) 
			{
				//create the 2 disk sets for recursion
				vector<double> new1, new2;
				double ratio = greedySplit(diskList, new1, new2);
				if (ratio < 0) return false;
				vector<pair<double, Point> > answer1, answer2;

				//create the 2 triangles
				double newScale = scale - sqrt(2.0)*r1;
				Point origin1(0.0, 0.0), origin2(sqrt(2.0)*r1, sqrt(2.0)*r1);

				//call the 2 recursive functions
				bool possible = right(new1, answer1, origin1, newScale, 0);
				possible = possible && right(new2, answer2, origin2, newScale, 0);

				if (!possible) return false;

				//place the first disk and copy the others into the output
				output.resize(0);
				output.push_back(make_pair(r1, Point(scale - r1/sqrt(2.0), r1/sqrt(2.0))));
				for (int i = 0 ; i < answer1.size() ; i ++)
				{
					output.push_back(answer1[i]);
				}
				for (int i = 0 ; i < answer2.size() ; i ++)
				{
					output.push_back(answer2[i]);
				}
			}

			//case 4
			else
			{
				//define all the angles
				double cos1 = (0.5*scale*scale + r1*r1 - r2*r2)/(sqrt(2.0)*r1*scale);
				double cos2 = (0.5*scale*scale + r2*r2 - r1*r1)/(sqrt(2.0)*r2*scale);
				double sin1 = sqrt(1 - cos1*cos1);
				double sin2 = sqrt(1 - cos2*cos2);
				
				//define the coordinates of the rectangle to recurse on 
				Point p1(sqrt(2.0)*r1*(cos1 + sin1), 0.0),
					p2(sqrt(2.0)*r1*(cos1 + sin1), sqrt(2.0)*r1*(cos1 - sin1)),
					p3(scale, sqrt(2.0)*r1*(cos1 - sin1)),
					p4(scale, 0.0);

				//create the diskList
				vector<double> newDiskList(diskList.size() - 2);
				for (int i = 2 ; i < diskList.size() ; i ++)
				{
					newDiskList[i-2] = diskList[i];
				}
				vector<pair<double, Point> > answer;

				//make the recursive call
				double width = scale - sqrt(2.0)*r1*(cos1 + sin1);
				double height = sqrt(2.0)*r1*(cos1 - sin1);
				bool possible;
				if (height/width <= 3)
				{
					possible = good_rectangle_cover(newDiskList, p1, p2, p3, p4, answer);
				}
				else
				{
					possible = bad_rectangle_cover(newDiskList, p1, p2, p3, p4, answer);
				}
				if (!possible) return false;

				//place the first 2 disks and then copy the remaining stuff
				output.resize(0);
				output.push_back(make_pair(r1, Point( r1*(sin1 + cos1)/sqrt(2.0) , r1*(cos1 - sin1)/sqrt(2.0) )));
				output.push_back(make_pair(r2, Point( r1*(sin1 + cos1)/sqrt(2.0) + (r1*cos1 + r2*cos2)/sqrt(2.0)
					 , r1*(cos1 - sin1)/sqrt(2.0) + (r1*cos1 + r2*cos2)/sqrt(2.0) )));
				for (int i = 0 ; i < answer.size() ; i ++)
				{
					output.push_back(answer[i]);
				}
			}
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
			output[i].second.y = newY;
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
