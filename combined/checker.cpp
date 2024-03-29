#include <iostream>
#include <vector>
#include <cmath>
#include "covering.h"

using namespace std;

double epsilon = 1e-9;
int numIter = 1000;

bool checkRectangle(vector<pair<double, Point> > & solution, Point bl, Point tl, Point tr, Point br)
{
	int n = solution.size();
	for (int i = 0 ; i <= numIter ; i ++)
	{
		for (int j = 0 ; j <= numIter ; j ++)
		{
			double x = bl.x + (br.x - bl.x)*(double(i)/numIter);
			double y = bl.y + (tl.y - bl.y)*(double(j)/numIter);
			bool covered = false;
			for (int k = 0 ; k < n && (!covered) ; k ++)
			{
				if ( (x - solution[k].second.x)*(x - solution[k].second.x) +  (y - solution[k].second.y)
					*(y - solution[k].second.y) <= solution[k].first*solution[k].first + epsilon)
				covered = true;
			}
			if (!covered)
			{
				cout<<"Point that is not covered: "<<x<<" "<<y<<endl;
				return false;
			}
		}
	}
	return true;
}


bool checkEq(vector<pair<double, Point> > & solution, Point point0, double scale)
{
	int n = solution.size();
	for (int i = 0 ; i <= numIter ; i ++)
	{
		double y = point0.y + double(i)*scale*sqrt(3.0)/(2.0*double(numIter));
		int xIter = numIter*(sqrt(3.0)*scale/2.0 - y)*2.0/sqrt(3.0)+1;
		for (int j = 0 ; j <= xIter ; j ++)
		{
			double x = point0.x + scale/2.0 + ((double(j) - double(xIter)/2.0)/double(xIter))*scale*((sqrt(3.0)*scale/2.0 - y)*2.0/(sqrt(3.0)*scale));
			bool covered = false;
			//cout<<x<<" "<<y<<" "<<i<<" "<<xIter<<" "<<j<<endl;
			for (int k = 0 ; k < n && (!covered); k ++)
			{
				if ( (x - solution[k].second.x)*(x - solution[k].second.x) + (y - solution[k].second.y)
				*(y - solution[k].second.y) <= solution[k].first*solution[k].first + epsilon )
					covered = true;
			}
			if (! covered)
			{
				cout<<"A point that is not covered is: ("<<x<<","<<y<<")\n";
				return false;
			}
		}
	}
	return true;
}


bool checkRight(vector<pair<double, Point> > & solution, Point point0, double scale)
{
	int n = solution.size();
	for (int i = 0 ; i <= numIter ; i ++)
	{
		double y = double(i)*scale/double(numIter);
		int xIter = int(double(numIter)*(scale - y)/scale) + 1;
		for (int j = 0 ; j <= xIter ; j ++)
		{
			double x = scale - double(j)*(scale - y)/double(xIter);
			bool covered = false;
			for (int k = 0 ; k < n ; k ++)
			{
				if ( (x - solution[k].second.x)*(x - solution[k].second.x) + (y - solution[k].second.y)
				*(y - solution[k].second.y) <= solution[k].first*solution[k].first + epsilon )
					covered = true;
			}

			if (!covered)
			{
				cout<<"A point that is not covered is: ("<<x<<","<<y<<")\n";
				return false;
			}
		}
	}
	return true;
}
