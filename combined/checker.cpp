#include <iostream>
#include <vector>
#include <cmath>
#include "covering.h"

using namespace std;

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
					*(y - solution[k].second.y) <= solution[k].first*solution[k].first)
				covered = true;
			}
			if (!covered)
			{
				cout<<x<<" "<<y<<endl;
				return false;
			}
		}
	}
	return true;
}
