#include <iostream>
#include <cmath>
#include <vector>
#include "covering.h"
#include <random>

using namespace std;

double pi = 3.141593;

int main()
{
	
	int numIter = 5000;
	
	default_random_engine generator;
	uniform_real_distribution<double> distribution(1.0, 3.0);
	
	for (int ii = 0 ; ii <= numIter ; ii ++)
	{
		double factor = (11.0/12.0)*(double(ii)/double(numIter));
		int count = 0;

		for(int i = 0 ; i < numIter ; i ++)
		{
			vector<double> diskList;

			double areaYet = 0.0;
			double lam = distribution(generator);
			uniform_real_distribution<double> diskDistribution(0.0, sqrt(1 + lam*lam)/2.0);
			double r = diskDistribution(generator);
			while (areaYet + r*r < factor*lam)
			{
				areaYet += r*r;
				diskList.push_back(r);
				r = diskDistribution(generator);
			}
			diskList.push_back(sqrt(factor*lam - areaYet + .00001));
			sort(diskList.begin(), diskList.end());
			reverse(diskList.begin(), diskList.end());

			vector<pair<double, Point> > answer;
			Point a(0,0), b(0,lam), c(1, lam), d(1, 0);
			//Point a(0,0);
			count += good_rectangle_cover(diskList, a, b, c, d, answer);
		}
		cout<<factor*pi<<" "<<double(count)/double(numIter)<<endl;
	}


	return 0;
}


