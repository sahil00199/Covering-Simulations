#include <iostream>
#include <cmath>
#include <vector>
#include "covering.h"
#include <random>

using namespace std;

int main()
{
	/*vector<double> disks;
	disks.push_back(0.51);
	disks.push_back(0.03);
	disks.push_back(0.03);
	vector<pair<double, Point> > answer;
	bool possible = coverEqTriangles(disks, answer);
	if (possible)
		for (int i = 0 ; i < answer.size() ; i ++)
		{
			cout<<answer[i].first<<"\t";
			answer[i].second.print();
			cout<<endl;
		}
	else
		cout<<"not possible as per the algorithm\n";


*/


	double delta = 0.2;

	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,sqrt(77.0)/6.0);

	for (double lam = 1.0 ; lam <= 10.0 ; lam += 0.005)
	{
		
		vector<double> diskList;

		double areaYet = 0.0;
		double r = distribution(generator);
		while (areaYet + r*r < 11.0*lam/12.0)
		{
			areaYet += r*r;
			diskList.push_back(r);
			r = distribution(generator);
		}
		diskList.push_back(sqrt(11.0*lam/12.0 - areaYet + .0001));
		sort(diskList.begin(), diskList.end());
		reverse(diskList.begin(), diskList.end());

		vector<pair<double, Point> > answer;
		//Point a(0,0), b(0,lam), c(1, lam), d(1, 0);
		Point a(0,0), b(0, 1), c(lam, 1), d(lam, 0);
		bool possible = bad_rectangle_cover_bounded_radii(diskList, a, b, c, d, answer);
		bool correct = checkRectangle(answer, a, b, c, d);
		if (!(possible && correct))
		{
			cout<<possible<<" "<<correct<<" ";
			cout<<lam<<endl;
			for (int i = 0 ; i < diskList.size() ; i ++)
			{
				cout<<diskList[i]<<" ";
			}
			cout<<endl;
		}
		cout<<lam<<endl;
	}


	return 0;
}


