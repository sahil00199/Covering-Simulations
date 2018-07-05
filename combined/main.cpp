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
	uniform_real_distribution<double> distribution(0.0,0.7);

	for (int ii = 0 ; ii < 5000 ; ii ++)
	{
		
		vector<double> diskList;

		double areaYet = 0.0;
		double r = distribution(generator);
		while (areaYet + r*r < 0.5)
		{
			areaYet += r*r;
			diskList.push_back(r);
			r = distribution(generator);
		}
		diskList.push_back(sqrt(0.5 - areaYet + .0001));
		sort(diskList.begin(), diskList.end());
		reverse(diskList.begin(), diskList.end());

		vector<pair<double, Point> > answer;
		//Point a(0,0), b(0,lam), c(1, lam), d(1, 0);
		Point a(0,0);
		bool possible = right(diskList, answer, a, 1.0, 0);
		bool correct = checkRight(answer, a, 1.0);
		if (!(possible && correct))
		{
			cout<<possible<<" "<<correct<<"\n";
			for (int i = 0 ; i < diskList.size() ; i ++)
			{
				cout<<diskList[i]<<" ";
			}
			cout<<endl;
		}
		if (ii%10 == 1) cout<<ii<<endl;
	}


	return 0;
}


