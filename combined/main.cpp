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
	double lam = 1.6;

	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,sqrt(lam*lam + 1)/2.0);

	cout<<distribution(generator)<<endl;

	vector<double> diskList(5);
	diskList[0] = 0.5;
	diskList[1] = 0.5;
	diskList[2] = 0.5;
	diskList[3] = 0.5;
	diskList[4] = 0.5;
	vector<pair<double, Point> > answer;
	Point a(0,0), b(0,1), c(lam, 1), d(lam, 0);
	bool possible = good_rectangle_cover(diskList, a, b, c, d, answer);
	if (!possible)
	{
		cout<<"Not possible\n";
	}
	else
	{
		cout<<"Possible\n";
		bool correct = checkRectangle(answer, a, b, c, d);
		if (correct) cout<<"Correct!\n";
		else cout<<"Incorrect\n";
	}

	return 0;
}
