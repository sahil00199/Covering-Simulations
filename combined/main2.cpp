#include <iostream>
#include "covering.h"
#include <vector>

using namespace std;

int main()
{
	double lam = 1.6
	vector<double> diskList;
	for (int i = 0 ; i < 40 ; i ++)
	{
		diskList.push_back(0.2);
	}
	for (int i = 0 ; i < 20 ; i ++)
	{
		diskList.push_back(0.15);
	}
	vector<pair<double, Point> > answer;
	Point a(0,0), b(0,1), c(lam, 1), d(lam, 0);
	bool possible = bad_rectangle_cover_bounded_radii(diskList, a, b, c, d, answer);
	if (!possible)
	{
		cout<<"Not possible\n";
	}
	else
	{
		for (int i = 0 ; i < answer.size() ; i++)
		{
			cout<<answer[i].first<<" ";
			answer[i].second.print();
			cout<<endl;
		}
	}



	return 0;
}