#include <iostream>
#include "covering.h"
#include <vector>

using namespace std;

int main()
{
	vector<double> diskList(4);
	diskList[0] = 0.496634;
	diskList[1] = 0.483233;
	diskList[2] = 0.138902;
	diskList[3] = 0.0254341;
	vector<pair<double, Point> > answer;
	Point a(0,0);
	bool possible = eq(diskList, answer, a, 1.0);
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