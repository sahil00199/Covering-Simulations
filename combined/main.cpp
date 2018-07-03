#include <iostream>
#include <cmath>
#include <vector>
#include "covering.h"

using namespace std;

int main()
{
	vector<double> disks;
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



	return 0;
}
