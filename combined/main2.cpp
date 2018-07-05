#include <iostream>
#include "covering.h"
#include <vector>

using namespace std;

int main()
{
	vector<double> diskList(8);
	diskList[0] = 0.348654;
	diskList[1] = 0.346328;
	diskList[2] = 0.314253;
	diskList[3] = 0.25955;
	diskList[4] = 0.178727;
	diskList[5] = 0.166943;
	diskList[6] = 0.144057;
	diskList[7] = 0.109136;
	vector<pair<double, Point> > answer;
	Point a(0,0);
	bool possible = right(diskList, answer, a, 1.0, 0);
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
	//checkRight(answer, a, 1.0);



	return 0;
}