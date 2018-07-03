#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double greedySplit(vector<double> & input, vector<double> & out1, vector<double> & out2)
{
	out1.resize(0);
	out2.resize(0);
	double area1 = 0.0, area2 = 0.0;
	for (int i = 0 ; i < input.size() ; i ++)
	{
		if (area1 <= area2)
		{
			out1.push_back(input[i]);
			area1 += input[i]*input[i];
		}
		else
		{
			out2.push_back(input[i]);
			area2 += input[i]*input[i];
		}
	}

	if (area2 > area1)
	{
		vector<double> temp = out1;
		out1 = out2;
		out2 = temp;
		return area2/area1;
	}

	return area1/area2;
}

int main()
{
	vector<double> input;
	int numIter = 3;
	for (int i = 0 ; i < numIter ; i ++)
	{
		input.push_back(numIter - i);
	}
	vector<double> a, b;
	double ratio = greedySplit(input, a, b);
	cout<<ratio<<endl;



	return 0;
}
