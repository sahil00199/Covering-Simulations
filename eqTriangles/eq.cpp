#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

struct Point
{	
	double x, y;

	void print()
	{
		cout<<"("<<x<<","<<y<<") ";

		return;
	}

	Point()
	{
		x = -1;
		y = -1;
	}

	Point (const Point & p)
	{
		x = p.x;
		y = p.y;
	}

	
};

double greedySplit(vector<double> & input, vector<double> & out1, vector<double> & out2)
{
	//assumes the input to be sorted, but that does not really matter
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


bool eq(vector<double> & diskList, vector<pair<double, Point> > & output, Point point0, double scale) 
//assuming sorted input
{
	
	output.resize(0);

	if (diskList.size() == 0 || (diskList.size() == 1 && diskList[0] < scale/sqrt(3.0)) ) 
	{
		return false;
	}
	
	/*for (int i = 0 ; i < diskList.size() ; i ++)
	{
		pair<double, Point> temp;
		temp.first = diskList[i];
		temp.second.x = -1;
		temp.second.y = -1;
		output.push_back(temp);
	}
	*/

	double r1 = diskList[0];
	
	if (r1 >= scale/sqrt(3.0)) //case 1
	{
		pair<double, Point> temp;
		output.resize(1);
		output[0].first = r1;
		output[0].second.x = 0.5*scale;
		output[0].second.y = (sqrt(3)/6.0)*scale;
		//output.push_back(temp);
	}

	else if(r1 >= 0.5) //case 2
	{
		//recursion
		vector<double> newList(diskList.size() - 1);
		vector<pair<double, Point> > newOutput;
		Point newOrigin;
		double x = (0.5*scale + sqrt(3.0)*sqrt(r1*r1 - 0.25*scale*scale));
		double newScale = scale - x;
		newOrigin.x = x/2.0;
		newOrigin.y = sqrt(3.0)*x/2.0;
		for (int i = 1 ; i < diskList.size() ; i ++)
		{
			newList[i-1] = diskList[i];
		}
		eq(newList, newOutput, newOrigin, newScale);

		//entering the first disk in the output
		pair<double, Point> temp;
		temp.first = r1;
		temp.second.x = 0.5*scale;
		temp.second.y = sqrt(r1*r1 - 0.25*scale*scale);
		//cout<<newOutput.size()<<endl;cout<<newOutput[0].first<<endl; newOutput[0].second.print(); cout<<endl;
		output.push_back(temp);

		//copy the recursed answers back to the final answer
		for (int i = 0 ; i < newOutput.size() ; i++)
		{
			output.push_back(newOutput[i]);
		}
	}

	for (int i = 0 ; i < output.size() ; i ++)
	{
		//output[i].first *= scale;
		output[i].second.x += point0.x;
		output[i].second.y += point0.y;
		//output[i].second.print(); cout<<endl;
	}

	return true;

}


int main()
{
	vector<double> disks;
	disks.push_back(0.51);
	disks.push_back(0.3);
	disks.push_back(0.3);
	vector<pair<double, Point> > answer;
	Point origin;
	origin.x = 0.0;
	origin.y = 0.0;
	eq(disks, answer, origin, 1.0);
	for (int i = 0 ; i < answer.size() ; i ++)
	{
		cout<<answer[i].first<<"\t";
		answer[i].second.print();
		cout<<endl;
	}



	return 0;
}

