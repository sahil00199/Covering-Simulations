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

	Point(double xx, double yy)
	{
		x = xx;
		y = yy;
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
		bool possible = eq(newList, newOutput, newOrigin, newScale);

		if (! possible) return false;

		//entering the first disk in the output
		pair<double, Point> temp;
		temp.first = r1;
		temp.second.x = 0.5*scale;
		temp.second.y = sqrt(r1*r1 - 0.25*scale*scale);
		output.push_back(temp);

		//copy the recursed answers back to the final answer
		for (int i = 0 ; i < newOutput.size() ; i++)
		{
			output.push_back(newOutput[i]);
		}
	}

	else if (r1/scale >= (3*sqrt(3.0) + sqrt(10.0))/17.0) //case 6
	{
		double height = (sqrt(3.0)/2.0)*(scale - sqrt(3.0)*r1);
		double width = scale/2.0 - sqrt(sqrt(3.0)*scale*r1 - 3.0*scale*scale/4.0);

		double r2 = diskList[1]; //we know that there are 2 disks, else it would have 
		//returned false earlier

		if (r2 < 0.1111) //case 6a (HALF-SPLIT to cover the 2 rectangles)
		{
			vector<double> new1, new2;
			vector<double>newDiskList(diskList.size() - 1);
			for (int i = 1 ; i < diskList.size() ; i ++)
			{
				newDiskList[i-1] = diskList[i];
			}

			//recursion
			double ratio = greedySplit(newDiskList, new1, new2);
			if (ratio < 0) return false;
			Point l1(0.0, 0.0), l2(0.0, height), l3(width, height), l4(width, 0.0);
			Point ri1(scale - width, 0.0), ri2(scale - width, height);
			Point ri3(scale, height), ri4(scale, 0.0);
			vector<pair<double, Point> > answer1, answer2;
			bool possible = good_rectangle_cover(new1, l1, l2, l3, l4, answer1);
			possible = possible && good_rectangle_cover(new2, ri1, ri2, ri3, ri4, answer2);

			if (!possible) return false;

			//now place the first disk and copy the outputs of the recursive calls
			output.resize(0);
			output.push_back(make_pair( r1 , Point(scale/2.0, sqrt(3.0)*scale/2.0 - r1) ));
			for (int i = 0 ; i < answer1.size() ; i ++)
			{
				output.push_back(answer1[i]);
			}
			for (int i = 0 ; i < answer2.size() ; i ++)
			{
				output.push_back(answer2[i]);
			}
		}

		else if (r2 < 0.47816*scale) //case 6b
		{
			//recursion
			vector<double> newDiskList(diskList.size() - 2);
			for (int i = 2 ; i < diskList.size() ; i ++)
			{
				newDiskList[i-2] = diskList[i];
			}
			Point ri1(scale - width, 0.0), ri2(scale - width, height);
			Point ri3(scale, height), ri4(scale, 0.0);
			vector<pair<double, Point> > answer;
			bool possible = good_rectangle_cover(newDiskList, ri1, ri2, ri3, ri4, answer);

			if (!possible) return false;

			//now place the first 2 disks and then copy the output of the recursive call
			output.resize(0);
			output.push_back(make_pair(r1, Point(scale/2.0, sqrt(3.0)*scale/2.0 - r1)));
			output.push_back(make_pair(r2, Point(r2, 0.0)));
			for (int i = 0 ; i < answer.size() ; i ++)
			{
				output.push_back(answer[i]);
			}
		}

		else
		{
			//recursion on the tiny eq triangle at the bottom right
			double eps = 0.5*scale - r2;
			double newScale = 7.5*eps;
			Point newOrigin(scale - newScale, 0.0);
			vector<double> newDiskList(diskList.size() - 2);
			for (int i = 2 ; i < diskList.size() ; i ++)
			{
				newDiskList[i-2] = diskList[i];
			}
			vector<pair<double, Point> > answer;
			bool possible = eq(newDiskList, answer, newOrigin, newScale);

			if (!possible) return false;

			//now place the firast 2 disks and then copy the output of the recursive call
			output.resize(0);
			output.push_back(make_pair(r1, Point(scale/2.0, sqrt(3.0)*scale/2.0 - r1)));
			output.push_back(make_pair(r2, Point(0.0, r1)));
			for (int i = 0 ; i < answer.size() ; i ++)
			{
				output.push_back(answer[i]);
			}
		}
	}

	else if (r1/scale >= 11.0/16.0 - sqrt(249.0/256.0 - 11.0*sqrt(3.0)/24.0))
	{
		double width = scale;
		double height = sqrt(3.0)*scale/2.0 - 2*r1;

		//recursion 
		vector<pair<double, Point> >answer;
		vector<double> newDiskList(diskList.size() - 1);
		for (int i = 1 ; i < diskList.size() ; i ++) //we know that there are at least 2 disks
		{
			newDiskList[i-1] = diskList[i];
		}
		Point bl(0.0, 0.0), tl(0.0, height), tr(width, height), br(width, 0.0);
		bool possible;

		//check if the rectangle is skewed or not
		if (width/height < 3) //not skewed (width is always more than the height, so not checking the reciprocal)
		{
			possible = good_rectangle_cover(newDiskList, bl, tl, tr, br, answer);
		}

		else
		{
			possible = bad_rectangle_cover(newDiskList, bl, tl, tr, br, answer);
		}

		if (!possible) return false;

		//place the first disk and copy the answer of the recursive call
		output.resize(0);
		output.push_back(make_pair(r1, Point(scale/2.0, sqrt(3.0)*scale/2.0 - r1)));
		for (int i = 0 ; i < answer.size() ; i ++)
		{
			output.push_back(answer[i]);
		}
	}

	else if (r1/scale >= (48.0 - 22.0*sqrt(3.0))/(39.0*sqrt(2.0)))
	{
		//first split the disks into 2 sets
		int index = 1;
		double lCritical = (6.0/sqrt(77.0))*(11.0/16.0 - sqrt(249.0/256.0 - 11.0*sqrt(3.0)/24.0))*scale;
		double area = 0.0;
		vector<double> new1;
		while (area < 11.0*scale*lCritical/12.0)
		{
			if (index == diskList.size()) return false; //there are not enough disks to cover!
			new1.push_back(diskList[index]);
			area += diskList[index]*diskList[index];
			index ++;
		}
		vector<double> new2;
		while (index < diskList.size())
		{
			new2.push_back(diskList[index]);
			index ++;
		}

		//now compute the 2 rectangles lo(wer) and up(per) to recurse on
		double lowerWidth = scale;
		double lowerHeight = 12*area/(11.0*scale);
		Point lo1(0.0,0.0), lo2(0.0, lowerHeight), lo3(lowerWidth, lowerHeight), lo4(lowerWidth, 0.0);
		vector<pair<double, Point> > answer1, answer2;

		double upperHeight = sqrt(3.0)*scale/2.0 - lowerHeight - 2*r1;
		double upperWidth = scale - 2.0*lowerHeight/sqrt(3.0);
		Point upperOrigin(lowerHeight/sqrt(3), lowerHeight);
		Point up1(upperOrigin.x, upperOrigin.y);
		Point up2(upperOrigin.x, upperOrigin.y + upperHeight);
		Point up3(upperOrigin.x + upperWidth, upperOrigin.y + upperHeight);
		Point up4(upperOrigin.x + upperWidth, upperOrigin.y);

		//make the 2 recursive calls
		bool possible = bad_rectangle_cover_bounded_radii(new1, lo1, lo2, lo3, lo4, answer1);
		double skew = upperWidth/upperHeight;
		if (1/3 <= skew && skew <= 3.0)
		{
			possible = possible && good_rectangle_cover(new2, up1, up2, up3, up4, answer2);
		}
		else
		{
			possible = possible && good_rectangle_cover(new2, up1, up2, up3, up4 ,answer2);
		}

		if (!possible) return false;

		//place the first disks and copy the results of the recurisve calls
		output.resize(0);
		output.push_back(make_pair(r1, Point(scale/2.0, sqrt(3.0)*scale/2.0 - r1)));
		for (int i = 0 ; i < answer1.size() ; i++)
		{
			output.push_back(answer1[i]);
		}
		for (int i = 0 ; i < answer2.size() ; i ++)
		{
			output.push_back(answer2[i]);
		}
	}


	else
	{
		double x = sqrt(2.0)*r1/(sqrt(3.0) - 1.0);
		double lam = sqrt(3.0) - 1.0;

		//compute the dimensions of the triangle and the rectangle on the right
		double newScale = scale - lam*x;
		double height = lam*x*sqrt(3.0)/2.0;
		double width = scale - x + lam*x/2.0;

		//form the 2 sets of disks
		double areaAddedYet = 0.0;
		double areaNeeded = 0.5*newScale*newScale;
		vector<double> new1;
		int index = 1;
		while (areaAddedYet < areaNeeded)
		{
			if (index == diskList.size()) return false; //there are not enough disks!
			new1.push_back(diskList[index]);
			areaAddedYet += (diskList[index]*diskList[index]);
			index ++;
		}
		vector<double> new2 ;
		while (index < diskList.size())
		{
			new2.push_back(diskList[index]);
			index ++;
		}

		//make the recursive calls
		Point newOrigin(lam*x/2.0, sqrt(3.0)*lam*x/2.0);
		vector<pair<double, Point> > answer1, answer2;
		bool possible = eq(new1, answer1, newOrigin, newScale);
		Point p1(x - lam*x/2.0, 0.0), p2(x - lam*x/2.0, sqrt(3.0)*lam*x/2.0),
		p3(scale, sqrt(3.0)*lam/2.0), p4(scale, 0.0);
		possible = possible && bad_rectangle_cover_bounded_radii(new2, p1, p2, p3, p4, answer2);

		if (!possible) return false;

		//place the largest disk and copy the output of the recursive calls
		output.resize(0);
		output.push_back(make_pair(r1, Point(x/2.0, sqrt(r1*r1 - x*x/4.0))));
		for (int i = 0 ; i < answer1.size() ; i ++)
		{
			output.push_back(answer1[i]);
		}
		for (int i = 0 ; i < answer2.size() ; i ++)
		{
			output.push_back(answer2[i]);
		}
	}



	for (int i = 0 ; i < output.size() ; i ++)
	{
		output[i].second.x += point0.x;
		output[i].second.y += point0.y;
	}

	return true;

}



bool coverEqTriangles(vector<double> & diskList, vector<pair<double, Point> > & output)
{
	Point p0;
	p0.x = 0.0;
	p0.y = 0.0;
	return eq(diskList, output, p0, 1.0);
}



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

