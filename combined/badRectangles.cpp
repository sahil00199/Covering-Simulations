#include <iostream>
#include <vector>
#include "covering.h"

using namespace std;

bool bad_rectangle_cover_bounded_radii(vector<double> & diskList, Point bl, Point tl, Point tr, Point br
	, vector<pair<double, Point> > & output)
{
	double width = br.x - bl.x;
	double height = tl.y - bl.y;

	if (width >= height)
	{
		double skew = width/height;
 		if (skew <= 3)
		{
			return good_rectangle_cover(diskList, bl, tl, tr, br, output);
		}

		else if (skew <= 10.0/3.0) //directly split into 2
		{
			vector<double> new1, new2;
			double ratio = greedySplit(diskList, new1, new2);

			if (ratio < 0) return false;

			double splitX = bl.x + (br.x - bl.x)*(1.0/(1.0 + 1.0/ratio));

			//recurse on the 2 good rectangles
			Point bMid(splitX, bl.y), tMid(splitX, tl.y);
			vector<pair<double, Point> > answer1, answer2;
			bool possible = good_rectangle_cover(new1, bl, tl, tMid, bMid, answer1);
			possible = possible && good_rectangle_cover(new2, bMid, tMid, tr, br, answer2);

			if (!possible) return false;

			output.resize(0);
			for (int i = 0 ; i < answer1.size() ; i ++)
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
			//create a new set and a good rectangle to go with it and recurse on the rest
			vector<double> newDiskSet;
			vector<double> remainingDiskSet;
			vector<pair<double, Point> > answer1, answer2;
			double areaYet = 0.0;
			int index = 0;
			while (areaYet < height*11.0/36.0)
			{
				if (index == diskList.size()) return false;//not enough disks
				newDiskSet.push_back(diskList[index]);
				areaYet += diskList[index]*diskList[index];
				index ++;
			}
			while (index < diskList.size())
			{
				remainingDiskSet.push_back(diskList[index]);
				index ++;
			}

			//compute the rectangle 
			double newWidth = (areaYet*12.0/11.0)/height;
			Point bMid(bl.x + newWidth, bl.y), tMid(tl.x + newWidth, tl.y);

			//make the recursive calls
			bool possible = good_rectangle_cover(newDiskSet, bl, tl, tMid, bMid, answer1);
			possible = possible && bad_rectangle_cover_bounded_radii(remainingDiskSet, bMid, tMid,
				tr, br, answer2);

			if (!possible) return false;

			output.resize(0);
			for (int i = 0 ; i < answer1.size() ; i ++)
			{
				output.push_back(answer1[i]);
			}
			for (int i = 0 ; i < answer2.size() ; i ++)
			{
				output.push_back(answer2[i]);
			}
		}
	}

	else
	{
		double skew = height/width;
 		if (skew <= 3)
		{
			return good_rectangle_cover(diskList, bl, tl, tr, br, output);
		}

		else if (skew <= 10.0/3.0) //directly split into 2
		{
			vector<double> new1, new2;
			double ratio = greedySplit(diskList, new1, new2);

			if (ratio < 0) return false;

			double splitY = bl.y + (tl.y - bl.y)*(1.0/(1.0 + 1.0/ratio));

			//recurse on the 2 good rectangles
			Point midL(bl.x, splitY), midR(, splitY);
			vector<pair<double, Point> > answer1, answer2;
			bool possible = good_rectangle_cover(new1, bl, midL, midR, br, answer1);
			possible = possible && good_rectangle_cover(new2, midL, tl, tr, midR, answer2);

			if (!possible) return false;

			output.resize(0);
			for (int i = 0 ; i < answer1.size() ; i ++)
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
			//create a new set and a good rectangle to go with it and recurse on the rest
			vector<double> newDiskSet;
			vector<double> remainingDiskSet;
			vector<pair<double, Point> > answer1, answer2;
			double areaYet = 0.0;
			int index = 0;
			while (areaYet < width*11.0/36.0)
			{
				if (index == diskList.size()) return false;//not enough disks
				newDiskSet.push_back(diskList[index]);
				areaYet += diskList[index]*diskList[index];
				index ++;
			}
			while (index < diskList.size())
			{
				remainingDiskSet.push_back(diskList[index]);
				index ++;
			}

			//compute the rectangle 
			double newHeight = (areaYet*12.0/11.0)/width;
			Point midL(bl.x, bl.y + newHeight), midR(tr.x, br.y + newHeight);

			//make the recursive calls
			bool possible = good_rectangle_cover(newDiskSet, bl, midL, midR, br, answer1);
			possible = possible && bad_rectangle_cover_bounded_radii(remainingDiskSet, midL, tl,
				tr, midR, answer2);

			if (!possible) return false;

			output.resize(0);
			for (int i = 0 ; i < answer1.size() ; i ++)
			{
				output.push_back(answer1[i]);
			}
			for (int i = 0 ; i < answer2.size() ; i ++)
			{
				output.push_back(answer2[i]);
			}
		}

	}


	return true;
}
