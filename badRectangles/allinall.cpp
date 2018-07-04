#define rep(i,a,b) for(int i=a;i<b;++i)
#define repr(i,a,b) for(int i=a,i > b;--i)
#define mm(lamb, tttt) memset(lamb, tttt, sizeof lamb)

#define null NULL
#define eps 0.000000001
#define mod 1000000007
#define PI 3.14159265358979323846
#define pb push_back
#define pf push_front
#define mp make_pair
#define fi first
#define se second
#define ALL(V) V.begin(), V.end()
#define sz(V) (ll)V.size()
#define _ <<" "<<

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stack>
#include <queue>
#include <deque>
#include <vector>
#include <iterator>
#include <bitset>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <limits.h>
#include <iomanip>
#include <cctype>
#include <numeric>
#include <complex>

using namespace std;

typedef long long ll;
typedef vector <int> vi;
typedef pair <double, double> ii;
typedef pair<int, pair<int,int> > iii;
typedef vector<ii> vii;


struct Point{
	double x,y;
	Point():x(-1),y(-1){}
	Point(double _x,double _y):x(_x),y(_y){}
};

bool deter(double r_1,double r_2,double x,double y){
	if(x >= y) swap(x,y);
	if(r_1 >= 0.5*sqrt(x*x+y*y)) return true;
	if(r_2 < 0.5*x+eps) return false;
	if(2*sqrt(r_1*r_1-0.25*x*x)+2*sqrt(r_2*r_2-0.25*x*x) >= y+eps) return true;
	return false;
}

// return the coordinates of the disk,assume the disk is sorted 
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

bool good_rectangle_cover(vector<double> & radii,Point lowleft,Point upleft,Point upright,Point lowright,vector<pair<double,Point> > & ans){
	double x = upleft.y - lowleft.y,y=lowright.x-lowleft.x;
	//cout<<"In rec lowleft "<<lowleft.x<<" "<<lowleft.y<<" (ver,hor)= "<<x<<" "<<y<<" number of disks "<<radii.size()<<endl;
	if(radii.empty()) return false;
	// basecases
	
	if(radii.size() == 1){
	//	cout<<"Base case 1 : side lengths are "<<x<<" "<<y<<" radii is "<<radii[0]<<" lowleft is "<<lowleft.x<<" "<<lowleft.y<<endl;
		double g_0 = 0.5*sqrt(x*x+y*y);
		if(radii[0] < g_0) {
			ans.pb(mp(radii[0],Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return false;
		}
		else {// #2.1
	//		cout<<"2.1 Placing circle of r="<<radii[0]<<" at x= "<<0.5*(lowright.x+lowleft.x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
			ans.pb(mp(radii[0],Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return true;
		}
	}
	
	double r_1 = radii[0],r_2=radii[1];
	
	if(radii.size() == 2){
	//	cout<<"Base case 2 : side lengths are "<<x<<" "<<y<<" radii are "<<r_1<<" "<<r_2<<" lowleft is "<<lowleft.x<<" "<<lowleft.y<<" "<<deter(r_1,r_2,x,y)<<endl;
		if(!deter(r_1,r_2,x,y)) return 0;
		double small = min(x,y),large = max(x,y);
		double g_0 = 0.5*sqrt(x*x+y*y);
		if(r_1 >= g_0) {//#3.1
			ans.pb(mp(radii[0],Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return true;
		} 
		else if(r_2 < 0.5*min(x,y)){//#3.2 // they cant possibly cover it
			if(x <= y) { 
				if(r_1 >= 0.5*x) ans.pb(mp(radii[0],Point(lowleft.x+sqrt(r_1*r_1-0.25*small*small),0.5*(upleft.y+lowleft.y))));
			}
			else if(r_1 >= 0.5*y) ans.pb(mp(radii[0],Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*small*small))));
			return false;
		}

		if (2.0*sqrt(r_1*r_1 - 0.25*(small*small))+2.0*sqrt(r_2*r_2 - 0.25*(small*small)) < large) {
			if(x <= y){
				ans.pb(mp(radii[0],Point(lowleft.x+sqrt(r_1*r_1-0.25*small*small),0.5*(upleft.y+lowleft.y))));
				ans.pb(mp(radii[1],Point(lowright.x-sqrt(r_2*r_2-0.25*small*small),0.5*(upleft.y+lowleft.y))));
			}
			else {
				ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-small*small*0.25))));
				ans.pb(mp(r_2,Point(0.5*(lowleft.x+lowright.x),upleft.y-sqrt(r_1*r_1-small*small*0.25))));
			}
			return false;
		}
		else{// there are two cases depending on the orientation
			if(x <= y){
				ans.pb(mp(radii[0],Point(lowleft.x+sqrt(r_1*r_1-0.25*small*small),0.5*(upleft.y+lowleft.y))));
				ans.pb(mp(radii[1],Point(lowright.x-sqrt(r_2*r_2-0.25*small*small),0.5*(upleft.y+lowleft.y))));
				return true;
			}
			else{
				ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-small*small*0.25))));
				ans.pb(mp(r_2,Point(0.5*(lowleft.x+lowright.x),upleft.y-sqrt(r_2*r_2-small*small*0.25))));
				return true;
			}
		}
	}
	//recursion
	
	
	if(x <= y){ // horizontal side is longer 
		// case wise
		double g_0 = 0.5*sqrt(x*x+y*y), g_1 = 0.5*(sqrt(x*x+(y-(x/6.0))*(y-(x/6.0)))), g_2 =  x*sqrt((484.0-11.0*(sqrt(1360)))/288.0);
		if(r_1 >= g_0) {
			Point center(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y));
			ans.pb(mp(r_1,center));
			return true;
		}
		else if(r_1 >= g_1){
			ans.pb(mp(r_1,Point(lowleft.x+sqrt(r_1*r_1-0.25*x*x),0.5*(upleft.y+lowleft.y))));
			double l_0 = 0.5*x*sqrt(1.0+1.0/36.0),l_1 = 0.5*x*sqrt(0.25+1.0/36.0),l_2 = 0.5*x*sqrt(1.0/16.0+1.0/36.0);
			if(r_2 >= l_0){
				ans.pb(mp(r_2,Point(lowleft.x+(y-x/12.0),lowright.y+x/2)));
				return true;
			}
			else if(r_2 >= l_1){
				// first deal with the boundary case
				if( y/x >= 2.999938 && (r_1 >= 1.5805*x &&  r_2 >= 0.49804*x)){//#4.3.1
					// scale everythig up and find epsilon to calculate z and then scale evryhing down 
					double scale = (1.0/x);
					r_1*=scale;r_2*=scale;x*=scale;y*=scale;
					double epsilon = 0.5*sqrt(10)-r_1;
					double z = (8*epsilon*sqrt(1.0+y*y));
					z/=scale;x/=scale;y/=scale;r_1/=scale;
					double len = y-2.0*sqrt(r_1*r_1-0.25*x*x);
					ans.pb(mp(r_2,Point(lowright.x - len/2.0,upright.y-sqrt(r_2*r_2-0.25*(len*len)))));
					radii.erase(radii.begin());radii.erase(radii.begin());
					return good_rectangle_cover(radii,Point(lowright.x-z/3.0,lowright.y),Point(lowright.x-z/3.0,lowright.y+z),Point(lowright.x,lowright.y+z),lowright,ans);	
				}
				else{ 
					double len = y-2.0*sqrt(r_1*r_1-0.25*x*x);
					double len1 = x - 2.0*sqrt(r_2*r_2-0.25*len*len);
					ans.pb(mp(r_2,Point(lowright.x - len/2.0,upright.y-sqrt(r_2*r_2-0.25*(len*len)))));
					radii.erase(radii.begin());radii.erase(radii.begin());
					if(len/len1 <= 3 && len1/len <=3 ) return good_rectangle_cover(radii,Point(lowright.x-len,lowright.y),Point(lowright.x-len,lowright.y+len1),Point(lowright.x,lowright.y+len1),lowright,ans);
					else if(len1 >= len) return good_rectangle_cover(radii,Point(lowright.x-len1/3.0,lowright.y),Point(lowright.x-len1/3.0,lowright.y+len1),Point(lowright.x,lowright.y+len1),lowright,ans);
					else return good_rectangle_cover(radii,Point(lowright.x-len,lowright.y),Point(lowright.x-len,lowright.y+len/3.0),Point(lowright.x,lowright.y+len/3.0),lowright,ans);
				}
			}
			else if(r_2 >= l_2){
				ans.pb(mp(r_2,Point(lowleft.x+(y-x/12.0),upleft.y-x/8.0)));
				radii.erase(radii.begin());
				radii.erase(radii.begin());
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				double ratio = greedySplit(radii,a1,a2);
				if(ratio > 2) cerr<<"ratio shouldnnt be bad than 2"<<endl;
			    bool f1 = good_rectangle_cover(a1,Point(lowright.x-x/6.0,lowright.y + 0.75*x - (ratio*0.75*x/(ratio+1.0))),Point(lowright.x-x/6.0,lowright.y + 0.75*x),Point(lowright.x,lowleft.y+0.75*x),Point(lowright.x,lowright.y + 0.75*x - (ratio*0.75*x/(ratio+1.0))),ans1);
				bool f2 = good_rectangle_cover(a2,Point(lowright.x-x/6.0,lowright.y),Point(lowright.x-x/6.0,lowright.y+(0.75*x/(ratio+1.0))),Point(lowright.x,lowright.y+(0.75*x/(ratio+1.0))),lowright,ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;	
			}
			else {
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				radii.erase(radii.begin());
				double ratio = greedySplit(radii,a1,a2);
				bool f1 = good_rectangle_cover(a1,Point(lowright.x-x/6.0,lowright.y),Point(lowright.x-x/6.0,lowright.y+x/2.0),Point(lowright.x,lowright.y+x/2.0),lowright,ans1);
				bool f2 = good_rectangle_cover(a2,Point(lowright.x-x/6.0,lowright.y+x/2.0),Point(lowright.x-x/6.0,upright.y),upright,Point(lowright.x,lowright.y+x/2.0),ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;
			}
		}
		else if(r_1 >= g_2){ 
	//		cout<<"8.8 place r= "<<r_1<<" x= "<<lowleft.x + sqrt(r_1*r_1-0.25*x*x)<<" y= "<<0.5*(lowleft.y+upleft.y)<<endl;
			ans.pb(mp(r_1,Point(lowleft.x + sqrt(r_1*r_1-0.25*x*x),0.5*(lowleft.y+upleft.y))));
			radii.erase(radii.begin());
			// two cases depending on whether the recatngle is good.
			double rem = y-2*sqrt(r_1*r_1 - 0.25*x*x);
			if(rem >= x/3.0){
				return good_rectangle_cover(radii,Point(lowright.x-rem,lowright.y),Point(upright.x-rem,upright.y),upright,lowright,ans);
			}
			else {
				return good_rectangle_cover(radii,Point(lowright.x - x/3.0,lowright.y),Point(lowright.x-x/3.0,upright.y),upright,lowright,ans);
			}
		}
		else{
			vector<double> a1,a2;
			vector<pair<double,Point> > ans1,ans2;
			double ratio = greedySplit(radii,a1,a2);
			if(ratio >= 3) cerr<<"Ratio shouldnt be bad than 3"<<endl;
			bool f1 = good_rectangle_cover(a2,lowleft,upleft,Point(upleft.x+(y)/(ratio+1.0),upleft.y),Point(lowleft.x+(y)/(ratio+1.0),lowleft.y),ans1);
			bool f2 = good_rectangle_cover(a1,Point(lowleft.x+(y)/(ratio+1.0),lowleft.y),Point(upleft.x+(y)/(ratio+1.0),upleft.y),upright,lowright,ans2);
			for(auto itr : ans1 ) ans.pb(itr);
			for(auto itr : ans2 ) ans.pb(itr);
			return f1 & f2;
		}
	}
	else{ // vertical side is longer 
		double g_0 = 0.5*sqrt(x*x+y*y), g_1 = 0.5*(sqrt(y*y+(x-(y/6.0))*(x-(y/6.0)))), g_2 =  y*sqrt((484.0-11.0*(sqrt(1360)))/288.0);
		if(r_1 >= g_0) {
			ans.pb(mp(r_1,Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return true;
		}
		else if(r_1 >= g_1){
			ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
			double l_0 = 0.5*y*sqrt(1.0+1.0/36.0),l_1 = 0.5*y*sqrt(0.25+1.0/36.0),l_2 = 0.5*y*sqrt(1.0/16.0+1.0/36.0);
			if(r_2 >= l_0) {
				ans.pb(mp(r_2,Point(0.5*(lowleft.x+lowright.x),upleft.y-y/12.0))); 
				return true;
			} 
			else if(r_2 >= l_1){
				if( x/y >= 2.999938 && (r_1 >= 1.5805*y &&  r_2 >= 0.49804*y)){
					// scale everythig up and fnd epsilon to calculate z and then scale evryhing down 
					double scale = (1.0/y);
					r_1*=scale;r_2*=scale;x*=scale;y*=scale;
					double epsilon = 0.5*sqrt(10)-r_1;
					double z = 2.0*x*x*(8*epsilon*sqrt(1.0+x*x));
					z/=scale;x/=scale;y/=scale;r_1/=scale;
					double len = x-2.0*sqrt(r_1*r_1-0.25*y*y);
					ans.pb(mp(r_2,Point(upleft.x+sqrt(r_2*r_2-0.25*len*len),upleft.y-0.5*len)));
					radii.erase(radii.begin());
					radii.erase(radii.begin());
					return good_rectangle_cover(radii,Point(upright.x-z,upright.y-z/3.0),Point(upright.x-z,upright.y),upright,Point(upright.x,upright.y-z/3.0),ans);	
				}
				else{
					double len = x-2.0*sqrt(r_1*r_1-0.25*y*y);
					double len1 = y - 2.0*sqrt(r_2*r_2-0.25*len*len);
					ans.pb(mp(r_2,Point(lowleft.x+sqrt(r_2*r_2-0.25*len*len),upleft.y-0.5*len)));
					radii.erase(radii.begin());radii.erase(radii.begin());
					if(len/len1 <= 3.0 && len1/len <=3.0) return good_rectangle_cover(radii,Point(upright.x-len1,upright.y-len),Point(upright.x-len1,upright.y),upright,Point(upright.x,upright.y-len),ans);
					else if(len1 > len) return good_rectangle_cover(radii,Point(upright.x-len1,upright.y-len1/3.0),Point(upright.x-len1,upright.y),upright,Point(upright.x,upright.y-len1/3.0),ans);
					else return good_rectangle_cover(radii,Point(upright.x-len/3.0,upright.y-len),Point(upright.x-len/3.0,upright.y),upright,Point(upright.x,upright.y-len),ans);
				}
			}
			else if(r_2 >= l_2){
				ans.pb(mp(r_2,Point(upleft.x+y/8.0,upleft.y-y/12.0)));
				radii.erase(radii.begin());
				radii.erase(radii.begin());
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				double ratio = greedySplit(radii,a1,a2);
				if(ratio >= 2) cerr<<"Ratio shouldnt be bad than 2"<<endl;
				bool f1 = good_rectangle_cover(a1,Point(upleft.x+y*0.25,upleft.y-y/6.0),Point(upleft.x+y*0.25,upleft.y),Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y),Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y-y/6.0),ans1);
				bool f2 = good_rectangle_cover(a1,Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y-y/6.0),Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y),upright,Point(upright.x,upright.y-y/6.0),ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;
			}
			else {
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				radii.erase(radii.begin());
				double ratio = greedySplit(radii,a1,a2);
				bool f1 = good_rectangle_cover(a1,Point(upleft.x,upleft.y-y/6.0),upleft,Point(upleft.x+y/2.0,upleft.y),Point(upleft.x+y/2.0,upleft.y-y/6.0),ans1);
				bool f2 = good_rectangle_cover(a2,Point(upleft.x+y/2.0,upleft.y-y/6.0),Point(upleft.x+y/2.0,upleft.y),upright,Point(upright.x,upright.y-y/6.0),ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;
			}
		}
		else if(r_1 >= g_2){
	//		cout<<"8.9 place r= "<<r_1<<" x= "<<0.5*(lowleft.x+lowright.x)<<" y= "<<lowleft.y+sqrt(r_1*r_1-0.25*y*y)<<" rem is "<<x-2.0*sqrt(r_1*r_1 - 0.25*y*y)<<endl;
			ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
			radii.erase(radii.begin());
			// two cases depending on whether the rectangle is good.
			double rem = x-2.0*sqrt(r_1*r_1 - 0.25*y*y);
			if(rem >= y/3.0) {
				return good_rectangle_cover(radii,Point(upleft.x,upleft.y-rem),upleft,upright,Point(upright.x,upright.y-rem),ans);
			}
			else {
				return good_rectangle_cover(radii,Point(upleft.x,upleft.y-y/3.0),upleft,upright,Point(upright.x,upright.y-y/3.0),ans);
			}
		}
		else {
			vector<double> a1,a2;
			vector<pair<double,Point> > ans1,ans2;
			double ratio = greedySplit(radii,a1,a2);
			if(ratio >= 3) cerr<<"Ratio shouldnt be bad than 3"<<endl;
			bool f1 = good_rectangle_cover(a1,lowleft,Point(lowleft.x,lowleft.y+ratio*x*(1.0/(ratio+1))),Point(lowright.x,lowright.y+ratio*x*(1.0/(ratio+1))),lowright,ans1);
			bool f2 = good_rectangle_cover(a2,Point(lowleft.x,lowleft.y+ratio*x*(1.0/(ratio+1))),upleft,upright,Point(lowright.x,lowright.y+ratio*x*(1.0/(ratio+1))),ans2);
			for(auto itr : ans1 ) ans.pb(itr);
			for(auto itr : ans2 ) ans.pb(itr);
			return f1 & f2;
		}
	}
	
	return false;
}

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
			Point midL(bl.x, splitY), midR(br.x, splitY);
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

bool bad_rectangle_cover(vector<double> & radii,Point lowleft,Point upleft,Point upright,Point lowright,vector<pair<double,Point> > & ans){
	//// basecases
	if(radii.empty()) return false;
	double y = lowright.x-lowleft.x,x=upright.y-lowright.y,r_1=radii[0];
	if(radii.size()==1) {
		if(r_1 >= 0.5*sqrt(x*x+y*y)) {ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),0.5*(upleft.y+lowleft.y)))); return true;}
		else return false;	
	}
	if(x <= y){
		double g_0 = x*(sqrt(77.0))/6.0;
		if(r_1 <= g_0) return bad_rectangle_cover_bounded_radii(radii,lowleft,upleft,upright,lowright, ans);
		ans.pb(mp(r_1,Point(lowleft.x+sqrt(r_1*r_1-0.25*x*x),0.5*(lowleft.y+upleft.y))));
		double rem = y - 2*sqrt(r_1*r_1-0.25*x*x);
		radii.erase(radii.begin());
		if(rem >= 3.0*x) return bad_rectangle_cover(radii,Point(lowright.x-rem,lowright.y),Point(upright.x-rem,upright.y),upright,lowright,ans);
		else if(rem >= x/3.0) return good_rectangle_cover(radii,Point(lowright.x-rem,lowright.y),Point(upright.x-rem,upright.y),upright,lowright,ans);
		else return bad_rectangle_cover(radii,Point(lowright.x-rem,lowright.y),Point(upright.x-rem,upright.y),upright,lowright,ans);
	}
	else {
		double g_0 = y*(sqrt(77))/6.0;
		if(r_1 <= g_0) return bad_rectangle_cover_bounded_radii(radii,lowleft,upleft,upright,lowright, ans);
		ans.pb(mp(r_1,Point(0.5*(lowright.x+lowleft.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
		double rem = x - 2*sqrt(r_1*r_1-0.25*y*y);
		radii.erase(radii.begin());
		if(rem >= 3.0*y) return bad_rectangle_cover(radii,Point(upleft.x,upleft.y-rem),upleft,upright,Point(upright.x,upright.y-rem),ans);
		else if(rem >= x/3.0) return good_rectangle_cover(radii,Point(upleft.x,upleft.y-rem),upleft,upright,Point(upright.x,upright.y-rem),ans);
		else return bad_rectangle_cover(radii,Point(upleft.x,upleft.y-rem),upleft,upright,Point(upright.x,upright.y-rem),ans);
	}
} 

bool checker(double x,double y,pair<double,Point> itr){
	double x1 = itr.se.x,y1=itr.se.y,r1=itr.fi;
	if(((x-x1)*(x-x1) +(y-y1)*(y-y1)) <= r1*r1 + eps) return true;
	else return false; 
}

int main(){
	cout<<"enter lambda"<<endl;
	double lambda;cin>>lambda;
	cout<<"Enter number of distinct disks"<<endl;
	int n;cin>>n;
	std::vector<double> C;
	cout<<"Enter theri frequency followed by radius"<<endl;
	double Area=0;
	rep(i,0,n){
		int f;double x;cin>>f>>x; rep(i,0,f) C.pb(x);Area+=f*PI*x*x;
	}
	cout<<Area<<" "<<PI*(lambda*lambda+2)/4.0<<endl;
	vector<pair<double,Point> > ans;
	sort(C.begin(),C.end());
	reverse(C.begin(),C.end());
	bool f = bad_rectangle_cover(C,Point(0,0),Point(0,1),Point(lambda,1),Point(lambda,0),ans);
	ofstream myfile;
	myfile.open("input");
	cout<<f<<endl;
	for(auto itr : ans){
		myfile<<itr.fi<<" ("<<itr.se.x<<","<<itr.se.y<<")"<<endl;
		//cout<<itr.fi<<" "<<itr.se.x<<" "<<itr.se.y<<endl;
	}
	myfile.close();
	
	if(f){
		for(double x = 0.0 ; x <=lambda ;x+=0.001){
			for(double y = 0.0 ; y <=1.0 ;y+=0.0001){
				bool f=0;
				for(auto itr : ans){
					if(checker(x,y,itr)) {f=1;break;}
				}
				if(!f) {cout<<"Found a point not covered (x,y)== "<<x<<" "<<y<<endl;
				double mindist = 10;
				for(auto itr : ans){
					double x1 = itr.se.x,y1=itr.se.y,r1=itr.fi;
					mindist = min(mindist,((x-x1)*(x-x1) +(y-y1)*(y-y1)) - r1*r1 );
				}
				cout<<mindist<<endl;
				return 0;}	
			}
		}
	}
}
