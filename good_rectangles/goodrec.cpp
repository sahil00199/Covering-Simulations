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
	Point(double _x,double _y):x(_x),y(_y){}
};



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

// mistake 1 :- find g_2 in terms of x

bool good_rectangle_cover(vector<double> & radii,Point lowleft,Point upleft,Point upright,Point lowright,vector<pair<double,Point> > & ans){
	if(radii.empty()) return false;
	// basecases
	double x = upleft.y - lowleft.y,y=lowright.x-lowleft.x;
	if(radii.size() == 1){
		cout<<"Base case 1 "<<x<<" "<<y<<" "<<radii[0]<<endl;
		double g_0 = 0.5*sqrt(x*x+y*y);
		if(radii[0] < g_0) {
			ans.pb(mp(radii[0],Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return false;
		}
		else {// #2.1
			cout<<"2.1 Placing circle of r="<<radii[0]<<" at x= "<<0.5*(lowright.x+lowleft.x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
			ans.pb(mp(radii[0],Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return true;
		}
	}
	
	double r_1 = radii[0],r_2=radii[1];
	
	if(radii.size() == 2){
		cout<<"Base case 1 "<<x<<" "<<y<<" "<<r_1<<" "<<r_2<<endl;
		double small = min(upleft.y - lowleft.y,lowright.x-lowleft.x),large = max(upleft.y - lowleft.y,lowright.x-lowleft.x);
		double g_0 = 0.5*sqrt(x*x+y*y);
		if(r_1 >= g_0) {//#3.1
			cout<<"3.1 Placing circle of r="<<radii[0]<<" at x= "<<0.5*(lowright.x+lowleft.x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<" "<<upleft.y<<" "<<lowleft.y<<endl;
			ans.pb(mp(radii[0],Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return true;
		} 
		else if(r_2 < 0.5){//#3.2
			cout<<"3.2 Placing circle of r="<<radii[0]<<" at x= "<<lowleft.x+sqrt(r_1*r_1-0.25*small*small)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
			ans.pb(mp(radii[0],Point(lowleft.x+sqrt(r_1*r_1-0.25*small*small),0.5*(upleft.y+lowleft.y))));
			return false;
		}

		if (2.0*sqrt(r_1*r_1 - 0.25*(small*small))+2.0*sqrt(r_2*r_2 - 0.25*(small*small)) < large) return false;
		else{// there are two cases depending on the orientation
			if(x <= y){//#3.3
				cout<<"3.3 Placing circle of r="<<radii[0]<<" at x= "<<lowleft.x+sqrt(r_1*r_1-0.25*small*small)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
				cout<<"3.3 Placing circle of r="<<radii[1]<<" at x= "<<lowright.x-sqrt(r_2*r_2-0.25*small*small)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
				ans.pb(mp(radii[0],Point(lowleft.x+sqrt(r_1*r_1-0.25*small*small),0.5*(upleft.y+lowleft.y))));
				ans.pb(mp(radii[1],Point(lowright.x-sqrt(r_2*r_2-0.25*small*small),0.5*(upleft.y+lowleft.y))));
				return true;
			}
			else{//#3.4
				cout<<"3.4 Placing circle of r="<<radii[0]<<" at x= "<<0.5*(lowleft.x+lowright.x)<<" y= "<<lowleft.y+2.0*sqrt(r_1*r_1-small*small*0.25)<<endl;
				cout<<"3.4 Placing circle of r="<<radii[1]<<" at x= "<<0.5*(lowleft.x+lowright.x)<<" y= "<<upleft.y-2.0*sqrt(r_1*r_1-small*small*0.25)<<endl;
				ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+2.0*sqrt(r_1*r_1-small*small*0.25))));
				ans.pb(mp(r_2,Point(0.5*(lowleft.x+lowright.x),upleft.y-2.0*sqrt(r_1*r_1-small*small*0.25))));
				return true;
			}
		}
	}
	//recursion
	
	
	if(x <= y){ // horizontal side is longer 
		// case wise
		double g_0 = 0.5*sqrt(x*x+y*y), g_1 = 0.5*(sqrt(x*x+(y-(1.0/6.0*x))*(y-(1.0/6.0*x)))), g_2 =  sqrt((484.0-11.0*(sqrt(1360)))/288.0);
		if(r_1 >= g_0) {//#4.1
			Point center(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y));
			cout<<"4.1 Placing circle of r="<<r_1<<" at x= "<<0.5*(lowright.x+lowleft.x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
			ans.pb(mp(r_1,center));
			return true;
		}
		else if(r_1 >= g_1){//#4.2
			cout<<"4.2 Placing circle of r="<<radii[0]<<" at x= "<<lowleft.x+sqrt(r_1*r_1-0.25*x*x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
			ans.pb(mp(radii[0],Point(lowleft.x+sqrt(r_1*r_1-0.25*x*x),0.5*(upleft.y+lowleft.y))));
			double l_0 = 0.5*x*sqrt(1.0+1.0/36.0),l_1 = 0.5*x*sqrt(0.25+1.0/36.0),l_2 = 0.5*x*sqrt(1.0/16.0+1.0/36.0);
			if(r_2 >= l_0){//#4.2.1
				cout<<"4.2.1 Placing circle of r="<<r_2<<" at x= "<<lowleft.x+(y-x/12.0)<<" y= "<<x/2<<endl;
				ans.pb(mp(r_2,Point(lowleft.x+(y-x/12.0),x/2)));
				return true;
			}
			else if(r_2 >= l_1){//#4.3
				// first deal with the boundar case
				if( y/x >= 2.999938 && (r_1 >= 1.5805*x &&  r_2 >= 0.49804*x)){//#4.3.1
					// scale everythig up and fnd epsilon to calculate z and then scale evryhing down 
					double scale = (1.0/x);
					r_1*=scale;r_2*=scale;x*=scale;y*=scale;
					double epsilon = 0.5*sqrt(10)-r_1;
					double z = 2.0*y*y*(8*epsilon*sqrt(1.0+y*y));
					z/=scale;x/=scale;y/=scale;r_1/=scale;
					cout<<"4.3 Placing circle of r="<<r_1<<" at x= "<<lowleft.x+sqrt(r_1*r_1-0.25*x*x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
					ans.pb(mp(r_1,Point(lowleft.x+sqrt(r_1*r_1-0.25*x*x),0.5*(upleft.y+lowleft.y))));
					double len = y-2.0*sqrt(r_1*r_1-0.25*x*x);
					cout<<"4.3 Placing circle of r="<<r_1<<" at x= "<<lowleft.x+sqrt(r_1*r_1-0.25*x*x)<<" y= "<<0.5*(upleft.y+lowleft.y)<<endl;
					ans.pb(mp(r_2,Point(lowright.x - len/2.0,upright.y-sqrt(r_2*r_2-0.25*(len*len)))));
					radii.erase(radii.begin());radii.erase(radii.begin());
					return good_rectangle_cover(radii,Point(lowright.x-z/3.0,lowright.y),Point(lowright.x-z/3.0,lowright.y+z),Point(lowright.x,lowright.y+z),lowright,ans);	
				}
				else{ //#4.3.2
					ans.pb(mp(r_1,Point(lowleft.x+sqrt(r_1*r_1-0.25*x*x),0.5*(upleft.y+lowleft.y))));
					double len = y-2.0*sqrt(r_1*r_1-0.25*x*x);
					double len1 = x - 2.0*sqrt(r_2*r_2-0.25*len*len);
					ans.pb(mp(r_2,Point(lowright.x - y/12.0,upright.y-sqrt(r_2*r_2-0.25*(len*len)))));
					radii.erase(radii.begin());radii.erase(radii.begin());
					return good_rectangle_cover(radii,Point(lowright.x-len,lowright.y),Point(lowright.x-len,lowright.y+len1),Point(lowright.x,lowright.y+len1),lowright,ans);
				}
			}
			else if(r_2 >= l_2){//#4.4
				ans.pb(mp(r_2,Point(lowleft.x+(y-x/12.0),upleft.y-x/8.0)));
				radii.erase(radii.begin());
				radii.erase(radii.begin());
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				double ratio = greedySplit(radii,a1,a2);
				if(ratio > 2) cerr<<"ratio shouldnnt be bad than 2"<<endl;
			    bool f1 = good_rectangle_cover(a1,Point(lowright.x-x/6.0,lowleft.y + 0.75*x - (ratio*0.75/(ratio+1))),Point(lowright.x-x/6.0,lowleft.y + 0.75*x),Point(lowright.x,lowleft.y+0.75*x),Point(lowright.x,lowleft.y + 0.75*x - (ratio*0.75/(ratio+1))),ans1);
				bool f2 = good_rectangle_cover(a2,Point(lowright.x-x/6.0,lowleft.y),Point(lowright.x-x/6.0,lowleft.y+(0.75/(ratio+1))),Point(lowright.x,lowleft.y+(0.75/(ratio+1))),Point(lowright.x,lowleft.y),ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;	
			}
			else {//#4.5
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				double ratio = greedySplit(radii,a1,a2);
				bool f1 = good_rectangle_cover(a1,Point(lowleft.x+y-x/6.0,lowleft.y),Point(lowright.x-x/6.0,lowleft.y+x/2.0),Point(lowright.x,lowright.y+x/2.0),lowright,ans1);
				bool f2 = good_rectangle_cover(a2,Point(lowright.x-x/6.0,lowleft.y+x/2.0),Point(lowright.x-x/6.0,upleft.y),upright,Point(lowright.x,lowright.y+x/2.0),ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;
			}
		}
		else if(r_1 >= g_2){
			ans.pb(mp(r_1,Point(lowleft.x + sqrt(r_1*r_1-0.25*x*x),0.5*(lowleft.y+upleft.y))));
			radii.erase(radii.begin());
			// two cases depending on whether the recatangle is good.
			double rem = y-2*sqrt(r_1*r_1 - 0.25*x*x);
			if(rem >= x/3.0){
				return good_rectangle_cover(radii,Point(lowright.x-rem,lowleft.y),Point(upright.x-rem,upright.y),Point(upright.x,upright.y),Point(lowright.x,lowright.y),ans);
			}
			else {
				return good_rectangle_cover(radii,Point(lowright.x - x/3.0,lowright.y),Point(lowright.x-x/3.0,upright.y),Point(upright.x,upright.y),Point(lowright.x,lowright.y),ans);
			}
		}
		else{
			vector<double> a1,a2;
			vector<pair<double,Point> > ans1,ans2;
			double ratio = greedySplit(radii,a1,a2);
			if(ratio >= 3) cerr<<"Ratio shouldnt be bad than 3"<<endl;
			bool f1 = good_rectangle_cover(a2,Point(lowleft.x,lowleft.y),Point(upleft.x,upleft.y),Point(upleft.x+(y)/(ratio+1.0),upleft.y),Point(lowleft.x+(y)/(ratio+1.0),lowleft.y),ans1);
			bool f2 = good_rectangle_cover(a1,Point(lowleft.x+(y)/(ratio+1.0),lowleft.y),Point(upleft.x+(y)/(ratio+1.0),upleft.y),Point(upright.x,upright.y),Point(lowright.x,lowright.y),ans2);
			for(auto itr : ans1 ) ans.pb(itr);
			for(auto itr : ans2 ) ans.pb(itr);
			return f1 & f2;
		}
	}
	else{ // vertical side is longer 
		double g_0 = 0.5*sqrt(x*x+y*y), g_1 = 0.5*(sqrt(y*y+(x-(1.0/6.0*y))*(x-(1.0/6.0*y)))), g_2 =  sqrt((484.0-11.0*(sqrt(1360)))/288.0);
		if(r_1 >= g_0) {
			ans.pb(mp(r_1,Point(0.5*(lowright.x+lowleft.x),0.5*(upleft.y+lowleft.y))));
			return true;
		}
		else if(r_1 >= g_1){
			ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
			double l_0 = 0.5*y*sqrt(1.0+1.0/36.0),l_1 = 0.5*y*sqrt(0.25+1.0/36.0),l_2 = 0.5*y*sqrt(1.0/16.0+1.0/36.0);
			if(r_2 >= l_0) {
				ans.pb(mp(r_2,Point(0.5*(lowleft.x+lowright.x),upleft.y-x/12.0)));
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
					ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
					double len = x-2.0*sqrt(r_1*r_1-0.25*y*y);
					ans.pb(mp(r_2,Point(lowleft.x+sqrt(r_2*r_2-0.25*len*len),upleft.y-0.5*len)));
					radii.erase(radii.begin());radii.erase(radii.begin());
					return good_rectangle_cover(radii,Point(upright.x-z,upright.y-z/3.0),Point(upright.x-z,upright.y),upright,Point(upright.x,upright.y-z/3.0),ans);	
				}
				else{
					ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
					double len = x-2.0*sqrt(r_1*r_1-0.25*y*y);
					double len1 = y - 2.0*sqrt(r_2*r_2-0.25*len*len);
					ans.pb(mp(r_2,Point(lowleft.x+sqrt(r_2*r_2-0.25*len*len),upleft.y-0.5*len)));
					radii.erase(radii.begin());radii.erase(radii.begin());
					return good_rectangle_cover(radii,Point(upright.x-len1,upright.y-len),Point(upright.x-len1,upright.y),upright,Point(upright.x,upright.y-len),ans);
				}
			}
			else if(r_2 >= l_2){
				ans.pb(mp(r_2,Point(lowleft.x+y/8.0,upleft.y-y/12.0)));
				radii.erase(radii.begin());
				radii.erase(radii.begin());
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				double ratio = greedySplit(radii,a1,a2);
				if(ratio >= 2) cerr<<"Ratio shouldnt be bad than 2"<<endl;
				bool f1 = good_rectangle_cover(a1,Point(lowleft.x+y*0.25,upleft.y-y/6.0),Point(upleft.x+y*0.25,upleft.y),Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y),Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y-y/6.0),ans1);
				bool f2 = good_rectangle_cover(a1,Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y-y/6.0),Point(upleft.x+ratio*(y*0.75/(1.0+ratio)),upleft.y),Point(upright.x,upright.y),Point(upright.x,upright.y-y/6.0),ans1);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;
			}
			else {
				vector<double> a1,a2;
				vector<pair<double,Point> > ans1,ans2;
				double ratio = greedySplit(radii,a1,a2);
				bool f1 = good_rectangle_cover(a1,Point(lowleft.x,upleft.y-y/6.0),Point(upleft.x,upleft.y),Point(upleft.x+y/2.0,upleft.y),Point(upleft.x+y/2.0,upleft.y-y/6.0),ans1);
				bool f2 = good_rectangle_cover(a2,Point(upleft.x+y/2.0,upleft.y-y/6.0),Point(upleft.x+y/2.0,upleft.y),Point(upright.x,upright.y),Point(upright.x,upright.y-y/6.0),ans2);
				for(auto itr : ans1 ) ans.pb(itr);
				for(auto itr : ans2 ) ans.pb(itr);
				return f1 & f2;
			}
		}
		else if(r_1 >= g_2){
			ans.pb(mp(r_1,Point(0.5*(lowleft.x+lowright.x),lowleft.y+sqrt(r_1*r_1-0.25*y*y))));
			radii.erase(radii.begin());
			// two cases depending on whether the recatangle is good.
			double rem = y-2*sqrt(r_1*r_1 - 0.25*y*y);
			if(rem >= y/3.0) {
				return good_rectangle_cover(radii,Point(lowleft.x,lowleft.y+2*sqrt(r_1*r_1 - 0.25*y*y)),Point(upleft.x,upleft.y),Point(upright.x,upright.y),Point(upright.x,lowleft.y+2*sqrt(r_1*r_1 - 0.25*y*y)),ans);
			}
			else {
				return good_rectangle_cover(radii,Point(lowleft.x,upleft.y-y/3.0),Point(upleft.x,upleft.y),Point(upright.x,upright.y),Point(upright.x,upright.y-y/3.0),ans);
			}
		}
		else {
			vector<double> a1,a2;
			vector<pair<double,Point> > ans1,ans2;
			double ratio = greedySplit(radii,a1,a2);
			if(ratio >= 3) cerr<<"Ratio shouldnt be bad than 3"<<endl;
			bool f1 = good_rectangle_cover(a1,Point(lowleft.x,lowleft.y),Point(lowleft.x,lowleft.y+ratio*x*(1.0/(ratio+1))),Point(lowright.x,lowright.y+ratio*x*(1.0/(ratio+1))),Point(lowright.x,lowright.y),ans1);
			bool f2 = good_rectangle_cover(a2,Point(lowleft.x,lowleft.y+ratio*x*(1.0/(ratio+1))),upleft,upright,Point(lowright.x,lowright.y+ratio*x*(1.0/(ratio+1))),ans2);
			for(auto itr : ans1 ) ans.pb(itr);
			for(auto itr : ans2 ) ans.pb(itr);
			return f1 & f2;
		}
	}
	
	return false;
}

bool wrapper_rc(vector<double> & radii,vector<pair<double,Point> > & ans,double lambda){
	return good_rectangle_cover(radii,Point(0,0),Point(0,1),Point(lambda,1),Point(lambda,0),ans);
}

int main(){
	vector<double>  C;
	vector<pair<double,Point> > ans;
	cout<<"enter number of disks "<<endl;
	int n;cin>>n;
	cout<<"enter there radii"<<endl;
	double area=0.0;

	rep(i,0,n){
		double x;cin>>x;C.pb(x);area+=PI*(x*x);
	}

	cout<<area<<" M "<<(11.0/12.0)*PI*3.0<<endl;
	if(area < (11.0/12.0)*PI*3.0) cout<<"WARNING"<<endl;
	
	sort(C.begin(),C.end());
	reverse(C.begin(),C.end());
	cout<<wrapper_rc(C,ans,3.0)<<endl;
	ofstream myfile;
	myfile.open("input");
	for(auto itr : ans){
		myfile<<itr.fi<<" ("<<itr.se.x<<","<<itr.se.y<<")"<<endl;
	}
	myfile.close();
}

