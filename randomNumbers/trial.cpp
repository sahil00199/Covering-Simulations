#include <iostream>
#include <random>

using namespace std;

int main()
{
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,sqrt(77.0)/6.0);
	for (int i = 0 ; i < 10 ; i ++)
		cout<<distribution(generator)<<endl;


	return 0;
}