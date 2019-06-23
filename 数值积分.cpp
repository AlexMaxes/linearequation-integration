#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <iomanip>
#include "Timer.h"

using namespace std;

double f1(double x) //f(x) = x * e^x
{
	return x * exp(x);
}

double f2(double x) //f(x) = e^x
{
	return exp(x);
}

double f3(double x) //f(x) = e^x * sqrt(1 - x^2)
{
	return exp(x) * sqrt(1 - x * x);
}

double f4(double x) //f(x) = sin(x)  
{
	return sin(x);
}

double f5(double x) //f(x) = sin(PI/4 * x + PI/4)  
{
	return sin(M_PI / 4 * x + M_PI / 4);
}

double gauss_ch1(double (*f)(double), int n)
{
	double res = 0;
	for (int k = 1; k <= n; k++)
		res += f(cos((2 * k - 1) * M_PI / (2 * n)));
	res = res * M_PI / n;
	return res;
}

double gauss_ch2(double (*f)(double), int n)
{
	double res=0;
	for (int k = 1; k <= n; k++)
		res += f(cos(k * M_PI / (n + 1))) * sin(k * M_PI / (n + 1)) * sin(k * M_PI / (n + 1));
	res = res * M_PI / (n + 1);
	return res;
}

double comp_trep(double (*f)(double), double a, double b)
{
	double T0 = (f(a) + f(b)) * (b - a) / 2, Tm;
	int m = 1, count = 0;
	while(true)
	{
		double Hm = (b - a) / pow(2, m), sum = 0;
		for(int i = 1; i <= pow(2, m - 1); i++)
		{
			sum += f(Hm * (2 * i - 1) + a);
			count++;
		}
		sum *= Hm;
		Tm = T0 / 2 + sum;
		if(fabs(Tm - T0) > 0.0000000001)
			T0 = Tm;
		else
			break;
		m++;
	}
	cout<<"f(x)运行次数："<<count<<endl<<"迭代次数："<<m<<endl; 
	return Tm;
}

double romberg(double (*f)(double), double a, double b)
{
	double T[4][10000];
	T[0][0] = (f(a) + f(b)) * (b - a) / 2;
	int m = 1, count = 0;
	while(true)
	{
		double sum = 0, Hm = (b - a) / pow(2, m);
		for( int i = 1; i <= pow(2, m - 1); i++)
		{
			sum += f(Hm * (2 * i - 1) + a);
			count++;
		}
		T[0][m] = T[0][m - 1] / 2 + sum * Hm;
		int min;
		if(m > 3)
			min = 3;
		else
			min = m;
		for(int j = 1; j <= min; j++)
			T[j][m - j] = (pow(4, j) * T[j - 1][m - j + 1] - T[j - 1][m - j]) / (pow(4, j) - 1);
		if(fabs((T[3][m - 3] - T[3][m - 4]) / T[3][m - 3]) < 0.0000000001)
			break;
		m++;
	}
	cout<<"f(x)运行次数："<<count<<endl<<"迭代次数："<<m<<endl; 
	return T[3][m - 3];
}

double gauss_leg_9(double (*f)(double)) 
{
	double Xi[9] = {0.9681602395, -0.9681602395, 0.8360311073, -0.8360311073, 0.6133714327, -0.6133714327, 0.3242534234, -0.3242534234, 0};
	double Ai[9] = {0.0812743884, 0.0812743884, 0.1806481607, 0.1806481607, 0.2606106964, 0.2606106964, 0.3123470770, 0.3123470770, 0.3302393550};
	double sum = 0;
	for (int i = 0; i < 9; i++)
	{	
		sum += Ai[i] * f(Xi[i]);
		//cout<<Ai[i]<<" "<<f(Xi[i])<<" "<<sum<<endl;
	}		
	return sum;
}


int main()
{
	Timer time;
	int n = 20;
	double res;
	cout<<"实验一："<<endl<<endl; 
	
	time.reset();
	time.start();
	res = gauss_ch1(f1,n);
	time.stop();
	cout<<"gauss_ch1  方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	time.reset();
	time.start();
	res = gauss_ch2(f2,n);
	time.stop();
	cout<<"gauss_ch2  方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	time.reset();
	time.start();
	res = comp_trep(f3, -1, 1);
	time.stop();
	cout<<"comp_trep  方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	time.reset();
	time.start();
	res = romberg(f3, -1, 1);
	time.stop();
	cout<<"romberg    方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	time.reset();
	time.start();
	res = gauss_leg_9(f3);
	time.stop();
	cout<<"gauss_leg_9方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	cout<<"实验二："<<endl<<endl;
	 
	time.reset();
	time.start();
	res = comp_trep(f4, 0, M_PI / 2);
	time.stop();
	cout<<"comp_trep  方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"误差为             ："<<setiosflags(ios::fixed)<<setprecision(20)<<fabs(1 - res)<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	time.reset();
	time.start();
	res = romberg(f4, 0, M_PI / 2);
	time.stop();
	cout<<"romberg    方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"误差为             ："<<setiosflags(ios::fixed)<<setprecision(20)<<fabs(1 - res)<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	
	time.reset();
	time.start();
	res = gauss_leg_9(f5) * M_PI / 4;
	time.stop();
	cout<<"gauss_leg_9方法结果："<<setiosflags(ios::fixed)<<setprecision(20)<<res<<endl;
	cout<<"误差为             ："<<setiosflags(ios::fixed)<<setprecision(20)<<fabs(1 - res)<<endl;
	cout<<"运行时间           ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	return 0;
}
