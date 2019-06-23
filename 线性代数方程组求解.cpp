#include <iostream>
#include <math.h>
#include <iomanip>
#include "Timer.h"

using namespace std;

bool lu(double* a, int* pivot, int n){
	double main_elem, t;      
	for (int i = 0; i < n - 1; i++){       
		main_elem = fabs(a[n*i + i]);         
		pivot[i]=i; 
		for (int j = i + 1; j < n; j++){        
			if (fabs(a[n*j  + i]) > main_elem) {            
				main_elem = fabs(a[n*j + i]);
				pivot[i] = j;
			}
		}
		if (!fabs(a[n*n - 1])) 
			return true;
		if(pivot[i] != i) {         
			for (int j = i; j < n; j++){           
				t = a[n * i  + j];
				a[n * i  + j] = a[n * pivot[i] + j];                 
				a[n * pivot[i] + j] = t;
			}
		}
		//U 
		for (int j = i + 1; j < n; j++)             
			a[n*j + i] = a[n * j + i] / a[n * i + i]; 
		for (int j = i + 1; j < n; j++)          
			for(int k = i + 1; k < n; k++)  
				a[n * j  + k] = a[n * j + k] - a[n * j + i] * a[n * i + k];     
	}
	//L 
	for (int i = 0; i < n - 2; i++)    
		for(int k = i + 1; k < n-1; k++) {          
			t = a[n * pivot[k] + i];
			a[n * pivot[k] + i] = a[k * n  + i];             
			a[k * n  + i] = t;         
		}
	if (!fabs(a[n * n - 1])) 
		return true;
	else
		return false;
}

bool guass(double const* lu, int const* p, double* b, int n) {
	double t;
	for (int i = 0; i < n; i++)
		if (!fabs(lu[i * n + i])) 
			return true;
	//计算Qb 
	for (int i = 0; i < n - 1; i++) {      
		t = b[p[i]];         
		b[p[i]] = b[i];         
		b[i] = t;
	}
	//Ly=b 
	for(int i = 0 ; i < n; i++)          
		for(int j = 0; j < i; j++) 
			b[i] = b[i] - lu[n * i + j] * b[j];       
	//Uy=x
	for (int i = n - 1; i >= 0; i--)  {    
		for (int j = n - 1; j > i; j--)
			b[i] = b[i] - lu[n * i + j] * b[j];          
		b[i] = b[i] / lu[n * i + i];
	}
	return false;
}

void qr(double* a, double* d, int n) {  
	double T, main_elem, t[n];
	for (int i = 0; i < n - 1; i++){    
		main_elem = T = 0;
		for (int j = i; j < n; j++)
			main_elem  += a[n*j  + i] * a[n*j  + i ];           
		a[n*i  + i ] > 0 ? main_elem = -sqrt(main_elem) : main_elem = sqrt(main_elem); 
		d[i] = main_elem;
		a[n * i  + i] -= main_elem;         
		for (int j = i; j <= n - 1; j++)
			T += a[n*j  + i] * a[n*j  + i];         
		T = sqrt(T);
		for (int j = i; j <= n - 1; j++)
			a[n*j + i] /= T;
		for (int k = i + 1; k < n; k++) {         
			for (int j = i; j < n; j++)  {           
				T = 0;
				for (int l = i; l < n; l++)
					T += a[n*j  + i] * a[n*l  + i] * a[n*l  + k];                 
				t[j] = a[j*n + k] - 2 * T;
			}
			for (int j = i; j < n; j++)
				a[j * n + k] = t[j];
		}
	}
	d[n - 1] = a[(n - 1) * n + n - 1];
}

bool hshld(double const *qr, double const *d, double *b, int n) {
	double T, t[n];  
	for (int i = 0; i < n; i++)
		if (!fabs(d[i])) 
			return true;
	//b = Qtb 
	for (int i = 0; i < n - 1; i++) {    
		for (int j = i; j < n; j++) {        
			T  = 0;
			for (int k = i; k < n; k++)
				T  += qr[k * n + i] * qr[j * n + i] * b[k];             
			t[j] = b[j] - 2 * T;
		}
		for (int j = i; j < n; j++)             
			b[j] = t[j];
	}
	//x = R\b
	for (int j = n - 1; j > -1; j--)  {   
		for (int k = n - 1; k != j; k--)
			b[j] -= b[k] * qr[j*n + k];         
		b[j] /= d[j];
	}
	return false;
}

double* build_h(int n) //Hilbert矩阵 
{
	double* res = new double[n * n];
	for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				res[n * i + j] = 1.0 / (i + j + 1);
	return res;
}

double* build_nn(int n) //n*n矩阵 
{
	double* res = new double[n * n];
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			if(i == j || j == n - 1)
				res[n * i + j] = 1.0;
			else if(j > i)
				res[n * i + j] = 0.0;
			else
				res[n * i + j] = -1.0;
	return res;
}

double* build_bh(double * h,int n) 
{
	double* res = new double[n]; 
	for(int i = 0; i < n; i++)
	{
		double sum = 0.0;
		for(int j = 0; j < n; j++)
			sum += h[n * i + j];
		res[i] = sum;
	}
	return res;
}

double* build_bn(int n)
{
	double* res = new double[n];
	for(int i =0; i < n; i++)
		res[i] = 2.0 - i;
	return res;
}

double* build_bn2(int n)
{
	double* res = new double[n];
	for(int i =0; i < n; i++)
		res[i] = 2.0 + i;
	return res;
}

void call_guass(int n, int id)
{
	Timer time;
	time.reset();
	time.start();
	cout<<"当n="<<n<<"时："<<endl<<endl;
	double *a,*b,*A,*B;
	if(id == 1)
	{
		a = build_h(n);
		A = build_h(n);
		b = build_bh(a,n);
		B = build_bh(a,n);
	}
	else if(id == 2) 
	{
		a = build_nn(n);
		A = build_nn(n);
		b = build_bn(n);
		B = build_bn(n);
	}
	else
	{
		a = build_nn(n);
		A = build_nn(n);
		b = build_bn2(n);
		B = build_bn2(n);
	}
	
	int pivot[n];
	if(lu(a,pivot,n))
		//cout<<"不能分解";
		return; 
	if(guass(a,pivot,b,n))
		return;
	for(int i = 0; i < n; i++)
	{
		cout<<setiosflags(ios::fixed)<<setprecision(6)<<b[i]<<" ";
		if((i + 1 ) % 5 == 0)
			cout<<endl; 
	}
	cout<<endl;
	time.stop();
	if(id != 3)
	{
		double sum = 0;
		for(int i =0; i < n; i++)
			sum += fabs(b[i] - 1);
		cout<<"绝对误差为："<<sum / n<<endl; 
		cout<<"相对误差为："<<(sum / n) / 1.0 * 100<<"%"<<endl;
	}
	else
	{
		double sum_r = 0,sum_x = 0;
		for(int i = 0; i < n; i++)
		{
			double r = 0;
			for(int j = 0; j < n; j++)
				r += A[i * n + j] * b[j];
			r -= B[i];
			sum_r += r * r;
			sum_x += b[i] * b[i];
		}
		cout<<"残量为："<<sqrt(sum_r)<<endl;
		cout<<"相对残量为："<<sqrt(sum_r) / sqrt(sum_x)<<endl;
	}	
	cout<<"运行时间 ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	cout<<"--------------------------------------"<<endl;
}

void call_hd(int n, int id)
{
	Timer time;
	time.reset();
	time.start();
	cout<<"当n="<<n<<"时："<<endl<<endl;
	double *a,*b,*A,*B;
	if(id == 1)
	{
		a = build_h(n);
		A = build_h(n);
		b = build_bh(a,n);
		B = build_bh(a,n);
	}
	else if(id == 2)
	{
		a = build_nn(n);
		A = build_nn(n);
		b = build_bn(n);
		B = build_bn(n);
	}
	else
	{
		a = build_nn(n);
		A = build_nn(n);
		b = build_bn2(n);
		B = build_bn2(n);
	}
	double d[n];
	qr(a,d,n);
	if(hshld(a,d,b,n))
		return;
	for(int i = 0; i < n; i++)
	{
		cout<<setiosflags(ios::fixed)<<setprecision(6)<<b[i]<<" ";
		if((i + 1 ) % 5 == 0)
			cout<<endl; 
	}
	cout<<endl;
	time.stop();
	if(id != 3)
	{
		double sum = 0;
		for(int i =0; i < n; i++)
			sum += fabs(b[i] - 1);
		cout<<"绝对误差为："<<sum / n<<endl; 
		cout<<"相对误差为："<<(sum / n) / 1.0 * 100<<"%"<<endl;
	}
	else
	{
		double sum_r = 0,sum_x = 0;
		for(int i = 0; i < n; i++)
		{
			double r = 0;
			for(int j = 0; j < n; j++)
				r += A[i * n + j] * b[j];
			r -= B[i];
			sum_r += r * r;
			sum_x += b[i] * b[i];
		}
		cout<<"残量为："<<sqrt(sum_r)<<endl;
		cout<<"相对残量为："<<sqrt(sum_r) / sqrt(sum_x)<<endl;
	}
	cout<<"运行时间 ："<<setiosflags(ios::fixed)<<setprecision(20)<<time.seconds()<<endl<<endl;
	cout<<"--------------------------------------"<<endl;
}

int main()
{
	int n;
	cout<<"实验一："<<endl<<endl; 
	cout<<"列选主元Gauss消去法：" <<endl<<endl; 
	call_guass(5,1);
	call_guass(10,1);
	call_guass(15,1);

	cout<<"Householder变换法：" <<endl<<endl; 
	call_hd(5,1);
	call_hd(10,1);
	call_hd(15,1);
	
	cout<<"实验二:"<<endl<<endl;
	cout<<"列选主元Gauss消去法：" <<endl<<endl; 
	call_guass(10,2);
	call_guass(30,2);
	call_guass(60,2);
	
	cout<<"Householder变换法：" <<endl<<endl; 
	call_hd(10,2);
	call_hd(30,2);
	call_hd(60,2);
	
	cout<<"实验三:"<<endl<<endl;
	cout<<"列选主元Gauss消去法：" <<endl<<endl; 
	call_guass(10,3);
	call_guass(30,3);
	call_guass(60,3);
	
	cout<<"Householder变换法：" <<endl<<endl; 
	call_hd(10,3);
	call_hd(30,3);
	call_hd(60,3);
	return 0;
}
