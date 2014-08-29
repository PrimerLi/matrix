#include <iostream>
#include <cstdlib>
#include <cmath>
#define PI acos(-1.0)
using namespace std;
class complex
{
private:
  double real;
  double imag;
public:
  complex();
  complex(double, double);
  complex(double);
  complex(const complex &);
  complex operator+ (const complex &) const;
  complex operator- (const complex &) const;
  complex operator* (const complex &) const;
  complex operator* (double) const;
  complex operator/ (const complex &) const;
  complex operator/ (double) const;
  const complex & operator= (const complex &); 
  const complex & operator= (double *);
  const complex & operator= (double);
  const complex & operator= (int);
  complex scale(double) const;
  void print() const;
  complex conjugate() const;
  double norm() const;
  double arg() const;
  double Real();
  double Imag();
  friend ostream & operator<< (ostream &os, const complex &);
  friend complex operator* (double, const complex &);
  friend complex operator/ (double, const complex &);
  friend double Re(const complex &);
  friend double Im(const complex &);
  friend double norm(const complex &);
  friend complex operator+ (double, const complex &);
  friend complex operator- (double, const complex &);
  friend complex exp(const complex & z);
};
complex::complex()
{
  real=0;
  imag=0;
}
complex::complex(double a, double b)
{
  real=a;
  imag=b;
}
complex::complex(double a)
{
    real=a;
    imag=0;
}
complex::complex(const complex &argument)
{
    real=argument.real;
    imag=argument.imag;
}
complex complex::operator+ (const complex & parameter) const
{
  return complex(real+parameter.real, imag+parameter.imag);
}
complex complex::operator- (const complex & parameter) const
{
  return complex(real-parameter.real, imag-parameter.imag);
}
complex complex::operator* (const complex & parameter) const
{
  complex temp;
  temp.real=real*parameter.real-imag*parameter.imag;
  temp.imag=real*parameter.imag+imag*parameter.real;
  return temp;
}
complex complex::operator* (double factor) const
{
  complex temp;
  temp.real=real*factor;
  temp.imag=imag*factor;
  return temp;
}
complex complex::operator/ (const complex & parameter) const
{
  if(parameter.real==0 && parameter.imag==0)
  {
    cout<<"You cannot divide by zero. "<<endl;
	exit(0);
  }
  else 
  {
    complex temp;
	temp.real=(real*parameter.real+imag*parameter.imag)/(parameter.real*parameter.real+parameter.imag*parameter.imag);
	temp.imag=(imag*parameter.real-real*parameter.imag)/(parameter.real*parameter.real+parameter.imag*parameter.imag);
	return temp;
  }
}
complex complex::operator/ (double a) const
{
  complex temp;
  if(a==0)
  {
    cout<<"Denominator cannot be zero. "<<endl;
	exit(0);
  }
  else
  {
    temp.real=this->real/a;
	temp.imag=this->imag/a;
	return temp;
  }
}
const complex & complex::operator= (const complex & parameter)
{
  real=parameter.real;
  imag=parameter.imag;
  return *this;
}
const complex & complex::operator= (double *array)
{
  this->real=array[0];
  this->imag=array[1];
  return *this;
}
const complex & complex::operator= (double a)
{
  this->real=a;
  this->imag=0;
  return *this;
}
const complex & complex::operator= (int a)
{
  this->real=a;
  this->imag=0;
  return *this;
}
complex complex::scale(double factor) const
{
  complex temp;
  temp.real=real*factor;
  temp.imag=imag*factor;
  return temp;
}
complex complex::conjugate() const
{
  complex temp;
  temp.real=real;
  temp.imag=-imag;
  return temp;
}
void complex::print() const
{
   if(imag==0)
   {
     cout<<real;
   }
   else if(real==0)
   {
     if(imag==1)
	 {
	   cout<<"i";
	 }
	 else if(imag==-1)
	 {
	   cout<<"-i";
	 }
	 else 
	 cout<<imag<<"i";
   }
   else 
   {
     if(imag==1)
	 cout<<real<<"+i";
	 else if(imag==-1)
	 cout<<real<<"-i";
	 else if(imag<0)
	 cout<<real<<"-"<<(-imag)<<"i";
	 else 
	 cout<<real<<"+"<<imag<<"i";
   }
}
double complex::norm() const
{
  double norm;
  norm=sqrt(real*real+imag*imag);
  return norm;
}
double complex::arg() const
{
  double arg;
  arg=atan(imag/real);
  arg=arg*180/PI;
  if(real>0 && imag>0)
  {
    return arg;
  }
  else if(real<0 && imag>0)
  {
    arg=arg+180;
	return arg;
  }
  else if(real<0 && imag<0)
  {
    arg=arg+180;
  }
  else if(real<0 && imag==0)
  {
    return 180;
  }
  else if(real==0 && imag==0)
  {
    cout<<"The argument for zero is uncertain. "<<endl;
	exit(0);
  }
  else if(real>0 && imag<0)
  {
    arg=arg+360;
	return arg;
  }
  else if(real==0 && imag>0)
  {
    return 90;
  }
  else if(real==0 && imag<0)
  {
    return 270;
  }
  else if(real>0 && imag==0)
  {
    return 0;
  }
}

double complex::Real()
{
  return this->real;
}
double complex::Imag()
{
  return this->imag;
}

ostream & operator<< (ostream &os, const complex & param)
{
   if(param.imag==0)
   {
     os<<param.real;
   }
   else if(param.real==0)
   {
     if(param.imag==1)
	 {
	   os<<"i";
	 }
	 else if(param.imag==-1)
	 {
	   os<<"-i";
	 }
	 else 
	 os<<param.imag<<"i";
   }
   else 
   {
     if(param.imag==1)
	 os<<param.real<<"+i";
	 else if(param.imag==-1)
	 os<<param.real<<"-i";
	 else if(param.imag<0)
	 os<<param.real<<"-"<<(-param.imag)<<"i";
	 else 
	 os<<param.real<<"+"<<param.imag<<"i";
   }
   return os;
}

complex operator* (double factor, const complex & Z)
{
  complex temp;
  temp.real=factor*Z.real;
  temp.imag=factor*Z.imag;
  return temp;
}

complex operator/ (double numerator, const complex & argument)
{
  complex temp;
  temp.real=numerator*(argument.real/(argument.real*argument.real+argument.imag*argument.imag));
  temp.imag=numerator*(-argument.imag/(argument.real*argument.real+argument.imag*argument.imag));
  return temp;
}
double Re(const complex &Z)
{
  return Z.real;
}
double Im(const complex &Z)
{
  return Z.imag;
}
double norm(const complex & parameter)
{
  return sqrt(parameter.real*parameter.real+parameter.imag*parameter.imag);
}
complex operator+ (double a, const complex &argument)
{
    complex temp(a+argument.real, argument.imag);
    return temp;
}
complex operator- (double a, const complex &argument)
{
    complex temp(a-argument.real, -argument.imag);
    return temp;
}

complex exp(const complex & z)
{
    double x = z.real;
    double y = z.imag;
    return complex(exp(x)*cos(y), exp(x)*sin(y));
}
