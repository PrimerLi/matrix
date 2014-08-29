#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

class Vector
{
private:
    int length;
    double *array;
public:
    Vector() {}
    Vector(int);
    Vector(const Vector &);
    void Initialize();
    void Set_Length(int);
    double Sum() const;
    double Product() const;
    Vector operator+ (const Vector &) const;
    Vector operator- (const Vector &) const;
    double operator* (const Vector &) const;
    Vector operator* (double) const;
    double operator[] (int) const;
    Vector operator/ (double) const;
    const Vector & operator= (const Vector &);
    const Vector & operator= (double *);
    double * Get_Element() const;
    int Get_Length() const;
    void Print() const;
    double Norm() const;
    ~Vector();
    friend ostream & operator<< (ostream &, const Vector &);
    friend Vector operator* (double, Vector &);
};

Vector::Vector(int n)
{
    length=n;
    array=new double [length];
    for(int i=0; i<length; i++)
        array[i]=0;
}

Vector::Vector(const Vector &argument)
{
    length=argument.length;
    for(int i=0; i<length; i++)
        array[i]=argument.array[i];
}

void Vector::Initialize()
{
    for(int i=0; i<length; i++)
    {
        cin>>array[i];
    }
}

void Vector::Set_Length(int n)
{
    length=n;
    array=new double [length];
    for(int i=0; i<length; i++)
        array[i]=0;
}

double Vector::Sum() const
{
    double s=0;
    for(int i=0; i<length; i++)
        s=s+array[i];
    return s;
}

double Vector::Product() const
{
    double p=1;
    for(int i=0; i<length; i++)
        p=p*array[i];
    return p;
}

Vector Vector::operator+ (const Vector &argument) const
{
    Vector temp(length);
    for(int i=0; i<length; i++)
    {
        temp.array[i]=array[i]+argument.array[i];
    }
    return temp;
}

Vector Vector::operator- (const Vector &argument) const
{
    Vector temp(length);
    for(int i=0; i<length; i++)
    {
        temp.array[i]=array[i]-argument.array[i];
    }
    return temp;
}

double Vector::operator* (const Vector &argument) const
{
    double product=0;
    for(int i=0; i<length; i++)
    {
        product+=array[i]*argument.array[i];
    }
    return product;
}

Vector Vector::operator* (double factor) const
{
    Vector temp(length);
    for(int i=0; i<length; i++)
        temp.array[i]=factor*array[i];
    return temp;
}

double Vector::operator[] (int i) const
{
    return array[i];
}

Vector Vector::operator/ (double factor) const
{
    Vector temp(length);
    for(int i=0; i<length; i++)
        temp.array[i]=array[i]/factor;
    return temp;
}

const Vector & Vector::operator= (const Vector &argument)
{
    if(this==&argument)
    {
        return *this;
    }
    length=argument.length;
    delete []array;
    array=new double [length];
    for(int i=0; i<length; i++)
        array[i]=argument.array[i];
    return *this;
}

const Vector & Vector::operator= (double *a)
{
    for(int i=0; i<length; i++)
        array[i]=a[i];
    return *this;
}

double * Vector::Get_Element() const
{
    double *a;
    a=new double [length];
    for(int i=0; i<length; i++)
        a[i]=array[i];
    return a;
}

int Vector::Get_Length() const
{
    return length;
}

void Vector::Print() const
{
    for(int i=0; i<length; i++)
        cout<<array[i]<<" ";
    cout<<endl;
}

double Vector::Norm() const
{
    double sum=0;
    for(int i=0; i<length; i++)
    {
        sum=sum+array[i]*array[i];
    }
    sum=sqrt(sum);
    return sum;
}

Vector::~Vector()
{
    delete []array;
}

ostream & operator<< (ostream & os, const Vector &argument)
{
    for(int i=0; i<argument.length; i++)
        os<<argument.array[i]<<" ";
    os<<endl;
    return os;
}

Vector operator* (double factor, Vector &argument)
{
    Vector temp(argument.length);
    for(int i=0; i<argument.length; i++)
    {
        temp.array[i]=factor*argument.array[i];
    }
    return temp;
}