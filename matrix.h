#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <assert.h>
#include "vector.h"
#include "complex.h"
using namespace std;

class Matrix
{
private:
    int dim;
    double **mat;
public:
    Matrix () {}
    Matrix(int);
    Matrix(const Matrix & argument);
    void Initialize();
    void Set_Dimension(int);
    void Identity();
    Matrix Transpose() const;
    Matrix operator+ (const Matrix &) const;
    Matrix operator- (const Matrix &) const;
    Matrix operator* (const Matrix &) const;
    Matrix operator* (double) const;
    double * operator* (double *) const;
    Vector operator* (const Vector &) const;
    double & operator() (int, int);
    const Matrix & operator= (const Matrix &);
    const Matrix & operator= (double **);
    double **Get_Element() const;
    void Print() const;
    void Print_mathematica() const;
    double Det() const;
    double Trace() const;
    Matrix Inverse() const;
    Matrix QRDecomposition() const;
    complex * Eigenvalues(int, double **) const;
    ~Matrix();
    friend ostream & operator<< (ostream &, const Matrix &);
    friend Matrix operator* (double, const Matrix &);
};

Matrix::Matrix(int a)
{
    dim=a;
    mat=new double * [dim];
    for(int i=0; i<dim; i++)
    {
        mat[i]=new double [dim];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            mat[i][j]=0;
        }
    }
}

Matrix::Matrix(const Matrix & argument)
{
    dim=argument.dim;
    mat=new double *[dim];
    for(int i=0; i<dim; i++)
    {
        mat[i]=new double [dim];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            mat[i][j]=argument.mat[i][j];
        }
    }
}

void Matrix::Initialize()
{
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            cin>>mat[i][j];
        }
    }
}

void Matrix::Set_Dimension(int n)
{
    dim=n;
    mat=new double *[dim];
    for(int i=0; i<dim; i++)
    {
        mat[i]=new double [dim];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            mat[i][j]=0;
        }
    }
}

void Matrix::Identity()
{
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            if(i==j)
                mat[i][j]=1;
            else
                mat[i][j]=0;
        }
    }
}

Matrix Matrix::Transpose() const
{
    Matrix temp(dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            temp.mat[i][j]=mat[j][i];
        }
    }
    return temp;
}

Matrix Matrix::operator+ (const Matrix &argument) const
{
    Matrix temp(dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            temp.mat[i][j]=mat[i][j]+argument.mat[i][j];
        }
    }
    return temp;
}

Matrix Matrix::operator- (const Matrix &argument) const
{
    Matrix temp(dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            temp.mat[i][j]=mat[i][j]-argument.mat[i][j];
        }
    }
    return temp;
}

Matrix Matrix::operator* (const Matrix &argument) const
{
    Matrix temp(dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            temp.mat[i][j]=0;
            for(int k=0; k<dim; k++)
            {
                temp.mat[i][j]+=mat[i][k]*argument.mat[k][j];
            }
        }
    }
    return temp;
}

Matrix Matrix::operator* (double factor) const
{
    Matrix temp(dim);
    for(int i=0; i<dim; i++)
        for(int j=0; j<dim; j++)
            temp.mat[i][j]=factor*mat[i][j];
    return temp;
}

double * Matrix::operator* (double *argument) const
{
    double *result;
    result=new double [dim];
    for(int i=0; i<dim; i++)
    {
        result[i]=0;
        for(int j=0; j<dim; j++)
        {
            result[i]+=mat[i][j]*argument[j];
        }
    }
    return result;
}

Vector Matrix::operator* (const Vector &argument) const
{
    assert(argument.Get_Length()==dim);
    double *temp;
    temp=argument.Get_Element();
    Vector Result(dim);
    double *result;
    result=new double [dim];
    for(int i=0; i<dim; i++)
    {
        result[i]=0;
        for(int j=0; j<dim; j++)
        {
            result[i]=result[i]+mat[i][j]*temp[j];
        }
    }
    Result=result;
    delete []temp;
    delete []result;
    return Result;
}

double & Matrix::operator() (int i, int j)
{
    return mat[i][j];
}

const Matrix & Matrix::operator= (const Matrix &argument)
{
    if(this==&argument)
    {
        return *this;
    }
    dim=argument.dim;
    for(int i=0; i<dim; i++)
        delete []mat[i];
    delete []mat;
    mat=new double *[dim];
    for(int i=0; i<dim; i++)
        mat[i]=new double [dim];
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            mat[i][j]=argument.mat[i][j];
        }
    }
    return *this;
}

const Matrix & Matrix::operator= (double **a)
{
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            mat[i][j]=a[i][j];
        }
    }
    return *this;
}

double Matrix::Det() const
{
    int count=0;
    double **a;
    a=new double *[dim];
    for(int i=0; i<dim; i++)
    {
        a[i]=new double [dim];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            a[i][j]=mat[i][j];
        }
    }

    int column=0;
    while(column<dim)
    {
        int nonzero=column;
        while(a[nonzero][column]==0)
        {
            nonzero++;
            if(nonzero==dim)
            {
                cout<<"This matrix is singular."<<endl;
                exit(0);
            }
        }
         
        if(nonzero!=column)
        {
	    count++;
            double temp;
            for(int j=0; j<dim; j++)
            {
                temp=a[column][j];
                a[column][j]=a[nonzero][j];
                a[nonzero][j]=temp;
            }
        }
        double factor=1;
        for(int i=column+1; i<dim; i++)
        {
            factor=-a[i][column]/a[column][column];
            for(int j=0; j<dim; j++)
            {
                a[i][j]=a[i][j]+factor*a[column][j];
            }
        }
        column++;
    }   
    double temp=1;
    for(int i=0; i<dim; i++)
    {
    	temp=temp*a[i][i];
    }
    if(count%2==0)
        return temp;
    else
        return -temp;
}
double ** Matrix::Get_Element() const
{
    double **temp;
    temp=new double *[dim];
    for(int i=0; i<dim; i++)
    {
        temp[i]=new double [dim];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            temp[i][j]=mat[i][j];
        }
    }
    return temp;
}
void Matrix::Print() const
{
    double small=0.000001;
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            if(fabs(mat[i][j])<small)
                cout<<0<<"  ";
            else
                cout<<mat[i][j]<<"  ";
        }
        cout<<endl;
    }
}

double Matrix::Trace() const
{
    double sum=0;
    for(int i=0; i<dim; i++)
        sum=sum+mat[i][i];
    return sum;
}

void Matrix::Print_mathematica() const
{
    double small=0.000001;
    cout<<"{";
    for(int i=0; i<dim; i++)
    {
        cout<<"{";
        for(int j=0; j<dim; j++)
        {
            if(j!=dim-1)
                cout<<mat[i][j]<<", ";
            else
                cout<<mat[i][j];
        }
        if(i!=dim-1)
            cout<<"}, ";
        else
            cout<<"}";
    }
    cout<<"}";
    cout<<endl;
}

Matrix Matrix::Inverse() const
{
    double small=0.000001;
    
    double **a;
    a=new double *[dim];
    for(int i=0; i<dim; i++)
    {
        a[i]=new double [dim];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            a[i][j]=mat[i][j];
        }
    }
    
    double **augment;
    int dim_a=2*dim;
    augment=new double *[dim];
    for(int i=0; i<dim; i++)
    {
        augment[i]=new double [dim_a];
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            augment[i][j]=a[i][j];
        }
    }
    for(int i=0; i<dim; i++)
    {
        for(int j=dim; j<dim_a; j++)
        {
            if(i==j-dim)
                augment[i][j]=1;
            else
                augment[i][j]=0;
        }
    }
    
    int column=0;
    while(column<dim)
    {
        int nonzero=column;
        while(augment[nonzero][column]==0)
        {
            nonzero++;
            if(nonzero==dim)
            {
                cout<<"This matrix is singular."<<endl;
                exit(0);
            }
        }
        
        /*for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim_a; j++)
            {
                if(fabs(augment[i][j])<small)
                    cout<<0<<"  ";
                else
                    cout<<augment[i][j]<<"  ";
            }
            cout<<endl;
        }
        cout<<endl;*/
         
        if(nonzero!=column)
        {
            double temp;
            for(int j=0; j<dim_a; j++)
            {
                temp=augment[column][j];
                augment[column][j]=augment[nonzero][j];
                augment[nonzero][j]=temp;
            }
        }
        double factor=1;
        for(int i=column+1; i<dim; i++)
        {
            factor=-augment[i][column]/augment[column][column];
            for(int j=0; j<dim_a; j++)
            {
                augment[i][j]=augment[i][j]+factor*augment[column][j];
            }
        }
        column++;
    }
    
    /*for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim_a; j++)
        {
            if(fabs(augment[i][j])<small)
                cout<<0<<"  ";
            else
                cout<<augment[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl;*/
    
    int row=dim-1;
    while(row>0)
    {
        for(int i=row-1; i>=0; i--)
        {
            double factor=-augment[i][row]/augment[row][row];
            for(int j=0; j<dim_a; j++)
            {
                augment[i][j]=augment[i][j]+augment[row][j]*factor;
            }
        }
        row--;
        /*for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim_a; j++)
            {
                if(fabs(augment[i][j])<small)
                    cout<<0<<"  ";
                else
                    cout<<augment[i][j]<<"  ";
            }
            cout<<endl;
        }
        cout<<endl;*/
    }
    
    row=0;
    while(row<dim)
    {
        double factor=(1/augment[row][row]);
        for(int j=row; j<dim_a; j++)
        {
            augment[row][j]=augment[row][j]*factor;
        }
        row++;
    }
    
    /*for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim_a; j++)
        {
            if(fabs(augment[i][j])<small)
                cout<<0<<"  ";
            else
                cout<<augment[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl;*/
    
    Matrix temp(dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            temp.mat[i][j]=augment[i][j+dim];
        }
    }
    for(int i=0; i<dim; i++)
    {
        delete []a[i];
        delete []augment[i];
    }
    delete []a;
    delete []augment;
    return temp;
}

Matrix Matrix::QRDecomposition() const
{
    Vector *a;
    Vector *u;
    Vector *e;
    a=new Vector[dim];
    u=new Vector[dim];
    e=new Vector[dim];
    for(int i=0; i<dim; i++)
    {
        a[i].Set_Length(dim);
        u[i].Set_Length(dim);
        e[i].Set_Length(dim);
    }
    
    int column=0;
    while(column<dim)
    {
        double *temp;
        temp=new double [dim];
        for(int i=0; i<dim; i++)
            temp[i]=mat[i][column];
        a[column]=temp;
        delete []temp;
        column++;
    }
    
    column=0;
    while(column<dim)
    {
        if(column==0)
        {
            u[column]=a[column];
            e[column]=u[column]/u[column].Norm();
        }
        else
        {
            Vector project(dim);
            for(int i=0; i<=column-1; i++)
            {
                project=project+(a[column]*e[i])*e[i];
            }
            u[column]=a[column]-project;
            e[column]=u[column]/u[column].Norm();
        }
        column++;
    }
    
    /*cout<<"Here I begin. "<<endl;
    for(int i=0; i<dim; i++)
    {
        cout<<e[i];
    }
    cout<<"Here I end. "<<endl;*/
    
    Matrix Q(dim);
    column=0;
    while(column<dim)
    {
        double *temp;
        temp=new double [dim];
        temp=e[column].Get_Element();
        for(int i=0; i<dim; i++)
        {
            Q.mat[i][column]=temp[i];
        }
        delete []temp;
        column++;
    }
    
    delete []a;
    delete []u;
    delete []e;
    return Q;
}

complex * Matrix::Eigenvalues(int Iteration_max, double **UMatrix) const
{
    complex *Eigen_Value_Vector;
    Eigen_Value_Vector=new complex [dim];
    for(int i=0; i<dim; i++)
    {
        Eigen_Value_Vector[i]=0;
    }
    
    Matrix *A;
    Matrix *Q;
    Matrix *R;
    A=new Matrix [Iteration_max+1];
    Q=new Matrix [Iteration_max];
    R=new Matrix [Iteration_max];
    for(int i=0; i<Iteration_max+1; i++)
        A[i].Set_Dimension(dim);
    for(int i=0; i<Iteration_max; i++)
    {
        Q[i].Set_Dimension(dim);
        R[i].Set_Dimension(dim);
    }
    
    int iteration=0;
    while(iteration<Iteration_max)
    {
        if(iteration==0)
        {
            A[0]=*this;
            Q[0]=A[0].QRDecomposition();
            R[0]=Q[0].Transpose()*A[0];
        }
        else
        {
            A[iteration]=R[iteration-1]*Q[iteration-1];
            Q[iteration]=A[iteration].QRDecomposition();
            R[iteration]=Q[iteration].Transpose()*A[iteration];
        }
        iteration++;
    }
    A[Iteration_max]=R[Iteration_max-1]*Q[Iteration_max-1];
    Matrix T(dim);
    T=A[Iteration_max];
    
    int flag=0;
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<i; j++)
        {
            if(i!=j && fabs(T.mat[i][j])>0.01)
            {
                flag++;
            }
        }
    }
    
    Matrix U(dim);
    U.Identity();
    for(int i=0; i<Iteration_max; i++)
        U=U*Q[i];
    
    double **temp_Umatrix;
    temp_Umatrix=new double *[dim];
    for(int i=0; i<dim; i++)
        temp_Umatrix[i]=new double [dim];
    temp_Umatrix=U.Get_Element();
    
    //cout<<"Similarity transformation matrix is "<<endl;
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            UMatrix[i][j]=temp_Umatrix[i][j];
            //cout<<UMatrix[i][j]<<"  ";
        }
        //cout<<endl;
    }
    
    for(int i=0; i<dim; i++)
        delete []temp_Umatrix[i];
    delete []temp_Umatrix;
    
    if(flag>0)
    {
        double small=0.001;
        cout<<"This matrix may contain complex eigenvalues. "<<endl;
        complex *complex_eigenvalues;
        complex_eigenvalues=new complex [dim];
        int row=0;
        int count=0;
        while(row<dim)
        {
            if(row<dim-1 && fabs(T.mat[row+1][row])<small)
            {
                //cout<<"One real eigenvalue. "<<endl;
                complex_eigenvalues[count]=T.mat[row][row];
                //cout<<complex_eigenvalues[count]<<endl;
                count++;
                row++;
            }
            else if(row<dim-1 && fabs(T.mat[row+1][row])>=small)
            {
                //cout<<"Two complex eigenvalues. "<<endl;
                complex e1;
                complex e2;
                double **block;
                block=new double *[2];
                for(int i=0; i<2; i++)
                    block[i]=new double [2];
                for(int i=0; i<2; i++)
                {
                    for(int j=0; j<2; j++)
                    {
                        block[i][j]=T.mat[row+i][row+j];
                    }
                }
                double trace=block[0][0]+block[1][1];
                double determinant=block[0][0]*block[1][1]-block[1][0]*block[0][1];
                double discriminant=trace*trace-4*determinant;
                if(discriminant>0.1*small)
                {
                    cout<<"Something may be wrong. The Eigenvalue should be complex. "<<endl;
                    exit(0);
                }
                else
                {
                    complex sqrt_delta(0, sqrt(-discriminant));
                    e1=0.5*(trace+sqrt_delta);
                    e2=e1.conjugate();
                }
                complex_eigenvalues[count]=e1;
                complex_eigenvalues[count+1]=e2;
                //cout<<complex_eigenvalues[count]<<"\t"<<complex_eigenvalues[count+1]<<endl;
                count=count+2;
                row=row+2;
                
                for(int i=0; i<2; i++)
                    delete []block[i];
                delete []block;
            }
            else if(row==dim-1)
            {
                //cout<<"Last eigenvalue. "<<endl;
                //cout<<T.mat[row][row]<<endl;
                complex_eigenvalues[row]=T.mat[row][row];
                break;
            }
        }
        //cout<<"Here are your complex eigenvalues: "<<endl;
        for(int i=0; i<dim; i++)
        {
            //cout<<complex_eigenvalues[i]<<"     ";
            Eigen_Value_Vector[i]=complex_eigenvalues[i];
        }
        //cout<<endl;
        /*complex sum;
        for(int i=0; i<dim; i++)
            sum=sum+complex_eigenvalues[i];
        complex product;
        product=1;
        for(int i=0; i<dim; i++)
            product=product*complex_eigenvalues[i];
        cout<<"The trace of the matrix is "<<sum<<endl;
        cout<<"The determinant of the matrix is "<<product<<endl;*/
        delete []complex_eigenvalues;
    }
    else if(flag==0)
    {
        /*cout<<"U=";
        U.Print_mathematica();
        U.Print();
        cout<<"T=";
        T.Print_mathematica();
        T.Print();*/
        Matrix Test(dim);
        Test=U*T*U.Transpose();
        //cout<<"Test: "<<endl;
        //Test.Print();
        Matrix Error(dim);
        Error=*this-Test;
        if((Error*Error.Transpose()).Trace()<0.001)
        {
            //cout<<"\nSuccess"<<endl;
            double *real_eigenvalues;
            real_eigenvalues=new double [dim];
            for(int i=0; i<dim; i++)
                real_eigenvalues[i]=T.mat[i][i];
            //cout<<"Here are the real eigenvalues: "<<endl;
            for(int i=0; i<dim; i++)
            {
                //cout<<real_eigenvalues[i]<<"    ";
                Eigen_Value_Vector[i]=real_eigenvalues[i];
            }
            //cout<<endl;
            /*double sum=0;
            for(int i=0; i<dim; i++)
                sum=sum+real_eigenvalues[i];
            double product=1;
            for(int i=0; i<dim; i++)
                product=product*real_eigenvalues[i];
            cout<<"The trace of the matrix is "<<sum<<endl;
            cout<<"The determinant of the matrix is "<<product<<endl;*/
            delete []real_eigenvalues;
        }
        else
        {
            cout<<"\nFail"<<endl;
            exit(0);
        }
    }
    
    delete []A;
    delete []Q;
    delete []R;
    
    return Eigen_Value_Vector;
}

Matrix::~Matrix()
{
    for(int i=0; i<dim; i++)
        delete []mat[i];
    delete []mat;
}

ostream & operator<< (ostream & os, const Matrix &argument)
{
    double small=0.00001;
    for(int i=0; i<argument.dim; i++)
    {
        for(int j=0; j<argument.dim; j++)
        {
            if(fabs(argument.mat[i][j])<small)
                os<<0<<"    ";
            else
                os<<argument.mat[i][j]<<"   ";
        }
        os<<endl;
    }
    return os;
}

Matrix operator* (double factor, const Matrix &argument)
{
    Matrix temp(argument.dim);
    for(int i=0; i<argument.dim; i++)
        for(int j=0; j<argument.dim; j++)
            temp.mat[i][j]=factor*argument.mat[i][j];
    return temp;
}
