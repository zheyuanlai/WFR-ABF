//------------------------------------------
//
//          Simple matrix class
//
//-------------------------------------------

#include "matrix.hpp"

//---------------------------------------------------
//             Creations, destructions
//---------------------------------------------------

//------------ creator ------------
Matrix::Matrix()
{
  row = 0;
  col = 0;
  coord = NULL;
}

//------------ destructor ------------
Matrix::~Matrix()
{
  if (coord != NULL) 
    delete[] coord;
  coord = NULL;
  row = 0;
  col = 0;
}

//------------ initialization ------------
void Matrix::set_size(int Row, int Col)
{
  row = Row;
  col = Col;
  if (coord != NULL) 
    delete[] coord;
  coord = new double[col*row];
}

Matrix::Matrix(int Row, int Col)
{
  row = Row;
  col = Col;
  coord = new double[Row*Col];
}

//------------ creation same matrice ------------
Matrix & Matrix::operator= (const Matrix & a)
{
  if (this != &a) 
    {
      row = a.row;
      col = a.col;
      if (coord != NULL) 
	delete[] coord;
      if (a.coord != NULL) 
	{
	  coord = new double[row*col];
	  for (int i=0;i<row*col;i++) 
	    coord[i] = a.coord[i];
	} else {
	  coord = NULL;
	}
    }
  return *this;
}

//------------ set to zero ------------
void Matrix::zeros()
{
  if (coord != NULL) 
    {
      for (int i = 0; i < col*row; i++) 
	coord[i] = 0.0;
    }
}

//------------ access an element ------------
double & Matrix::operator() (int a,int b)
{
  if (coord != NULL) 
    {
      return coord[a*col+b];
    } else {
      cerr << " ### ERROR ### requiring the component of an empty matrix! " << endl;
      exit(1);
    }
}

//---------------------------------------------------
//        Overloaded operations
//---------------------------------------------------

//------------ addition ------------
Matrix Matrix::operator+ (const Matrix & a)
{
  Matrix ret;
  if ( (row != a.row) && (col != a.col) )
    {
      cerr << " ### ERROR ### adding matrices of differents sizes! " << endl;
      exit(1);
    } else {
      if ((coord != NULL) && (a.coord != NULL)) 
	{
	  ret.row = a.row;
	  ret.col = a.col;
	  ret.coord = new double[col*row];
	  for (int i = 0; i < col*row; i++) 
	    ret.coord[i] = coord[i] + a.coord[i];
	}
    }
  return ret;
}

//------------ substraction ------------
Matrix Matrix::operator- (const Matrix & a)
{
  Matrix ret;
  if ( (row != a.row) && (col != a.col) )
    {
      cerr << " ### ERROR ### substracting matrices of differents sizes! " << endl;
      exit(1);
    } else {
      if ((coord != NULL) && (a.coord != NULL)) 
	{
	  ret.row = a.row;
	  ret.col = a.col;
	  ret.coord = new double[col*row];
	  for (int i = 0; i < col*row; i++) 
	    ret.coord[i] = coord[i] - a.coord[i];
	}
    }
  return ret;
}

//------------ multiplication by a constant ------------
Matrix Matrix::operator* (const double c)
{
  Matrix ret;
  if (coord != NULL) 
    {
      ret.row = row;
      ret.col = col;
      ret.coord = new double[col*row];
      for (int i = 0; i < col*row; i++) 
	ret.coord[i] = coord[i]*c;
    }
  return ret;
}

//------------ multiplication by matrix a ------------
Matrix Matrix::operator*(const Matrix & a) 
{
  Matrix ret;
  if(a.row != col) 
    return ret; 
  ret.set_size(row,a.col);
  ret.zeros();
  for(int i = 0; i < row; i++)
    for(int j = 0; j < a.col; j++)
      for(int k = 0; k < col; k++)
	ret(i,j) = ret(i,j) + coord[i*col+k]*a.coord[k*a.col+j]; 
  return ret; 
}

