#include "matrix.hpp"
//=============================================================================

// members of matrix
void matrix::AllocMem(void)
{
	this->ptr = new double*[this->n];
	#ifdef _DEBUG
	if(!ptr) throw "matrix::AllocMem(): Error!";
	#endif
	for(unsigned int i=0; i<this->n; i++){
		this->ptr[i] = new double[this->m];
		#ifdef _DEBUG
		if(!ptr[i]) throw "matrix::AllocMem(): Error!";
    #endif
	}
}

void matrix::FreeMem(void)
{
	for(unsigned int i=0; i<this->n; i++) delete[] this->ptr[i];
	delete[] this->ptr;
}

 matrix::matrix(const matrix &mt) // copy constructor
{
	this->n=mt.n; this->m=mt.m; this->AllocMem();
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) this->ptr[i][j] = mt.ptr[i][j];
}

matrix matrix::operator=(const matrix& mt)
{
	if(this->ptr == mt.ptr) return *this;
	this->FreeMem();
	this->n=mt.n; this->m=mt.m;
	this->AllocMem();
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) this->ptr[i][j] = mt.ptr[i][j];
	return *this;

}

matrix matrix::operator+(const matrix &mt) const
{
	#ifdef _DEBUG
	if(this->n!=mt.n || this->m!=mt.m) throw " matrix::operator+ : Error! ";
	#endif
	matrix rez(*this);
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) rez[i][j] += mt[i][j];
	return rez;
}

matrix matrix::operator-(const matrix &mt) const
{
	#ifdef _DEBUG
	if(this->n!=mt.n || this->m!=mt.m) throw " matrix::operator- : Error! ";
	#endif
	matrix rez(*this);
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) rez[i][j] -= mt[i][j];
	return rez;
}

matrix matrix::operator*(double x) const
{
	matrix rez(this->n, this->m);
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) rez[i][j] = x*this->ptr[i][j];
	return rez;
}

// real matrix multiplication
matrix matrix::operator*(const matrix &mt) const
{
	#ifdef _DEBUG
	if(this->m!=mt.n) throw " matrix::operator* : Error! ";
	#endif
	matrix rez(this->n, mt.m);
	for(unsigned int i=0; i<rez.n; i++)
		for(unsigned int j=0; j<rez.m; j++){
			rez[i][j] = 0;
			for(unsigned int k=0; k<this->m; k++) rez[i][j] += this->ptr[i][k]*mt[k][j];
		}
	return rez;
}

matrix matrix::operator*=(double x)
{
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) this->ptr[i][j] *= x;
	return *this;
}

matrix matrix::operator/(double x) const
{
	#ifdef _DEBUG
	if(x+1 == 1) throw " matrix::operator/ : Error! ";
	#endif
	matrix rez(this->n, this->m);
	for(unsigned int i=0; i<rez.n; i++)
		for(unsigned int j=0; j<rez.m; j++) rez[i][j] = this->ptr[i][j]/x;
	return rez;
}

bool matrix::operator==(const matrix &mt) const
{
	if(this->n != mt.n || this->m != mt.m) return false;
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)
			if(this->ptr[i][j] != mt[i][j]) return false;
	return true;
}

void matrix::Fill(double x)
{
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) this->ptr[i][j] = x;
}

void matrix::Identity(void)
{
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) this->ptr[i][j] = (i==j);
}

// friends of "matrix"
double norm(const matrix& mt)
{
	double s = 0;
	for(unsigned int i=0; i<Str(mt); i++)
		for(unsigned int j=0; j<Col(mt); j++) s += SQR(mt[i][j]);
	return s;
}

double Max(const matrix& mt) // searching of the square of maximum
{
	double mx = 0;
	for(unsigned int i=0; i<Str(mt); i++){
		for(unsigned int j=0; j<Col(mt); j++){
			double bf = SQR(mt[i][j]);
			if(mx<bf) mx = bf;
		}
	}
	return mx;
}

std::ofstream& operator<<(std::ofstream& out, const matrix& m)
{
		unsigned int Imax = Str(m), Jmax = Col(m);
		for(unsigned int i=0; i<Imax; i++){
				for(unsigned int j=0; j<Jmax; j++){
						out << (m[i][j]);
					if(((i+1!=Imax) || (j+1!=Jmax))) out << " ";
				}
		}
	return out;
}
//------------------------------------------------------------------------------

// members of matrixC
void matrixC::AllocMem(void)
{
	this->ptr = new complex*[this->n];
	#ifdef _DEBUG
	if(!ptr) throw "matrixC::AllocMem(): Error!";
	#endif
	for(unsigned int i=0; i<this->n; i++){
		this->ptr[i] = new complex[this->m];
		#ifdef _DEBUG
		if(!ptr[i]) throw "matrixC::AllocMem(): Error!";
		#endif
	}
}

void matrixC::FreeMem(void)
{
	for(unsigned int i=0; i<this->n; i++) delete[] this->ptr[i];
	delete[] this->ptr;
}

 matrixC::matrixC(const matrixC &mt)
{
	this->n=mt.n; this->m=mt.m; this->AllocMem();
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) this->ptr[i][j] = mt[i][j];
}

matrixC matrixC::operator=(const matrixC& mt)
{
	if(this->ptr == mt.ptr) return *this;
	this->FreeMem();
	this->n=mt.n; this->m=mt.m;
	this->AllocMem();
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) this->ptr[i][j] = mt.ptr[i][j];
	return *this;
}

matrixC matrixC::operator+(const matrixC& mt) const
{
	#ifdef _DEBUG
	if(this->n!=mt.n || this->m!=mt.m) throw " matrixC::operator+ : Error! ";
	#endif
	matrixC m(mt);
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) m[i][j] += this->ptr[i][j];
	return m;
}

matrixC matrixC::operator-(const matrixC &mt) const
{
	#ifdef _DEBUG
	if(this->n!=mt.n || this->m!=mt.m) throw " matrixC::operator- : Error!. ";
	#endif
	matrixC rez(*this);
	for(unsigned int i=0; i<mt.n; i++)
		for(unsigned int j=0; j<mt.m; j++) rez[i][j] -= mt[i][j];
	return rez;
}

void matrixC::Fill(const complex& z)
{
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) this->ptr[i][j] = z;
}

void matrixC::Identity(void)
{
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) this->ptr[i][j] = (i==j);
}

matrixC matrixC::operator*(const complex& z) const
{
	matrixC rez(this->n, this->m);
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) rez[i][j] = this->ptr[i][j]*z;
	return rez;
}

// complex matrix multiplication
matrixC matrixC::operator*(const matrixC &mt) const
{
	#ifdef _DEBUG
	if(this->m!=mt.n) throw " matrixC::operator* : Error! ";
	#endif
	matrixC rez(this->n, mt.m);
	for(unsigned int i=0; i<rez.n; i++) // this loop is by elements if result matrix
		for(unsigned int j=0; j<rez.m; j++){
			rez[i][j] = 0;
			for(unsigned int k=0; k<this->m; k++)
				rez[i][j] += this->ptr[i][k]*mt[k][j];
		}
	return rez;
}

matrixC matrixC::operator/(const complex &z) const
{
	#ifdef _DEBUG
	if(z+1 == 1) throw " matrixC::operator/ : Error! ";
	#endif
	matrixC rez(this->n, this->m);
	for(unsigned int i=0; i<rez.n; i++)
		for(unsigned int j=0; j<rez.m; j++) rez[i][j] = this->ptr[i][j]/z;
	return rez;
}

matrixC matrixC::operator/=(const complex& z)
{
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++) this->ptr[i][j] /= z;
	return *this;
}

bool matrixC::operator==(const matrixC &mt) const
{
	if(this->n != mt.n || this->m != mt.m) return false;
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)
			if(this->ptr[i][j] != mt[i][j]) return false;
	return true;
}

// friends of "matrixC"
double norm(const matrixC& mt)
{
	double s = 0;
	for(unsigned int i=0; i<Str(mt); i++)
		for(unsigned int j=0; j<Col(mt); j++) s += norm(mt[i][j]);
	return s;
}

double Max(const matrixC& mt)
{
	double mx = 0;
	for(unsigned int i=0; i<Str(mt); i++){ // searching for square of maximal element
		for(unsigned int j=0; j<Col(mt); j++){
			double bf = norm(mt[i][j]);
			if(mx<bf) mx = bf;
		}
	}
	return mx;
}

std::ofstream& operator<<(std::ofstream& out, const matrixC& m)
{
	out.precision(10);
	for(unsigned int i=0; i<Str(m); i++){
		for(unsigned int j=0; j<Col(m); j++)
			out<<"{"<<real(m[i][j])<<", "<<imag(m[i][j])<<"} ";
		out<<std::endl;
	}
	return out;
}
//-----------------------------------------------------------------------------

