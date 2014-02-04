#include "PhysMtr.hpp"

//==============================================================================

// members of "Arr2D"
void Arr2D::AllocMem(double*** _ptr)
{
	this->ptr = new double**[this->N];
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++){
		this->ptr[i] = new double*[this->M];
        for(unsigned int j=0; j<this->M; j++){
			this->ptr[i][j] = new double[nm];
			for(unsigned int k=0; k<nm; k++)
				this->ptr[i][j][k] = (_ptr == NULL ? 0 : _ptr[i][j][k]);
		}
	}
}

void Arr2D::FreeMem(void)
{
	for(unsigned int i=0; i<this->N; i++){
		for(unsigned int j=0; j<this->M; j++) delete[] this->ptr[i][j];
		delete[] this->ptr[i];
	}
	delete[] this->ptr;
}

 Arr2D::Arr2D(const Arr2D& Arr)
{
	this->N=Arr.N; this->M=Arr.M; this->n=Arr.n; this->m=Arr.m;
	this->AllocMem();
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++) this->ptr[i][j][k] = Arr.ptr[i][j][k];
}

matrix Arr2D::operator() (unsigned int _N, unsigned int _M) const
{
	#ifdef _DEBUG
	if(_N>=this->N || _M>=this->M)
		throw " Arr2D::operator() : Incorrect access! ";
	#endif
	double* const pt = this->ptr[_N][_M];
	matrix mt(this->n, this->m);
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)	mt[i][j] = pt[i*n+j];
	return mt;
}

Arr2D Arr2D::operator=(const Arr2D &Arr)
{
	if(this->ptr == Arr.ptr) return *this;
	this->FreeMem();
	this->N=Arr.N; this->M=Arr.M; this->n=Arr.n; this->m=Arr.m;
	this->AllocMem();
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++){
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++) this->ptr[i][j][k] = Arr.ptr[i][j][k];
	}
	return *this;
}

Arr2D Arr2D::operator+=(const Arr2D &arr)
{
	#ifdef _DEBUG
	if(arr.n!=this->n || arr.m!=this->m || arr.N!=this->N || arr.M!=this->M)
		throw " Arr2D::operator+= : Incorrect assignment! ";
	#endif
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++)	this->ptr[i][j][k] += arr.ptr[i][j][k];
	return *this;
}

void Arr2D::ClearArr(void) const
{
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++)	this->ptr[i][j][k] = 0;
}

void Arr2D::insert(unsigned int _N, unsigned int _M, const matrix& mt)
{
	#ifdef _DEBUG
	if(Str(mt)!=this->n || Col(mt)!=this->m || _N>=this->N || _M>=this->M)
		throw " Arr2D::insert() : Incorrect insertion! ";
	#endif
	double* const pt = this->ptr[_N][_M];
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)	pt[i*this->n+j] += mt[i][j];
}

void Arr2D::replace(unsigned int _N, unsigned int _M, const matrix& mt)
{
	#ifdef _DEBUG
	if(Str(mt)!=this->n || Col(mt)!=this->m || _N>=this->N || _M>=this->M)
		throw " Arr2D::replace() : Incorrect insertion! ";
	#endif
	double* const pt = this->ptr[_N][_M];
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)	pt[i*this->n+j] = mt[i][j];
}


// friends of "Arr2D"
double Max(const Arr2D& Arr)
{
	double mx = 0;
	for(unsigned int i=0; i<Arr.n; i++)
		for(unsigned int j=0; j<Arr.m; j++)
			for(unsigned int ii=0; ii<StrArr(Arr); ii++)
				for(unsigned int jj=0; jj<ColArr(Arr); jj++){
					double bf = Arr.ptr[ii][jj][i*Arr.n+j];
					if(bf > mx) mx = bf;
				}
	return mx;
}

matrix SumArr(const Arr2D& Arr)
{
	matrix mt(Arr.n, Arr.m);
	for(unsigned int i=0; i<Arr.n; i++)
		for(unsigned int j=0; j<Arr.m; j++){
			mt[i][j] = 0;
			for(unsigned int ii=0; ii<StrArr(Arr); ii++)
				for(unsigned int jj=0; jj<ColArr(Arr); jj++)
					mt[i][j] += Arr.ptr[ii][jj][i*Arr.n+j];
		}
	return mt;
}

Arr2D Arr2D::operator*(double x)
{
	Arr2D tmp(this->N,this->M,this->n,this->m);

	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++)
				tmp.ptr[i][j][k] = this->ptr[i][j][k]*x;
	return tmp;
}

//==============================================================================
// members of "Arr2DC"
void Arr2DC::AllocMem(complex*** _ptr)
{
	this->ptr = new complex**[this->N];
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++){
		this->ptr[i] = new complex*[this->M];
		for(unsigned int j=0; j<this->M; j++){
			this->ptr[i][j] = new complex[nm];
			for(unsigned int k=0; k<nm; k++)
				this->ptr[i][j][k] = (_ptr == NULL ? complex(0.0) : _ptr[i][j][k]);
		}
	}
}

void Arr2DC::FreeMem(void)
{
	for(unsigned int i=0; i<this->N; i++){
		for(unsigned int j=0; j<this->M; j++) delete[] this->ptr[i][j];
			delete[] this->ptr[i];
	}
	delete[] this->ptr;
}

 Arr2DC::Arr2DC(const Arr2DC& Arr)
{
	this->N=Arr.N; this->M=Arr.M; this->n=Arr.n; this->m=Arr.m;
	this->AllocMem();
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++) this->ptr[i][j][k] = Arr.ptr[i][j][k];
}

matrixC Arr2DC::operator() (unsigned int _N, unsigned int _M) const
{
	#ifdef _DEBUG
	if(_N>=this->N || _M>=this->M)
		throw " Arr2DC::operator() : Incorrect access! ";
	#endif
	complex* const pt = this->ptr[_N][_M];
	matrixC mt(this->n, this->m);
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)	mt[i][j] = pt[i*n+j];
	return mt;
}

Arr2DC Arr2DC::operator=(const Arr2DC &Arr)
{
	if(this->ptr == Arr.ptr) return *this;
	this->FreeMem();
	this->N=Arr.N; this->M=Arr.M; this->n=Arr.n; this->m=Arr.m;
	this->AllocMem();
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++){
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++) this->ptr[i][j][k] = Arr.ptr[i][j][k];
	}
	return *this;
}

Arr2DC Arr2DC::operator+=(const Arr2DC &arr)
{
	#ifdef _DEBUG
	if(arr.n!=this->n || arr.m!=this->m || arr.N!=this->N || arr.M!=this->M)
		throw " Arr2DC::operator+= : Incorrect assignment! ";
	#endif
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++)	this->ptr[i][j][k] += arr.ptr[i][j][k];
	return *this;
}

void Arr2DC::ClearArr(void) const
{
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++)	this->ptr[i][j][k] = 0;
}

void Arr2DC::insert(unsigned int _N, unsigned int _M, const matrixC& mt)
{
	#ifdef _DEBUG
	if(Str(mt)!=this->n || Col(mt)!=this->m || _N>=this->N || _M>=this->M)
		throw " Arr2DC::insert() : Incorrect insertion! ";
	#endif
	complex* const pt = this->ptr[_N][_M];
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)	pt[i*this->n+j] += mt[i][j];
}

void Arr2DC::replace(unsigned int _N, unsigned int _M, const matrixC& mt)
{
	#ifdef _DEBUG
	if(Str(mt)!=this->n || Col(mt)!=this->m || _N>=this->N || _M>=this->M)
		throw " Arr2DC::replace() : Incorrect insertion! ";
	#endif
	complex* const pt = this->ptr[_N][_M];
	for(unsigned int i=0; i<this->n; i++)
		for(unsigned int j=0; j<this->m; j++)	pt[i*this->n+j] = mt[i][j];
}


// friends of "Arr2DC"
matrixC SumArr(const Arr2DC& Arr)
{
	matrixC mt(Arr.n, Arr.m);
	for(unsigned int i=0; i<Arr.n; i++)
		for(unsigned int j=0; j<Arr.m; j++){
			mt[i][j] = 0;
			for(unsigned int ii=0; ii<StrArr(Arr); ii++)
				for(unsigned int jj=0; jj<ColArr(Arr); jj++)
					 mt[i][j] += Arr.ptr[ii][jj][i*Arr.n+j];
		}
	return mt;
}

Arr2DC Arr2DC::operator/=(double x)
{
	for(unsigned int i=0, nm = this->n*this->m; i<this->N; i++)
		for(unsigned int j=0; j<this->M; j++)
			for(unsigned int k=0; k<nm; k++)	this->ptr[i][j][k] /= x;
	return *this;
}
//------------------------------------------------------------------------------

