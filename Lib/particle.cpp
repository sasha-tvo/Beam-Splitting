#include "particle.hpp"

#define Sqrt3_2 0.86602540378443864676372317075294  // sqrt(3)/2
#define Sqrt3  1.7320508075688772935274463415059  // sqrt(3)
#define Sqrt6  2.4494897427831780981972840747059  // sqrt(6)
#define Sqrt2  1.4142135623730950488016887242097  // sqrt(2)
#define Sqrt1_3 0.57735026918962576450914878050196  // 1/sqrt(3)

//------------------------------------------------------------------------------
// the initialization of prism vertices
void Prism::SetVertices(void)
{
	const double bf1 = 0.5*this->w, bf2 = Sqrt3_2*this->w;

	this->p[3].x  = this->p[9].x  = -(this->p[0].x = this->p[6].x = this->w);
	this->p[3].y  = this->p[9].y  =   this->p[0].y = this->p[6].y = 0;
	this->p[6].z  = this->p[9].z  = -(this->p[0].z = this->p[3].z = this->h);

	this->p[2].x  = this->p[8].x  = -(this->p[1].x = this->p[7].x = bf1);
	this->p[2].y  = this->p[8].y  =   this->p[1].y = this->p[7].y = bf2;
	this->p[7].z  = this->p[8].z  = -(this->p[1].z = this->p[2].z = this->h);

	this->p[4].x  = this->p[10].x = -(this->p[5].x = this->p[11].x = bf1);
	this->p[4].y  = this->p[10].y =   this->p[5].y = this->p[11].y = -bf2;
	this->p[10].z = this->p[11].z = -(this->p[4].z = this->p[5].z  = this->h);
}

// the initialization of prism facets
void Prism::SetFacets(void)
{
	int* q;
	// the basal facets
	q = this->Gr[0]; q[0]=0; q[1]=1; q[2]=2; q[3]=3; q[4]=4; q[5]=5; q[6]=0;
	q = this->Gr[7]; q[0]=6; q[1]=11; q[2]=10; q[3]=9; q[4]=8; q[5]=7; q[6]=6;
	// the rectangular facets
	for(int i=1; i<=6; i++) {
		q=this->Gr[i]; q[0]=q[4]=i-1; q[1]=i+5; q[3]=i%6; q[2]=q[3]+6; q[5]=q[6]=-1;
	}
}
//==============================================================================
void Pyramid::SetVertices(void)
{
	const double bf1 = 0.5*this->w, bf2 = Sqrt3_2*this->w;

	this->p[0].x=this->p[0].y=0; this->p[0].z=this->h; // top of the pyramid

	this->p[1].x = this->w;	 this->p[1].y = 0;
	this->p[2].x = bf1;		 this->p[2].y = bf2;
	this->p[3].x = -bf1;	 this->p[3].y = bf2;
	this->p[4].x = -this->w; this->p[4].y = 0;
	this->p[5].x = -bf1;	 this->p[5].y = -bf2;
	this->p[6].x = bf1;		 this->p[6].y = -bf2;

	this->p[1].z = this->p[2].z = this->p[3].z =
	this->p[4].z = this->p[5].z = this->p[6].z = -this->h;
}
// the initialization of the pyramid facets
void Pyramid::SetFacets(void)
{
	int* q; // the bottom
	q = this->Gr[0]; q[0]=1; q[1]=6; q[2]=5; q[3]=4; q[4]=3; q[5]=2; q[6]=1;
	for(int i=1; i<=6; i++) { // the tapered facets
		q=this->Gr[i]; q[0]=q[3]=0; q[1]=i; q[2]=i%6+1; q[4]=q[5]=q[6]=-1;
	} 
}
//==============================================================================
// members of class "Bullet":
// the initialization of bullet vertices
void Bullet::SetVertices(void)
{
	const double bf1 = 0.5*this->w, bf2 = Sqrt3_2*this->w;

	this->p[3].x  = this->p[9].x  = -(this->p[0].x = this->p[6].x = this->w);
	this->p[3].y  = this->p[9].y  =   this->p[0].y = this->p[6].y = 0;
	this->p[6].z  = this->p[9].z  = -(this->p[0].z = this->p[3].z = this->h);

	this->p[2].x  = this->p[8].x  = -(this->p[1].x = this->p[7].x = bf1);
	this->p[2].y  = this->p[8].y  =   this->p[1].y = this->p[7].y = bf2;
	this->p[7].z  = this->p[8].z  = -(this->p[1].z = this->p[2].z = this->h);

	this->p[4].x  = this->p[10].x = -(this->p[5].x = this->p[11].x = bf1);
	this->p[4].y  = this->p[10].y =   this->p[5].y = this->p[11].y = -bf2;
	this->p[10].z = this->p[11].z = -(this->p[4].z = this->p[5].z  = this->h);
	// top of pyramid
	this->p[12].x = this->p[12].y = 0; this->p[12].z = -this->h-this->hp;
}
// the initialization of bullet facets
void Bullet::SetFacets(void)
{
	// basal facet (bottom)
	int* q=this->Gr[0]; q[0]=0; q[1]=1; q[2]=2; q[3]=3; q[4]=4; q[5]=5; q[6]=0;
	for(int i=1; i<=12; i++) { // ..and all other facets
		const int ij = i%6;
		q=this->Gr[i]; q[5]=q[6]=-1;
		(i<=6) ? (q[0]=q[4]=i-1, q[1]=i+5, q[3]=ij, q[2]=q[3]+6) : // prism
			   (q[0]=q[3]=12, q[1]=ij+6, q[2]=i-1, q[4]=-1);   // pyramid
	}
}
//==============================================================================
// members of class "TaperedPrism"
// the initialization of the prism vertices
void TaperedPrism::SetVertices(void)
{
	const double b1 = 0.5*this->w, b2 = Sqrt3_2*this->w,
				 b3 = 0.5*this->wp, b4 = Sqrt3_2*this->wp;

	this->p[0].x = this->wp; 	this->p[0].y = 0;		this->p[0].z = this->h;
	this->p[1].x = b3;			this->p[1].y = b4;		this->p[1].z = this->h;
	this->p[2].x = -b3;			this->p[2].y = b4;		this->p[2].z = this->h;
	this->p[3].x = -this->wp;	this->p[3].y = 0;		this->p[3].z = this->h;
	this->p[4].x = -b3;			this->p[4].y = -b4;		this->p[4].z = this->h;
	this->p[5].x = b3;			this->p[5].y = -b4;		this->p[5].z = this->h;

	this->p[6].x =  this->w;	this->p[6].y =  0;		this->p[6].z = this->hp;
	this->p[7].x =  b1;			this->p[7].y =  b2;		this->p[7].z = this->hp;
	this->p[8].x = -b1;			this->p[8].y =  b2;		this->p[8].z = this->hp;
	this->p[9].x = -this->w; 	this->p[9].y =  0;		this->p[9].z = this->hp;
	this->p[10].x = -b1;		this->p[10].y = -b2;	this->p[10].z = this->hp;
	this->p[11].x = b1;			this->p[11].y = -b2;	this->p[11].z = this->hp;

	this->p[12].x = this->w;	this->p[12].y = 0;		this->p[12].z = -this->hp;
	this->p[13].x = b1;			this->p[13].y = b2;		this->p[13].z = -this->hp;
	this->p[14].x = -b1;		this->p[14].y = b2;		this->p[14].z = -this->hp;
	this->p[15].x = -this->w; 	this->p[15].y = 0;		this->p[15].z = -this->hp;
	this->p[16].x = -b1;		this->p[16].y = -b2;	this->p[16].z = -this->hp;
	this->p[17].x = b1;			this->p[17].y = -b2;	this->p[17].z = -this->hp;

	this->p[18].x = this->wp;	this->p[18].y = 0;		this->p[18].z = -this->h;
	this->p[19].x = b3;			this->p[19].y = b4;		this->p[19].z = -this->h;
	this->p[20].x = -b3;		this->p[20].y = b4;		this->p[20].z = -this->h;
	this->p[21].x = -this->wp; 	this->p[21].y = 0;		this->p[21].z = -this->h;
	this->p[22].x = -b3;		this->p[22].y = -b4;	this->p[22].z = -this->h;
	this->p[23].x = b3;			this->p[23].y = -b4;	this->p[23].z = -this->h;
}
// the initialization of the prism facets
void TaperedPrism::SetFacets(void)
{
	int* q;  // the basal facets:
	q=this->Gr[0]; q[0]=0; q[1]=1; q[2]=2; q[3]=3; q[4]=4; q[5]=5; q[6]=0;
	q=this->Gr[19]; q[0]=18; q[1]=23; q[2]=22; q[3]=21; q[4]=20; q[5]=19; q[6]=18;
	for(int i=1; i<=18; i++) {  // the rectangular facets
		q=this->Gr[i];
		q[0]=q[4]=i-1; q[1]=i+5; q[3]=(i%6) ? i : i-6; q[2]=q[3]+6; q[5]=q[6]=-1;
	} //*/
}
//==============================================================================
// members of class "CavityPrism":
// the initialization of prism vertices
void CavityPrism::SetVertices(void)
{
	const double b1 = 0.5*this->w, b2 = Sqrt3_2*this->w;

	this->p[3].x  = this->p[9].x  = -(this->p[0].x = this->p[6].x = this->w);
	this->p[3].y  = this->p[9].y  =   this->p[0].y = this->p[6].y = 0;
	this->p[6].z  = this->p[9].z  = -(this->p[0].z = this->p[3].z = this->h);

	this->p[2].x  = this->p[8].x  = -(this->p[1].x = this->p[7].x = b1);
	this->p[2].y  = this->p[8].y  =   this->p[1].y = this->p[7].y = -b2;
	this->p[7].z  = this->p[8].z  = -(this->p[1].z = this->p[2].z = this->h);

	this->p[4].x  = this->p[10].x = -(this->p[5].x = this->p[11].x = b1);
	this->p[4].y  = this->p[10].y =   this->p[5].y = this->p[11].y = b2;
	this->p[10].z = this->p[11].z = -(this->p[4].z = this->p[5].z = this->h);

	this->p[12].x = this->p[12].y = 0; this->p[12].z = this->c-this->h;
}
// the initialization of bullet facets
void CavityPrism::SetFacets(void)
{
	int* q;
	q = this->Gr[0];  q[0]=0;  q[1]=1;  q[2]=2;  q[3]=3;  q[4]=4; q[5]=5; q[6]=0;
	q = this->Gr[1];  q[0]=0;  q[1]=6;  q[2]=7;  q[3]=1;  q[4]=0; q[5]=q[6]=-1;
	q = this->Gr[2];  q[0]=1;  q[1]=7;  q[2]=8;  q[3]=2;  q[4]=1; q[5]=q[6]=-1;
	q = this->Gr[3];  q[0]=8;  q[1]=9;  q[2]=3;  q[3]=2;  q[4]=8; q[5]=q[6]=-1;
	q = this->Gr[4];  q[0]=3;  q[1]=9;  q[2]=10; q[3]=4;  q[4]=3; q[5]=q[6]=-1;
	q = this->Gr[5];  q[0]=4;  q[1]=10; q[2]=11; q[3]=5;  q[4]=4; q[5]=q[6]=-1;
	q = this->Gr[6];  q[0]=0;  q[1]=5;  q[2]=11; q[3]=6;  q[4]=0; q[5]=q[6]=-1;
	q = this->Gr[7];  q[0]=6;  q[1]=11; q[2]=12; q[3]=6;  q[4]=q[5]=q[6]=-1;
	q = this->Gr[8];  q[0]=11; q[1]=10; q[2]=12; q[3]=11; q[4]=q[5]=q[6]=-1;
	q = this->Gr[9];  q[0]=10; q[1]=9;  q[2]=12; q[3]=10; q[4]=q[5]=q[6]=-1;
	q = this->Gr[10]; q[0]=9;  q[1]=8;  q[2]=12; q[3]=9;  q[4]=q[5]=q[6]=-1;
	q = this->Gr[11]; q[0]=8;  q[1]=7;  q[2]=12; q[3]=8;  q[4]=q[5]=q[6]=-1;
	q = this->Gr[12]; q[0]=7;  q[1]=6;  q[2]=12; q[3]=7;  q[4]=q[5]=q[6]=-1;
}
//==============================================================================

// members of class "Cup"
void Cup::SetVertices(void)
{
	double bf1 = 0.5*this->wp, bf2 = Sqrt3_2*this->wp, bf3 = -this->hc;

	this->p[0].x = this->p[6].x = this->wp;
	this->p[0].y = this->p[6].y = 0;
	this->p[0].z = this->hp;  this->p[6].z = 0;

	this->p[1].x = this->p[7].x = bf1;
	this->p[1].y = this->p[7].y = bf2;
	this->p[1].z = this->hp;  this->p[7].z = 0;

	this->p[2].x = this->p[8].x = -bf1;
	this->p[2].y = this->p[8].y = bf2;
	this->p[2].z = this->hp;  this->p[8].z = 0;

	this->p[3].x = this->p[9].x = -this->wp;
	this->p[3].y = this->p[9].y = 0;
	this->p[3].z = this->hp;  this->p[9].z = 0;

	this->p[4].x = this->p[10].x = -bf1;
	this->p[4].y = this->p[10].y = -bf2;
	this->p[4].z = this->hp;  this->p[10].z = 0;

	this->p[5].x = this->p[11].x = bf1;
	this->p[5].y = this->p[11].y = -bf2;
	this->p[5].z = this->hp;  this->p[11].z = 0;

	bf1 = 0.5*this->wc; bf2 = Sqrt3_2*this->wc;

	this->p[12].x = this->wc;	this->p[12].y = 0;		this->p[12].z = bf3;
	this->p[13].x = bf1;		this->p[13].y = bf2;	this->p[13].z = bf3;
	this->p[14].x = -bf1;		this->p[14].y = bf2;	this->p[14].z = bf3;
	this->p[15].x = -this->wc;	this->p[15].y = 0;		this->p[15].z = bf3;
	this->p[16].x = -bf1;		this->p[16].y = -bf2;	this->p[16].z = bf3;
	this->p[17].x = bf1;		this->p[17].y = -bf2;	this->p[17].z = bf3;
}
// the initialization of cup facets
void Cup::SetFacets(void)
{
	int* q; // basal facets
	q=this->Gr[0]; q[0]=0; q[1]=1; q[2]=2; q[3]=3; q[4]=4; q[5]=5;  q[6]=0;
	q=this->Gr[13]; q[0]=12; q[1]=17; q[2]=16; q[3]=15; q[4]=14; q[5]=13; q[6]=12;
	for(int i=1; i<=12; i++) {  // the rectangular facets
		q=this->Gr[i];
		q[0]=q[4]=i-1; q[1]=i+5; q[3]=(i%6) ? i : i-6; q[2]=q[3]+6; q[5]=q[6]=-1;
	} //*/
}
//==============================================================================
void TrianglePrism::SetVertices(void)
{
	/*
	const double bf1 = 0.5*this->w, bf2 = Sqrt3_2*this->w;
	this->p[0].x = this->p[3].x = 0;
	this->p[0].y = this->p[3].y = -this->w;
	this->p[3].z = this->p[5].z =  this->p[4].z = -(this->p[0].z = this->p[1].z = this->p[2].z = this->h);
	this->p[2].x = this->p[5].x = -(this->p[1].x =  this->p[4].x = bf2);
	this->p[2].y = this->p[4].y =  this->p[1].y =  this->p[5].y = bf1;
	*/
	const double bf1 = 0.5*this->w, bf2=Sqrt3_2*this->w;
	this->p[0].x = this->p[3].x = this->w;
	this->p[0].y = this->p[3].y = 0;
	this->p[3].z = this->p[4].z =   this->p[5].z = -(this->p[0].z = this->p[1].z = this->p[2].z = this->h);
	this->p[1].x = this->p[4].x =   this->p[2].x =   this->p[5].x = bf1;
	this->p[2].y = this->p[5].y = -(this->p[1].y =   this->p[4].y = bf2);
}
//==============================================================================
void TrianglePrism::SetFacets(void)
{
	int* q;
	// the basal facets
	q = this->Gr[0]; q[0]=0; q[1]=1; q[2]=2; q[3]=0; q[4]=-1;
	q = this->Gr[4]; q[0]=3; q[1]=5; q[2]=4; q[3]=3; q[4]=-1;
	// the rectangular facets
	for(int i=1; i<=3; i++) {
		q=this->Gr[i]; q[0]=q[4]=i-1; q[1]=i+2; q[3]=i%3; q[2]=q[3]+3;
	}
}
//==============================================================================
void Parallelepiped::SetVertices(void)
{
	this->p[1].x = this->p[5].x = -(this->p[0].x = this->p[4].x = this->w/2.);
	this->p[2].x = this->p[6].x = -(this->p[3].x = this->p[7].x = this->w/2.);

	this->p[3].y = this->p[7].y = -(this->p[0].y = this->p[4].y = this->wp/2.);
	this->p[2].y = this->p[6].y = -(this->p[1].y = this->p[5].y = this->wp/2.);

	  this->p[4].z = this->p[5].z = this->p[6].z = this->p[7].z =
	-(this->p[0].z = this->p[1].z = this->p[2].z = this->p[3].z = this->h);
}
//==============================================================================
void Parallelepiped::SetFacets(void)
{
	int* q;
	// the basal facets
	q = this->Gr[0]; q[0]=0; q[1]=1; q[2]=2; q[3]=3; q[4]=0;
	q = this->Gr[5]; q[0]=4; q[1]=7; q[2]=6; q[3]=5; q[4]=4;
	// the rectangular facets
	for(int i=1; i<=4; i++) {
		q=this->Gr[i]; q[0]=q[4]=i-1; q[1]=i+3; q[3]=i%4; q[2]=q[3]+4;
	}
}
//==============================================================================
void Tetrahedron::SetVertices(void)
{
	double	bf1 = this->w*Sqrt3/6., bf2 = this->w/2.,
			bf3 = this->w*Sqrt6/12., bf4 = this->w*Sqrt3/3.;

	this->p[0].x = bf1;		this->p[0].y = bf2;		this->p[0].z = -bf3;
	this->p[1].x = -bf4;	this->p[1].y = 0;		this->p[1].z = -bf3;
	this->p[2].x = bf1;		this->p[2].y = -bf2;	this->p[2].z = -bf3;
	this->p[3].x = 0;		this->p[3].y = 0;		this->p[3].z = this->w*Sqrt6/4.;
}
//==============================================================================
void Tetrahedron::SetFacets(void)
{
	int* q;
	// the basal facets
	q = this->Gr[0]; q[0]=0; q[1]=2; q[2]=1; q[3]=0;
	// the triangular facets
	for(int i=1; i<=3; i++) {
		q=this->Gr[i]; q[0]=q[3]=i-1; q[1]=i%3; q[2]=3;
	}
}
//==============================================================================
void Dodecahedron::SetVertices(void)
{
	const double bf1 = Sqrt3_2*this->w, bf2 = 0.5*this->w;

	  this->p[6].x = this->p[18].x = this->p[9].y  = this->p[21].y =
	-(this->p[0].x = this->p[12].x = this->p[3].y  = this->p[15].y = this->w);
	  this->p[6].y = this->p[18].y = this->p[0].y  = this->p[12].y =
	  this->p[3].x = this->p[9].x  = this->p[15].x = this->p[21].x = 0;
	 this->p[12].z = this->p[18].z = this->p[15].z = this->p[21].z =
	-(this->p[0].z = this->p[6].z  = this->p[3].z  = this->p[9].z  = this->h);

	   this->p[1].x  = this->p[11].x = this->p[13].x = this->p[23].x =
	 -(this->p[5].x  = this->p[7].x  = this->p[17].x = this->p[19].x = -bf1);
	   this->p[1].y  = this->p[13].y = this->p[5].y  = this->p[17].y =
	 -(this->p[11].y = this->p[23].y = this->p[7].y  = this->p[19].y = -bf2);
	   this->p[1].z  = this->p[11].z = this->p[5].z  = this->p[7].z  =
	 -(this->p[13].z = this->p[23].z = this->p[17].z = this->p[19].z = -this->h);

	   this->p[2].x  = this->p[10].x = this->p[14].x = this->p[22].x =
	 -(this->p[4].x  = this->p[8].x  = this->p[16].x = this->p[20].x = -bf2);
	   this->p[2].y  = this->p[4].y  = this->p[14].y = this->p[16].y =
	 -(this->p[10].y = this->p[22].y = this->p[8].y  = this->p[20].y = -bf1);
	   this->p[2].z  = this->p[10].z = this->p[4].z  = this->p[8].z  =
	 -(this->p[14].z = this->p[22].z = this->p[16].z = this->p[20].z = -this->h);
}
//==============================================================================
void Dodecahedron::SetFacets(void)
{
	int* q;
	// the basal facets
	q = this->Gr[0];	q[0]=0; q[1]=1; q[2]=2; q[3]=3; q[4]=4; q[5]=5;
						q[6]=6; q[7]=7; q[8]=8; q[9]=9; q[10]=10; q[11]=11; q[12] = 0;

	q = this->Gr[13];	q[0]=12; q[1]=13; q[2]=14; q[3]=15; q[4]=16; q[5]=17; q[6]=18;
						q[7]=19; q[8]=20; q[9]=21; q[10]=22; q[11]=23; q[12] = 12;

	// the rectangular facets
	for(int i=1; i<=12; i++) {
		q=this->Gr[i]; q[0]=q[4]=i-1; q[1]=i+11; q[3]=i%12; q[2]=q[3]+12;
		q[5]=q[6]=q[7]=q[8]=q[9]=q[10]=q[11]=q[12]=-1;
	}
}
//============================================================================== 
