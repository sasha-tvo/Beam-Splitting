#include <QCoreApplication>
#include <QDir>
#include <QtCore/QtGlobal>

#include <conio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <list>
#include <time.h>
#include <math.h>

#include "..\Lib\trajectory.hpp"
#include "..\Lib\PhysMtr.hpp"
#include "..\Lib\Mueller.hpp"

#include "..\Lib\particle.hpp"



using namespace std;
//==============================================================================
Crystal* Body = NULL;						///< Crystal particle
unsigned int 	KoP,						///< Kind of particle
				AoP56,						///< Flag, 1 means angle of tip 56 deg, 0 means sizes of tip defined in data file
				GammaNumber,				///< Number of steps for Gamma rotation
				BettaNumber,				///< Number of steps for Betta rotation
				ThetaNumber,				///< Total number of bins in output file
				Itr,						///< Maximal number of internal reflection for every beam
				EDF;						///< EDF = 1 means to extract the delta-function in forward direction, EDF = 0 otherwise
complex 		_RefI(0,0);					///< Refraction index
unsigned int 	Sorting=0,					///< Sorting = 0 if all trajectories are taking into account, otherwise Sorting is equal to number of trajectories to calculate
				_NoF,						///< Number of facets of the crystal
				**Face;						///< An array for masking the trajectories which are out of interest
double			Radius,						///< Radius of the particle				
				Halh_Height,				///< Half height of the particle
				TipRadius,					///< Radius of the tip (if AoP56 = 0)
				TipHeight,					///< Height of the tip (if AoP56 = 0)
				SizeBin,					///< The size of the bin for Theta angle (radians)
				P = 0;						///< Probability distribution for Betta angle 
Arr2D 			mxd(0,0,0,0);				///< An array of output Mueller matrixes
matrix 			back(4,4),					///< Mueller matrix in backward direction 
				forw(4,4);					///< Mueller matrix in forward direction 
list<Chain> mask;							///< List of trajectories to take into account
Point3D			k(0,0,1),					///< Direction on incident wave
				Ey(0,1,0);					///< Basis for polarization characteristic of light
//==============================================================================
///< handler for the emitted beams
void  Handler(Beam& bm);
//==============================================================================

/// Reads the parameters from data file
int ReadFile(char* name, double* params, unsigned int n);
/// Fill in the \b mask and \b **Face from data file 
void MaskAppend(char s[], unsigned int n);
/// Shows the title
void ShowTitle(void);
/// Shows current time
void ShowCurrentTime(void);
/// Deletes \b **Face structure
void DelFace(void);


/// Main()
int main(int argc, char* argv[])
{
	QCoreApplication a(argc, argv); ///< QT needs it
	ShowTitle();
	cout << "\nLoading settings... ";

	// Go to current directory
	QString QCurDir;
	QCurDir=QCoreApplication::applicationDirPath();
	QDir dir;
	dir.cd(QCurDir);
	QDir::setCurrent(dir.absolutePath());

	unsigned int NumberOfParameters = 13;	///< Number of lines in data files, except trajectories
	double params[NumberOfParameters];		///< array of input data	
	// read parameters from data file
	try
	{
		if(ReadFile((char*)"params.dat", params, NumberOfParameters)) {
			cout << "\nError! Incorrect input file. Press any key for exit.";
			getch(); return 1;
		}
	}
	catch(const char* s) {
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}
	//---------------------------------------------------------------------------
	KoP = 			params[0];
	AoP56 = 		params[1];
	Radius = 		params[2];
	Halh_Height = 	params[3];
	TipRadius = 	params[4];
	TipHeight = 	params[5];
	_RefI =			complex(params[6],params[7]);
	GammaNumber = 	params[8];
	BettaNumber = 	params[9];
	ThetaNumber = 	params[10];
	Itr = 			params[11];
	EDF = 			params[12];
	//----------------------------------------------------------------------------
	double  
		NormAng = 			sqrt(3.0)/(2.0*tan(0.48869219055841228153863341517681)), 
		NormGammaAngle = 	0,
		NormBettaAngle =  	0;
		SizeBin = 			M_PI/ThetaNumber; // the size of the bin (radians)
	//----------------------------------------------------------------------------
	double hp;
	switch(KoP) { // choosing the kind of the particle
		case 0: // the hexagonal prizm
			Body= new Prism(_RefI, Radius, Halh_Height,Itr,k,Ey);
			NormGammaAngle = 	M_PI/(3.0*GammaNumber);
			NormBettaAngle =  	M_PI/(2.0*BettaNumber);
		break;
		case 1: // the hexagonal bullet
			hp = AoP56 ? NormAng*Radius : TipHeight;
			Body = new Bullet(_RefI, Radius, Halh_Height, hp,Itr,k,Ey);
			NormGammaAngle = 	M_PI/(3.0*GammaNumber);
			NormBettaAngle =  	M_PI/(BettaNumber);
		break;
		case 2: // the hexagonal pyramid
			hp = AoP56 ? NormAng*Radius : Halh_Height;
			Body = new Pyramid(_RefI, Radius, hp,Itr,k,Ey);
			NormGammaAngle = 	M_PI/(3.0*GammaNumber);
			NormBettaAngle =  	M_PI/(BettaNumber);			
		break;
		case 3: // the hexagonal tapered prizm
			hp = AoP56 ? NormAng*(Radius-TipRadius) : TipHeight;
			Body = new TaperedPrism(_RefI,
									Radius, Halh_Height, TipRadius, Halh_Height-hp,Itr,k,Ey);
			NormGammaAngle = 	M_PI/(3.0*GammaNumber);
			NormBettaAngle =  	M_PI/(BettaNumber);
		break;
		case 4: // cup
			hp = AoP56 ? NormAng*(Radius-TipRadius) : TipHeight;
			Body = new Cup(_RefI,
							Radius, Halh_Height, TipRadius, hp,Itr,k,Ey);
			NormGammaAngle = 	M_PI/(3.0*GammaNumber);
			NormBettaAngle =  	M_PI/(BettaNumber);
		break;
		
	}

	back.Fill(0); forw.Fill(0);
	mxd = Arr2D(1, ThetaNumber+1, 4, 4);
	mxd.ClearArr();

	Body->Phase() = false;
	const double NumOrient = GammaNumber*BettaNumber;
	//----------------------------------------------------------------------------
	#ifdef _Windows
		gotoxy(1,6);  ShowCurrentTime();   // current time (Win32 API)
	#endif
    clock_t t = clock();
    cout << BettaNumber;
	//----------------------------------------------------------------------------
	double s = 0, tetta, gamma;
	try {
		for(unsigned int i=0; i<BettaNumber; i++) {
			if(!i) cout << "\n" << 0;
			tetta = (i+0.5)*NormBettaAngle;
			for(unsigned int j=0; j<GammaNumber; j++) {
				gamma = (j+0.5)*NormGammaAngle;
				Body->ChangePosition(tetta, gamma,0.0);
				P = sin(tetta);
				s += P*Body->FTforConvexCrystal(Handler);
				if(!(j%100)) cout<<'.';
			}
			if(!KoP || KoP==3)
				cout << "\n" << (1-cos(tetta));
			else
				cout << "\n" << (1-cos(tetta))/2.0;
		}
	}
	catch(char* s) {
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}
	t = clock()-t;

	cout << "\nTotal time of calculation = " << t/CLK_TCK << " seconds";
	char text1[] = "M11 m22 m33 m44",
	text2[] = "tetta M11 m12 m21 m22 m33 m34 m43 m44";
	//----------------------------------------------------------------------------
	double b[3], f[3];
	b[0] = back[0][0]; b[1] = (back[1][1]-back[2][2])/2.; b[2] = back[3][3];
	f[0] = forw[0][0]; f[1] = (forw[1][1]+forw[2][2])/2.; f[2] = forw[3][3];
	//----------------------------------------------------------------------------
	double D_tot = b[0]+f[0];
	for(int j=0; j<=ThetaNumber; j++)
		D_tot += mxd(0,j,0,0);
	//----------------------------------------------------------------------------
	if(EDF) {
		ofstream bck("back.dat", ios::out), frw("forward.dat", ios::out);
		for(unsigned int j=0;j<3;j++) {
			f[j] /= D_tot; b[j] /= D_tot;
		}
		frw << text1; bck << text1;
		if(f[0]<=DBL_EPSILON)
			frw << "\n0 0 0 0";
		else
			frw << "\n" << f[0] << " " << f[1]/f[0] << " " <<  f[1]/f[0] << " " << f[2]/f[0];
		if(b[0]<=DBL_EPSILON)
			bck << "\n0 0 0 0";
		else
			bck << "\n" << b[0] << " " << b[1]/b[0] << " " << -b[1]/b[0] << " " << b[2]/b[0];
		bck.close();
		frw.close();
	}
	else {
		mxd(0,ThetaNumber,0,0) += f[0]; mxd(0,0,0,0) += b[0];
		mxd(0,ThetaNumber,1,1) += f[1]; mxd(0,0,1,1) += b[1];
		mxd(0,ThetaNumber,2,2) += f[1]; mxd(0,0,2,2) -= b[1];
		mxd(0,ThetaNumber,3,3) += f[2]; mxd(0,0,3,3) += b[2];
	}
	//----------------------------------------------------------------------------
	D_tot /= 2.0*NumOrient/M_PI;
	const double NRM = 4.0*NumOrient;
	ofstream M("M.dat", ios::out), out("out.dat", ios::out);
	M << text2;
	for(int j=ThetaNumber;j>=0;j--) {
		double sn;
		M << '\n' << 180.0/ThetaNumber*(ThetaNumber-j);
		sn = (j==0 || j==ThetaNumber) ? 1-cos(SizeBin/2.0) : (cos((j-0.5)*SizeBin)-cos((j+0.5)*SizeBin));
		matrix bf = mxd(0,j);
		if(bf[0][0] <= DBL_EPSILON)
			M << " 0 0 0 0 0 0 0 0";
		else
			M<<' '<<bf[0][0]/(NRM*sn*D_tot)<<' '
			<<bf[0][1]/bf[0][0]<<' '<<bf[1][0]/bf[0][0]<<' '<<bf[1][1]/bf[0][0]<<' '
			<<bf[2][2]/bf[0][0]<<' '<<bf[2][3]/bf[0][0]<<' '<<bf[3][2]/bf[0][0]<<' '<<bf[3][3]/bf[0][0];
	}
	M.close();
	//----------------------------------------------------------------------------
	// information for log-file
	out << "\nTotal time of calculation = " << t/CLK_TCK << " seconds";
	out << "\nTotal number of insident photons = " << s;
	out << "\nTotal number of body orientation = " << NumOrient;
	out << "\nMean incident energy per orientation = " << s/(NumOrient*m_2pi)*M_PI/2.0;
	out << "\nMean particle cross section = " << D_tot;
	out << "\nMean scattered energy per orientation = " << D_tot/m_2pi;
	out << "\nScattering albedo = " << D_tot*(2*NumOrient/M_PI)/s;
	out.close();
	//----------------------------------------------------------------------------
	// some information for user
	cout << "\nMean incident energy = " << s/(NumOrient*m_2pi)*M_PI/2.0
	<< "\nMean scattered energy = " << D_tot/m_2pi
	<< "\nScattering albedo = " << D_tot*(2.0*NumOrient/M_PI)/s
	<< "\nAll done. Please, press any key.";
	delete Body;
	if (Sorting)
		DelFace();
	return 0;
}
//==============================================================================


void  Handler(Beam& bm)
{
	double Area = P*CrossSection(bm);
	matrix bf = Mueller(bm());
	//----------------------------------------------------------------------------
	if(bm.r.z >= 1-DBL_EPSILON) {
		back += Area*bf;
	}
	else 
		if(bm.r.z <= DBL_EPSILON-1) {
			forw += Area*bf;
		}
		else {
			const unsigned int ZenAng = round(acos(bm.r.z)/SizeBin);
			double tmp = SQR(bm.r.y);
			if(tmp > DBL_EPSILON) {
				tmp = acos(bm.r.x/sqrt(SQR(bm.r.x)+tmp));
				if(bm.r.y<0) tmp = m_2pi-tmp;
					tmp *= 2;
				double cs = cos(tmp), sn = sin(tmp);
				RightRotateMueller(bf, cs, sn);
			}
			mxd.insert(0, ZenAng, Area*bf);
		}
}
//==============================================================================

const int size = 256;

int ReadFile(char* name, double* params, unsigned int n)
{
	char buf[size]=""; //буфер
	ifstream in(name, ios::in); //входной файл
	for(unsigned int i=0; i<n; i++) {
		if(in.eof()) return 1;
		in.getline(buf, size);
		params[i] = strtod(buf, NULL);
		if (params[i]<0) throw "All values must be positive in data file";
	}
	//----------------------------------------------------------------------------
	switch(int(params[0])) {
		case 0: _NoF =  8; break;
		case 1: _NoF = 13; break;
		case 2: _NoF =  7; break;
		case 3: _NoF = 20; break;
		case 4: _NoF = 14; break;
		case 5: _NoF = 6; break;
		case 6: _NoF = 4; break;
		case 7: _NoF = 5; break;
	}

	if(!in.eof()) {
	in.getline(buf, size);
	Sorting = strtod(buf, NULL);
	if(Sorting>0) {
		Face = new unsigned int*[int(Itr)+1];
		for(unsigned int i=0;i<Itr+1; i++)
			Face[i] = new unsigned int[_NoF];

		for(unsigned int i=0;i<Itr+1;i++)
			for(unsigned int j=0;j<_NoF;j++)
				Face[i][j] = 0;

		for(unsigned int j=0;j<Sorting;j++) {
			if(in.eof()) return 1;
			in.getline(buf, size);
			MaskAppend(buf,size);
		}
	}
	}
	in.close();
	return 0;
}


void MaskAppend(char s[], unsigned int n)
{
	list<unsigned int> ch;
	unsigned int intern_numb = 0; 
	char *buf,*end;
	end=s;
	do
	{
		buf=end;
		int facet_numb=strtol(buf,&end,10); 
		if (strlen(buf)!=strlen(end))
		{
			if ((facet_numb>_NoF)||(facet_numb<0)) throw "Error! Incorrect parameters of trajectories in data file";
			ch.push_back(facet_numb);
			Face[intern_numb][facet_numb] = 1; 
			intern_numb++; 
		}
	}
	while (strlen(buf)!=strlen(end));    
	if (ch.size()==0) throw "Error! There was not enough trajectories in data file";
	mask.push_back(ch); 
}
//==============================================================================

void DelFace(void)
{
    for(unsigned int i=0;i<_NoF; i++)
        delete[] Face[i];
    delete[] Face;
}

void ShowTitle(void)
 {
  cout << "*************************************************************\
         \n Light Scattering by Nonspherical Particles.                 \
         \n (c)Group of Wave Dispersion Theory,                    \
         \n    Institute of Atmospheric Optics RAS, Tomsk, Russia, 2014 \
         \n*************************************************************";
 }

//==============================================================================

void ShowCurrentTime(void)
 {
/*  using namespace std;
  SYSTEMTIME Time;
  GetLocalTime(&Time);
  cout << "Numerical calculations started at: " <<
          Time.wHour << ':' << Time.wMinute << ':' << Time.wSecond;*/
 }
