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

//Axillary functions
#include "..\Lib\trajectory.hpp"
#include "..\Lib\PhysMtr.hpp"
#include "..\Lib\Mueller.hpp"


//The beam-splitting
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
				EDF,						///< EDF = 1 means to extract the delta-function in forward direction, EDF = 0 otherwise
				isNotRandom;				///< if 1 particle rotates by discret steps, if 0 particle rotates with random function
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
///< Handler for the emitted beams
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
	srand (time(NULL));
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
	// Read parameters from data file
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
	//Assign the parameters from data file.
	//---------------------------------------------------------------------------
	KoP = 			params[0];
	AoP56 = 		params[1];
	Radius = 		params[2];
	Halh_Height = 	params[3];
	TipRadius = 	params[4];
	TipHeight = 	params[5];
	_RefI =			complex(params[6],0);
	GammaNumber = 	params[7];
	BettaNumber = 	params[8];
	ThetaNumber = 	params[9];
	Itr = 			params[10];
	EDF = 			params[11];
	isNotRandom =	params[12];

	//----------------------------------------------------------------------------
	double  
		NormAng = 			sqrt(3.0)/(2.0*tan(0.48869219055841228153863341517681)), 
		NormGammaAngle =	0,
		NormBettaAngle =	0;
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

	// the arrays for exact backscattering and forwardscattering Mueller matrices
	back.Fill(0); forw.Fill(0);

	// the array for the Mueller matrices dependent from zenith angle
	mxd = Arr2D(1, ThetaNumber+1, 4, 4);
	mxd.ClearArr();

	Body->Phase() = false;
	const double NumOrient = GammaNumber*BettaNumber;


	// The main loop
	//----------------------------------------------------------------------------

	clock_t t = clock();
	cout << "\nPercent done";
	double s = 0, tetta, gamma;
	try {
		for(unsigned int i=0; i<BettaNumber; i++) {
			if(!i) cout << "\n" << 0;
			if (isNotRandom)
				tetta = (i+0.5)*NormBettaAngle;
			for(unsigned int j=0; j<GammaNumber; j++) {
				if (isNotRandom)
					gamma = (j+0.5)*NormGammaAngle;
				else
				{
					gamma=2.0*M_PI*((double)(rand()) / ((double)RAND_MAX+1.0));
					tetta=acos(((double)(rand()) / ((double)RAND_MAX))*2.0-1.0);
				}
				Body->ChangePosition(tetta, gamma,0.0);
				//P = sin(tetta);
				P=1.0;
				s += P*Body->FTforConvexCrystal(Handler);
				if(!(j%100)) cout<<'.';
			}
			if(!KoP || KoP==3)
				cout << "\n" << (i+0.5)*NormBettaAngle/M_PI*200;
			else
				cout << "\n" << (i+0.5)*NormBettaAngle/M_PI*100;
		}
	}
	catch(char* s) {
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}
	t = clock()-t;
	// End of the main loop
	//----------------------------------------------------------------------------

	cout << "\nTotal time of calculation = " << t/CLK_TCK << " seconds";

	//Analytical averaging over alpha angle
	//----------------------------------------------------------------------------
	double b[3], f[3];	
	b[0] = back[0][0]; b[1] = (back[1][1]-back[2][2])/2.; b[2] = back[3][3];
	f[0] = forw[0][0]; f[1] = (forw[1][1]+forw[2][2])/2.; f[2] = forw[3][3];

	//Integrating
	//----------------------------------------------------------------------------
	double D_tot = b[0]+f[0];
	for(int j=0; j<=ThetaNumber; j++)
		D_tot += mxd(0,j,0,0);	

	//Normalizing coefficient
	const double NRM = 1.0/(4.0*M_PI);

	// Extracting the forward and backward peak in a separate file if needed
	//----------------------------------------------------------------------------
	if(EDF) {
		ofstream bck("back.dat", ios::out), frw("forward.dat", ios::out);
		frw << "M11 m22 m33 m44";
		bck << "M11 m22 m33 m44";
		if(f[0]<=DBL_EPSILON)
			frw << "\n0 0 0 0";
		else
			frw << "\n" << f[0]/(NRM*(1-cos(SizeBin/2.0))*D_tot) << " " << f[1]/f[0] << " " <<  f[1]/f[0] << " " << f[2]/f[0];
		if(b[0]<=DBL_EPSILON)
			bck << "\n0 0 0 0";
		else
			bck << "\n" << b[0]/(NRM*(1-cos(SizeBin/2.0))*D_tot) << " " << b[1]/b[0] << " " << -b[1]/b[0] << " " << b[2]/b[0];
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

	//Output the matrices
	ofstream M("M.dat", ios::out), out("out.dat", ios::out);
	M <<  "tetta M11 m12 m21 m22 m33 m34 m43 m44";
	for(int j=ThetaNumber;j>=0;j--) {
		double sn;
		//Special case in first and last step
		M << '\n' << 180.0/ThetaNumber*(ThetaNumber-j)+(j==0 ?-0.25*180.0/ThetaNumber:0)+(j==ThetaNumber ?0.25*180.0/ThetaNumber:0);
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
	// Information for log-file
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

//The function that takes outgoing beams into account
void  Handler(Beam& bm)
{
	uint  pn = 0, szP = SizeP(bm);
	uint vi=0;
	// If there is a special list of the beams to take into account, check if the beam "bm" is on the list
	if (Sorting)
	{		
		bool flag = false;
		list<Chain>::const_iterator c = mask.begin();
		for (;c!=mask.end();c++)
		{
			list<uint >::const_iterator it = c->Ch.begin();
			if (szP!=c->sz) continue;
			list<uint >::const_iterator fs = bm.BeginP();
			for (;it!=c->Ch.end() && (*it)==(*fs);it++, fs++);
			if (it==c->Ch.end())
			{
				flag = true;
				break;
			}
		}
		if (!flag) return;
	}


	double Area = P*CrossSection(bm);
	matrix bf = Mueller(bm());
	//----------------------------------------------------------------------------
	// Collect the beam in array
	const unsigned int ZenAng = round(acos(bm.r.z)/SizeBin);
	if (ZenAng==0)
	{
		back += Area*bf;
	}
	else 		
		if (ZenAng==ThetaNumber)
		{
			forw += Area*bf;
		}
		else {
			// Rotate the Mueller matrix of the beam to appropriate coordinate system

			double tmp = SQR(bm.r.y);
			if(tmp > DBL_EPSILON) {
				tmp = acos(bm.r.x/sqrt(SQR(bm.r.x)+tmp));
				if(bm.r.y<0) tmp = m_2pi-tmp;
				tmp *= -2.0;
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
	char buf[size]=""; //temp buffer
	ifstream in(name, ios::in);
	for(unsigned int i=0; i<n; i++) {
		if(in.eof()) return 1;
		in.getline(buf, size);
		params[i] = strtod(buf, NULL);
		if (params[i]<0) throw "All values must be positive in data file";
	}
	//Choose the _NoF (Number of facets) according to the chosen kind of particle
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
	Itr = params[10];

	if(!in.eof()) {
	in.getline(buf, size);
	Sorting = strtod(buf, NULL);
	// If there is a special list of beams
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
			MaskAppend(buf,size); //Add the beam into special list
		}
	}
	}
	in.close();
	return 0;
}


//This function interprets the text line into beam name and bushed this name into special list
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
	for(unsigned int i=0;i<Itr+1; i++)
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


