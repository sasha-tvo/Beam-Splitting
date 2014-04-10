#include <QCoreApplication>
#include <QDir>
#include <QTime>
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

Crystal* Body = NULL;						///< Crystal particle
unsigned int	KoP,						///< Kind of particle
				AoP56,						///< Flag, 1 means angle of tip 56 deg, 0 means sizes of tip defined in data file
				GammaNumber,				///< Number of steps for Gamma rotation
				BettaNumber,				///< Number of steps for Betta rotation
				Itr,						///< Maximal number of internal reflection for every beam
				F_Mt;
complex			_RefI(0,0);					///< Refraction index
unsigned int	Sorting=0,					///< Sorting = 0 if all trajectories are taking into account, otherwise Sorting is equal to number of trajectories to calculate
				_NoF,						///< Number of facets of the crystal
				**Face,						///< An array for masking the trajectories which are out of interest
				betta_i,
				gamma_j;
double			Radius,						///< Radius of the particle				
				Halh_Height,				///< Half height of the particle
				TipRadius,					///< Radius of the tip (if AoP56 = 0)
				TipHeight,					///< Height of the tip (if AoP56 = 0)
				P = 0.0,					///< Probability distribution for Betta angle
				cos_angle,
				lm;
Arr2D			mxd(0,0,0,0);				///< An array of output Mueller matrixes
matrix			back(4,4),					///< Mueller matrix in backward direction
				forw(4,4);					///< Mueller matrix in forward direction 
list<Chain>	mask,							///< List of trajectories to take into account
			Lbm;
vector<matrix>	nc,
				dif;
vector<double>	en;
Point3D			k(0,0,1),					///< Direction on incident wave
				Ey(0,1,0);					///< Basis for polarization characteristic of light

///< handler for the emitted beams
void Handler(Beam& bm);

/// Reads the parameters from data file
int ReadFile(char* name, double* params, unsigned int n);

/// Fill in the \b mask and \b **Face from data file 
void MaskAppend(char s[], unsigned int n);

/// Shows the title
void ShowTitle(void);

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
	// read parameters from data file
	try
	{
		if(ReadFile((char*)"params.dat", params, NumberOfParameters))
		{
			cout << "\nError! Incorrect input file. Press any key for exit.";
			getch(); return 1;
		}
	}
	catch(const char* s)
	{
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}
	//---------------------------------------------------------------------------
	KoP =			params[0];
	AoP56 =			params[1];
	Radius =		params[2];
	Halh_Height =	params[3];
	TipRadius =		params[4];
	TipHeight =		params[5];
	_RefI =			complex(params[6],0.0);
	cos_angle =		cos(params[7]*M_PI/180.0);
	GammaNumber =	params[8];
	BettaNumber =	params[9];
	lm =			params[10];
	F_Mt =			params[12];
	//----------------------------------------------------------------------------
	double
		NormAng =			sqrt(3.0)/(2.0*tan(0.48869219055841228153863341517681)),
		NormGammaAngle =	0.0,
		NormBettaAngle =	0.0;
	//----------------------------------------------------------------------------
	double hp;
	switch(KoP)
	{ // choosing the kind of the particle
		case 0: // the hexagonal prizm
			Body= new Prism(_RefI, Radius, Halh_Height, Itr, k, Ey);
			NormGammaAngle =	M_PI/(3.0*GammaNumber);
			NormBettaAngle =	M_PI/(2.0*BettaNumber);
		break;
		case 1: // the hexagonal bullet
			hp = AoP56 ? NormAng*Radius : TipHeight;
			Body = new Bullet(_RefI, Radius, Halh_Height, hp, Itr, k, Ey);
			NormGammaAngle =	M_PI/(3.0*GammaNumber);
			NormBettaAngle =	M_PI/(BettaNumber);
		break;
		case 2: // the hexagonal pyramid
			hp = AoP56 ? NormAng*Radius : Halh_Height;
			Body = new Pyramid(_RefI, Radius, hp, Itr, k, Ey);
			NormGammaAngle =	M_PI/(3.0*GammaNumber);
			NormBettaAngle =	M_PI/(BettaNumber);
		break;
		case 3: // the hexagonal tapered prizm
			hp = AoP56 ? NormAng*(Radius-TipRadius) : TipHeight;
			Body = new TaperedPrism(_RefI, Radius, Halh_Height, TipRadius,
									Halh_Height-hp, Itr, k, Ey);
			NormGammaAngle =	M_PI/(3.0*GammaNumber);
			NormBettaAngle =	M_PI/(BettaNumber);
		break;
		case 4: // cup
			hp = AoP56 ? NormAng*(Radius-TipRadius) : TipHeight;
			Body = new Cup(_RefI, Radius, Halh_Height, TipRadius, hp, Itr, k,
						   Ey);
			NormGammaAngle =	M_PI/(3.0*GammaNumber);
			NormBettaAngle =	M_PI/(BettaNumber);
		break;
	}

	nc.clear();
	dif.clear();
	en.clear();

	Body->Phase() = false;
	const double NumOrient = GammaNumber*BettaNumber;

	if(Sorting)
	{
		for(unsigned int i=0; i<Itr+1; i++)
			for(unsigned int j=0; j<_NoF; j++)
				Body->Facets[i][j]=Face[i][j];
		DelFace();
	}
	//----------------------------------------------------------------------------
	QTime time = QTime::currentTime();
	QString str = time.toString("hh_mm_ss");
	cout << "\nNumerical calculations started at: " << (time.toString("hh:mm:ss")).toStdString();
	clock_t t = clock();

	//----------------------------------------------------------------------------
	string fRes = "Data_"+(time.toString("hh_mm_ss")).toStdString();
	dir.mkdir(QString::fromStdString(fRes));
	dir.cd(QString::fromStdString(fRes));
	QDir::setCurrent(dir.absolutePath());

	double s = 0, betta, gamma;
	try
	{
		for(betta_i=0; betta_i<BettaNumber; betta_i++)
		{
			if(!betta_i)
				cout << "\n" << 0;
			betta = (betta_i+0.5)*NormBettaAngle;
			for(gamma_j=0; gamma_j<GammaNumber; gamma_j++)
			{
				gamma = (gamma_j+0.5)*NormGammaAngle;
				Body->ChangePosition(betta, gamma, 0.0);
				P = sin(betta);
				s += P*Body->FTforConvexCrystal(Handler);
				if(!(gamma_j%10)) cout<<'.';
			}
			unsigned int size = Lbm.size();
			ofstream f("log.dat", ios::app);
			f << endl << betta_i << " " << betta/M_PI*180.0 << " " << size;
			f.close();
			if(!KoP || KoP==3)
				cout << "\n" << (1-cos(betta)) << " (" << size << ") ";
			else
				cout << "\n" << (1-cos(betta))/2.0 << " (" << size << ") ";
		}
	}
	catch(char* s)
	{
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}

	const double NRM =2.0*NumOrient/M_PI;
	//----------------------------------------------------------------------------
	if(F_Mt)
	{
		ofstream res("res.dat", ios::out);
		res << "trajectory\tenergy\n";
		list<Chain>::const_iterator c = Lbm.begin();
		for(unsigned int is=0; c!=Lbm.end(); c++, is++)
		{
			double sum = 0.0;
			string _tr = "", tr = "", name1, name2;
			for(list<unsigned int>::const_iterator it = c->Ch.begin(); it!=c->Ch.end(); it++)
			{
				_tr += "_"+to_string(*it);
				tr += " "+to_string(*it);
			}
			name1 = "nr"+_tr+".dat";
			name2 = "dif"+_tr+".dat";
			ofstream f1(name1.c_str(), ios::out), f2(name2.c_str(), ios::out);
			f1 << name1 << endl << BettaNumber << " " << GammaNumber;
			f2 << name2 << endl << BettaNumber << " " << GammaNumber;
			for(unsigned int ib=0; ib<BettaNumber; ib++)
			{
				f1 << endl; f2 << endl;
				for(unsigned int ig=0; ig<GammaNumber; ig++)
				{
					double v_nc = (nc[is])[ib][ig], v_dif = (dif[is])[ib][ig]/NRM; // нормировка
					f1 << " " << v_nc; f2 << " " << v_dif;
					sum += v_dif;
				}
			}
			f1.close(); f2.close();
			res << (tr+"//") << '\t' << sum << endl;
		}
		res.close();
	}
	else
	{
		ofstream res("res.dat", ios::out);
		res << "trajectory\tenergy\n";
		list<Chain>::const_iterator c = Lbm.begin();
		for(unsigned int is=0; c!=Lbm.end(); c++, is++)
		{
			string tr = "";
			for(list<unsigned int>::const_iterator it = c->Ch.begin(); it!=c->Ch.end(); it++)
				tr += " "+to_string(*it);
			double v_dif = en[is]/NRM;
			res << (tr+"//") << '\t' << v_dif << endl;
		}
		res.close();
	}
	//--------------------------------------------------------------------------
	delete Body;
	t = clock()-t;
	cout << "\nTotal time of calculation = " << t/CLK_TCK << " seconds";
	cout << "\nAll done. Please, press any key.";
	return 0;
}
//==============================================================================


void Handler(Beam& bm)
{
	list<unsigned int> facets;
	for(list<unsigned int>::const_iterator b = bm.BeginP(); b!=bm.EndP(); b++)
		facets.push_back(*b);

	if(Sorting)
	{
		bool flag = false;
		for(list<Chain>::const_iterator c = mask.begin(); c!=mask.end(); c++)
		{
			list<unsigned int>::const_iterator it = c->Ch.begin(), fs = facets.begin();
			if(SizeP(bm)!=c->Ch.size())
				continue;
			for(;it!=c->Ch.end() && (*it)==(*fs); it++, fs++);
			if(it==c->Ch.end())
			{
				flag = true;
				break;
			}
		}
		if(!flag) return;
	}

	if((bm.r*k)<cos_angle)
	{
		matrix m = Mueller(bm());

		unsigned int nst = 0; // номер в списке
		list<Chain>::const_iterator c = Lbm.begin();
		for(; c!=Lbm.end(); c++, nst++)
		{
			list<unsigned int>::const_iterator it = c->Ch.begin(), fs = facets.begin();
			if(SizeP(bm)!=c->Ch.size())
				continue;
			for(; it!=c->Ch.end() && (*it)==(*fs); it++, fs++);
			if(it==c->Ch.end())
				break;
		}
		if(c==Lbm.end())
		{
			Lbm.push_back(facets);
			if(F_Mt)
			{
				matrix temp(BettaNumber, GammaNumber);
				temp.Fill(0.0);
				nc.push_back(temp);
				dif.push_back(temp);
			}
			else
				en.push_back(0.0);
		}

		std::list<Point3D>::const_iterator pt = bm.Begin();
		Point3D rn = k,
				rn_pr = rn-bm.r*(bm.r*rn),		// проекци€ на плоскость перпендикул€рную пучку.
				cur = (*pt)-bm.r*(bm.r*(*pt));	// проецируем первую вершину на плоскость перпендикул€рную пучку,

		double	curr = cur*rn_pr,				// вычисл€ем скол€рное произведение проекции нормали и проекции точки
				min_ro_n = curr,
				max_ro_n = curr;				// полагаем максимуму и минимум равным первому рассто€нию

		for(pt++; pt!=bm.End(); pt++)
		{
			Point3D cur=(*pt)-bm.r*(bm.r*(*pt));	// проецируем вершину на плоскость перпендикул€рную пучку,
			double curr=cur*rn_pr;					// вычисл€ем скал€рное произведение проекции нормали и проекции точки
			if(curr>max_ro_n)
				max_ro_n=curr;
			if(curr<min_ro_n)
				min_ro_n=curr;
		}
		max_ro_n = max_ro_n-min_ro_n;
		double	NumberOfRing = max_ro_n/lm,
				area = AreaOfBeam(bm),
				x = M_PI*NumberOfRing,
				fn = m[0][0]*SQR(area)/lm/lm;

		if(x>2.2)
			fn *= 2.6/x/x/x;
		else
			fn *= 0.0688734*x*x*x-0.3142876*x*x+0.0192811*x+0.9998942;

		if(F_Mt)
		{
			(nc[nst])[betta_i][gamma_j] = NumberOfRing;
			(dif[nst])[betta_i][gamma_j] = fn*P;
		}
		else
			en[nst] += fn*P;
	}
}
//==============================================================================


int ReadFile(char* name, double* params, unsigned int n)
{
	const int size = 256;
	char buf[size]=""; //буфер
	ifstream in(name, ios::in); //входной файл
	for(unsigned int i=0; i<n; i++)
	{
		if(in.eof())
			return 1;
		in.getline(buf, size);
		params[i] = strtod(buf, NULL);
		if(params[i]<0)
			throw "All values must be positive in data file";
	}

	switch(int(params[0]))
	{
		case 0: _NoF =  8; break;
		case 1: _NoF = 13; break;
		case 2: _NoF =  7; break;
		case 3: _NoF = 20; break;
		case 4: _NoF = 14; break;
		case 5: _NoF = 6;  break;
		case 6: _NoF = 4;  break;
		case 7: _NoF = 5;  break;
	}

	Itr = params[11];

	if(!in.eof())
	{
		in.getline(buf, size);
		Sorting = strtod(buf, NULL);
		if(Sorting>0)
		{
			Face = new unsigned int*[int(Itr)+1];
			for(unsigned int i=0; i<Itr+1; i++)
				Face[i] = new unsigned int[_NoF];

			for(unsigned int i=0; i<Itr+1; i++)
				for(unsigned int j=0; j<_NoF; j++)
					Face[i][j] = 0;

			for(unsigned int j=0; j<Sorting; j++)
			{
				if(in.eof())
					return 1;
				in.getline(buf, size);
				MaskAppend(buf,size);
			}
		}
	}
	in.close();
	return 0;
}
//==============================================================================


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
		if(strlen(buf)!=strlen(end))
		{
			if((facet_numb>_NoF)||(facet_numb<0))
				throw "Error! Incorrect parameters of trajectories in data file";
			ch.push_back(facet_numb);
			Face[intern_numb][facet_numb] = 1;
			intern_numb++;
		}
	}
	while(strlen(buf)!=strlen(end));
	if(ch.size()==0)
		throw "Error! There was not enough trajectories in data file";
	mask.push_back(ch); 
}
//==============================================================================


void DelFace(void)
{
	for(unsigned int i=0;i<_NoF; i++)
		delete[] Face[i];
	delete[] Face;
}
//==============================================================================


void ShowTitle(void)
{
	cout << "*************************************************************\
			\nLight Scattering by Nonspherical Particles.                 \
			\n(c)Group of Wave Dispersion Theory,                         \
			\n   Institute of Atmospheric Optics RAS, Tomsk, Russia, 2014 \
			\n*************************************************************\n";
}
//==============================================================================
