#include <QCoreApplication>
#include <QDir>
#include <QTime>
#include <QDate>
//#include <QtCore/QtGlobal>

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
				Def_Angle,					// Angle of deformation
				P = 0.0,					///< Probability distribution for Betta angle
				cos_angle,
				lm,
				T;							///< Max tilt angle
Arr2D			mxd(0,0,0,0);				///< An array of output Mueller matrixes

list<Chain>	mask;							///< List of trajectories to take into account
vector<Chain> Lbm;
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
void MaskAppend(char s[]);

/// Shows the title
void ShowTitle(void);

/// Deletes \b **Face structure
void DelFace(void);

//------------------------------------------------------------------------------
int NumMaxEnergy(vector<double> a)
{
	uint sz = a.size();
	if(sz==0)
		return -1; // список пустой
	int max_i = 0;
	double max_val = a[0];
	for(uint i=1; i<sz; i++)
	{
		if(a[i]>max_val)
		{
			max_i = i;
			max_val = a[i];
		}
	}
	return max_i;
}

int PoiskVzTrajectory(Chain A, vector<Chain> B)
{
	uint sza = A.Size(), szb = B.size();
	if(szb==0)
		return -1;  // список пустой
	int num = -1;   // номер в списке
	for(uint i=0; i<szb; i++)
	{
		if(sza!=B[i].Size())
			continue;
		else
		{
			list<uint>::const_iterator a = A.Begin(), b = B[i].End();
			b--;
			for(; a!=A.End() && (*a)==(*b); a++, b--);
			if(a==A.End())
			{
				num = i;
				break;
			}
		}
	}
	return num;
}

//------------------------------------------------------------------------------

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

	unsigned int NumberOfParameters = 14,	///< Number of lines in data files, except trajectories
				 GammaMaxAngle,			    // ћаксимальный угол поворота по углу гамма (дл€ хаотич. ориентации)
				 BettaMaxAngle;				// ћаксимальный угол поворота по углу бетта	(дл€ хаотич. ориентации)
	double params[NumberOfParameters];		///< array of input data
	// read parameters from data file
	try
	{
		if(ReadFile((char*)"params.dat", params, NumberOfParameters))
		{
			cout << "\nError! Incorrect input file. Press any key for exit. (1)";
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
	Def_Angle =		params[1];
	Radius =		params[2];
	Halh_Height =	params[3];
	_RefI =			complex(params[6],0.0);
	cos_angle =		cos(params[7]*M_PI/180.0);
	GammaNumber =	params[8];
	T =				params[9];
	BettaNumber =	params[10];
	lm =			params[11];
	F_Mt =			params[13]; // ѕараметр указывает записывать ли "картины" номеров колец и
								// дифракционного вклада (дл€ более детального исследовани€)

	double	NormGammaAngle =	0.0,
			NormBettaAngle =	0.0,
			Def_Angle_lim = 90.0,
			Rad = M_PI/180.0;

	switch(KoP) {
	  case 0: break;
	  case 1: Def_Angle_lim = atan(Radius/(2.0*Halh_Height))/Rad; break;
	  case 2: Def_Angle_lim = atan(Halh_Height/(2.0*Radius))/Rad; break;
	  case 3: Def_Angle_lim = atan(Halh_Height/(2.0*Radius))/Rad; break;
	  case 4: Def_Angle_lim = atan(Radius/(2.0*Halh_Height))/Rad; break;
	  default: {
		 cout << "\nError! Incorrect input file. Press any key for exit.(2)";
		 getch(); return 1;
	  }
	 }
	 if(Def_Angle>Def_Angle_lim) {
	  cout << "\nError! Incorrect input file. Press any key for exit.(3)";
	  getch(); return 1;
	 }

	if(fabs(Def_Angle)<FLT_EPSILON)
		KoP = 0;
	//----------------------------------------------------------------------------
	switch(KoP)	{ // choosing the kind of the particle
	 case 0: { // the perfect hexagonal prizm
		Body = new Prism(_RefI, Radius, Halh_Height, Itr, k, Ey);
		GammaMaxAngle = 60;
		BettaMaxAngle = 90;
		break;
	 }
	 case 1: {
		GammaMaxAngle = 60;
		BettaMaxAngle = 180;
		Body = new DeformatedPrism(_RefI, Radius, Halh_Height, KoP, Def_Angle, Itr, k, Ey);
		break;
	 }
	 case 2: {
	   GammaMaxAngle = 360;
	   BettaMaxAngle = 180;
	   Body = new DeformatedPrism(_RefI, Radius, Halh_Height, KoP, Def_Angle, Itr, k, Ey);
	   break;
	 }
	 case 3: {
	   GammaMaxAngle = 360;
	   BettaMaxAngle = 90;
	   Body = new DeformatedPrism(_RefI, Radius, Halh_Height, KoP, Def_Angle, Itr, k, Ey);
	   break;
	 }
	 case 4: {
	   GammaMaxAngle = 360;
	   BettaMaxAngle = 180;
	   Body = new DeformatedPrism(_RefI, Radius, Halh_Height, KoP, Def_Angle, Itr, k, Ey);
	   break;
	 }
   }

	NormGammaAngle = (double)GammaMaxAngle*Rad/(double)GammaNumber;
	NormBettaAngle = (double)BettaMaxAngle*Rad/(double)(BettaNumber-1);


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
	QDate date = QDate::currentDate();
	QTime time = QTime::currentTime();
	QString str = time.toString("hh_mm_ss");
	string time_start = "\nNumerical calculations started at: "+(time.toString("hh:mm:ss")).toStdString();
	cout << time_start;
	clock_t t = clock();

	QString fRes = "Data_"+date.toString("ddMMyy")+"_"+(time.toString("hhmmss"));
	dir.mkdir(fRes);
	dir.cd(fRes);
	QDir::setCurrent(dir.absolutePath());

	QString f_From = QCurDir+"/params.dat",
			f_To = "copy_params.dat";
	QFile::copy(f_From, f_To);

	double s = 0, betta, gamma;
	ofstream f("log.dat", ios::out);
	f << (QFileInfo(QCoreApplication::applicationFilePath()).fileName()).toStdString();
	f << time_start;
	f.close();

	uint gamma_25p = floor(0.25*GammaNumber);

	//--------------------------------------------------------------------------
	double p_Betta[BettaNumber], dcos[BettaNumber];
	if(!KoP || KoP==3) {
	 for (uint i=0; i<BettaNumber; i++)
	 {
		double betta = (double)i*NormBettaAngle;
		p_Betta[i] = 1.0;
		if (!i)
			dcos[i] = 1.0-cos(0.5*NormBettaAngle);
		else
			if(i==(BettaNumber-1))
				dcos[i] = cos(betta-0.5*NormBettaAngle)-cos(betta);
			else
				dcos[i] = cos(betta-0.5*NormBettaAngle)-cos(betta+0.5*NormBettaAngle);
	 }
	}
	else {
	 for (uint i=0; i<BettaNumber; i++)
	 {
	   double betta = (double)i*NormBettaAngle;
	   p_Betta[i] = 0.5;
	   if(!i || i==(BettaNumber-1))
		dcos[i] = 1.0-cos(0.5*NormBettaAngle);
	   else
		dcos[i] = cos(betta-0.5*NormBettaAngle)-cos(betta+0.5*NormBettaAngle);
	 }
	}

	if (fabs(T)>FLT_EPSILON)
	{
		NormBettaAngle = T/double(BettaNumber)*M_PI/180.0;
		double norm1 = 0.0, norm2 = 0.0, sigma = T/2.0*M_PI/180.;
		for (uint i=0; i<BettaNumber; i++)
		{
			double betta = (i+0.5)*NormBettaAngle, h;
			if (!i)
				dcos[i] = 1-cos(0.5*NormBettaAngle), h = 0.5*NormBettaAngle;
			else
				if(i==(BettaNumber-1))
					dcos[i] = cos(betta-0.5*NormBettaAngle)-cos(betta), h = 0.5*NormBettaAngle;
				else
					dcos[i] = cos(betta-0.5*NormBettaAngle)-cos(betta+0.5*NormBettaAngle), h = NormBettaAngle;
			p_Betta[i] = exp(-betta*betta/(2.0*sigma*sigma));
			norm1 += p_Betta[i]*dcos[i];
			norm2 += p_Betta[i]*h;
		}

		for (uint i=0; i<BettaNumber; i++)
			p_Betta[i] /= norm1;

		ofstream f("probability.dat", ios::out);
		f << "betta p_Betta dcos";
		for (uint i=0; i<BettaNumber; i++)
		{
			double betta = (i+0.5)*NormBettaAngle;
			f << endl << betta/M_PI*180. << " " << p_Betta[i] << " " << dcos[i];
		}
		f << "\nnorm = " << 1.0/(sqrt(2.0*M_PI)*sigma*180.0/M_PI)*norm2*180.0/M_PI;
		f.close();

	}
	//--------------------------------------------------------------------------
	try
	{
		for(betta_i=0; betta_i<BettaNumber; betta_i++)
		{
			if(!betta_i)
				cout << "\n" << 0;
			betta = (double)betta_i*NormBettaAngle;
			P = p_Betta[betta_i]*dcos[betta_i];
			for(gamma_j=0; gamma_j<GammaNumber; gamma_j++)
			{
				gamma = (gamma_j+0.5)*NormGammaAngle;
				Body->ChangePosition(betta, gamma, 0.0);
				s += P*Body->FTforConvexCrystal(Handler);
				if(!(gamma_j%gamma_25p)) cout<<'.';
			}
			unsigned int size = Lbm.size();
			ofstream f("log.dat", ios::app);
			f << endl << betta_i << " " << betta/Rad << " " << size;
			f.close();
			cout << endl << 100.0*(betta_i+1)/(double)BettaNumber;
		}
	}
	catch(char* s)
	{
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}

	const double NRM = GammaNumber;
	//----------------------------------------------------------------------------
	if(F_Mt)
	{
		en.clear();
		ofstream res("res.dat", ios::out);
		res << "trajectory\tenergy\n";
		vector<Chain>::const_iterator c = Lbm.begin();
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
			f1 << name1 << endl << BettaNumber << " " << GammaNumber << " " << params[9];
			f2 << name2 << endl << BettaNumber << " " << GammaNumber << " " << params[9];
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
			en.push_back(sum);
		}
		res.close();
	}
	else
	{
		ofstream res("res.dat", ios::out);
		res << "trajectory\tenergy\n";
		vector<Chain>::const_iterator c = Lbm.begin();
		for(unsigned int is=0; c!=Lbm.end(); c++, is++)
		{
			string tr = "";
			for(list<unsigned int>::const_iterator it = c->Ch.begin(); it!=c->Ch.end(); it++)
				tr += " "+to_string(*it);
			en[is] /= NRM;
			res << (tr+"//") << '\t' << en[is] << endl;
		}
		res.close();

	}
	//--------------------------------------------------------------------------
	ofstream sort_res("sort.dat", ios::out);
	sort_res << "trajectory\tenergy\tnum" << endl;
	int i = 0;
	do
	{
		int imax = NumMaxEnergy(en);
		if(imax==-1)
			break;
		Chain A = Lbm[imax];
		// записываем траекторию с max-вкладом в файл
		for(list<unsigned int>::const_iterator it = A.Begin(); it!=A.End(); it++)
			sort_res << " " << to_string(*it);
		sort_res << "//\t" << en[imax] << "\t" << i << endl;
		vector<Chain>::iterator _l = Lbm.begin();
		vector<double>::iterator _e = en.begin();
		advance(_l,imax); Lbm.erase(_l);   // удал€ем из списков траекторию с max-вкладом
		advance(_e,imax); en.erase(_e);
		int ivz = PoiskVzTrajectory(A, Lbm); // ищем номер взаимной траектории
		if(ivz>=0)
		{
			Chain B = Lbm[ivz];
			// записываем траекторию с max-вкладом в файл
			for(list<unsigned int>::const_iterator it = B.Begin(); it!=B.End(); it++)
				sort_res << " " << to_string(*it);
			sort_res << "//\t" << en[ivz] << "\t" << i << endl;
			_l = Lbm.begin();
			_e = en.begin();
			advance(_l,ivz); Lbm.erase(_l);   // удал€ем из списков взаимную траекторию
			advance(_e,ivz); en.erase(_e);
		}
		i++;
	}
	while(!Lbm.empty());
	sort_res.close();
	//--------------------------------------------------------------------------
	delete Body;
	t = (clock()-t)/CLK_TCK;

	ofstream g("log.dat", ios::app);
	g << "\nTotal time of calculation = " << t << " seconds";
	g << "\ns = " << s/NRM;
	g.close();


	cout << "\nTotal time of calculation = " << t << " seconds";
	cout << "\nAll done. Please, press any key.";
	return 0;
}
//==============================================================================


void Handler(Beam& bm)
{
	// отсеиваем траектории с направлением рассе€ни€ вне окрестности, заданной в файле параметров
	if((bm.r*k)<cos_angle)
		return;

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

	double M11 = 0.5*norm(bm());

	unsigned int nst = 0; // номер в списке
	vector<Chain>::iterator c = Lbm.begin();
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
			area = CrossSection(bm),
			x = M_PI*NumberOfRing,
			fn = M11*SQR(area)/lm/lm;

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

	_NoF =  8;

	Itr = params[12];

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
				MaskAppend(buf);
			}
		}
	}
	in.close();
	return 0;
}
//==============================================================================


void MaskAppend(char s[])
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
	for(unsigned int i=0;i<Itr+1; i++)
		delete[] Face[i];
	delete Face;
}
//==============================================================================


void ShowTitle(void)
{
	cout << "*************************************************************\
			\nLight Scattering by Nonspherical Particles.                 \
			\n(c)Group of Wave Dispersion Theory,                         \
			\n   Institute of Atmospheric Optics RAS, Tomsk, Russia, 2015 \
			\n*************************************************************\n";
}
//==============================================================================
