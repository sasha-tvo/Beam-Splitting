#include <QCoreApplication>
#include <QDir>

//---------------------------------------------------------------------------
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
//typedef unsigned int uint;

using namespace std;
//==============================================================================
Crystal* Body = NULL;						///< Crystal particle
uint			KoP,						///< Kind of particle
				AoP56,						///< Flag, 1 means angle of tip 56 deg, 0 means sizes of tip defined in data file
				ThetaNumber,				///< Total number of steps by Theta in output file
				PhiNumber,					///< Total number of steps by Phi in output file
				Itr,						///< Maximal number of internal reflection for every beam
				EDF;						///< EDF = 1 means to extract the delta-function in forward direction, EDF = 0 otherwise
int				GammaNumber,				///< Number of steps for Gamma rotation
				BettaNumber;				///< Number of steps for Betta rotation
complex 		_RefI(0,0);					///< Refraction index
uint			NumberOfTrajectory=0,		///< NumberOfTrajectory = 0 if all trajectories are taking into account, otherwise NumberOfTrajectory is equal to number of trajectories to calculate
				_NoF,						///< Number of facets of the crystal
				**Face;						///< An array for masking the trajectories which are out of interest
double			Radius,						///< Radius of the particle
				Halh_Height,				///< Half height of the particle
				TipRadius,					///< Radius of the tip (if AoP56 = 0)
				TipHeight,					///< Height of the tip (if AoP56 = 0)
				SizeBin,					///< The size of the bin for Theta angle (radians)
				betta_min_Rad,				///< Initial Betta angle
				betta_max_Rad,				///< Final Betta angle
				gamma_cnt,					///< Basis point to gamma rotating
				GammaLimmExtendedCoeff,		///< Coeffitient by which Gamma Limm will be extended
				lm,							///< Wave Length
				ConusGrad,					///< Conus in backscattering (grad)
				ConusRad,					///< Conus in backscattering (rad)
				P = 0;						///< Probability distribution for Betta angle
Arr2D 			mxd(0,0,0,0);				///< An array of output Mueller matrixes
matrix 			back(4,4),					///< Mueller matrix in backward direction
				forw(4,4);					///< Mueller matrix in forward direction
list<Chain>		mask;						///< List of trajectories to take into account
list<Chain>		AllowedMask;				///< Mask of allowed trajectories, noly this trajectories can be in mask
Point3D			k(0,0,1),					///< Direction on incident wave
				Ey(0,1,0);					///< Basis for polarization characteristic of light
bool			perpend_diff;				///< If true - diffraction will be calculated with Shifted screen, perpendicular to propagation direction; if false - diffraction will be calculated with incline screen
//==============================================================================
const uint		NumberOfParameters = 16;		///< Number of lines in data files, except trajectories
double			params[NumberOfParameters];		///< Array of input data
double			NormGammaAngle = 	0,			///< Normalize coefficient for Gamma
				NormBettaAngle =  	0;			///< Normalize coefficient for Betta
const double	Rad = M_PI/180.;				///< Grad to Rad coefficient
double			df = 0.0,						///< dPhi
				dt = 0.0;						///< dTheta
const uint		AllowedTrajectoryNumber = 28,	///< Number of trajectories that can be calculated
				AllowedMuellerMatrixNumber = 16;	///< Maximal allowed number of muller matrixes

/// The list of trajectories that can be calculated
const string	allowed_beams_mask[AllowedTrajectoryNumber]
												= {"0", "0 7 0", "3", "3 6 3", "2 6 7 4", "2 7 6 4", "4 6 7 2", "4 7 6 2", "3 6 7 3", "3 7 6 3", "0 6 7 0",
												"0 7 6 0", "3 7 0 6 7 3", "3 7 6 0 7 3", "0 6 3 7 6 0", "0 6 7 3 6 0", "2 6 7 3", "2 7 6 3", "4 6 7 3",
												"4 7 6 3", "2 1 6 7 4", "4 7 6 1 2", "2 7 6 5 4", "4 5 6 7 2", "2 7 0 6 7 4", "2 7 6 0 7 4", "4 7 0 6 7 2", "4 7 6 0 7 2"};
/// The mask of allowed muller matrixes
int				Mueller_mask[AllowedMuellerMatrixNumber]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
uint			n_tr,							// номер траектории, для которой вычисляется gamma_lim
				Lim;							// номер текущей области из массива sort_lim
double			NumberOfRing,					// номер дифракционного кольца для пучка с номером n_tr
				bm_tetta,						// tetta_координата для пучка с номером n_tr
				lim[AllowedMuellerMatrixNumber]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}, sort_lim[AllowedMuellerMatrixNumber];
vector<Arr2DC>	Jones_temp;


///< handler for the emitted beams
void Handler(Beam& bm);

void HandlerTemp(Beam& bm);
//==============================================================================

/// Reads the parameters from data file
int ReadFile(char* name, double* params, uint n);

/// Fill in the \b mask and \b **Face from data file
void MaskAppend(const char *s);

/// Fill in the \b AllowedMask from list of trajectories we can take into account
void AllowedMaskAppend(const char *s);

/// Shows the title
void ShowTitle(void);

/// Shows current time
void ShowCurrentTime(void);

/// Deletes \b **Face structure
void DelFace(void);

/// Find out the GammaLim for every trajectory we take into account
uint FillLim(double ConusRad, double betta, double sc_coef, const char *name);


//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	srand (time(NULL));
	QCoreApplication a(argc, argv);
	ShowTitle();
	cout << "\nLoading settings... ";

	QString QCurDir;
	QCurDir=QCoreApplication::applicationDirPath();
	QDir dir;
	dir.cd(QCurDir);
	QDir::setCurrent(dir.absolutePath());

	try
	{
		if (ReadFile((char*)"params.dat", params, NumberOfParameters))
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

	KoP =					0; //Only hexagonal
	AoP56 =					0; // No Tip
	Halh_Height =			params[0]/2.0;
	Radius =				params[1]/2.0;
	_RefI =					complex(params[2],0);
	betta_min_Rad =			params[3]*Rad;
	betta_max_Rad =			params[4]*Rad;
	BettaNumber =			params[5];
	GammaLimmExtendedCoeff = params[6];
	gamma_cnt =				params[7]*Rad;
	GammaNumber =			params[8];	
	lm =					params[11];
	ConusGrad =				params[12];
	ThetaNumber =			params[13];
	PhiNumber =				params[14];
	perpend_diff =			(params[15]>0.5?true:false);

	if (Halh_Height<=0 || Radius<=0)
	{
		cout << "\nError! Incorrect input data. Press any key for exit.";
		getch(); return 1;
	}

	cout << "OK";
	//----------------------------------------------------------------------------

	#ifdef _Windows
		gotoxy(1,7);  ShowCurrentTime();   // current time (Win32 API)
	#endif
	clock_t tm = clock();

	switch(KoP)
	{ // choosing the kind of the particle
		case 0: // the hexagonal prizm
			Body= new Prism(_RefI, Radius, Halh_Height,Itr,k,Ey);
			NormGammaAngle = 	M_PI/(3.0*GammaNumber);
			NormBettaAngle =  	M_PI/(2.0*BettaNumber);
		break;
		default:
			cout << "Wring kind of particle!";
			return 1;
	}
	//----------------------------------------------------------------------------
	Body->Phase() = true;

	double 	dBettaRad = 0.0, dGammaRad = 1.0, ConusRad = ConusGrad*Rad;//, Tgamma;
	if (BettaNumber) dBettaRad = (betta_max_Rad-betta_min_Rad)/(double)BettaNumber;
	if (ThetaNumber) dt = ConusRad/(double)ThetaNumber;
	if (PhiNumber) df = m_2pi/(double)(PhiNumber+1);

	string  folder_name = "Data", fRes = folder_name,
	str = to_string(ConusGrad)+" "+to_string(PhiNumber)+" "+to_string(ThetaNumber);

	dir.mkdir(QString::fromStdString(fRes));
	dir.cd(QString::fromStdString(fRes));
	QDir::setCurrent(dir.absolutePath());

	Arr2D M_(PhiNumber+1,ThetaNumber+1,4,4);
	M_.ClearArr();



	vector<Arr2D> M;
	M.clear();
	for (uint q=0;q<AllowedMuellerMatrixNumber; q++)
		M.push_back(M_);

	uint orn = 0;
	double dcos_sum = 0.0;
	cout <<endl;
	//----------------------------------------------------------------------------
	AllowedMask.clear();
	for (uint i=0; i<AllowedTrajectoryNumber; i++)
		AllowedMaskAppend(allowed_beams_mask[i].c_str());

	if (NumberOfTrajectory)
	{
		//Run over all trajectories in params.dat
		list<Chain>::const_iterator cm = mask.begin();
		for (;cm!=mask.end();cm++)
		{
			uint  pn = 0, sm = cm->Size();
			//Run over list of allowed trajectories
			for (list<Chain>::const_iterator ct = AllowedMask.begin();ct!=AllowedMask.end();ct++, pn++)
			{
				uint  st = ct->Size();
				// check size
				if (sm!=st) continue;

				list<unsigned int>::const_iterator im = cm->Ch.begin(), it = ct->Ch.begin();
				// run over all equal elements
				for (;im!=cm->Ch.end() && (*im)==(*it);im++, it++);
				// if ran all elements, trajectories are equal
				if (im==cm->Ch.end())
				{
					uint  vi;
					string tr;
					if (!pn) { vi = 0; tr = "0"; }
					else if (pn==1) { vi = 1; tr = "070"; }
					else if (pn==2) { vi = 2; tr = "3"; }
					else if (pn==3) { vi = 3; tr = "363"; }
					else if (pn>3 && pn<8) { vi = 4; tr = "2674_2764_4672_4762"; }
					else if (pn==8 || pn==9) { vi = 5; tr = "3673_3763"; }
					else if (pn==10 || pn==11) { vi = 6; tr = "0670_0760"; }
					else if (pn==12 || pn==13) { vi = 7; tr = "370673_376073"; }
					else if (pn==14 || pn==15) { vi = 8; tr = "063760_067360"; }
					else if (pn==16 || pn==17) { vi = 9; tr = "2673_2763"; }
					else if (pn==18 || pn==19) { vi = 10; tr = "4673_4763"; }
					else if (pn==20) { vi =  11; tr = "21674"; }
					else if (pn==21) { vi =  12; tr = "47612"; }
					else if (pn==22) { vi =  13; tr = "27654"; }
					else if (pn==23) { vi =  14; tr = "45672"; }
					else { vi =  15;  tr = "2674XL"; }
					Mueller_mask[vi] = 1;
					ofstream out(("ga_"+tr+".dat").c_str(), ios::out);
					out << to_string(BettaNumber) << " " << to_string(ThetaNumber) << " " << to_string(dBettaRad/Rad) << " " << to_string(dt/Rad);
					out.close();
					break;
				}
			}
		}
	}
	else
	{
		Mueller_mask[0] = 1;
		ofstream out("ga_all.dat", ios::out);
		out << to_string(BettaNumber) << " " << to_string(ThetaNumber) << " " << to_string(dBettaRad/Rad) << " " << to_string(dt/Rad);
		out.close();
	}


	//----------------------------------------------------------------------------
	ofstream out("out.dat", ios::out);
	out.close();	
	double dgamma_sum = 0.0, norm = (GammaNumber>0 ? 1.0/(M_PI/3.0) : 1.0);
	try
	{
		for (int Betta_i=0; Betta_i<=BettaNumber; Betta_i++)
		{
			if (!Betta_i)
				cout << endl << "0% ";
			double bettaRad = betta_min_Rad+(double)Betta_i*dBettaRad,
				dcos;
			if (BettaNumber)
			{
				if (!Betta_i)
					dcos = (cos(bettaRad)-cos(bettaRad+dBettaRad))/2.0;
				else if (Betta_i==BettaNumber)
					dcos = (cos(bettaRad-dBettaRad)-cos(bettaRad))/2.0;
				else
					dcos = (cos(bettaRad-dBettaRad)-cos(bettaRad+dBettaRad))/2.0;
			}
			else
				dcos = 1.0;


			dcos_sum += dcos;
			uint  NumLim=1;
			sort_lim[0] = M_PI;


			if (NumberOfTrajectory)
				NumLim = FillLim(ConusRad, bettaRad, GammaLimmExtendedCoeff, "out.dat");

			if (!NumLim) continue;
			dgamma_sum = 0.0;
			//-------------------------------------------------------------------------
			if (!GammaNumber)
				NumLim = 1;
			for (Lim=0;Lim<NumLim;Lim++)
			{
				cout << ". ";
				if (GammaNumber)
					dGammaRad = sort_lim[Lim]/(double)GammaNumber;
				int Gamma_Number_to_skip=0;
				if (Lim)
					Gamma_Number_to_skip = 1+floor(sort_lim[Lim-1]/dGammaRad); // GammaNumber - число "внутренних" шагов в цикле, которые нужно пропустить
				for (int Gamma_j=-GammaNumber; Gamma_j<=GammaNumber; Gamma_j++)
				{
					if (Lim && abs(Gamma_j)<Gamma_Number_to_skip) continue; // пропускаем внутренние шаги
					double gamma = gamma_cnt+Gamma_j*dGammaRad,
						   Pgamma = dGammaRad;
					if (Lim && abs(Gamma_j) == Gamma_Number_to_skip)
						Pgamma = ((double)Gamma_Number_to_skip+0.5)*dGammaRad-sort_lim[Lim-1];
					if (GammaNumber && (abs(Gamma_j) == GammaNumber))
						Pgamma -= 0.5*dGammaRad;
					dgamma_sum += Pgamma;
					Body->ChangePosition(bettaRad, gamma, 0);
					Jones_temp.clear();
					Arr2DC tmp_(PhiNumber+1,ThetaNumber+1,2,2);
					tmp_.ClearArr();
					for (uint q=0;q<AllowedMuellerMatrixNumber; q++)
						Jones_temp.push_back(tmp_);
					Body->FTforConvexCrystal(Handler);
					for (uint q=0; q<AllowedMuellerMatrixNumber; q++)
					{
						if (Mueller_mask[q]==1)
						{
							for (uint  j_tt=0; j_tt<=ThetaNumber; j_tt++)
								for (uint  i_fi=0; i_fi<=PhiNumber; i_fi++)
								{
									matrix m_tmp = Mueller(Jones_temp[q](i_fi,j_tt));
									M[q].insert(i_fi,j_tt,dcos*Pgamma*norm*m_tmp);
								}
						}
					}
				}
				//cout << endl << to_string(Lim);
			}
			orn++;
			cout <<"\r                        \r";
			cout << (orn*100./(BettaNumber+1))<<"%"<<" ";
		}
	}
	catch(char* s)
	{
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
	}

	// вывод в файл
	string name_f;
	if (GammaNumber)
		name_f = "random_";
	else
		name_f = "b_"+to_string(betta_min_Rad/Rad)+"_";

	if (NumberOfTrajectory)
		for (uint q=0; q<AllowedMuellerMatrixNumber; q++)
	{
		string tr;
		if (!q) tr = "0";
		else if (q==1) tr = "070";
		else if (q==2) tr = "3";
		else if (q==3) tr = "363";
		else if (q==4) tr = "2674_2764_4672_4762";
		else if (q==5) tr = "3673_3763";
		else if (q==6) tr = "0670_0760";
		else if (q==7) tr = "370673_376073";
		else if (q==8) tr = "063760_067360";
		else if (q==9) tr = "2673_2763";
		else if (q==10) tr = "4673_4763";
		else if (q==11) tr = "21674";
		else if (q==12) tr = "47612";
		else if (q==13) tr = "27654";
		else if (q==14) tr = "45672";
		else if (q==15) tr = "2674XL";
		if (Mueller_mask[q]==1)
		{

			ofstream f((name_f+"_"+tr+".dat").c_str(), ios::out),
					res(("ga_"+tr+".dat").c_str(), ios::app);
			f << (to_string(ConusGrad)+" "+to_string(PhiNumber)+" "+to_string(ThetaNumber)).c_str();
			res << endl<< "tetta M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44";
			f.precision(10);
			res.precision(10);


			matrix sum(4,4);
			for (uint j_tt=0; j_tt<=ThetaNumber; j_tt++)
			{
				sum.Fill(0);
				double tt = (double)j_tt*dt/Rad;
				for (uint i_fi=0; i_fi<=PhiNumber; i_fi++)
				{
					double fi = -((double)i_fi)*df;
					matrix m = M[q](i_fi,j_tt), L(4,4);
					f << endl << tt << " " << -fi/Rad << " "; f << m;
					L[0][0] = 1.0; L[0][1] = 0.0;			L[0][2] = 0.0;			L[0][3] = 0.0;
					L[1][0] = 0.0; L[1][1] = cos(2.0*fi);	L[1][2] = sin(2.0*fi);	L[1][3] = 0.0;
					L[2][0] = 0.0; L[2][1] =-sin(2.0*fi);	L[2][2] = cos(2.0*fi);	L[2][3] = 0.0;
					L[3][0] = 0.0; L[3][1] = 0.0;			L[3][2] = 0.0;			L[3][3] = 1.0;
					if (!j_tt)
						sum += L*m*L;
					else
						sum += m*L;
				}
				sum /= ((double)PhiNumber+1.0);
				res << endl << tt << " ";
				res << sum;
			}
			f.close();
			res << endl;
			res.close();
		}
	}
	else
	{
		ofstream f((name_f+"_all.dat").c_str(), ios::out),
				res("ga_all.dat", ios::app);
		f << (to_string(ConusGrad)+" "+to_string(PhiNumber)+" "+to_string(ThetaNumber)).c_str();
		res << endl<< "tetta M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44";
		f.precision(10);
		res.precision(10);


		matrix sum(4,4);
		for (uint j_tt=0; j_tt<=ThetaNumber; j_tt++)
		{
			sum.Fill(0);
			double tt = (double)j_tt*dt/Rad;
			for (uint i_fi=0; i_fi<=PhiNumber; i_fi++)
			{
				double fi = -((double)i_fi)*df;
				matrix m = M[0](i_fi,j_tt), L(4,4);
				f << endl << tt << " " << -fi/Rad << " "; f << m;
				L[0][0] = 1.0; L[0][1] = 0.0;			L[0][2] = 0.0;			L[0][3] = 0.0;
				L[1][0] = 0.0; L[1][1] = cos(2.0*fi);	L[1][2] = sin(2.0*fi);	L[1][3] = 0.0;
				L[2][0] = 0.0; L[2][1] =-sin(2.0*fi);	L[2][2] = cos(2.0*fi);	L[2][3] = 0.0;
				L[3][0] = 0.0; L[3][1] = 0.0;			L[3][2] = 0.0;			L[3][3] = 1.0;
				if (!j_tt)
					sum += L*m*L;
				else
					sum += m*L;
			}
			sum /= ((double)PhiNumber+1.0);
			res << endl << tt << " ";
			res << sum;
		}
		f.close();
		res << endl;
		res.close();
	}

	//----------------------------------------------------------------------------
	ofstream _out("out.dat", ios::out);
	tm = clock()-tm;
	_out << "\n\nTotal time of calculation = " << tm/CLK_TCK << " seconds";
	_out.close();
	//----------------------------------------------------------------------------
	cout << "\n\nTotal time of calculation = " << tm/CLK_TCK << " seconds";
	cout << "\nAll done. Please, press any key.";
	getch();

	return 0;
}


void Handler(Beam& bm)
{
	uint  pn = 0, szP = SizeP(bm);
	if (NumberOfTrajectory)
	{
		// поиск - нужно ли учитывать пучок с данной траекторией?
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

		// поиск - есть ли такая траектория в "фиксированном" списке?
		flag = false;
		c = AllowedMask.begin();
		for (;c!=AllowedMask.end();c++, pn++)
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

	//double ctetta = bm.r*k;

	//if (ctetta < 0.17364817766693034885171662676931) return;

	//----------------------------------------------------------------------------
	uint  vi;
	if (pn<=3) vi = pn;
	else if (pn>3 && pn<8) vi = 4;
	else if (pn==8 || pn==9) vi = 5;
	else if (pn==10 || pn==11) vi = 6;
	else if (pn==12 || pn==13) vi = 7;
	else if (pn==14 || pn==15) vi = 8;
	else if (pn==16 || pn==17) vi = 9;
	else if (pn==18 || pn==19) vi = 10;
	else if (pn==20) vi = 11;
	else if (pn==21) vi = 12;
	else if (pn==22) vi = 13;
	else if (pn==23) vi = 14;
	else vi = 15;

	if (NumberOfTrajectory)
		if (sort_lim[Lim]>lim[vi]) return;
	//----------------------------------------------------------------------------

	bm.F = bm.e;
	bm.T = bm.F%bm.r;
	bm.T /=length(bm.T); // базис выходящего пучка

	Point3D r0 = CenterOfBeam(bm);

	double lng_proj0 = bm.lng+r0*bm.r;

	if (perpend_diff)
		bm.SetCoefficients_abcd(bm.r, bm.T, bm.F, r0);
	else
	{
		Point3D Nx, Ny;
		const double tmp = sqrt(SQR(bm.N.x)+SQR(bm.N.y));
		if (fabs(bm.N.z)>1-DBL_EPSILON)
		{
			Nx=Point3D(0,-bm.N.z,0);
			Ny=Point3D(1,0,0);
		}
		else
		{
			Nx=Point3D(bm.N.y/tmp,-bm.N.x/tmp,0);
			Ny=bm.N%Nx;
			//Ny /= length(Ny);
		}
		bm.SetCoefficients_abcd(bm.N, Nx, Ny, r0);
	}
	for (uint  i_fi=0; i_fi<=PhiNumber; i_fi++)
		for (uint  j_tt=0; j_tt<=ThetaNumber; j_tt++)
		{

			double f = (double)i_fi*df,
				t = (double)j_tt*dt,
				cf = cos(f),
				sf = sin(f),
				ct = cos(t),
				st = sin(t);
			Point3D vr(st*cf,st*sf,ct),
				vf,
				vt;

			double ctetta = bm.r*vr;

			if (ctetta > 0.17364817766693034885171662676931)
			{
				if (!j_tt)
					vf = -Ey;
				else
					vf = Point3D(-sf,cf,0);
				vt = vf%vr;
				vt /= length(vt);

				matrixC Jn_rot(2,2);
				complex fn(0,0);
				if (perpend_diff)
				{
					Jn_rot[0][0]=-bm.F*vf; Jn_rot[0][1]= bm.T*vf;
					Jn_rot[1][0]= bm.F*vt; Jn_rot[1][1]=-bm.T*vt;
					fn = bm.DiffractionShiftedPr(vr, lm);
					//fn = bm.DiffractionShifted(vr, lm);
				}
				else
				{
					Jn_rot[0][0]=-(bm.N%bm.T)*vf; Jn_rot[0][1]=-(bm.N%bm.F)*vf;
					Jn_rot[1][0]= (bm.N%bm.T)*vt; Jn_rot[1][1]= (bm.N%bm.F)*vt;
					//fn = bm.DiffractionInclinePr(vr, lm);
					fn = bm.DiffractionIncline(vr, lm);
				}
				matrixC fn_jn = exp_im(m_2pi*(lng_proj0-vr*r0)/lm)*bm();
				Jones_temp[vi].insert(i_fi,j_tt,fn*Jn_rot*fn_jn);
				//Jones_temp[vi].insert(i_fi,j_tt,fn*Jn_rot*bm());
			}
		}
}


//------------------------------------------------------------------------------
void HandlerTemp(Beam& bm )
{
	//----------------------------------------------------------------------------
	uint  pn, szP = SizeP(bm);
	if (NumberOfTrajectory)
	{
		// поиск - нужно ли учитывать пучок с данной траекторией?
		bool flag = false;
		list<Chain>::const_iterator c = mask.begin();
		for (pn=0;c!=mask.end();c++, pn++)
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

		// поиск - есть ли такая траектория в "фиксированном" списке?
		pn = 0;
		flag = false;
		c = AllowedMask.begin();
		for (;c!=AllowedMask.end();c++, pn++)
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
		if (pn!=n_tr) return;
	}
	//----------------------------------------------------------------------------

	bm_tetta = bm.Spherical().tetta;

	std::list<Point3D>::const_iterator pt = bm.Begin();
	Point3D rn = k,
	rn_pr = rn-bm.r*(bm.r*rn),     // проекция на плоскость перпендикулярную пучку.
	cur = (*pt)-bm.r*(bm.r*(*pt)); // проецируем первую вершину на плоскость перпендикулярную пучку,

	float curr = cur*rn_pr, // вычисляем сколярное произведение проекции нормали и проекции точки
	min_ro_n = curr,
	max_ro_n = curr;  // пологаем максимуму и минимум равным первому расстоянию

	for (pt++; pt!=bm.End(); pt++)
	{
		Point3D cur=(*pt)-bm.r*(bm.r*(*pt)); // проецируем вершину на плоскость перпендикулярную пучку,
		float curr=cur*rn_pr; // вычисляем скалярное произведение проекции нормали и проекции точки
		if (curr>max_ro_n) max_ro_n=curr;
		if (curr<min_ro_n) min_ro_n=curr;
	}
	max_ro_n = max_ro_n-min_ro_n;
	NumberOfRing = max_ro_n/lm;
}

const int size = 256;

int ReadFile(char* name, double* params, uint n)
{
	char buf[size]=""; //буфер
	ifstream in(name, ios::in); //входной файл
	for (uint i=0; i<n; i++)
	{
		if (in.eof()) return 1;
		in.getline(buf, size);
		params[i] = strtod(buf, NULL);
		if (params[i]<0) throw "All values must be positive in data file";
	}
	//----------------------------------------------------------------------------
	switch(0)
	{ // only for hexagonal
	case 0: _NoF =  8; break;
	case 1: _NoF = 13; break;
	case 2: _NoF =  7; break;
	case 3: _NoF = 20; break;
	case 4: _NoF = 14; break;
	case 5: _NoF = 6; break;
	case 6: _NoF = 4; break;
	case 7: _NoF = 5; break;
	}
	Itr =					params[9];
	if (!in.eof())
	{
		in.getline(buf, size);
		NumberOfTrajectory = strtod(buf, NULL);
		if (NumberOfTrajectory>0)
		{
			Face = new uint*[int(Itr)+1];
			for (uint i=0;i<Itr+1; i++)
				Face[i] = new uint[_NoF];

			for (uint i=0;i<Itr+1;i++)
				for (uint j=0;j<_NoF;j++)
					Face[i][j] = 0;

			for (uint j=0;j<NumberOfTrajectory;j++)
			{
				if (in.eof()) return 1;
				in.getline(buf, size);
				MaskAppend(buf);
			}
		}
	}
	in.close();
	return 0;
}


void MaskAppend(const char *s)
{
	list<uint> ch;
	uint intern_numb = 0;
	char *buf,*end;
	end=(char*)s;
	do
	{
		buf=end;
		int facet_numb=strtol(buf,&end,10);
		if (strlen(buf)!=strlen(end))
		{
			if ((facet_numb>(int)_NoF)||(facet_numb<0)) throw "Error! Incorrect parameters of trajectories in data file";
			ch.push_back(facet_numb);
			Face[intern_numb][facet_numb] = 1;
			intern_numb++;
		}
	}
	while (strlen(buf)!=strlen(end));
	if (ch.size()==0) throw "Error! There was not enough trajectories in data file";
	mask.push_back(ch);
}

//---------------------------------------------------------------------------
void AllowedMaskAppend(const char *s)
{
	list<uint> ch;
	char *buf,*end;
	end= (char*)s;
	do
	{
		buf=end;
		int facet_numb=strtol(buf,&end,10);
		if (strlen(buf)!=strlen(end))
		{
			if ((facet_numb>(int)_NoF)||(facet_numb<0)) throw "Error! Incorrect parameters of trajectories in data file";
			ch.push_back(facet_numb);
		}
	}
	while (strlen(buf)!=strlen(end));
	if (ch.size()==0) throw "Error! There was not enough trajectories in data file";
	AllowedMask.push_back(ch);
}


//==============================================================================

void DelFace(void)
{
	for (uint i=0;i<_NoF; i++)
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
	/*
	using namespace std;
	SYSTEMTIME Time;
	GetLocalTime(&Time);
	cout << "Numerical calculations started at: " <<
		  Time.wHour << ':' << Time.wMinute << ':' << Time.wSecond;
	*/
}

//---------------------------------------------------------------------------
uint FillLim(double ConusRad, double betta, double sc_coef, const char* name)
{
	ofstream out(name, ios::app);
	out << endl << endl << "betta = " << betta/Rad;
	uint  nl = 0;
	for (uint ns=0;ns<AllowedMuellerMatrixNumber;ns++)
	{
		if (!Mueller_mask[ns]) continue;
		switch (ns)
		{
		case 0: n_tr = 0; break;
		case 1: n_tr = 1; break;
		case 2: n_tr = 2; break;
		case 3: n_tr = 3; break;
		case 4: n_tr = 4; break;
		case 5: n_tr = 8; break;
		case 6: n_tr = 10; break;
		case 7: n_tr = 12; break;
		case 8: n_tr = 14; break;
		case 9: n_tr = 16; break;
		case 10: n_tr = 18; break;
		case 11: n_tr = 20; break;
		case 12: n_tr = 21; break;
		case 13: n_tr = 22; break;
		case 14: n_tr = 23; break;
		case 15: n_tr = 24; break;
		}
		if (ns<2)
		{
			out << endl << '\t' << ns << " " << 30;
			lim[ns] = M_PI/6.0; nl++; continue;
		}
		float gamma_ConusRad, gamma_start=0.0, gamma_end=0.0;
		do
		{
			// ищем gamma_ConusRad - предельный аксиальный угол, когда пучок выходит за конус

			// подбираем предельный угол gamma_end
			Body->ChangePosition(betta, M_PI/6.0+gamma_end,0.0);
			NumberOfRing = -1.0;
			Body->FTforConvexCrystal(HandlerTemp);
			gamma_end += 5.0*Rad;
		} while((bm_tetta<ConusRad)&&(gamma_end<88.0*Rad));  // <88.0???

		if (NumberOfRing<0.0)
		{
			lim[ns] = -1.0; continue; // нет такой траектории
		}

		if (bm_tetta<ConusRad)    //считаем весь диапазон
			gamma_end=30.0*Rad;

		//методом половинного деления ищем gamma_lim
		//с точностью gamma_pres процентов
		float gamma_pres=0.005;// *100 процентов
		while((gamma_end-gamma_start)/gamma_end>gamma_pres)
		{
			gamma_ConusRad=(gamma_end+gamma_start)/2.0;
			Body->ChangePosition(betta, M_PI/6.0+gamma_ConusRad,0.0); //поворот кристалла
			NumberOfRing = -1.0;
			Body->FTforConvexCrystal(HandlerTemp);
			if (bm_tetta>ConusRad)
				gamma_end = gamma_ConusRad;
			else
				gamma_start = gamma_ConusRad;
		}
		gamma_ConusRad = (gamma_end+gamma_start)/2.0;

		if (gamma_ConusRad>M_PI/6.0)
		{
			out << endl << '\t' << ns << " " << 30;
			lim[ns] = M_PI/6.0; nl++; continue;
		}
		else
		{
			if (gamma_ConusRad<=FLT_EPSILON)
			{
				continue;                           // ???
			}
			else
			{
				float MaxNRing = 5.0,                            // номер макс кольца, от которого будем учитывать вклад
				gamma_lim;                                 // предельный аксиальный угол для MaxNRing колец
				gamma_start=0.0; gamma_end=0.0;                  // задаем мин и макc углы для поиска
				do
				{
					// подбираем предельный угол gamma_end, когда выходим за MaxNRing кольцо
					Body->ChangePosition(betta, M_PI/6.0+gamma_end,0.0);
					NumberOfRing = -1.0;
					Body->FTforConvexCrystal(HandlerTemp);
					gamma_end += 5.0*Rad;
				} while((NumberOfRing<MaxNRing)&&(gamma_end<88.0*Rad)); // <88.0???

				if (NumberOfRing<MaxNRing) //считаем весь диапазон
					gamma_end=30.0*Rad;

				//методом половинного деления ищем gamma_lim
				//с точностью gamma_pres процентов
				float gamma_pres=0.005;// *100 процентов
				while((gamma_end-gamma_start)/gamma_end>gamma_pres)
				{
					gamma_lim=(gamma_end+gamma_start)/2.0;
					Body->ChangePosition(betta, M_PI/6.0+gamma_lim,0.0); //поворот кристалла
					NumberOfRing=-1.0;
					Body->FTforConvexCrystal(HandlerTemp);
					if (NumberOfRing>MaxNRing)
						gamma_end = gamma_lim;
					else
						gamma_start = gamma_lim;
				}
				gamma_lim = (gamma_end+gamma_start)/2.0;

				if (gamma_lim>M_PI/6.0) gamma_lim=M_PI/6.0;

				double Tgamma = gamma_ConusRad+gamma_lim*sc_coef;
				if (Tgamma>M_PI/6.0) Tgamma=M_PI/6.0;
				lim[ns] = Tgamma; nl++;
				out << endl << '\t' << ns << " " << Tgamma/Rad;
			}
		}
	}

	double copy_lim[AllowedMuellerMatrixNumber];
	for (uint i=0; i<AllowedMuellerMatrixNumber; copy_lim[i]=lim[i], sort_lim[i]=0.0, i++);
	for (int i=nl-1; i>=0; i--)
	{
		double max = -1.0;
		int mx=0;
		for (uint j=0; j<AllowedMuellerMatrixNumber; j++)
			if (copy_lim[j]>max) { max = copy_lim[j]; mx = j; }
		sort_lim[i] = max;
		copy_lim[mx] = -1.0;
	}
	out.close();
	return nl;
	}
