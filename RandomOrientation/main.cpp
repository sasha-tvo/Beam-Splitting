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
				GammaLimmExtendedCoeff,		///< Coeffitient by which Gamma Limm will be extended
				lm,							///< Wave Length
				ConusGrad,					///< Conus in backscattering (grad)
				ConusRad,					///< Conus in backscattering (rad)
				P = 0;						///< Probability distribution for Betta angle
Arr2D 			mxd(0,0,0,0);				///< An array of output Mueller matrixes
matrix 			back(4,4),					///< Mueller matrix in backward direction
				forw(4,4);					///< Mueller matrix in forward direction
list<Chain>		mask;						///< List of trajectories to take into account
Point3D			k(0,0,1),					///< Direction on incident wave
				Ey(0,1,0);					///< Basis for polarization characteristic of light
bool			perpend_diff;				///< If true - diffraction will be calculated with Shifted screen, perpendicular to propagation direction; if false - diffraction will be calculated with incline screen
//==============================================================================

struct segment
{
	float A,B,Gm;
};

struct Gamma_segments
{

	float betta;
	vector<bool> group_mask;
	int segment_count;
	vector<segment> segm;
};


const uint		NumberOfParameters = 16;		///< Number of lines in data files, except trajectories
double			params[NumberOfParameters];		///< Array of input data
double			NormGammaAngle = 	0,			///< Normalize coefficient for Gamma
				NormBettaAngle =  	0;			///< Normalize coefficient for Betta
const double	Rad = M_PI/180.;				///< Grad to Rad coefficient
double			df = 0.0,						///< dPhi
				dt = 0.0;						///< dTheta
uint			MuellerMatrixNumber=0;
Gamma_segments	CurrentSegment,PriviusSegment,NextSegment;

vector<string>			names_of_groups;
vector<matrixC>			Jones_temp;
list<Gamma_segments>	Map;



///< handler for the emitted beams
void Handler(Beam& bm);

/// Reads the parameters from data file
int ReadFile(char* name, double* params, uint n);

/// Fill in the \b mask and \b **Face from data file
int MaskAppend(const char *s);

/// Build Names Of Groups in \b names_of_groups by mask
void BuildNamesOfGroups(void);

/// Shows the title
void ShowTitle(void);

/// Shows current time
void ShowCurrentTime(void);

/// Deletes \b **Face structure
void DelFace(void);


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
	Halh_Height =			params[0]/1.0;
	Radius =				params[1]/1.0;
	_RefI =					complex(params[2],0);
	betta_min_Rad =			params[3]*Rad;
	betta_max_Rad =			params[4]*Rad;
	BettaNumber =			params[5];
	GammaLimmExtendedCoeff = params[6];
	GammaNumber =			params[8];
	Itr =					params[9];
	lm =					params[11];
	ConusGrad =				params[12];
	ThetaNumber =			params[13];
	PhiNumber =				params[14];
	perpend_diff =			(params[15]>0.5?true:false);

	if (Halh_Height<=0 || Radius<=0 || NumberOfTrajectory==0)
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
			cout << "Wrong kind of particle!";
			return 1;
	}
	//----------------------------------------------------------------------------
	Body->Phase() = true;

	double 	gamma_cnt = M_PI/6.0, dBettaRad = 0.0, dGammaRad = 1.0, ConusRad = ConusGrad*Rad;//, Tgamma;
	if (BettaNumber) dBettaRad = (betta_max_Rad-betta_min_Rad)/(double)BettaNumber;
	if (ThetaNumber) dt = ConusRad/(double)ThetaNumber;
	if (PhiNumber) df = m_2pi/(double)(PhiNumber+1);

	string  folder_name = "Data", fRes = folder_name,
		str = to_string(ConusGrad)+" "+to_string(PhiNumber)+" "+to_string(ThetaNumber);

	dir.mkdir(QString::fromStdString(fRes));
	dir.cd(QString::fromStdString(fRes));
	QDir::setCurrent(dir.absolutePath());

	matrix M_(4,4);
	M_.Fill(0);
	vector<matrix> M;
	M.clear();
	for (uint q=0;q<=MuellerMatrixNumber; q++)
		M.push_back(M_);

	uint orn = 0;
	cout <<endl;
	for (uint vi=0;vi<=MuellerMatrixNumber;vi++)
	{
		string tr=names_of_groups[vi];
		ofstream out(("ga_"+tr+".dat").c_str(), ios::out);
		out << "betta M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44"<<endl;
		out << to_string(BettaNumber) << " " << to_string(ThetaNumber) << " " << to_string(dBettaRad/Rad) << " " << to_string(dt/Rad)<<endl;
		out.close();
	}

	{
		ofstream res("ga_all.dat", ios::out);
		res << "betta M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44"<<endl;
		//res << to_string(BettaNumber) << " " << to_string(ThetaNumber) << " " << to_string(dBettaRad/Rad) << " " << to_string(dt/Rad)<<endl;
		res.close();
	}

	list<Gamma_segments>::iterator p=Map.begin(),pp=Map.begin(),pn=Map.begin();
	pn++;
	//CurrentSegment=*p;
	//PriviusSegment=*pp;
	//NextSegment=*pn;

	ofstream out("out.dat", ios::out);
	out.close();	
	double norm = 1.0/(M_PI/3.0);
	//double
	try
	{
		for (int Betta_i=0; Betta_i<=BettaNumber; Betta_i++)
		{
			if (!Betta_i)
				cout << endl << "0% ";
			double bettaRad = betta_min_Rad+(double)Betta_i*dBettaRad;
			// do while pp.betta <bettaRad < pn.betta
			for (;(pn->betta<bettaRad) && (pn!=Map.end());pn++,pp++);
			if (abs(bettaRad-pn->betta)<abs(bettaRad-pp->betta))
				CurrentSegment=*pn;
			else
				CurrentSegment=*pp;
			cout <<endl << pp->betta<< " "<<CurrentSegment.betta << " " <<bettaRad<<" "<<pn->betta;
			cout <<" "<< CurrentSegment.segment_count;

			for (int gs=0; gs<CurrentSegment.segment_count;gs++)
			{
				double A,B,Gm,Gm_curr;
				int steps=0;
				A=CurrentSegment.segm[gs].A;
				B=CurrentSegment.segm[gs].B;
				Gm=CurrentSegment.segm[gs].Gm;
				steps=(B-A)/Gm;
				if (steps==0) steps=1;
				Gm_curr= (B-A)/steps;
				cout << ". ";
				double gamma=A+Gm_curr*0.5;
				double Pgamma = Gm_curr;
				for (int st=0;  st<steps; st++)
				{					
					//cout << " "<<gamma/Rad;
					Body->ChangePosition(bettaRad, gamma, 0);
					Jones_temp.clear();
					matrixC tmp_(2,2);
					tmp_.Fill(complex(0,0));
					for (uint q=0;q<=MuellerMatrixNumber; q++)
						Jones_temp.push_back(tmp_);
					Body->FTforConvexCrystal(Handler);
					for (uint q=0; q<=MuellerMatrixNumber; q++)
					{
						matrix m_tmp = Mueller(Jones_temp[q]);
						M[q]+=Pgamma*norm*m_tmp;
					}
					gamma = gamma+Gm_curr;
				}

			}
			orn++;
			cout <<endl;
			cout <<"\r                        \r";
			cout << (orn*100./(BettaNumber+1))<<"%"<<" ";

			//------------------
			matrix m_tot(4,4);
			m_tot.Fill(0);

			for (uint q=0; q<=MuellerMatrixNumber; q++)
			{
				string tr;
				tr=names_of_groups[q];
				if (tr!="")
				{
					ofstream res(("ga_"+tr+".dat").c_str(), ios::app);
					res.precision(10);
					matrix m = M[q];
					M[q].Fill(0);

					double tt1,tt2;
					tt1=m[1][1];
					tt2=m[2][2];
					m[1][1]=(tt1-tt2)*0.5;
					m[2][2]=(tt2-tt1)*0.5;
					m[0][1] = 0.0; m[0][2] = 0.0;	m[0][3] = 0.0;
					m[1][0] = 0.0; m[1][2] = 0.0;	m[1][3] = 0.0;
					m[2][0] = 0.0; m[2][1] = 0.0;	m[2][3] = 0.0;
					m[3][0] = 0.0; m[3][1] = 0.0;	m[3][2] = 0.0;
					res <<" "<< Betta_i << " ";
					res << m;
					res << endl;
					m_tot+=m;
					res.close();
				}
			}
			ofstream res("ga_all.dat", ios::app);
			res.precision(10);
			res <<" "<< bettaRad/Rad << " ";
			res << m_tot;
			res << endl;
			res.close();

			//---------------
		}
	}
	catch(const char* s)
	{
		cout << endl << s << "\nPress any key.";
		getch(); return 1;
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
	uint szP = SizeP(bm);
	uint vi=0;
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
				vi=c->group;
				break;
			}
		}
		if (!flag) return;

	}
	else
		throw " No one trajectory!";
	if (!CurrentSegment.group_mask[vi]) return;

	double ctetta = bm.r*k;

	if (ctetta < 0.17364817766693034885171662676931) return;

	bm.F = bm.e;
	bm.T = bm.F%bm.r;
	bm.T /=length(bm.T); // базис выходящего пучка

	Point3D r0 = CenterOfBeam(bm);

	double lng_proj0 = bm.lng+r0*bm.r;

	Point3D vr(0,0,1),vt,vf;
	vf=-Ey;

	vt = -Ey%vr;
	vt /= length(vt);

	matrixC Jn_rot(2,2);
	complex fn(0,0);
	Point3D n0=bm.r;
	Point3D ns=bm.N;
	Point3D nr=vr;
	switch(perpend_diff)
	{
	case 0:
		{
			Jn_rot[0][0]= bm.F*Ey; Jn_rot[0][1]= -bm.T*Ey;
			Jn_rot[1][0]= bm.F*vt; Jn_rot[1][1]= -bm.T*vt;
			fn = bm.DiffractionShifted(vr, lm);
		}; break;
	case 1:
		{
			Jn_rot[0][0]= (bm.N%bm.T)*Ey; Jn_rot[0][1]= (bm.N%bm.F)*Ey;
			Jn_rot[1][0]= (bm.N%bm.T)*vt; Jn_rot[1][1]= (bm.N%bm.F)*vt;
			fn = bm.DiffractionIncline(vr, lm);
		}; break;
	case 2:
		{
			Jn_rot[0][0]=((nr%(ns%bm.T))-nr%(nr%(n0%(ns%bm.T))))*vt/2.0; Jn_rot[0][1]=((nr%(ns%bm.F))-nr%(nr%(ns%(n0%bm.F))))*vt/2.0;
			Jn_rot[1][0]=((nr%(ns%bm.T))-nr%(nr%(n0%(ns%bm.T))))*vf/2.0; Jn_rot[1][1]=((nr%(ns%bm.F))-nr%(nr%(ns%(n0%bm.F))))*vf/2.0;
			fn = bm.DiffractionIncline(vr, lm);
		}; break;
	default:  cout<< "Error in diffraction type!!!";
	}
	matrixC fn_jn = exp_im(m_2pi*(lng_proj0-vr*r0)/lm)*bm();
	Jones_temp[vi]+=fn*Jn_rot*fn_jn;
}


int ReadFile(char* name, double* params, uint n)
{
	const int size = 256;
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
	in.close();

	ifstream map_f("map.dat", ios::in); //входной файл для карты рассчетов

	map_f.getline(buf, size);
	NumberOfTrajectory=strtod(buf, NULL);

	if (NumberOfTrajectory>0)
	{
		Face = new uint*[int(params[9])+1];
		for (uint i=0;i<params[9]+1; i++)
			Face[i] = new uint[_NoF];

		for (uint i=0;i<params[9]+1;i++)
			for (uint j=0;j<_NoF;j++)
				Face[i][j] = 0;
		vector<int> tr_to_group;

		for (uint j=0;j<NumberOfTrajectory;j++)
		{
			if (map_f.eof()) return 1;
			map_f.getline(buf, size);
			tr_to_group.push_back(MaskAppend(buf));
		}
		BuildNamesOfGroups();

		//-------------------
		int bt_c,gm_c;
		float bt_max_rad,gm_min_rad,gm_max_rad,dbt_rad,dgm_rad;
		map_f >>bt_c>>gm_c>>bt_max_rad>>gm_min_rad>>gm_max_rad;
		bt_max_rad*=Rad;
		gm_min_rad*=Rad;
		gm_max_rad*=Rad;
		dbt_rad=bt_max_rad/((float)bt_c);
		dgm_rad=(gm_max_rad-gm_min_rad)/(float)gm_c;
		for (uint i=0;i<bt_c;i++)
		{
			int t;
			segment seg;
			Gamma_segments gm_seg;
			map_f >>t;
			gm_seg.betta=((float)t+0.5)*dbt_rad;
			for (uint j=0; j<MuellerMatrixNumber;j++)
			{
				gm_seg.group_mask.push_back(false);
			}
			for (uint j=0; j<NumberOfTrajectory;j++)
			{
				map_f>>t;
				if (t>0) gm_seg.group_mask[tr_to_group[j]]=true;
			}
			int seg_count;
			map_f >>seg_count;
			gm_seg.segment_count=0;
			for (uint j=0;j<seg_count;j++)
			{
				float A,B,Gm;
				map_f >>A >>B>>Gm;
				if (Gm<gm_c+50)
				{
					seg.A=A*dgm_rad+gm_min_rad;
					seg.B=B*dgm_rad+gm_min_rad;
					seg.Gm=Gm*dgm_rad/params[8];
					gm_seg.segm.push_back(seg);
					gm_seg.segment_count++;
				}
			}
			Map.push_back(gm_seg);
		}
/*
		vector<Gamma_segments>::iterator p=Map.begin();
		for (;p!=Map.end();p++)
		{
			Gamma_segments gm_seg=*p;
			cout <<endl<< gm_seg.betta << " ";
			for (uint j=0; j<NumberOfTrajectory;j++)
			{
				cout << gm_seg.group_mask[tr_to_group[j]] << " ";
			}
			for (uint j=0;j<gm_seg.segment_count;j++)
				cout<< gm_seg.segm[j].A<<" " << gm_seg.segm[j].B<<" "<< gm_seg.segm[j].Gm<<" ";

			cout <<endl;
		}
		throw "Done";
*/




		//-------------------



	}

	map_f.close();
	return 0;
}


int MaskAppend(const char *s)
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
	buf++;
	int facet_group=strtol(buf,&end,10);
	if (strlen(buf)==strlen(end)) throw "Error! there is no group for a trajectory";
	if (facet_group>(int)MuellerMatrixNumber) MuellerMatrixNumber=facet_group;
	if (ch.size()==0) throw "Error! There was not enough symbols in a trajectory in data file";

	mask.push_back(Chain(ch,facet_group));
	return facet_group;
}

void BuildNamesOfGroups()
{
	for (uint i=0;i<=MuellerMatrixNumber;i++)
	{
		string name="";
		list<Chain>::const_iterator c = mask.begin();
		for ( ;c!=mask.end();c++)
			if (c->group==i)
			{
				list<uint>::const_iterator n = c->Ch.begin();
				for (;n!=c->Ch.end();n++)
					name+=to_string(*n)+'_';
				name+='-';
			}
		names_of_groups.push_back(name);
	}
}


//==============================================================================

void DelFace(void)
{
	for (uint i=0;i<Itr+1; i++)
		delete[] Face[i];
	delete[] Face;
}

void ShowTitle(void)
{
	cout << "*************************************************************\
		 \n Light Scattering by Nonspherical Particles.                 \
		 \n (c)Group of Wave Dispersion Theory,                    \
		 \n    Institute of Atmospheric Optics RAS, Tomsk, Russia, 2014 \
		 \n GNU GPL v.2 http://www.gnu.org/licenses/old-licenses/gpl-2.0.html \
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

