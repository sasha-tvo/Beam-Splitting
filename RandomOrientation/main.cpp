#include <QCoreApplication>
#include <QDir>

//---------------------------------------------------------------------------
//#include <fstream>
//#include <stdlib.h>
//#include <conio.h>  // input, output: getch
#include <time.h>  // time functions: clock
//#include <list>   // standart list structure
//#include <string.h>

//#include <math.h>

//==============================================================================
#include "..\Lib\Geometry.hpp"
#include "..\Lib\particle.hpp"
#include "..\Lib\PhysMtr.hpp"
#include "..\Lib\Prizm\Mueller.hpp"
//#include "..\Lib\ConsoleBackward.hpp"
#include "..\Lib\trajectory.hpp"

using namespace std;
//==============================================================================
Crystal* Body = NULL;

const int NumberOfParameters = 16;  // количество строк в файле "params.dat"
unsigned int  Itr, Sorting = 0, TypeOfOrientation = 0, TypeOfCalculation = 0, NumFi, NumTetta;
double Rad = M_PI/180., lm, df = 0.0, dt = 0.0, area_min;
list<Chain> mask;
Point3D k(0,0,1.0), Ey(0,1.0,0);

list<Chain> TrBeams;  //маска пучков
const unsigned int  NumTr = 28,  // число траекторий в фиксированном списке
              NumSum = 16; // макс число матриц Мюллера на выходе программы
const string tr_beams[NumTr] = {"0", "0 7 0", "3", "3 6 3", "2 6 7 4", "2 7 6 4", "4 6 7 2", "4 7 6 2", "3 6 7 3", "3 7 6 3", "0 6 7 0",
                                    "0 7 6 0", "3 7 0 6 7 3", "3 7 6 0 7 3", "0 6 3 7 6 0", "0 6 7 3 6 0", "2 6 7 3", "2 7 6 3", "4 6 7 3",
                                    "4 7 6 3", "2 1 6 7 4", "4 7 6 1 2", "2 7 6 5 4", "4 5 6 7 2", "2 7 0 6 7 4", "2 7 6 0 7 4", "4 7 0 6 7 2", "4 7 6 0 7 2"};
int maska[NumSum]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    n_tr, // номер траектории, для которой вычисляется gamma_lim
    Lim;  // номер текущей области из массива sort_lim
double EPS_MODIF,
       NumberOfRing, // номер дифракционного кольца для пучка с номером n_tr
       bm_tetta,     // tetta_координата для пучка с номером n_tr
       lim[NumSum]={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0}, sort_lim[NumSum];
vector<Arr2DC> j_tmp;
bool fl_pl;

//==============================================================================
// handler of the emitted rays

void   Handler(Beam& bm, double, unsigned int )
 {
  unsigned int  pn = 0, szP = SizeP(bm);
  if(Sorting) {
   // поиск - нужно ли учитывать пучок с данной траекторией?
   bool flag = false;
   list<Chain>::const_iterator c = mask.begin();
   for(;c!=mask.end();c++) {
     list<unsigned int >::const_iterator it = c->Ch.begin();
     if(szP!=c->sz) continue;
     list<unsigned int >::const_iterator fs = bm.BeginP();
     for(;it!=c->Ch.end() && (*it)==(*fs);it++, fs++);
     if(it==c->Ch.end()) { flag = true; break; }
   }
   if(!flag) return;

   // поиск - есть ли такая траектория в "фиксированном" списке?
   flag = false;
   c = TrBeams.begin();
   for(;c!=TrBeams.end();c++, pn++) {
     list<unsigned int >::const_iterator it = c->Ch.begin();
     if(szP!=c->sz) continue;
     list<unsigned int >::const_iterator fs = bm.BeginP();
     for(;it!=c->Ch.end() && (*it)==(*fs);it++, fs++);
     if(it==c->Ch.end()) { flag = true; break; }
   }
   if(!flag) return;
  }

  double ctetta = bm.r*k;

  if(ctetta < 0.17364817766693034885171662676931) return;

  //----------------------------------------------------------------------------
  unsigned int  vi;
  if(pn<=3) vi = pn;
  else if(pn>3 && pn<8) vi = 4;
  else if(pn==8 || pn==9) vi = 5;
  else if(pn==10 || pn==11) vi = 6;
  else if(pn==12 || pn==13) vi = 7;
  else if(pn==14 || pn==15) vi = 8;
  else if(pn==16 || pn==17) vi = 9;
  else if(pn==18 || pn==19) vi = 10;
  else if(pn==20) vi = 11;
  else if(pn==21) vi = 12;
  else if(pn==22) vi = 13;
  else if(pn==23) vi = 14;
  else vi = 15;

  if(sort_lim[Lim]>lim[vi]) return;
  //----------------------------------------------------------------------------

  bm.F = bm.e;
  bm.T = bm.F%bm.r; bm.T /=length(bm.T); // базис выходящего пучка

  Point3D r0 = CenterOfBeam(bm);

  double lng_proj0 = bm.lng+r0*bm.r;

  if(fl_pl) bm.InputPr(bm.r, bm.T, bm.F, r0);
  else {
   Point3D Nx, Ny;
   const double tmp = sqrt(SQR(bm.N.x)+SQR(bm.N.y));
   if(fabs(bm.N.z)>1-DBL_EPSILON) {
    Nx=Point3D(0,-bm.N.z,0);
    Ny=Point3D(1,0,0);
   }
   else {
    Nx=Point3D(bm.N.y/tmp,-bm.N.x/tmp,0);
    Ny=bm.N%Nx; Ny /= length(Ny);
   }
   bm.InputPr(bm.N, Nx, Ny, r0);
  }
  for(unsigned int  i_fi=0; i_fi<=NumFi; i_fi++)
   for(unsigned int  j_tt=0; j_tt<=NumTetta; j_tt++) {
     double f = (double)i_fi*df, t = (double)j_tt*dt,
            cf = cos(f), sf = sin(f), ct = cos(t), st = sin(t);
     Point3D vr(st*cf,st*sf,ct), vf, vt;
     if(!j_tt)
      vf = -Ey;
     else
      vf = Point3D(-sf,cf,0);
     vt = vf%vr; vt /= length(vt);

     matrixC Jn_rot(2,2);
     complex fn(0,0);
     if(fl_pl) {
     Jn_rot[0][0]=-bm.F*vf; Jn_rot[0][1]= bm.T*vf;
     Jn_rot[1][0]= bm.F*vt; Jn_rot[1][1]=-bm.T*vt;
     fn = bm.DiffractionShiftedPr(vr, lm);
     }
     else {
      Jn_rot[0][0]=-(bm.N%bm.T)*vf; Jn_rot[0][1]=-(bm.N%bm.F)*vf;
      Jn_rot[1][0]= (bm.N%bm.T)*vt; Jn_rot[1][1]= (bm.N%bm.F)*vt;
      fn = bm.Diffraction_inclinePr(vr, lm);
     }




     matrixC fn_jn = exp_im(m_2pi*(lng_proj0-vr*r0)/lm)*bm();
     j_tmp[vi].insert(i_fi,j_tt,fn*Jn_rot*fn_jn);
  }
 }


//------------------------------------------------------------------------------
void   HandlerTemp(Beam& bm, double, unsigned int ) {
  //----------------------------------------------------------------------------
  unsigned int  pn, szP = SizeP(bm);
  if(Sorting) {
   // поиск - нужно ли учитывать пучок с данной траекторией?
   bool flag = false;
   list<Chain>::const_iterator c = mask.begin();
   for(pn=0;c!=mask.end();c++, pn++) {
     list<unsigned int >::const_iterator it = c->Ch.begin();
     if(szP!=c->sz) continue;
     list<unsigned int >::const_iterator fs = bm.BeginP();
     for(;it!=c->Ch.end() && (*it)==(*fs);it++, fs++);
     if(it==c->Ch.end()) { flag = true; break; }
   }
   if(!flag) return;

   // поиск - есть ли такая траектория в "фиксированном" списке?
   pn = 0;
   flag = false;
   c = TrBeams.begin();
   for(;c!=TrBeams.end();c++, pn++) {
     list<unsigned int >::const_iterator it = c->Ch.begin();
     if(szP!=c->sz) continue;
     list<unsigned int >::const_iterator fs = bm.BeginP();
     for(;it!=c->Ch.end() && (*it)==(*fs);it++, fs++);
     if(it==c->Ch.end()) { flag = true; break; }
   }
   if(!flag) return;
   if(pn!=n_tr) return;
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

  for(pt++; pt!=bm.End(); pt++) {
   Point3D cur=(*pt)-bm.r*(bm.r*(*pt)); // проецируем вершину на плоскость перпендикулярную пучку,
   float curr=cur*rn_pr; // вычисляем скалярное произведение проекции нормали и проекции точки
   if(curr>max_ro_n) max_ro_n=curr;
   if(curr<min_ro_n) min_ro_n=curr;
 }
 max_ro_n = max_ro_n-min_ro_n;
 NumberOfRing = max_ro_n/lm;
}
//---------------------------------------------------------------------------
void AppendTr(char s[])
{
 list<unsigned int > ch;
 int sz = 0;
 char sr[2]={" "};

 for(;s[sz+1]!=0;sz++);
 sz++;
 if(!sz) return;
 unsigned int  i=0,j=1;
 for(;j<sz; j++) {
   if(s[j]!=sr[0]) continue;
   char st[128]="";
   Sub(s,i,j-i, st);
   unsigned int  fst = strtod(st, NULL);
   ch.push_back(fst);
   for(;s[j+1]==sr[0]; j++);
   i=j+1;
   j=i;
 }
 char st[128]="";
 Sub(s,i,j-i+1,st);
 unsigned int  fst = strtod(st, NULL);
 ch.push_back(fst);
 TrBeams.push_back(ch);
}

//---------------------------------------------------------------------------
unsigned int    FillLim(double cone, double betta, double sc_coef, char* name)
{
 ofstream out(name, ios::app);
  out << endl << endl << "betta = " << betta/Rad;
 unsigned int  nl = 0;
 for(int ns=0;ns<NumSum;ns++) {
   if(!maska[ns]) continue;
   switch (ns) {
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
   if(ns<2) {
    out << endl << '\t' << ns << " " << 30;
    lim[ns] = M_PI/6.0; nl++; continue;
   }
   float gamma_cone, gamma_start=0.0, gamma_end=0.0;
   do {
    // ищем gamma_cone - предельный аксиальный угол, когда пучок выходит за конус

    // подбираем предельный угол gamma_end
    Body->ChangePosition(betta, M_PI/6.0+gamma_end,0.0);
    NumberOfRing = -1.0;
    Body->FTforConvexCrystal(HandlerTemp);
    gamma_end += 5.0*Rad;
   } while((bm_tetta<cone)&&(gamma_end<88.0*Rad));  // <88.0???

   if(NumberOfRing<0.0) {
    lim[ns] = -1.0; continue; // нет такой траектории
   }

   if(bm_tetta<cone)    //считаем весь диапазон
    gamma_end=30.0*Rad;

   //методом половинного деления ищем gamma_lim
   //с точностью gamma_pres процентов
   float gamma_pres=0.005;// *100 процентов
   while((gamma_end-gamma_start)/gamma_end>gamma_pres) {
    gamma_cone=(gamma_end+gamma_start)/2.0;
    Body->ChangePosition(betta, M_PI/6.0+gamma_cone,0.0); //поворот кристалла
    NumberOfRing = -1.0;
    Body->FTforConvexCrystal(HandlerTemp);
    if(bm_tetta>cone) gamma_end = gamma_cone;
    else gamma_start = gamma_cone;
   }
   gamma_cone = (gamma_end+gamma_start)/2.0;

   if(gamma_cone>M_PI/6.0) {
    out << endl << '\t' << ns << " " << 30;
    lim[ns] = M_PI/6.0; nl++; continue;
   }
   else if(gamma_cone<=FLT_EPSILON) {
    continue;                           // ???
   }
   else {
    float MaxNRing = 5.0,                            // номер макс кольца, от которого будем учитывать вклад
          gamma_lim;                                 // предельный аксиальный угол для MaxNRing колец
    gamma_start=0.0; gamma_end=0.0;                  // задаем мин и макc углы для поиска
    do {
     // подбираем предельный угол gamma_end, когда выходим за MaxNRing кольцо
     Body->ChangePosition(betta, M_PI/6.0+gamma_end,0.0);
     NumberOfRing = -1.0;
     Body->FTforConvexCrystal(HandlerTemp);
     gamma_end += 5.0*Rad;
    } while((NumberOfRing<MaxNRing)&&(gamma_end<88.0*Rad)); // <88.0???

    if(NumberOfRing<MaxNRing) //считаем весь диапазон
     gamma_end=30.0*Rad;

    //методом половинного деления ищем gamma_lim
    //с точностью gamma_pres процентов
    float gamma_pres=0.005;// *100 процентов
    while((gamma_end-gamma_start)/gamma_end>gamma_pres) {
     gamma_lim=(gamma_end+gamma_start)/2.0;
     Body->ChangePosition(betta, M_PI/6.0+gamma_lim,0.0); //поворот кристалла
     NumberOfRing=-1.0;
     Body->FTforConvexCrystal(HandlerTemp);
     if(NumberOfRing>MaxNRing) gamma_end = gamma_lim;
     else gamma_start = gamma_lim;
    }
    gamma_lim = (gamma_end+gamma_start)/2.0;

    if(gamma_lim>M_PI/6.0) gamma_lim=M_PI/6.0;

    double Tgamma = gamma_cone+gamma_lim*sc_coef;
    if(Tgamma>M_PI/6.0) Tgamma=M_PI/6.0;
    lim[ns] = Tgamma; nl++;
    out << endl << '\t' << ns << " " << Tgamma/Rad;
   }
 }

 double copy_lim[NumSum];
 for(int i=0; i<NumSum; copy_lim[i]=lim[i], sort_lim[i]=0.0, i++);
 for(int i=nl-1; i>=0; i--) {
   double max = -1.0;
   int mx;
   for(int j=0; j<NumSum; j++)
     if(copy_lim[j]>max) { max = copy_lim[j]; mx = j; }
   sort_lim[i] = max;
   copy_lim[mx] = -1.0;
 }
 out.close();
 return nl;
}

//---------------------------------------------------------------------------

#pragma argsused
int main(int argc, char* argv[])
{
  QCoreApplication a(argc, argv);

  clock_t t = clock();
  ShowTitle();
  cout << "\nLoading settings... ";

  QString QCurDir;
  QCurDir=QCoreApplication::applicationDirPath();
  QDir dir;
  dir.cd(QCurDir);
  QDir::setCurrent(dir.absolutePath());

  double params[NumberOfParameters];
  if(ReadFile("params.dat", params, NumberOfParameters)) {
   cout << "\nError! Incorrect input file. Press any key for exit.";
   getch(); return 1;
  }

  if(params[0]<=0 || params[1]<=0 || Sorting==0) {
   cout << "\nError! Incorrect input data. Press any key for exit.";
   getch(); return 1;
  }

  cout << "OK";
  //----------------------------------------------------------------------------

  #ifdef _Windows
  gotoxy(1,7);  ShowCurrentTime();   // current time (Win32 API)
  #endif
  clock_t tm = clock();
  //----------------------------------------------------------------------------
  Body = new Prizm(complex(params[2],0.0), params[1]/2.0, params[0]/2.0);
  Body->Phase() = true;

  int N_bt = params[5], N_gm = params[8], N_or = (N_bt+1)*(2*N_gm+1);
  double betta_min = params[3]*Rad, betta_max = params[4]*Rad,
         gamma_cnt = M_PI/6.0, dbetta = 0.0, dgamma = 1.0, cone = params[12]*Rad, Tgamma;

  if(N_bt) dbetta = (betta_max-betta_min)/N_bt;

  Itr = params[9];
  EPS_MODIF = pow(0.1,params[10]);
  lm = params[11];

  NumTetta = params[13];
  NumFi = params[14];
  if(NumTetta) dt = cone/(double)NumTetta;
  if(NumFi) df = m_2pi/(double)(NumFi+1);

  if(fabs(params[15])>FLT_EPSILON) fl_pl = false;
  else fl_pl = true;

  string  folder_name = "Data", fRes = folder_name,
          str = to_string(params[12])+" "+to_string(NumFi)+" "+to_string(NumTetta);

  dir.mkdir(QString::fromStdString(fRes));
  dir.cd(QString::fromStdString(fRes));
  QDir::setCurrent(dir.absolutePath());

  vector<Arr2D> M;
  unsigned int  orn = 0;
  double dcos_sum = 0.0;
  area_min = 2.5*lm*params[1];
  //cout << endl << folder_name.c_str();
  cout <<endl;
  //----------------------------------------------------------------------------
  TrBeams.clear();
  for(unsigned int  i=0; i<NumTr; i++)
      AppendTr(const_cast <char*> (tr_beams[i].c_str()));

  list<Chain>::const_iterator cm = mask.begin();
  for(;cm!=mask.end();cm++) {
    unsigned int  pn = 0, sm = cm->sz;
    for(list<Chain>::const_iterator ct = TrBeams.begin();ct!=TrBeams.end();ct++, pn++) {
      unsigned int  st = ct->sz;
      if(sm!=st) continue;
      list<unsigned int >::const_iterator im = cm->Ch.begin(), it = ct->Ch.begin();
      for(;im!=cm->Ch.end() && (*im)==(*it);im++, it++);
      if(im==cm->Ch.end()) {
       unsigned int  vi;
       string tr;
       if(!pn) { vi = 0; tr = "0"; }
       else if(pn==1) { vi = 1; tr = "070"; }
       else if(pn==2) { vi = 2; tr = "3"; }
       else if(pn==3) { vi = 3; tr = "363"; }
       else if(pn>3 && pn<8) { vi = 4; tr = "2674_2764_4672_4762"; }
       else if(pn==8 || pn==9) { vi = 5; tr = "3673_3763"; }
       else if(pn==10 || pn==11) { vi = 6; tr = "0670_0760"; }
       else if(pn==12 || pn==13) { vi = 7; tr = "370673_376073"; }
       else if(pn==14 || pn==15) { vi = 8; tr = "063760_067360"; }
       else if(pn==16 || pn==17) { vi = 9; tr = "2673_2763"; }
       else if(pn==18 || pn==19) { vi = 10; tr = "4673_4763"; }
       else if(pn==20) { vi =  11; tr = "21674"; }
       else if(pn==21) { vi =  12; tr = "47612"; }
       else if(pn==22) { vi =  13; tr = "27654"; }
       else if(pn==23) { vi =  14; tr = "45672"; }
       else { vi =  15;  tr = "2674XL"; }
       maska[vi] = 1;
       ofstream out(("ga_"+tr+".dat").c_str(), ios::out);
       out << to_string(N_bt) << " " << to_string(NumTetta) << " " << to_string(dbetta/Rad) << " " << to_string(dt/Rad);
       out.close();
       break;
      }
    }
  }

  //----------------------------------------------------------------------------
  ofstream out("out.dat", ios::out);
  out.close();
  double sum_dgamma = 0.0, norm = (N_gm>0) ? 1.0/(M_PI/3.0) : 1.0;
  try {
   for(int Betta_i=0; Betta_i<=N_bt; Betta_i++) {
     if(!Betta_i) cout << endl << "0% ";
     double betta = betta_min+(double)Betta_i*dbetta;
     unsigned int  NumLim = FillLim(cone,betta,params[6],"out.dat");
     if(!NumLim) continue;
     sum_dgamma = 0.0;
     M.clear();
     Arr2D M_(NumFi+1,NumTetta+1,4,4); M_.ClearArr();
     for(int q=0;q<NumSum; q++)
       M.push_back(M_);
     //-------------------------------------------------------------------------
     for(Lim=0;Lim<NumLim;Lim++) {
       cout << ". ";
       if(N_gm) dgamma = sort_lim[Lim]/N_gm;
       int n_gm;
       if(Lim) n_gm = floor(sort_lim[Lim-1]/dgamma); // n_gm - число "внутренних" шагов в цикле, которые нужно пропустить
       for(int Gamma_j=-N_gm; Gamma_j<=N_gm; Gamma_j++) {
         if(Lim && abs(Gamma_j)<n_gm) continue; // пропускаем внутренние шаги
         double gamma = gamma_cnt+Gamma_j*dgamma,
                Pgamma = dgamma;
         if(Lim && abs(Gamma_j) == n_gm) Pgamma = (n_gm+0.5)*dgamma-sort_lim[Lim-1];
         if(N_gm && (abs(Gamma_j) == N_gm)) Pgamma -= 0.5*dgamma;
         sum_dgamma += Pgamma;
         Body->ChangePosition(betta, gamma, 0);
         j_tmp.clear();
         Arr2DC tmp_(NumFi+1,NumTetta+1,2,2); tmp_.ClearArr();
         for(int q=0;q<NumSum; q++)
           j_tmp.push_back(tmp_);
         double s = Body->FTforConvexCrystal(Handler);
         for(int q=0; q<NumSum; q++) {
           if(maska[q]==1) {
            for(unsigned int  j_tt=0; j_tt<=NumTetta; j_tt++)
              for(unsigned int  i_fi=0; i_fi<=NumFi; i_fi++) {
                matrix m_tmp = Mueller(j_tmp[q](i_fi,j_tt));
                M[q].insert(i_fi,j_tt,Pgamma*norm*m_tmp);
              }
           }
         }
     }
     //cout << endl << to_string(Lim);
   }
   // вывод в файл
   string name_f = "b_"+to_string(betta/Rad);
   for(int q=0; q<NumSum; q++) {
     string tr;
     if(!q) tr = "0";
     else if(q==1) tr = "070";
     else if(q==2) tr = "3";
     else if(q==3) tr = "363";
     else if(q==4) tr = "2674_2764_4672_4762";
     else if(q==5) tr = "3673_3763";
     else if(q==6) tr = "0670_0760";
     else if(q==7) tr = "370673_376073";
     else if(q==8) tr = "063760_067360";
     else if(q==9) tr = "2673_2763";
     else if(q==10) tr = "4673_4763";
     else if(q==11) tr = "21674";
     else if(q==12) tr = "47612";
     else if(q==13) tr = "27654";
     else if(q==14) tr = "45672";
     else if(q==15) tr = "2674XL";
     if(maska[q]==1) {

     ofstream f((name_f+"_"+tr+".dat").c_str(), ios::out), res(("ga_"+tr+".dat").c_str(), ios::app);
     f << str.c_str();

     matrix sum(4,4);
     for(int j_tt=0; j_tt<=NumTetta; j_tt++) {
       sum.Fill(0);
       double tt = (double)j_tt*dt/Rad;
       for(int i_fi=0; i_fi<=NumFi; i_fi++) {
         double fi = -((double)i_fi)*df;
         matrix m = M[q](i_fi,j_tt), L(4,4);
         f << endl << tt << " " << -fi/Rad << " "; f << m;
         L[0][0] = 1.0; L[0][1] = 0.0; L[0][2] = 0.0; L[0][3] = 0.0;
         L[1][0] = 0.0; L[1][1] = cos(2.0*fi); L[1][2] = sin(2.0*fi); L[1][3] = 0.0;
         L[2][0] = 0.0; L[2][1] =-sin(2.0*fi); L[2][2] = cos(2.0*fi); L[2][3] = 0.0;
         L[3][0] = 0.0; L[3][1] = 0.0; L[3][2] = 0.0; L[3][3] = 1.0;
         if(!j_tt) sum += L*m*L;
         else sum += m*L;
       }
       sum /= (NumFi+1.0);
       res << endl << betta/Rad << " " << tt << " ";;
       res << sum;
     }
     f.close();
     res << endl;
     res.close();
    }
   }
    orn++;
    cout <<"\r                        \r";
    cout << (orn*100./(N_bt+1))<<"%"<<" ";
  }
  }
  catch(char* s) {
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

  return 0;//a.exec();;
}
//---------------------------------------------------------------------------
