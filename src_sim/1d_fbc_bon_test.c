#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_sf_zeta.h>

//#include "ccalloc.h"

//=============================================================================

#define L 2048  
#define ti 21.0013
//12.362277115590198 
//
//8.4350701763286864
//8.4347423072790413   
//10.8421,7.3471,5.51954,4.36308
#define sg 0.1
#define MCini L
#define MCd 20
#define mea 50000
#define MCtot MCini+mea*MCd
//#define ti 2.5
#define tf ti
#define dt -0.2
#define bn L*(L-1)/2
#define pi 4*atan(1.0)
#define k_m pi/(L+1)

gsl_rng *gslran;

void alias_tables (int *a, double *bp, double Jtot, double *dist, int *bn1,int *bn2)
{
  int cal,coun,al,i;
  int rn[bn];
  double p[bn],re[bn]; 
//  char ttt[20];
  for (i = 0; i < bn; ++i) {
   p[i] = dist[bn2[i]-bn1[i]-1]/Jtot;
   bp[i] = p[i]*L*(L-1)/2;
//   rn[i] = 0;
//   re[i] = 0;
  } 
  cal = -1;
  coun = bn-1;
  for (i = 0; i < bn; ++i) {
    if ( bp[i] >= 1) {
      cal += 1;
      re[cal] = bp[i];
      rn[cal] = i;
    }
  }
  for (i = bn-1; i > -1; --i) {
    if ( bp[i] < 1 ) {
      cal += 1;
      coun -= 1;
      re[cal] = bp[i];
      rn[cal] = i;
    }
  }
  al = bn-1;
  while (al > -1 && coun > -1) {
    a[rn[al]] = rn[coun];
    re[coun] = re[coun]-1+re[al];
    al -= 1;
    if (re[coun] < 1) {
      coun -= 1;
    }
  }
  for (i = 0; i < bn; ++i) bp[rn[i]] = re[i];
//  for (i=0;i<bn;++i){ 
//    printf("i= %d a= %d bp= %1.4f \n",i,a[i],bp[i]);
//  }
//  scanf("%20[^\n]",ttt);
}

//=============================================================================

void lattice (int *Spin)
{
  int i;
  for (i = 0; i < L; ++i) {
    double s = gsl_rng_uniform(gslran);
    if (s < 0.5) {
      Spin[i] = -1;
    } else {
      Spin[i] = 1;
    }
  }
}
//=============================================================================
double G_Fourier (int *Spin,int m)
{
  int i;
  double G,sum_sin;
  sum_sin = 0;
  for (i = 0; i < L; ++i) {
    sum_sin += sin(k_m*m*i)*Spin[i];
  }
  G = (sum_sin*sum_sin)/(L);

  return G;
}
//=============================================================================
double gr_uncor (int *Spin, int m)
{
  int i,a;
  double gr;
  gr = 0;
  
  for (i = 0; i < L-m; ++i) { 
       a = i+m;
    gr += Spin[i]*Spin[a];
  }
  gr = gr/(L-m);
  return gr;
}
//=============================================================================
double g_Lhalf (int *Spin, int m)
{
  int i;
  double gr;
  gr = 0;
  for (i = 0; i < m; ++i) gr += Spin[i]*Spin[i+m];
      gr = gr/(L);
  return gr;
}
//=============================================================================
double av_mag (int *Spin, double mag, int n)
{
  double s_mag;
   if ( (int) mag == 0) {
     s_mag = Spin[n];
   } else {
     s_mag = Spin[n]*fabs(mag)/mag;
   }
  return s_mag;
}
//=============================================================================
double av_mag_k (int *Spin, double *av_mag, int n)
{
  int i;
  double f_mag,sum_sin;
   sum_sin = 0;
  for (i = 0; i < L; ++i) {
    sum_sin += sin(k_m*n*i)*av_mag[i];
  }
  f_mag = sqrt(sum_sin*sum_sin);
  return f_mag;
}
//=============================================================================
double modes (int *Spin, int m)
{
  int i;
  double chi,mag_sin;
   mag_sin = 0;
   for (i = 0; i < L; ++i) {
     mag_sin += sin(k_m*m*(i+1.0) )*Spin[i];
   }
      chi = (mag_sin*mag_sin)/(L);
  return chi;
}
//=============================================================================
int findroot (int i, int *clu)
{
int r,s;
  r = s = i;
  while (clu[r] >= 0) {
    clu[s] = clu[r];
    s = r;
    r = clu[r];
  }
  return r;
}

//=============================================================================

void todo (int *Spin, double *bp, int *a, double beta, double Jtot, int *clu, int *bn1, int *bn2, double *dist)
{
  int i,j,m,imc,k,summu,ms,cs,I1,I2,R1,R2,R;
  double lam,mag,ene,S_B,S_C;
    double sum,sume,summ,sume2,summ2,summ4,sumgr0,sumcgrL2,sumugrL2,sumx_i,sumx_i0;
  double summ_i,summ_i2;
  double sumCgrL2[L],sumUgrL2[L],cgr[L],ugr[L];
  double gr[L],G[L],G_r[L],G_u[L],s_mag[L],s_mag_k[L],ss_mag[L],ss_mag_k0[L],ss_mag_kall[L],sum_chi_k[L],sum_m_k[L],sumG_r[L],sumG_u[L];
  double chi_k[L],chi_i[L],chi_i0[L],m_k[L],m_i[L],data[L],save_s_mag_k[L],ss[L];
  sum = sume = summ = sume2 = summ2 = summ4 = sumgr0 = sumugrL2 = sumcgrL2 = summ_i = summ_i2 = sumx_i = 0;
for (i = 0; i < L; ++i) {
   gr[i] = 0; G[i] = 0; G_u[i] = 0; G_r[i] = 0;
   ss[i] = 0;
   sumCgrL2[i] = 0;
   sumUgrL2[i] = 0;
   sumG_r[i] = 0; sumG_u[i] = 0;
   sum_m_k[i] = 0; sum_chi_k[i] = 0;
  }
  
  FILE *gr_file;
//  for (i = 0; i < L; ++i) sumMS[i] = 0;
  for (imc = 0; imc < MCtot; ++imc)  {
    for (k = 0; k < L; ++k) clu[k] = -1;
    summu = 0;
    lam=2*Jtot*beta;
    int K = gsl_ran_poisson(gslran,lam);
    for (ms = 0; ms < K; ++ms) {
      cs = gsl_rng_uniform_int(gslran,bn);
      I1 = bn1[cs];
      I2 = bn2[cs];
      if (gsl_rng_uniform(gslran) > bp[cs]) {
        cs = a[cs];
        I1 = bn1[cs];
        I2 = bn2[cs];
      }
      if (Spin[I1]*Spin[I2] == 1) {
        summu += 1;
        R1 = findroot (I1,clu);
        R2 = findroot (I2,clu);
        if (R1 != R2) {
          if (clu[R1] >= clu[R2]) {
            clu[R2] = clu[R2] + clu[R1];
            clu[R1] = R2;
          } else {
            clu[R1] = clu[R1] + clu[R2];
            clu[R2] = R1;
          }
        }
      }
    }  
    for (k = 0; k < L; ++k) {
      if (clu[k] < 0) {
        if (gsl_rng_uniform(gslran) < 0.5) {
          Spin[k] = -Spin[k];
        }
      }
    }
    for (k = 0; k < L; ++k) {
      if (clu[k] >= 0) {
        R = clu[k];
        while (clu[R] >= 0) R = clu[R];
        Spin[k] = Spin[R];
      }
    }
    if ((imc > MCini) && (MCd*(imc/MCd)) == imc) {
      mag = 0;
      for (i = 0; i < L; ++i) mag += Spin[i];
         // FOURIER CORR FUNCT
       for (m = 0; m < L; m += 1) ugr[m] = gr_uncor(Spin,m);
    // for (m = 0; m < 3; ++m) G[m] = G_Fourier(Spin,m);
    // for (m = 0; m < 3; ++m) sumg[m] += G[m];
      // K MODES
      for (m = 0; m < 3; ++m) chi_k[m] = modes(Spin,m);
      for (m = 0; m < 3; ++m) m_k[m] = sqrt(chi_k[m]/(L));
     
        sum += 1.;
        sume += summu;
        summ += fabs(mag);
        summ2 += mag*mag;
        summ4 += mag*mag*mag*mag;
//        S_B = 0.5*(Spin[L/2-1]+Spin[L/2]);
//        S_C = 0.5*(Spin[0]+Spin[L-1]);;
      
      for (m = 0; m < 3; ++m) sum_chi_k[m] += chi_k[m];
      for (m = 0; m < 3; ++m) sum_m_k[m] += m_k[m];
     // MAG PROFILE
      
      //for (i = 0; i < L; ++i) s_mag[i] += av_mag(Spin,mag,i);
     
      for (m = 0; m < L; m += 1){
	 sumG_r[m] += G_r[m];
	 sumG_u[m] += G_u[m];
	 //sumCgrL2[m] += cgr[m]; 
        sumUgrL2[m] += ugr[m];
      }
      
      
    }
  }
//  double E = sume/(sum*L);
  double E = (Jtot-sume/(sum*beta))/L;
  double M = summ/(sum*L);
  double M2 = summ2/(sum*L*L);
  double M4 = summ4/(sum*L*L*L*L);
  double X = L*(M2-M*M);
  double Q = M2*M2/M4;
 // double xi = L/(2*pi)*pow((double)(abs(sumg[1]/sumg[2])-1.),1./sg);
  for (i = 0; i < L; ++i) s_mag[i] = s_mag[i]/(sum);
  double av_sum_m_k = 0;
       for (m = 0; m < L; m += 1) {
	sumCgrL2[m] = sumUgrL2[m]/(sum)-M*M; 
        sumUgrL2[m] = sumUgrL2[m]/(sum);
      }
  
 gr_file = fopen("./GR_2048_01_FBC_TC_s.dat","w+");
  printf("L= %d T= %2.4f E= %1.4f M= %1.8f X= %3.4f \n",L,1/beta,E,M,X);
  // printf("L= %d T= %2.4f E= %1.4f M= %1.8f X= %3.4f X_mo= %3.4f X_b= %3.4f X_i= %3.4f X_bi= %3.4f Gc[1]= %1.8f Gu[1]= %1.8f Gc[L/2]= %1.8f Gu[L/2]= %1.8f  \n",L,1/beta,E,M,X,X_mo,X_b,X_i,X_bi,sumCgrL2[1],sumUgrL2[1],sumcGrL2,sumuGrL2);
 // for (m = 0; m < L/2+1; ++m)
     for (m = 0; m < L; m += 1) fprintf(gr_file,"m= %1.8f %d Gc= %1.8f Gcf= %1.8f Gu= %1.8f Guf= %1.8f GuLd2= %1.8f \n",(double) m/(L),m,sumCgrL2[m],sumG_r[m],sumUgrL2[m],sumG_u[m],sqrt(L)*sumUgrL2[m]);
  fclose(gr_file);
}

//=============================================================================

int main()
{
  int tde,i,it,j1,j2;
  double beta,T,Jtot;
  double pow(double x, double y);
  int bcoun = 0;
  int bn1[bn],bn2[bn];
  int k = 0;
  //for (j1 = 0; j1 < bn; ++j1) {
  //    bn1[j1] = 0;
  //    bn2[j1] = 0;
 // }
  for (j1 = bcoun; j1 < L; ++j1) {
    for (j2 = bcoun+1; j2 < L; ++j2) {
      bn1[k] = j1;
      bn2[k] = j2;
      k += 1;
    }
    bcoun += 1;
  }
  int a[bn],Spin[L],clu[L];
  double dist[bn],bp[bn];
  gslran = gsl_rng_alloc(gsl_rng_mt19937);
  srand(time(NULL));
  int rseed = rand() % 1234567890 + 1;
  gsl_rng_set(gslran, rseed);
  
  tde = (int)((ti-tf)/dt)+1;
  Jtot=0;
  for (i = 0; i < bn; ++i) a[i] = -1;
  for (i = 0; i < L-1; ++i) {
    dist[i] = pow((double) i+1, -(1+sg));
//    double hz1 =  gsl_sf_hzeta((double)1+sg,(double)i/L);
//    double hz2 =  gsl_sf_hzeta((double)1+sg,(double)(L-i)/L);
//    dist [i] = (hz1+hz2)*pow((double)L,-(1+sg));
    Jtot += (L-i-1)*dist[i];
//    Jtot += dist[i];
  }
//  Jtot = Jtot*L/2;
//   
  alias_tables (a,bp,Jtot,dist,bn1,bn2);
//  printf("%3.4f \n",Jtot);
  for (it =0; it < tde; ++it) {
//    beta = 0.0476161;
//    T=1./beta;
    T=ti-it*dt;
    beta=1/T;
    lattice (Spin);
//   clock_t start = clock();
    todo (Spin,bp,a,beta,Jtot,clu,bn1,bn2,dist);
//    double u = (((double)clock() - start) / CLOCKS_PER_SEC);
//  printf("Time elapsed for L = %d: %5.9f s\n", L,u);
  }
  return L;
}
