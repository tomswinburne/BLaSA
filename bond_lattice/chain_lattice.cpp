#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <random>

/*
  ex = exp(-al*(b-a0))
  dex = -al*ex, ddex = al*al*ex

  E = D*( ex*ex -2.0*ex + 1.0)
  F = 2.0*D*al*ex*(ex-1.0)
  C = = -dF = 2.0*D*al*al*(2.0*ex-1)*ex
  C = 0 => ex = 0.5
  ex=1 => C = 2*D*al*al
*/
#define LOGTWO 0.6931471805599453
#define PI 3.141592654
#define RTTWO 1.4142135623730951

double denergy(double d, double al, double a0, double b) {
  /*
  d,al,a0 = Morse parameters
  b = bond magnitude
  */
  return 2.0*d*al*exp(-al*(b-a0))*(1.0-exp(-al*(b-a0))); // -> 2.0*D*AL*AL*(b-a0)
};

double energy(double d, double al, double a0, double b) {
  return d*(1.0-exp(-al*(b-a0)))*(1.0-exp(-al*(b-a0))); // -> D*AL*AL*(b-a0)^2
};



/*
d,al,a0 = Morse parameters
a = lattice vector for bond
btm = transverse magnitude
blm = longitudinal magnitude
*/

double energy(double d, double al, double rt, double a0, double blm, double btm) {
  // longitudinal energy
  double el = d * ( 1.0 - exp(-al*(blm-a0)) ) * ( 1.0 - exp(-al*(blm-a0)) );
  // transverse energy
  double et = rt*d*al*al*btm*btm;
  return el + et;
};

double energy_long(double d, double al, double rt, double a0, double blm) {
  // longitudinal energy
  return d * ( 1.0 - exp(-al*(blm-a0)) ) * ( 1.0 - exp(-al*(blm-a0)) );
};

double energy_tran(double d, double al, double rt, double a0, double btm) {
  // transverse energy
  return rt*d*al*al*btm*btm;
};


double force(double d, double al, double rt, double a0, double blm, double btm, double *f) {
  // force longitudinal
  f[0] = 2.0 * d * al * exp(-al*(blm-a0)) * (exp(-al*(blm-a0))-1.0); // -> -2.0*d*al*al*db
  // force traditional
  f[1] = -2.0 * rt * d * al * al;

  return sqrt(f[0]*f[0]+f[1]*f[1]*btm*btm);
};


extern "C" {

  void simpos3D(double a0, double am, int N, double AL, double D, double RT, double T,
    double min_r, double max_r, int bins, int *H, double *V, int *CH,
    double *r, unsigned seed, int tsample, int ttherm, int rank, bool fcc,
    int CorrelationType, bool CentralForce) {
    if(rank==0) {
        std::cout<<"MORSE LATTICE SYMM:\n";
        std::cout<<"FCC: "<<fcc<<"\nCentralForce: "<<CentralForce<<"\n";
        std::cout<<"kT/[Energy]: "<<T<<"\na/a0:"<<am<<"\n";
        std::cout<<"CorrelationType: "<<CorrelationType<<std::endl;
    }
    std::default_random_engine generator(seed);
    std::normal_distribution<double> eta(0.0,1.0);
    int i,t,j,ni,nnc=6;
    if(fcc) nnc = 12;
    int niv[3*nnc];
    double nnv[3*nnc],dp;
    for(int k=0;k<nnc*3;k++) {
      niv[k]=0;
      nnv[k]=0.0;
    }

    int blc[nnc],blp[4*nnc]; // for correlations


    double a = a0*am; // target lattice constant

    if(!fcc) {
      // simple cubic- only nn interactions
      // 100
      niv[3*0+0] = +1;  niv[3*0+1] = +0;   niv[3*0+2] = +0;
      nnv[3*0+0] = a;  nnv[3*0+1] = 0.0;  nnv[3*0+2] = 0.0;

      // -100
      niv[3*1+0] = -1;  niv[3*1+1] = +0;  niv[3*1+2] = +0;
      nnv[3*1+0] = -a; nnv[3*1+1] = 0.0; nnv[3*1+2] = 0.0;

      // 010
      niv[3*2+0] = +0;  niv[3*2+1] = +1;   niv[3*2+2] = +0;
      nnv[3*2+0] = 0.0; nnv[3*2+1] = a;   nnv[3*2+2] = 0.0;

      // -010
      niv[3*3+0] = +0;  niv[3*3+1] = -1;  niv[3*3+2] = +0;
      nnv[3*3+0] = 0.0; nnv[3*3+1] = -a; nnv[3*3+2] = 0.0;

      // 001
      niv[3*4+0] = +0;  niv[3*4+1] = +0;  niv[3*4+2] = +1;
      nnv[3*4+0] = 0.0; nnv[3*4+1] = 0.0; nnv[3*4+2] = a;

      // -001
      niv[3*5+0] = +0;  niv[3*5+1] = +0;  niv[3*5+2] = -1;
      nnv[3*5+0] = 0.0; nnv[3*5+1] = 0.0; nnv[3*5+2] = -a;
    } else {
      a *= sqrt(2.0); // a -> a'
      // fcc x,y,z =  a'/2[110],a'/2[011],a'/2[101] , a' = a*sqrt(2)

      //-----
      // [110] == 100
      niv[3*0+0] = +1;      niv[3*0+1] = +0;     niv[3*0+2] = +0;
      nnv[3*0+0] =  a/2.0; nnv[3*0+1] = a/2.0; nnv[3*0+2] = 0.0;
      // [-1-10] == -100
      niv[3*1+0] = -1;      niv[3*1+1] = +0;      niv[3*1+2] = +0;
      nnv[3*1+0] = -a/2.0; nnv[3*1+1] = -a/2.0; nnv[3*1+2] = 0.0;
      // [1-10] = [101]-[011] == 0-11
      niv[3*2+0] = 0;       niv[3*2+1] = -1;      niv[3*2+2] = +1;
      nnv[3*2+0] =  a/2.0; nnv[3*2+1] = -a/2.0; nnv[3*2+2] = 0.0;
      // [-110] == 01-1
      niv[3*3+0] = 0;       niv[3*3+1] = 1;      niv[3*3+2] = -1;
      nnv[3*3+0] = -a/2.0; nnv[3*3+1] = a/2.0; nnv[3*3+2] = 0.0;

      //-----
      // [011] == 010
      niv[3*4+0] = 0;    niv[3*4+1] = +1;     niv[3*4+2] = +0;
      nnv[3*4+0] =  0.0; nnv[3*4+1] = a/2.0; nnv[3*4+2] = a/2.0;
      // [0-1-1] == 0-10
      niv[3*5+0] = 0;   niv[3*5+1] = -1;      niv[3*5+2] = +0;
      nnv[3*5+0] = 0.0; nnv[3*5+1] = -a/2.0; nnv[3*5+2] = -a/2.0;
      // [0-11] = [101]-[110] = -101
      niv[3*6+0] = -1;   niv[3*6+1] = 0;       niv[3*6+2] = 1;
      nnv[3*6+0] =  0.0; nnv[3*6+1] = -a/2.0; nnv[3*6+2] = a/2.0;
      // [01-1] == 10-1
      niv[3*7+0] = 1;   niv[3*7+1] = +0;     niv[3*7+2] = -1;
      nnv[3*7+0] = 0.0; nnv[3*7+1] = a/2.0; nnv[3*7+2] = -a/2.0;

      //-----
      // [101] == 001
      niv[3*8+0] = 0;       niv[3*8+1] = +0;  niv[3*8+2] = +1;
      nnv[3*8+0] =  a/2.0; nnv[3*8+1] = 0.0; nnv[3*8+2] = a/2.0;
      // [-10-1] == 00-1
      niv[3*9+0] = 0;       niv[3*9+1] = +0;  niv[3*9+2] = -1;
      nnv[3*9+0] = -a/2.0; nnv[3*9+1] = 0.0; nnv[3*9+2] = -a/2.0;
      // [10-1] = [110] - [011] == 1-10
      niv[3*10+0] = +1;      niv[3*10+1] = -1;  niv[3*10+2] = +0;
      nnv[3*10+0] =  a/2.0; nnv[3*10+1] = 0.0; nnv[3*10+2] = -a/2.0;
      // [-101] = -[110] + [011] == -110
      niv[3*11+0] = -1;      niv[3*11+1] = +1;  niv[3*11+2] = +0;
      nnv[3*11+0] = -a/2.0; nnv[3*11+1] = 0.0; nnv[3*11+2] = a/2.0;

      a /= sqrt(2.0);  // a' -> a
    }

    // could hard code this, but find all pairs with dp == a^2/2
    if(CorrelationType>1) {
      for(i=0;i<nnc;i++) {
        t=0;
        for(j=0;j<nnc;j++) {
          dp = nnv[3*i+0]*nnv[3*j+0]/a/a;
          dp += nnv[3*i+1]*nnv[3*j+1]/a/a;
          dp += nnv[3*i+2]*nnv[3*j+2]/a/a;
          if(dp>0.49 && dp<0.51) {
            blp[4*i+t] = j;
            t++;
          }
          if(t==4) break;
        }
      }
    }

    double dx[N*N*N*3],fv[N*N*N*3];

    double dr  = 0.5*(max_r-min_r)/bins;

    for(i=0;i<3*N*N*N;i++) dx[i] = 0.0;

    for(i=0;i<bins;i++) {
      for(j=0;j<4;j++) H[i+j*bins] = 0;
      if(CorrelationType>0) // 1: bl_0,bt_0, 2: bl_0,bt_0,bl_0,bl_6, bl_0,bl_1
        for(j=0;j<bins*(1+2*(CorrelationType-1));j++)
          CH[i*bins+j] = 0;
      r[i] = min_r+i*(max_r-min_r)/bins;
      V[i] = energy(D,AL,a0,r[i]+dr);
    }

    for(i=0;i<6;i++) V[bins+i]=0.0;
    //if(rank==0)
    //  std::cout<<"V_bond(<b>)/[Energy]: "<<energy(D,AL,RT,a0,a,0.0)<<std::endl;

    int ib,ib2,ib3,ib4,correction_count=0,ix,iy,iz,nix,niy,niz;
    //double theta,phi;
    double ig = 0.001,blm,btm,bmd,cc;
    double iG = sqrt(2.0*T*ig);
    double mf=0.0,bv[3],bt[3],abv[3],nv[3],flt[2],bm,fm,et,el,blmax=0.0,blmin=2.0;

    for(t=0; t<ttherm+tsample; t++) {
      mf=0.0;
      for(i=0;i<N*N*N;i++) {
        // index i = ix + N*iy + N*N*iz

        ix = i%N;
        iy = (i/N)%N;
        iz = i/(N*N);


        fv[3*i+0] = 0.0;
        fv[3*i+1] = 0.0;
        fv[3*i+2] = 0.0;

        for(j=0;j<nnc;j++) {
          if(CorrelationType>1) blc[j] = -1;

          // bv = x[ni] - x[i] + b0
          nix = (ix+niv[3*j+0]+N)%N;
          //bc[0] = (ix+niv[3*j+0])==nix;

          niy = (iy+niv[3*j+1]+N)%N;
          //bc[1] = (iy+niv[3*j+1])==niy;

          niz = (iz+niv[3*j+2]+N)%N;
          //bc[2] = (iz+niv[3*j+2])==niz;

          ni = nix + N*niy + N*N*niz;

          // this is x_ni - x_i - a_0
          bv[0] = dx[3*ni+0]-dx[3*i+0];
          bv[1] = dx[3*ni+1]-dx[3*i+1];
          bv[2] = dx[3*ni+2]-dx[3*i+2];
          bmd = sqrt(bv[0]*bv[0]+bv[1]*bv[1]+bv[2]*bv[2]);

          // this is unit a_0, i.e. nv = a_0 / |a_0|, where |a_0| = am
          nv[0] = nnv[3*j+0]/a;
          nv[1] = nnv[3*j+1]/a;
          nv[2] = nnv[3*j+2]/a;

          // this is x_ni - x_i
          abv[0] = bv[0] + nnv[3*j+0];
          abv[1] = bv[1] + nnv[3*j+1];
          abv[2] = bv[2] + nnv[3*j+2];
          bm = sqrt(abv[0]*abv[0]+abv[1]*abv[1]+abv[2]*abv[2]);

          // longitudinal magnitude
          blm = abv[0]*nv[0]+abv[1]*nv[1]+abv[2]*nv[2];
          blmax = std::max(blm,blmax);
          blmin = std::min(blm,blmin);

          // transverse vector
          bt[0] = abv[0]-blm*nv[0];
          bt[1] = abv[1]-blm*nv[1];
          bt[2] = abv[2]-blm*nv[2];

          // transverse magnitude
          btm = sqrt(bt[0]*bt[0]+bt[1]*bt[1]+bt[2]*bt[2]);

          // histogram bin indices
          ib =  int(bins*((bm -min_r)/(max_r-min_r)));
          ib2 = int(bins*((blm-min_r)/(max_r-min_r)));
          ib3 = int(bins*((bmd-min_r)/(max_r-min_r)));
          ib4 = int(bins*(btm/(max_r-min_r))); // transverse, no offset

          if(t>=ttherm && ib>=0 && ib<bins) H[ib]++; // mag, double counting
          if(t>=ttherm && ib2>=0 && ib2<bins) {
            H[ib2+bins]++; // proj, double counting
            if(CorrelationType>1) blc[j] = ib2;
            if(CorrelationType>0 && ib4>=0 && ib4<bins) CH[ib2*bins+ib4]++;
          }
          if(t>=ttherm && ib3>=0 && ib3<bins) H[ib3+2*bins]++; // diff, double counting
          if(t>=ttherm && ib4>=0 && ib4<bins) H[ib4+3*bins]++; // trans, double counting

          /*
          abv = x_n - x_i
           -d/dx_i V(|x_n-x_i|) = dV * (x_n-x_i)/|x_n-x_i| -> 2.0*D*AL*AL*(x_n-x_i)
          */
          if(CentralForce) {
            el = energy(D,AL,a0,bm)-energy(D,AL,a0,a);
            et = RT*D*AL*AL*btm*btm;
            fm = denergy(D,AL,a0,bm);
            fv[3*i+0] += abv[0] / bm * fm - 2.0*RT*D*AL*AL*bt[0];
            fv[3*i+1] += abv[1] / bm * fm - 2.0*RT*D*AL*AL*bt[1];
            fv[3*i+2] += abv[2] / bm * fm - 2.0*RT*D*AL*AL*bt[2];
            fm = sqrt(fm*fm+RT*D*AL*AL*et*4.0);
          } else {
            el = energy_long(D,AL,RT,a0,blm)-energy_long(D,AL,RT,a0,a);
            et = energy_tran(D,AL,RT,a0,btm);
            fm = force(D,AL,RT,a0,blm,btm,flt);
            fv[3*i+0] += -nv[0]*flt[0] - bt[0]*flt[1];
            fv[3*i+1] += -nv[1]*flt[0] - bt[1]*flt[1];
            fv[3*i+2] += -nv[2]*flt[0] - bt[2]*flt[1];
          }
          if(t>=ttherm/2 && t<ttherm) {
            V[bins+0] += el / (ttherm/2) / (1.0*N*N*N) / (2.0) ;
            V[bins+1] += et / (ttherm/2) / (1.0*N*N*N) / (2.0) ;
          } else if(t>=ttherm) {
            V[bins+3] += el / tsample / (1.0*N*N*N) / (2.0);
            V[bins+4] += et / tsample / (1.0*N*N*N) / (2.0);
          }
          mf = std::max(mf,fabs(fm));
        }
        if(CorrelationType>1) {
          // l=m, l=m+7 correlation, i.e. innerproduct -> -a^2
          for(j=0;j<nnc/2;j++) if(blc[2*j]>=0 && blc[2*j+1]>=0) {
            CH[bins*blc[2*j]+blc[2*j+1]+bins*bins]++;
            CH[bins*blc[2*j+1]+blc[2*j]+bins*bins]++;
          }
          // all nnv whose inner product -> a^2/2
          for(j=0;j<nnc;j++) for(ib4=0;ib4<4;ib4++)
            if(blc[j]>=0 && blc[blp[4*j+ib4]]>=0) {
              CH[blc[blp[4*j+ib4]]+bins*blc[j]+2*bins*bins]++;
              CH[bins*blc[blp[4*j+ib4]]+blc[j]+2*bins*bins]++;
            }
         }
      }
      if(ig*mf>0.01*a) correction_count++;
      for(i=0;i<N*N*N;i++) for(j=0;j<3;j++) {
        dx[3*i+j] += ig*fv[3*i+j] + iG*eta(generator);
        if(t>=ttherm/2 && t<ttherm) V[bins+2] += -fv[3*i+j]*dx[3*i+j]/3.0/(N*N*N*ttherm/2);
        else if(t>=ttherm) V[bins+5] += -fv[3*i+j]*dx[3*i+j]/3.0/(N*N*N*tsample);
      }
      // if(rank==0 && (t-ttherm)%1000==0) std::cout<<t/1000<<"k/"<<(tsample+ttherm)/1000<<"k"<<" tsample="<<(t>ttherm)<<std::endl;
    }
    cc=0; for(i=0;i<bins;i++) cc += H[i];
    if(rank==0) {
      std::cout<<"Histgram Captured "<<int(cc)<<"/"<<tsample*N*N*N*nnc<<" Bond Values ; "<<blmin<<"< b_ll <"<<blmax<<std::endl;
      if(correction_count>0) {
        std::cout<<"Possible Large Force Counts: "<<correction_count<<"/"<<ttherm+tsample;
        std::cout<<"= "<<int(100.*correction_count/double(ttherm+tsample))<<"%, "<<std::endl;
      }
    }
  };
}
