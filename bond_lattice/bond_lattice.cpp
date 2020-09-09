#include "bond_lattice.hpp"
extern "C" {

  void simpos3D(double a0, double am, int N, double AL, double D, double RT, double T,
    double min_r, double max_r, int bins, int *H, double *V, int *CH,
    double *r, unsigned seed, int tsample, int ttherm, int rank, int LatticeType,
    int CorrelationType, bool CentralForce) {
    if(rank==0) {
        std::cout<<"MORSE LATTICE SYMM:\n";
        std::cout<<"LatticeType: "<<LatticeType<<"\nCentralForce: "<<CentralForce<<"\n";
        std::cout<<"kT/[Energy]: "<<T<<"\na/a0:"<<am<<"\n";
        std::cout<<"CorrelationType: "<<CorrelationType<<std::endl;
    }
    std::default_random_engine generator(seed);
    std::normal_distribution<double> eta(0.0,1.0);

    double a = a0*am; // target lattice constant

    int i,t,j,ni,nnc,npb,nnc1,nnc2,npb1,npb2,noverlaps;
    double *overlaps;


    if(LatticeType==1) { // bcc 111
      // 8 neighbors in 1st shell
      nnc1 = 8;
      // 6 neighbors in 2nd shell => 0 as not counted
      nnc2 = 0;
      // I am looking for 2 dot products: +/-1/3
      noverlaps = 2;
      overlaps = new double[noverlaps];
      overlaps[0]=+1.0/3.0;
      overlaps[1]=-1.0/3.0;
      // There are three (of 8/2-1=3) bonds that give each for the first set
      npb1 = 3;
      // There are three (of 6/2=3) bonds that give each for the 2nd set
      npb2 = 0;
    } else if(LatticeType==2) { // bcc 111 + 100 shells
      // 8 neighbors in 1st shell
      nnc1 = 8;
      // 6 neighbors in 2nd shell
      nnc2 = 6;
      // I am looking for 4 dot products: +/-1/3 and +/-1/rt(3)
      noverlaps = 4;
      overlaps = new double[noverlaps];
      overlaps[0]=+1.0/3.0;
      overlaps[1]=-1.0/3.0;
      overlaps[2]=+1.0/sqrt(3.0);
      overlaps[3]=-1.0/sqrt(3.0);

      // There are three (of 8/2-1=3) bonds that give each for the first set
      npb1 = 3;
      // There are three (of 6/2=3) bonds that give each for the 2nd set
      npb2 = 3;
    } else { // fcc 110
      // 12 neighbors in 1st shell
      nnc1 = 12;
      // I am looking for 2 dot products: +/-1/2
      noverlaps = 2; // number of inner products
      overlaps = new double[noverlaps];
      overlaps[0]=+1.0/2.0;
      overlaps[1]=-1.0/2.0;
      // There are four (of 12/2-1=5) bonds that give each value => 2*4 = 8 bond pair indices for each bond
      npb1 = 4;
      // only 1st shell
      nnc2 = 0;
      npb2 = 0;
    }

    nnc = nnc1 + nnc2; // total number of bonds
    npb = npb1 + npb2; // total number of projected bonds

    int niv[3*nnc]; // 3-index shift from a site to its neighbor
    for(int k=0;k<nnc*3;k++) niv[k]=0;

    double nnv[3*nnc];  // bond vectors
    for(int k=0;k<nnc*3;k++) nnv[k]=0.0;

    double nm[nnc]; // bond vector magnitudes (modified below)
    for(int k=0;k<nnc;k++) nm[k]  = a;

    int bli[nnc], bti[nnc]; // bond long/trans index
    for(int k=0;k<nnc;k++) {
      bli[k] = 0;
      bti[k] = 0;
    }

    int oli[npb*noverlaps*nnc]; // overlap index
    for(int k=0;k<npb*noverlaps*nnc;k++) oli[k]=0;

    double dp[nnc]; // overlap vector
    for(int k=0;k<nnc;k++) dp[k]=0.0;

    // define lattice
    if(LatticeType>0) {
      a *= 2.0/sqrt(3.);
      // [-111]/2 = 100
      niv[3*0+0] = +1;  niv[3*0+1] = +0;   niv[3*0+2] = +0;
      nnv[3*0+0] = -a/2.0;  nnv[3*0+1] = +a/2.0;  nnv[3*0+2] = +a/2.0;

      // [1-1-1]/2 = -100
      niv[3*1+0] = -1;      niv[3*1+1] = +0;      niv[3*1+2] = +0;
      nnv[3*1+0] = +a/2.0;  nnv[3*1+1] = -a/2.0;  nnv[3*1+2] = -a/2.0;

      // [1-11]/2 = 010
      niv[3*2+0] = +0;     niv[3*2+1] = +1;       niv[3*2+2] = +0;
      nnv[3*2+0] = +a/2.0; nnv[3*2+1] = -a/2.0;   nnv[3*2+2] = +a/2.0;

      // [-11-1]/2 = 0-10
      niv[3*3+0] = +0;      niv[3*3+1] = -1;      niv[3*3+2] = +0;
      nnv[3*3+0] = -a/2.0;  nnv[3*3+1] = +a/2.0;  nnv[3*3+2] = -a/2.0;

      // [11-1]/2 = 001
      niv[3*4+0] = +0;      niv[3*4+1] = +0;      niv[3*4+2] = +1;
      nnv[3*4+0] = +a/2.0;  nnv[3*4+1] = +a/2.0;  nnv[3*4+2] = -a/2.0;

      // [-1-11]/2 = -001
      niv[3*5+0] = +0;      niv[3*5+1] = +0;      niv[3*5+2] = -1;
      nnv[3*5+0] = -a/2.0;  nnv[3*5+1] = -a/2.0;  nnv[3*5+2] = +a/2.0;

      // [111]/2 = -111 + 1-11 + 11-1 = 1 1 1
      niv[3*6+0] = +1;      niv[3*6+1] = +1;      niv[3*6+2] = +1;
      nnv[3*6+0] = +a/2.0;  nnv[3*6+1] = +a/2.0;  nnv[3*6+2] = +a/2.0;

      // -[111]/2 = 1-1-1 + -11-1 + -1-11 = -1 -1 -1
      niv[3*7+0] = -1;      niv[3*7+1] = -1;      niv[3*7+2] = -1;
      nnv[3*7+0] = -a/2.0;  nnv[3*7+1] = -a/2.0;  nnv[3*7+2] = -a/2.0;


      if(LatticeType==2) {
        // [100] = [1-1-1]/2 + [111]/2 = 0 1 1
        niv[3*8+0] = +0;  niv[3*8+1] =  +1;  niv[3*8+2] =  +1;
        nnv[3*8+0] = +a;  nnv[3*8+1] = 0.0;  nnv[3*8+2] = 0.0;
        nm[8] = a;

        // [-100] = -[1-1-1]/2 - [111]/2 = 0 -1 -1
        niv[3*9+0] = +0;  niv[3*9+1] =  -1;  niv[3*9+2] =  -1;
        nnv[3*9+0] = -a;  nnv[3*9+1] = 0.0;  nnv[3*9+2] = 0.0;
        nm[9] = a;

        // [010] = [-11-1]/2 + [111]/2 = 1 0 1
        niv[3*10+0] =  +1;  niv[3*10+1] =  +0;  niv[3*10+2] =  +1;
        nnv[3*10+0] = 0.0;  nnv[3*10+1] =  +a;  nnv[3*10+2] = 0.0;
        nm[10] = a;

        // [-100] = -[-11-1]/2 - [111]/2 = -1 0 -1
        niv[3*11+0] = -1;   niv[3*11+1] =  +0;  niv[3*11+2] =  -1;
        nnv[3*11+0] = 0.0;  nnv[3*11+1] = -a;   nnv[3*11+2] = 0.0;
        nm[11] = a;

        // [001] = [-1-11]/2 + [111]/2 = 1 1 0
        niv[3*12+0] =  +1;  niv[3*12+1] =  +1;  niv[3*12+2] =  +0;
        nnv[3*12+0] = 0.0;  nnv[3*12+1] =  0.0;  nnv[3*12+2] = +a;
        nm[12] = a;

        // [00-1] = -[1-1-1]/2 - [111]/2 = -1 -1 0
        niv[3*13+0] = -1;   niv[3*13+1] =  -1;  niv[3*13+2] =  +0;
        nnv[3*13+0] = 0.0;  nnv[3*13+1] = 0.0;   nnv[3*13+2] = -a;
        nm[13] = a;
      }
      a *= sqrt(3.)/2.0;

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

    // could hard code this, but find all pairs with dp == overlap
    // noverlaps == 2*nshells
    for(i=0;i<nnc;i++) {
      for(j=0;j<nnc;j++) {
        dp[j] = nnv[3*i+0]*nnv[3*j+0]/nm[i]/nm[j];
        dp[j] += nnv[3*i+1]*nnv[3*j+1]/nm[i]/nm[j];
        dp[j] += nnv[3*i+2]*nnv[3*j+2]/nm[i]/nm[j];
      }
      // npb values for each overlap, npb*noverlaps values for each bond

      for(int k=0;k<noverlaps;k++) {
        t=0;
        for(j=0;j<nnc;j++) if(std::fabs(dp[j]-overlaps[k])<0.01) {
          oli[noverlaps*npb*i+t+k*npb1]=j;
          t++;
        }
      }
    }


    double dx[N*N*N*3],fv[N*N*N*3];

    double dr  = 0.5*(max_r-min_r)/bins;

    for(i=0;i<3*N*N*N;i++) dx[i] = 0.0;

    for(i=0;i<bins;i++) {
      for(j=0;j<4;j++) H[i+j*bins] = 0;

      if(CorrelationType==1) for(j=0;j<4*bins;j++) CH[i*bins+j] = 0;

      r[i] = min_r+i*(max_r-min_r)/bins;
      V[i] = energy(D,AL,nm[0]/am,r[i]+dr);
    }

    for(i=0;i<6;i++) V[bins+i]=0.0;

    int ib,ib2,ib3,ib4,ib_shift,correction_count=0,ix,iy,iz,nix,niy,niz;
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
          // record indices of long and trans bonds for correlation analysis.
          // Initialize to -1, only use later if >=0
          if(CorrelationType>0) {
            bli[j] = -1;
            bti[j] = -1;
          }

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

          // this is unit a_0, i.e. nv = a_0 / |a_0|, where |a_0| = am
          nv[0] = nnv[3*j+0]/nm[j];
          nv[1] = nnv[3*j+1]/nm[j];
          nv[2] = nnv[3*j+2]/nm[j];

          // this is x_ni - x_i
          abv[0] = bv[0] + nnv[3*j+0];
          abv[1] = bv[1] + nnv[3*j+1];
          abv[2] = bv[2] + nnv[3*j+2];

          bmd = sqrt(bv[0]*bv[0]+bv[1]*bv[1]+bv[2]*bv[2]); // |b-a|
          bm = sqrt(abv[0]*abv[0]+abv[1]*abv[1]+abv[2]*abv[2]); // |b|
          blm = abv[0]*nv[0]+abv[1]*nv[1]+abv[2]*nv[2]; // b.a

          blmax = std::max(blm,blmax);
          blmin = std::min(blm,blmin);

          // transverse vector
          bt[0] = abv[0]-blm*nv[0];
          bt[1] = abv[1]-blm*nv[1];
          bt[2] = abv[2]-blm*nv[2];
          btm = sqrt(bt[0]*bt[0]+bt[1]*bt[1]+bt[2]*bt[2]); // |b x a|

          // histogram bin indices
          ib =  int(bins*((bm -min_r)/(max_r-min_r))); // magnitude
          ib2 = int(bins*((blm-min_r)/(max_r-min_r))); // longitudinal
          ib3 = int(bins*((bmd-min_r)/(max_r-min_r))); // difference
          ib4 = int(bins*(btm/(max_r-min_r))); // transverse, no offset

          ib_shift=4*bins*int(j>=nnc-nnc2); // if j>=nnc-nnc2 we are in 2nd shell

          // multiple bond types
          if(t>=ttherm && ib>=0 && ib<bins)   H[ib +0*bins + ib_shift]++; // mag, double counting
          if(t>=ttherm && ib2>=0 && ib2<bins) H[ib2+1*bins + ib_shift]++; // proj, double counting
          if(t>=ttherm && ib3>=0 && ib3<bins) H[ib3+2*bins + ib_shift]++; // diff, double counting
          if(t>=ttherm && ib4>=0 && ib4<bins) H[ib4+3*bins + ib_shift]++; // trans, double counting

          // Correlation Data
          if(t>=ttherm && ib2>=0 && ib2<bins) bli[j] = ib2;
          if(t>=ttherm && ib4>=0 && ib4<bins) bti[j] = ib4;

          // force data
          if(CentralForce) {
            el = energy(D,AL,nm[j]/am,bm)-energy(D,AL,nm[j]/am,nm[j]);
            et = RT*D*AL*AL*btm*btm;
            fm = denergy(D,AL,nm[j]/am,bm);
            fv[3*i+0] += abv[0] / bm * fm - 2.0*RT*D*AL*AL*bt[0];
            fv[3*i+1] += abv[1] / bm * fm - 2.0*RT*D*AL*AL*bt[1];
            fv[3*i+2] += abv[2] / bm * fm - 2.0*RT*D*AL*AL*bt[2];
            fm = sqrt(fm*fm+RT*D*AL*AL*et*4.0);
          } else {
            el = energy_long(D,AL,RT,nm[j]/am,blm)-energy_long(D,AL,RT,nm[j]/am,nm[j]);
            et = energy_tran(D,AL,RT,nm[j]/am,btm);
            fm = force(D,AL,RT,nm[j]/am,blm,btm,flt);
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

        if(CorrelationType>0 && t>=ttherm) {
          // bli : bond longitudinal count
          // bti : bond transverse count
          // npb : size of overlap dataset for each bond (i.e. 4 for fcc)
          // oli : overlap index

          // transverse, 1st shell
          for(j=0;j<nnc1;j++) {
            if(bli[2*j]<0 || bti[2*j]<0) continue;
            CH[bli[j]*bins+bti[j]]++;
          }

          // Opposite, 1st shell
          for(j=0;j<nnc1/2;j++) {
            if(bli[2*j]<0 || bli[2*j+1]<0) continue;
            CH[bins*bli[2*j]+bli[2*j+1]+bins*bins]++;
            //CH[bins*bli[2*j+1]+bli[2*j]+bins*bins]++;
          }

          // 2 overlaps in 1st shell, long-long
          for(j=0;j<nnc1;j++) for(int k=0;k<npb1;k++) { // 1st shell
            ib = bli[j];
            ib2 = bli[ oli[noverlaps*npb*j+k] ]; // index of over bond is oli[noverlaps*npb*j+k]
            if(ib<0 || ib2 <0) continue;
            CH[ ib*bins + ib2 + 2*bins*bins]++;
            //CH[ ib2*bins + ib + 2*bins*bins]++;
          }

          for(j=0;j<nnc1;j++) for(int k=npb1;k<2*npb1;k++)  { // 1st shell
            ib = bli[j];
            ib2 = bli[ oli[noverlaps*npb*j+k] ]; // index of over bond is oli[noverlaps*npb*j+k]
            if(ib<0 || ib2 <0) continue;
            CH[ ib*bins + ib2 + 3*bins*bins]++;
            //CH[ ib2*bins + ib + 3*bins*bins]++;
          }

          // 2 overlaps in 1st/2nd shell, long-long
          for(j=nnc1;j<nnc;j++) for(int k=2*npb1;k<2*npb1+npb2;k++) { // 1st shell
            ib = bli[j];
            ib2 = bli[ oli[noverlaps*npb*j+k] ]; // index of over bond is oli[noverlaps*npb*j+k]
            if(ib<0 || ib2 <0) continue;
            CH[ ib*bins + ib2 + 4*bins*bins]++;
            //CH[ ib2*bins + ib + 4*bins*bins]++;
          }
          for(j=nnc1;j<nnc;j++) for(int k=2*npb1+npb2;k<2*npb1+2*npb2;k++) { // 1st shell
            ib = bli[j];
            ib2 = bli[ oli[noverlaps*npb*j+k] ]; // index of over bond is oli[noverlaps*npb*j+k]

            if(ib<0 || ib2 <0) continue;
            CH[ ib*bins + ib2 + 5*bins*bins]++;
            //CH[ ib2*bins + ib + 5*bins*bins]++;
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
