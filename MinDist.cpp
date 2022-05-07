//g++-5 MinDist.cpp -w `root-config --cflags` -I$FEDRA_ROOT/include -L$FEDRA_ROOT/lib -lEIO -lDataConversion -lEdb -lEbase -lEdr -lScan -lAlignment -lEmath -lEphys -lvt -lDataConversion `root-config --libs` -o minDist
//./minDist

#include<stdio.h>
#include<EdbDataSet.h>
#include <vector>
#include <cmath>

using namespace std;

main(int argc, char *argv[])
{
  EdbDataProc *dproc = new EdbDataProc("lnk.def");

  dproc->InitVolume( 100, "nseg>=3&&abs(t.eTX)<0.4&&abs(t.eTY)<0.4");

  EdbPVRec *pvr = dproc->PVR();

  int ntrk = pvr->Ntracks();

  TObjArray *trObjArr = new TObjArray;
  for (int itrk = 0; itrk < ntrk; itrk++)
  {
    EdbTrackP *t = pvr->GetTrack(itrk);
    trObjArr->Add(t);
  }

  EdbTrackP *refTrack = new EdbTrackP(8);
  EdbTrackP *secTrack = new EdbTrackP(8);

  int ifl = 0;
  double p1[3], p2[3], v1[3], v2[3], vt[3];
  double sum, sig;
  double slopy1,slopy2,slopz1,slopz2;
  double chi_tr1,chi_tr2;
  vtx[5]=0.0;
  vtx[6]=0.0;
  vtx[3]=0.0;
  vtx[2]=0.0;
  vtx[1]=0.0;

  int trObjArrSize = trObjArr->GetEntriesFast();

  for(int i = 0; i < trObjArrSize; i++)
  {
    refTrack = (EdbTrackP*)(trObjArr->At(i));

    p1[1]  = refTrack->X();
    p1[2]  = refTrack->Y();
    p1[0]  = refTrack->Z();

    slopy1 = refTrack->TX();
    slopz1 = refTrack->TY();

    v1[0]  = 1./sqrt(1.+slopy1*slopy1+slopz1*slopz1);
    v1[1]  = slopy1*v1[0];
    v1[2]  = slopz1*v1[0];

    for(int j = 0; j < trObjArrSize; j++)
    {
      secTrack = (EdbTrackP*)(trObjArr->At(j));

      p2[1]  = secTrack->X();
      p2[2]  = secTrack->Y();
      p2[0]  = secTrack->Z();

      slopy2 = secTrack->TX();
      slopz2 = secTrack->TY();

      v2[0]  = 1./sqrt(1.+slopy2*slopy2+slopz2*slopz2);
      v2[1]  = slopy2*v2[0];
      v2[2]  = slopz2*v2[0];

      vt[0]  = v1[1]*v2[2]-v1[2]*v2[1];
      vt[1]  = v1[2]*v2[0]-v1[0]*v2[2];
      vt[2]  = v1[0]*v2[1]-v1[1]*v2[0];
    }

    double vtabs  = sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);

    if (vtabs == 0) 
    {
      ifl = -1;
      return(-1);
    }

    vt[0]  = vt[0]/vtabs;
    vt[1]  = vt[1]/vtabs; 
    vt[2]  = vt[2]/vtabs;

    double vtv1   = p1[0]*vt[0]+p1[1]*vt[1]+p1[2]*vt[2];
    double vtv2   = p2[0]*vt[0]+p2[1]*vt[1]+p2[2]*vt[2];
    double v1v2   = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

    sig  = (vtv2-vtv1)/2.;

    if (v1v2 == 1) 
    { 
      ifl = -1;
      return(-1);
    }

    sum = 0;  
    for(int i = 0; i < 3; i++) { sum = sum + (p2[i]-p1[i])*(v1[i]-v2[i]*v1v2);}

    double faca = sum/(1.-v1v2*v1v2);
    
    for(int j=0;j<3;j++) { vtx[j] = p1[j] + faca*v1[j] + sig*vt[j];}
    
    vtx[3] = sqrt(sig*sig); 
    vtx[4]= 0.0;
    vtx[5]= vtx[0]-p1[0];
    vtx[6]= vtx[0]-p2[0];

    if(vtx[3] < 2.0 && vtx[5] > -1580.0 && vtx[6] > -1580.0 && vtx[5] < 100.0 && vtx[6] < 100.0 && sqrt((vtx[5]-vtx[6])*(vtx[5]-vtx[6])) < 800.)
    {
      if((chi_tr1 < 1.0 && chi_tr2 < 1.0) || (Two_seg[0] > 3 && Two_seg[1] > 3)) { vtx[4]=1.0; }
    }

  }

  /*
  //////////////////
  //---> SB TRACK///
  //////////////////
  p1[1]  = TSB->myseg[TSB->nseg-1]->x;
  p1[2]  = TSB->myseg[TSB->nseg-1]->y;
  p1[0]  = TSB->myseg[TSB->nseg-1]->z;
  slopy1 = TSB->ax;
  slopz1 = TSB->ay;
  chi_tr1= TSB->chi2/(double)TSB->ndf;
  v1[0]  = 1./sqrt(1.+slopy1*slopy1+slopz1*slopz1);
  v1[1]  = slopy1*v1[0];
  v1[2]  = slopz1*v1[0];
  //////////////////////
  //---> SECOND TRACK///
  //////////////////////
  p2[1]  = TTR->myseg[TTR->nseg-1]->x;
  p2[2]  = TTR->myseg[TTR->nseg-1]->y;
  p2[0]  = TTR->myseg[TTR->nseg-1]->z;
  chi_tr2= TTR->chi2/(double)TTR->ndf;
  slopy2 = TTR->ax;
  slopz2 = TTR->ay;
  v2[0]  = 1./sqrt(1.+slopy2*slopy2+slopz2*slopz2);
  v2[1]  = slopy2*v2[0];
  v2[2]  = slopz2*v2[0];
  /////////////////////////
  //Book Track parameters//
  /////////////////////////
  Two_sl[0]= slopy1; 
  Two_sl[1]= slopz1;
  Two_sl[2]= slopy2;
  Two_sl[3]= slopz2;
  Two_seg[0]= TSB->nseg;
  Two_seg[1]= TTR->nseg;
  ///////////////////////////////////////////////////
  //---> VECTOR WHICH IS TRANSVERSE TO BOTH TRACKS///
  ///////////////////////////////////////////////////
  vt[0]  = v1[1]*v2[2]-v1[2]*v2[1];
  vt[1]  = v1[2]*v2[0]-v1[0]*v2[2];
  vt[2]  = v1[0]*v2[1]-v1[1]*v2[0];
  double vtabs  = sqrt(vt[0]*vt[0]+vt[1]*vt[1]+vt[2]*vt[2]);
  if (vtabs==0.) {
    ifl = -1;
    return(-1);
  }
  vt[0]  = vt[0]/vtabs;
  vt[1]  = vt[1]/vtabs; 
  vt[2]  = vt[2]/vtabs;
  //////////////////////////////////////////////
  //---> MINIMAL HALF DISTANCE BETWEEN TRACKS///
  //////////////////////////////////////////////
  double vtv1   = p1[0]*vt[0]+p1[1]*vt[1]+p1[2]*vt[2];
  double vtv2   = p2[0]*vt[0]+p2[1]*vt[1]+p2[2]*vt[2];
  sig  = (vtv2-vtv1)/2.;
  //
  double v1v2   = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  if (v1v2==1) { 
    ifl = -1;
    return(-1);
  }
  //
  //start from track 1
  //
  sum = 0; 
  for(int i=0;i<3;i++) { 
    sum = sum + (p2[i]-p1[i])*(v1[i]-v2[i]*v1v2);
  }
  double faca = sum/(1.-v1v2*v1v2);
  for(int j=0;j<3;j++){   
    vtx[j] = p1[j] + faca*v1[j] + sig*vt[j];
  }
  vtx[3] = sqrt(sig*sig); 
  vtx[4]= 0.0;
  vtx[5]= vtx[0]-p1[0];
  vtx[6]= vtx[0]-p2[0];
  //
  if(vtx[3]<2.0&&vtx[5]>-1580.0&&vtx[6]>-1580.0&&vtx[5]<100.0&&vtx[6]<100.0&&sqrt((vtx[5]-vtx[6])*(vtx[5]-vtx[6]))<800.){
    if((chi_tr1<1.0&&chi_tr2<1.0)||(Two_seg[0]>3&&Two_seg[1]>3)){
      vtx[4]=1.0;
    }
  }
  */
  return 0;
}