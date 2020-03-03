#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <iostream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "Parton_Shower_Lib.hh"
#include "jet_fragmentation.hh"

#define NUMEVENTSNoteWorking 1000
#define NUMEVENTSSave 50000
#define InvertQgluon -2.0
#define InvertQquark -0.5
#define SetNumInverts 10000000
#define NumTries 10000
#define TheNumEbins 800

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

 int main( int argc, char *argv[]){

  long int RANDOMSEED = atoi(argv[1]);
  long int NUMEVENTS = atoi(argv[2]);
  double Q = atof(argv[3]);
  double Rmax = atof(argv[4]);
  long int Flavor = atoi(argv[5]);
  long int ProgramMode = atoi(argv[6]);

  long int i,j;

  char filenameInit[1000];
  char filenameLeadingJetSpectra[1000];
  char filenameEventWideSpectra[1000];
  char filenameEventWideSpectraLogBin[1000];
  clock_t begin, end;
  double time_spent;
  double **LeadingJetSpectra;
  double **EventWideSpectra;
  double **EventWideSpectraLogBin;

  begin = clock();
  printf("\n");
  printf("Random seed: %ld.\n",RANDOMSEED);

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //Initialize Random number generator
  gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);//gsl_rng_taus
  gsl_rng_set(r,RANDOMSEED);

  //Initialize values for shower
  long int NumRadii;
  double ActualJetRadii[100], ActualParams[10];
  double *params, *JetRadii;
  JetRadii = &ActualJetRadii[0];
  params = &ActualParams[0];

  /*
    We note the vales in params:
    *(params) = Minimum Energy Fraction or Zcut;
    *(params+1) = CA;
    *(params+2) = NF;
    *(params+3) = CF;
    *(params+4) = Jet Energy;
    *(params+5) = Jet Initial Opening Angle;
    *(params+6) = Minimum Mass Scale;
    */


  sprintf(filenameInit,"Initialize_Parton_Shower.txt");
  LoadShowerParams(filenameInit, params, &NumRadii, JetRadii);

  *(params+4) = Q;
  *(params+5) = Rmax;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sprintf(filenameLeadingJetSpectra, "LeadingJetSpectra_rseed%ld_Q%f_Rmax%f_ProgramMode%ld_Flavor%ld_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameEventWideSpectra, "InclusiveSpectra_rseed%ld_Q%f_Rmax%f_ProgramMode%ld_Flavor%ld_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);
  sprintf(filenameEventWideSpectraLogBin, "InclusiveSpectraLogBin_rseed%ld_Q%f_Rmax%f_ProgramMode%ld_Flavor%ld_Angle.txt",
	  RANDOMSEED,Q,Rmax,ProgramMode,Flavor);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////

  PSEmissionsList emissions, DaughterEmissions, EmissionsWithinJet;
  double t,AveT,RND,Z;
  double CurrentWTAaxis[2];

  PSemissions_init(&emissions);
  PSemissions_init(&DaughterEmissions);
  PSemissions_init(&EmissionsWithinJet);

  //This wipes out all the EmissionsList's, reseeds emissions with a parton pointed at the z axis with energy fraction 1, and flavor specified
  ReInitialize(&emissions, &DaughterEmissions, CurrentWTAaxis, Flavor, &t);

  AveT = 0;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Build Splitting Function Inversion
  printf("Computing Inversion of Splitting Functions.\n");
  double *InvertPg,*InvertPq;
  long int NumInverts = SetNumInverts;
  InvertPg = (double *) malloc( (NumInverts+1) * sizeof(double ));
  InvertPq = (double *) malloc( (NumInverts+1) * sizeof(double ));

  BuildSplitFunctionInversion(InvertPg, InvertPq, NumInverts, params, r);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Initialize Histograms

  LeadingJetSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  EventWideSpectra = (double **) malloc( (NumRadii+1) * sizeof(double *));
  EventWideSpectraLogBin = (double **) malloc( (NumRadii+1) * sizeof(double *));
  for(j=0;j<(NumRadii+1);j++){
    LeadingJetSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    EventWideSpectra[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
    EventWideSpectraLogBin[j] = (double *) malloc( ( TheNumEbins+1) * sizeof(double ));
  };

  for(j=0; j<(NumRadii+1); j++){
    for(i=0; i<(TheNumEbins+1); i++){
      LeadingJetSpectra[j][i]=0;
      EventWideSpectra[j][i]=0;
      EventWideSpectraLogBin[j][i]=0;
    };
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////


  i=0;
  while(i<NUMEVENTS){
      t=0;
      //Generate An Event
      JetFragmentation(&emissions, &DaughterEmissions, CurrentWTAaxis, JetRadii,
		       NumRadii, params, TheNumEbins, LeadingJetSpectra, EventWideSpectra,
		       EventWideSpectraLogBin, i+1, &t, InvertPg, InvertPq, NumInverts, ProgramMode, r);
      AveT = AveT + t;
      //This wipes out all the EmissionsList's, reseeds emissions with a parton pointed at the z axis with energy fraction 1, and flavor specified
      ReInitialize(&emissions, &DaughterEmissions, CurrentWTAaxis, Flavor, &t);

      if(0==((i+1)%NUMEVENTSNoteWorking)){
	printf("Working on Event %ld\n" , (i+1) );
      };
      if(0==((i+1)%NUMEVENTSSave)){
	  end = clock();
	  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
          printf("Saving Histograms at Event %ld\n",(i+1));
	  printf("Average MC time: %f\n",AveT/( (double)(i+1) ));
	  printf("Time spent so far seconds: %f\n",time_spent);
	  printf("Time spent so far minutes: %f\n",time_spent/60);
	  printf("Time spent so far hours: %f\n",time_spent/(60*60));
	  printf("Time spent per event in seconds: %f\n",time_spent/(i+1));
	  //SAVE EVENTS SO FAR


	  WriteToDiskHistgramsFRAG(filenameLeadingJetSpectra, (i+1), NumRadii,
				   TheNumEbins, LeadingJetSpectra, JetRadii, MinZBook,
				   Flavor, params, 0);
	  WriteToDiskHistgramsFRAG(filenameEventWideSpectra, (i+1), NumRadii,
				   TheNumEbins, EventWideSpectra, JetRadii, MinZBook,
				   Flavor, params, 0);
	  WriteToDiskHistgramsFRAG(filenameEventWideSpectraLogBin, (i+1), NumRadii,
				   TheNumEbins, EventWideSpectraLogBin, JetRadii, MinZBook,
				   Flavor, params, 1);

      };
    i = i + 1;
    };//WHILE




  for(j=0;j<(NumRadii+1);j++){
    free(LeadingJetSpectra[j]);
    free(EventWideSpectra[j]);
    free(EventWideSpectraLogBin[j]);
  };
  free(LeadingJetSpectra);
  free(EventWideSpectra);
  free(EventWideSpectraLogBin);

  PSemissions_free(&emissions);
  PSemissions_free(&DaughterEmissions);

  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("\n");
  printf("Time spent seconds: %f\n",time_spent);
  printf("Time spent minutes: %f\n",time_spent/60);
  printf("Time spent hours: %f\n",time_spent/(60*60));
  printf("Time spent per event in seconds: %f\n",time_spent/NUMEVENTS);

  printf("\n");
  return 0;
}

