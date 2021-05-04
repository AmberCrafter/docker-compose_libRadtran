
/*----------------------------------------------------------------------
Full solution of the multi-layer twostream equation (solar, thermal),
where each vertical layer is divided into three regions ("Tripleclouds"):
i)   one region of optically thick cloud ("ck"),
ii)  one region of optically thin cloud  ("cn"),
iii) one cloud-free region ("f");
The maximum-random overlap is assumed for regions of optically thick cloud
as well as for cloudy regions as a whole (i.e., ck + cn regions).
Author: Nina Črnivec, nina.crnivec@physik.uni-muenchen.de 
----------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 * $Id: twomaxrnd3C.c 2623 2011-12-23 10:52:38Z bernhard.mayer $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "equation.h"
#include "twostrebe.h"
#include "twomaxrnd.h"
#include "twomaxrnd3C.h"
#include "solver.h"
#include "bandec.h"

#define NFLUX 6  // number of diffuse fluxes defined at each level (Edn_ck, Edn_cn, Edn_f, Eup_ck, Eup_cn, Eup_f)

static int calc_pdn_transport_coefficients_twomaxrnd3C (int     nlev,
                                                        double* cf,
                                                        double* cf_ck,
                                                        double* cf_cn,
                                                        double* pdn_CKtCK,
                                                        double* pdn_CKtCN,
                                                        double* pdn_CKtF,
                                                        double* pdn_CNtCK,
                                                        double* pdn_CNtCN,
                                                        double* pdn_CNtF,
                                                        double* pdn_FtCK,
                                                        double* pdn_FtCN,
                                                        double* pdn_FtF);

static int calc_pup_transport_coefficients_twomaxrnd3C (int     nlev,
                                                        double* cf,
                                                        double* cf_ck,
                                                        double* cf_cn,
                                                        double* pup_CKtCK,
                                                        double* pup_CKtCN,
                                                        double* pup_CKtF,
                                                        double* pup_CNtCK,
                                                        double* pup_CNtCN,
                                                        double* pup_CNtF,
                                                        double* pup_FtCK,
                                                        double* pup_FtCN,
                                                        double* pup_FtF);

static int buildMatrixA_twomaxrnd3C (int      nlev,
                                     double   Ag,
                                     double*  ar_a11_ck,
                                     double*  ar_a11_cn,
                                     double*  ar_a11_f,
                                     double*  ar_a12_ck,
                                     double*  ar_a12_cn,
                                     double*  ar_a12_f,
                                     double*  pdn_CKtCK,
                                     double*  pdn_CKtCN,
                                     double*  pdn_CKtF,
                                     double*  pdn_CNtCK,
                                     double*  pdn_CNtCN,
                                     double*  pdn_CNtF,
                                     double*  pdn_FtCK,
                                     double*  pdn_FtCN,
                                     double*  pdn_FtF,
                                     double*  pup_CKtCK,
                                     double*  pup_CKtCN,
                                     double*  pup_CKtF,
                                     double*  pup_CNtCK,
                                     double*  pup_CNtCN,
                                     double*  pup_CNtF,
                                     double*  pup_FtCK,
                                     double*  pup_FtCN,
                                     double*  pup_FtF,
                                     double** matrixA);

static int buildVectorBsol_twomaxrnd3C (int     nlev,
                                        double  Ag,
                                        double  mu0,
                                        double* ar_a13_ck,
                                        double* ar_a13_cn,
                                        double* ar_a13_f,
                                        double* ar_a23_ck,
                                        double* ar_a23_cn,
                                        double* ar_a23_f,
                                        double* ar_S_ck,
                                        double* ar_S_cn,
                                        double* ar_S_f,
                                        double* pdn_CKtCK,
                                        double* pdn_CKtCN,
                                        double* pdn_CKtF,
                                        double* pdn_CNtCK,
                                        double* pdn_CNtCN,
                                        double* pdn_CNtF,
                                        double* pdn_FtCK,
                                        double* pdn_FtCN,
                                        double* pdn_FtF,
                                        double* vectB);

static int buildVectorBthe_twomaxrnd3C (int     nlev,
                                        double  Ag,
                                        double  Bg,
                                        double* ar_theComp1_ck,
                                        double* ar_theComp1_cn,
                                        double* ar_theComp1_f,
                                        double* ar_theComp2_ck,
                                        double* ar_theComp2_cn,
                                        double* ar_theComp2_f,
                                        double* cf,
                                        double* cf_ck,
                                        double* cf_cn,
                                        double* vectB);

static void freeMemory_twomaxrnd3C (int      nlev,
                                    double*  ar_a11_ck,
                                    double*  ar_a11_cn,
                                    double*  ar_a11_f,
                                    double*  ar_a12_ck,
                                    double*  ar_a12_cn,
                                    double*  ar_a12_f,
                                    double*  ar_a13_ck,
                                    double*  ar_a13_cn,
                                    double*  ar_a13_f,
                                    double*  ar_a23_ck,
                                    double*  ar_a23_cn,
                                    double*  ar_a23_f,
                                    double*  ar_a33_ck,
                                    double*  ar_a33_cn,
                                    double*  ar_a33_f,
                                    double*  pdn_CKtCK,
                                    double*  pdn_CKtCN,
                                    double*  pdn_CKtF,
                                    double*  pdn_CNtCK,
                                    double*  pdn_CNtCN,
                                    double*  pdn_CNtF,
                                    double*  pdn_FtCK,
                                    double*  pdn_FtCN,
                                    double*  pdn_FtF,
                                    double*  pup_CKtCK,
                                    double*  pup_CKtCN,
                                    double*  pup_CKtF,
                                    double*  pup_CNtCK,
                                    double*  pup_CNtCN,
                                    double*  pup_CNtF,
                                    double*  pup_FtCK,
                                    double*  pup_FtCN,
                                    double*  pup_FtF,
                                    double*  bb_sol,
                                    double*  bb_the,
                                    double*  ar_theComp1_ck,
                                    double*  ar_theComp1_cn,
                                    double*  ar_theComp1_f,
                                    double*  ar_theComp2_ck,
                                    double*  ar_theComp2_cn,
                                    double*  ar_theComp2_f,
                                    double*  S_ck,
                                    double*  S_cn,
                                    double*  S_f,
                                    double*  Edir_ck,
                                    double*  Edir_cn,
                                    double*  Edir_f,
                                    double*  Eup_ck,
                                    double*  Eup_cn,
                                    double*  Eup_f,
                                    double*  Edn_ck,
                                    double*  Edn_cn,
                                    double*  Edn_f,
                                    double*  bb,
                                    double*  xx,
                                    double** AA);

static int twostream_maxrand3C (double*  dtau_org_ck,
                                double*  omega0_org_ck,
                                double*  g_org_ck,  // "_ck" for optically thick cloudy region
                                double*  dtau_org_cn,
                                double*  omega0_org_cn,
                                double*  g_org_cn,  // "_cn" for optically thin cloudy region
                                double*  dtau_org_f,
                                double*  omega0_org_f,
                                double*  g_org_f,  // "_f"  for cloud-free region
                                double*  cf,
                                double   scale_cf,
                                int      nlev,
                                double   S0,
                                double   mu0,
                                double   Ag,
                                double   Bg,
                                double*  B,
                                int      delta,
                                int      flagSolar,
                                int      flagThermal,
                                double** Edir,
                                double** Edn,
                                double** Eup,
                                double** Lavg);

int twomaxrnd3C (float* dtau_cldk,
                 float* omega0_cldk,
                 float* g_cldk,
                 float* dtau_cldn,
                 float* omega0_cldn,
                 float* g_cldn,
                 float* dtau_clr,
                 float* omega0_clr,
                 float* g_clr,
                 float* cf,
                 float  scale_cf,
                 int    nlev,
                 double S0,
                 double mu0,
                 double Ag,
                 int    planck,
                 int    delta,
                 int    nzout,
                 float* zd,
                 float* temper,
                 float  btemp,
                 float  wvnmlo,
                 float  wvnmhi,
                 float* zout,
                 float* fldn,
                 float* flup,
                 float* fldir,
                 float* uavg) {
  int ilev   = 0;
  int ilyr   = 0;
  int lu     = 0;
  int status = 0;

  double* dtau_clr_org_d   = calloc (nlev - 1, sizeof (double));
  double* omega0_clr_org_d = calloc (nlev - 1, sizeof (double));
  double* g_clr_org_d      = calloc (nlev - 1, sizeof (double));

  double* dtau_cldk_org_d   = calloc (nlev - 1, sizeof (double));
  double* omega0_cldk_org_d = calloc (nlev - 1, sizeof (double));
  double* g_cldk_org_d      = calloc (nlev - 1, sizeof (double));

  double* dtau_cldn_org_d   = calloc (nlev - 1, sizeof (double));
  double* omega0_cldn_org_d = calloc (nlev - 1, sizeof (double));
  double* g_cldn_org_d      = calloc (nlev - 1, sizeof (double));

  double* cf_d = calloc (nlev - 1, sizeof (double));
  double* B    = calloc (nlev, sizeof (double));
  double *Edir = NULL, *Edn = NULL, *Eup = NULL, *Lavg = NULL;

  double taumax = 100.0;
  double Bg     = 0;  // Used only for thermal RT; unused (0.0) for solar RT;
  float  plkavg = 0;

  int flagSolar = 0;

  if (S0 > 0.0)
    flagSolar = 1;

  if (planck) {
    for (ilev = 0; ilev < nlev; ilev++) {  // level temperatures and Planck functions
      F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temper[ilev], &plkavg);
      B[ilev] = plkavg;
    }

    // surface temperature and Planck function
    F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &btemp, &plkavg);
    Bg = plkavg;
  }

  // copy float arrays to double arrays
  for (ilyr = 0; ilyr < nlev - 1; ilyr++) {
    dtau_cldk_org_d[ilyr] = dtau_cldk[ilyr];
    dtau_cldn_org_d[ilyr] = dtau_cldn[ilyr];
    dtau_clr_org_d[ilyr]  = dtau_clr[ilyr];

    // restrict layer optical thickness to 100
    if (dtau_cldk_org_d[ilyr] > taumax)
      dtau_cldk_org_d[ilyr] = taumax;

    if (dtau_cldn_org_d[ilyr] > taumax)
      dtau_cldn_org_d[ilyr] = taumax;

    if (dtau_clr_org_d[ilyr] > taumax)
      dtau_clr_org_d[ilyr] = taumax;

    omega0_cldk_org_d[ilyr] = omega0_cldk[ilyr];
    g_cldk_org_d[ilyr]      = g_cldk[ilyr];

    omega0_cldn_org_d[ilyr] = omega0_cldn[ilyr];
    g_cldn_org_d[ilyr]      = g_cldn[ilyr];

    omega0_clr_org_d[ilyr] = omega0_clr[ilyr];
    g_clr_org_d[ilyr]      = g_clr[ilyr];

    cf_d[ilyr] = cf[ilyr];
  }  //e-for

  // Call twostream_maxrand3C code:
  status = twostream_maxrand3C (dtau_cldk_org_d,
                                omega0_cldk_org_d,
                                g_cldk_org_d,
                                dtau_cldn_org_d,
                                omega0_cldn_org_d,
                                g_cldn_org_d,
                                dtau_clr_org_d,
                                omega0_clr_org_d,
                                g_clr_org_d,
                                cf_d,
                                scale_cf,
                                nlev,
                                S0,
                                mu0,
                                Ag,
                                Bg,
                                B,
                                delta,
                                flagSolar,
                                planck,
                                &Edir,
                                &Edn,
                                &Eup,
                                &Lavg);

  if (status != 0) {
    fprintf (stderr, "Error %d returned by twostream_maxrand3C()\n", status);
    return status;
  }

  // Copy results to final fields:
  for (ilev = 0; ilev < nlev; ilev++) {
    for (lu = 0; lu < nzout; lu++) {
      if (zout[lu] == zd[ilev]) {
        fldn[lu]  = Edn[ilev];
        flup[lu]  = Eup[ilev];
        fldir[lu] = Edir[ilev];
        uavg[lu]  = Lavg[ilev];
      }  //e-if
    }    //e-for
  }      //e-for

  // Free memory:
  free (Edir);
  free (Edn);
  free (Eup);
  free (Lavg);
  free (dtau_clr_org_d);
  free (omega0_clr_org_d);
  free (g_clr_org_d);
  free (dtau_cldk_org_d);
  free (omega0_cldk_org_d);
  free (g_cldk_org_d);
  free (dtau_cldn_org_d);
  free (omega0_cldn_org_d);
  free (g_cldn_org_d);
  free (B);

  return 0;
}  //e-twomaxrnd3C

//=============================
// FUNCTION twostream_maxrand3C
//=============================

/*
OUTPUT:
Edir = direct irradiance;
Edn  = diffuse downward irradiance;
Eup  = diffuse upward irradiance;
Lavg = radiance (set to NaN);
*/

static int twostream_maxrand3C (double*  dtau_org_ck,
                                double*  omega0_org_ck,
                                double*  g_org_ck,  // "_ck" for optically thick cloudy region
                                double*  dtau_org_cn,
                                double*  omega0_org_cn,
                                double*  g_org_cn,  // "_cn" for optically thin cloudy region
                                double*  dtau_org_f,
                                double*  omega0_org_f,
                                double*  g_org_f,  // "_f"  for cloud-free region
                                double*  cf,
                                double   scale_cf,
                                int      nlev,
                                double   S0,
                                double   mu0,
                                double   Ag,
                                double   Bg,
                                double*  B,
                                int      delta,
                                int      flagSolar,
                                int      flagThermal,
                                double** Edir,
                                double** Edn,
                                double** Eup,
                                double** Lavg) {
  int nlyr    = nlev - 1;
  int ilev    = 0;
  int ilyr    = 0;
  int iStatus = 0;

  double dtau_ck;
  double omega0_ck;
  double g_ck;

  double dtau_cn;
  double omega0_cn;
  double g_cn;

  double dtau_f;
  double omega0_f;
  double g_f;

  /*
  Eddington coefficients:
  a11 = transmission coefficient for diffuse radiation;
  a12 = reflection coefficient for diffuse radiation; 
  a13 = reflection coefficient for the primary scattered parallel solar radiation;
  a23 = transmission coefficient for the primary scattered parallel solar radiation;
  a33 = transmission coefficient for the direct parallel solar radiation;
  */

  double a11_ck;
  double a12_ck;
  double a13_ck;
  double a23_ck;
  double a33_ck;

  double a11_cn;
  double a12_cn;
  double a13_cn;
  double a23_cn;
  double a33_cn;

  double a11_f;
  double a12_f;
  double a13_f;
  double a23_f;
  double a33_f;

  double* ar_a11_ck = calloc (nlyr, sizeof (double));
  double* ar_a12_ck = calloc (nlyr, sizeof (double));
  double* ar_a13_ck = calloc (nlyr, sizeof (double));
  double* ar_a23_ck = calloc (nlyr, sizeof (double));
  double* ar_a33_ck = calloc (nlyr, sizeof (double));

  double* ar_a11_cn = calloc (nlyr, sizeof (double));
  double* ar_a12_cn = calloc (nlyr, sizeof (double));
  double* ar_a13_cn = calloc (nlyr, sizeof (double));
  double* ar_a23_cn = calloc (nlyr, sizeof (double));
  double* ar_a33_cn = calloc (nlyr, sizeof (double));

  double* ar_a11_f = calloc (nlyr, sizeof (double));
  double* ar_a12_f = calloc (nlyr, sizeof (double));
  double* ar_a13_f = calloc (nlyr, sizeof (double));
  double* ar_a23_f = calloc (nlyr, sizeof (double));
  double* ar_a33_f = calloc (nlyr, sizeof (double));

  // Components of vector B in the thermal spectral range:
  double theComp1_ck;
  double theComp1_cn;
  double theComp1_f;

  double theComp2_ck;
  double theComp2_cn;
  double theComp2_f;

  double* ar_theComp1_ck = NULL;
  double* ar_theComp1_cn = NULL;
  double* ar_theComp1_f  = NULL;

  double* ar_theComp2_ck = NULL;
  double* ar_theComp2_cn = NULL;
  double* ar_theComp2_f  = NULL;

  // Transport coefficients:
  double* pdn_CKtCK = calloc (nlyr, sizeof (double));
  double* pdn_CKtCN = calloc (nlyr, sizeof (double));
  double* pdn_CKtF  = calloc (nlyr, sizeof (double));

  double* pdn_CNtCK = calloc (nlyr, sizeof (double));
  double* pdn_CNtCN = calloc (nlyr, sizeof (double));
  double* pdn_CNtF  = calloc (nlyr, sizeof (double));

  double* pdn_FtCK = calloc (nlyr, sizeof (double));
  double* pdn_FtCN = calloc (nlyr, sizeof (double));
  double* pdn_FtF  = calloc (nlyr, sizeof (double));

  double* pup_CKtCK = calloc (nlyr, sizeof (double));
  double* pup_CKtCN = calloc (nlyr, sizeof (double));
  double* pup_CKtF  = calloc (nlyr, sizeof (double));

  double* pup_CNtCK = calloc (nlyr, sizeof (double));
  double* pup_CNtCN = calloc (nlyr, sizeof (double));
  double* pup_CNtF  = calloc (nlyr, sizeof (double));

  double* pup_FtCK = calloc (nlyr, sizeof (double));
  double* pup_FtCN = calloc (nlyr, sizeof (double));
  double* pup_FtF  = calloc (nlyr, sizeof (double));

  double** AA     = NULL;
  double*  bb_sol = NULL;
  double*  bb_the = NULL;
  double*  bb;
  double*  xx;  // result vector

  double* S_ck = NULL;
  double* S_cn = NULL;
  double* S_f  = NULL;

  double* Edir_ck = NULL;
  double* Edir_cn = NULL;
  double* Edir_f  = NULL;

  double* Eup_ck;
  double* Eup_cn;
  double* Eup_f;

  double* Edn_ck;
  double* Edn_cn;
  double* Edn_f;

  double* cf_ck = calloc (nlyr, sizeof (double));
  double* cf_cn = calloc (nlyr, sizeof (double));

  AA = calloc (NFLUX * nlev, sizeof (double*));

  for (ilev = 0; ilev < NFLUX * nlev; ilev++) {
    if ((AA[ilev] = calloc (NFLUX * nlev, sizeof (double))) == NULL) {
      fprintf (stderr, "Error allocating memory for AA[%d]\n", ilev);
      return -1;
    }  //e-if
  }    //e-for

  if (flagSolar) {
    bb_sol = calloc (NFLUX * nlev, sizeof (double));
  }  //e-if

  if (flagThermal) {
    ar_theComp1_ck = calloc (nlev, sizeof (double));
    ar_theComp1_cn = calloc (nlev, sizeof (double));
    ar_theComp1_f  = calloc (nlev, sizeof (double));

    ar_theComp2_ck = calloc (nlev, sizeof (double));
    ar_theComp2_cn = calloc (nlev, sizeof (double));
    ar_theComp2_f  = calloc (nlev, sizeof (double));

    bb_the = calloc (NFLUX * nlev, sizeof (double));
  }  //e-if

  bb = calloc (NFLUX * nlev, sizeof (double));
  xx = calloc (NFLUX * nlev, sizeof (double));

  S_ck = calloc (nlev, sizeof (double));
  S_cn = calloc (nlev, sizeof (double));
  S_f  = calloc (nlev, sizeof (double));

  Edir_ck = calloc (nlev, sizeof (double));
  Edir_cn = calloc (nlev, sizeof (double));
  Edir_f  = calloc (nlev, sizeof (double));

  Eup_ck = calloc (nlev, sizeof (double));
  Eup_cn = calloc (nlev, sizeof (double));
  Eup_f  = calloc (nlev, sizeof (double));

  Edn_ck = calloc (nlev, sizeof (double));
  Edn_cn = calloc (nlev, sizeof (double));
  Edn_f  = calloc (nlev, sizeof (double));

  // Final total irradiances (sum: ck + cn + f):
  *Edir = calloc (nlev, sizeof (double));
  *Edn  = calloc (nlev, sizeof (double));
  *Eup  = calloc (nlev, sizeof (double));
  *Lavg = calloc (nlev, sizeof (double));

  // Check if scale_cf is set properly (0.0 <= scale_cf >= 1.0):
  if (scale_cf < 0.0 || scale_cf > 1.0) {
    fprintf (stderr, "Error - invalid value of input parameter: scale_cf should be between 0.0 and 1.0 !!!\n");
    return -2;
  }  //e-if

  // Calculate cf_ck and cf_cn:
  for (ilyr = 0; ilyr < nlyr; ilyr++) {
    cf_ck[ilyr] = scale_cf * cf[ilyr];
    cf_cn[ilyr] = (1.0 - scale_cf) * cf[ilyr];
  }  //e-for

  // At the moment it is only possible to calculate solar OR thermal RT, but not both simultaneously;
  if ((flagSolar && flagThermal) || (!flagSolar && !flagThermal)) {
    fprintf (stderr, "Error - invalid input parameters - use flagSolar and flagThermal alternatingly\n");
    fprintf (stderr, "flagSolar = %d , flagThermal = %d\n", flagSolar, flagThermal);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -3;
  }  //e-if

  // Calculate vertical profiles of transport coefficients pdn_... from vertical profiles of partial cloudiness:
  iStatus = calc_pdn_transport_coefficients_twomaxrnd3C (nlev,
                                                         cf,
                                                         cf_ck,
                                                         cf_cn,
                                                         pdn_CKtCK,
                                                         pdn_CKtCN,
                                                         pdn_CKtF,
                                                         pdn_CNtCK,
                                                         pdn_CNtCN,
                                                         pdn_CNtF,
                                                         pdn_FtCK,
                                                         pdn_FtCN,
                                                         pdn_FtF);
  if (iStatus != 0) {
    fprintf (stderr, "Error calculating vertical profiles of pdn transport coefficients; ERROR=%d \n", iStatus);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -4;
  }  //e-if

  // Calculate vertical profiles of transport coefficients pup_... from vertical profiles of partial cloudiness:
  iStatus = calc_pup_transport_coefficients_twomaxrnd3C (nlev,
                                                         cf,
                                                         cf_ck,
                                                         cf_cn,
                                                         pup_CKtCK,
                                                         pup_CKtCN,
                                                         pup_CKtF,
                                                         pup_CNtCK,
                                                         pup_CNtCN,
                                                         pup_CNtF,
                                                         pup_FtCK,
                                                         pup_FtCN,
                                                         pup_FtF);
  if (iStatus != 0) {
    fprintf (stderr, "Error calculating vertical profiles of pup transport coefficients; ERROR=%d \n", iStatus);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -5;
  }  //e-if

  // Calculate vertical profiles of Eddington coefficients
  // and vertical profiles of thermal components:
  for (ilyr = 0; ilyr < nlyr; ilyr++) {

    // Delta scaling of optical properties of cloudy regions:
    if (delta) {
      delta_scale_hg (dtau_org_ck[ilyr], omega0_org_ck[ilyr], g_org_ck[ilyr], &dtau_ck, &omega0_ck, &g_ck);
      delta_scale_hg (dtau_org_cn[ilyr], omega0_org_cn[ilyr], g_org_cn[ilyr], &dtau_cn, &omega0_cn, &g_cn);
    } else {
      dtau_ck   = dtau_org_ck[ilyr];
      omega0_ck = omega0_org_ck[ilyr];
      g_ck      = g_org_ck[ilyr];

      dtau_cn   = dtau_org_cn[ilyr];
      omega0_cn = omega0_org_cn[ilyr];
      g_cn      = g_org_cn[ilyr];
    }  //e-if

    // Optical properties of cloud-free regions:
    dtau_f   = dtau_org_f[ilyr];
    omega0_f = omega0_org_f[ilyr];
    g_f      = g_org_f[ilyr];

    // omega0 should not be 1 (avoiding singularity problem), restrict omega0 to 0.999999:
    if (omega0_ck > 0.999999)
      omega0_ck = 0.999999;
    if (omega0_cn > 0.999999)
      omega0_cn = 0.999999;
    if (omega0_f > 0.999999)
      omega0_f = 0.999999;

    eddington_coeffc (dtau_ck, g_ck, omega0_ck, mu0, &a11_ck, &a12_ck, &a13_ck, &a23_ck, &a33_ck);
    eddington_coeffc (dtau_cn, g_cn, omega0_cn, mu0, &a11_cn, &a12_cn, &a13_cn, &a23_cn, &a33_cn);
    eddington_coeffc (dtau_f, g_f, omega0_f, mu0, &a11_f, &a12_f, &a13_f, &a23_f, &a33_f);

    ar_a11_ck[ilyr] = a11_ck;
    ar_a11_cn[ilyr] = a11_cn;
    ar_a11_f[ilyr]  = a11_f;

    ar_a12_ck[ilyr] = a12_ck;
    ar_a12_cn[ilyr] = a12_cn;
    ar_a12_f[ilyr]  = a12_f;

    ar_a13_ck[ilyr] = a13_ck;
    ar_a13_cn[ilyr] = a13_cn;
    ar_a13_f[ilyr]  = a13_f;

    ar_a23_ck[ilyr] = a23_ck;
    ar_a23_cn[ilyr] = a23_cn;
    ar_a23_f[ilyr]  = a23_f;

    ar_a33_ck[ilyr] = a33_ck;
    ar_a33_cn[ilyr] = a33_cn;
    ar_a33_f[ilyr]  = a33_f;

    if (flagThermal) {
      calcThermalComponents (ilyr, B, dtau_ck, omega0_ck, g_ck, &theComp1_ck, &theComp2_ck);
      calcThermalComponents (ilyr, B, dtau_cn, omega0_cn, g_cn, &theComp1_cn, &theComp2_cn);
      calcThermalComponents (ilyr, B, dtau_f, omega0_f, g_f, &theComp1_f, &theComp2_f);

      ar_theComp1_ck[ilyr] = theComp1_ck;
      ar_theComp1_cn[ilyr] = theComp1_cn;
      ar_theComp1_f[ilyr]  = theComp1_f;

      ar_theComp2_ck[ilyr] = theComp2_ck;
      ar_theComp2_cn[ilyr] = theComp2_cn;
      ar_theComp2_f[ilyr]  = theComp2_f;
    }  //e-if
  }    //e-for

  // Initialize vectors S_ck, S_cn and S_f:
  // S[0]= S0 = S_ck[0] + S_cn[0] + S_f[0] = cf[0]*S0 + (1.0-cf[0])*S0 = cf_ck[0]*S0 + cf_cn[0]*S0 + (1.0-cf[0])*S0
  if (flagSolar) {
    S_ck[0] = cf_ck[0] * S0;
    S_cn[0] = cf_cn[0] * S0;
    S_f[0]  = (1.0 - cf[0]) * S0;

    for (ilev = 1; ilev < nlev; ilev++) {
      S_ck[ilev] = ar_a33_ck[ilev - 1] * (pdn_FtCK[ilev - 1] * S_f[ilev - 1] + pdn_CKtCK[ilev - 1] * S_ck[ilev - 1] +
                                          pdn_CNtCK[ilev - 1] * S_cn[ilev - 1]);
      S_cn[ilev] = ar_a33_cn[ilev - 1] * (pdn_FtCN[ilev - 1] * S_f[ilev - 1] + pdn_CKtCN[ilev - 1] * S_ck[ilev - 1] +
                                          pdn_CNtCN[ilev - 1] * S_cn[ilev - 1]);
      S_f[ilev]  = ar_a33_f[ilev - 1] *
                  (pdn_FtF[ilev - 1] * S_f[ilev - 1] + pdn_CKtF[ilev - 1] * S_ck[ilev - 1] + pdn_CNtF[ilev - 1] * S_cn[ilev - 1]);
    }  //e-for
  }    //e-if

  // Equation system has the following form: xx = AA*xx + bb;

  // Build matrix AA:
  iStatus = buildMatrixA_twomaxrnd3C (nlev,
                                      Ag,
                                      ar_a11_ck,
                                      ar_a11_cn,
                                      ar_a11_f,
                                      ar_a12_ck,
                                      ar_a12_cn,
                                      ar_a12_f,
                                      pdn_CKtCK,
                                      pdn_CKtCN,
                                      pdn_CKtF,
                                      pdn_CNtCK,
                                      pdn_CNtCN,
                                      pdn_CNtF,
                                      pdn_FtCK,
                                      pdn_FtCN,
                                      pdn_FtF,
                                      pup_CKtCK,
                                      pup_CKtCN,
                                      pup_CKtF,
                                      pup_CNtCK,
                                      pup_CNtCN,
                                      pup_CNtF,
                                      pup_FtCK,
                                      pup_FtCN,
                                      pup_FtF,
                                      AA);
  if (iStatus != 0) {
    fprintf (stderr, "buildMatrixA ERROR=%d \n", iStatus);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -6;
  }  //e-if

  // Make matrix A2 = AA - IdentityMatrix and save it to the same matrix AA,
  // since the equation system in the form Ax=b must be passed to solve_gauss;
  // xx = AA*xx + bb;
  // (AA-II)*xx = -bb;

  iStatus = makeMatrixA2 (nlev, NFLUX, AA);

  if (iStatus != 0) {
    fprintf (stderr, "makeMatrixA2 ERROR=%d \n", iStatus);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -7;
  }  //e-if

  // Build vector bb_sol:
  if (flagSolar) {
    iStatus = buildVectorBsol_twomaxrnd3C (nlev,
                                           Ag,
                                           mu0,
                                           ar_a13_ck,
                                           ar_a13_cn,
                                           ar_a13_f,
                                           ar_a23_ck,
                                           ar_a23_cn,
                                           ar_a23_f,
                                           S_ck,
                                           S_cn,
                                           S_f,
                                           pdn_CKtCK,
                                           pdn_CKtCN,
                                           pdn_CKtF,
                                           pdn_CNtCK,
                                           pdn_CNtCN,
                                           pdn_CNtF,
                                           pdn_FtCK,
                                           pdn_FtCN,
                                           pdn_FtF,
                                           bb_sol);
    if (iStatus != 0) {
      fprintf (stderr, "buildVectorBsol ERROR=%d \n", iStatus);
      freeMemory_twomaxrnd3C (nlev,
                              ar_a11_ck,
                              ar_a11_cn,
                              ar_a11_f,
                              ar_a12_ck,
                              ar_a12_cn,
                              ar_a12_f,
                              ar_a13_ck,
                              ar_a13_cn,
                              ar_a13_f,
                              ar_a23_ck,
                              ar_a23_cn,
                              ar_a23_f,
                              ar_a33_ck,
                              ar_a33_cn,
                              ar_a33_f,
                              pdn_CKtCK,
                              pdn_CKtCN,
                              pdn_CKtF,
                              pdn_CNtCK,
                              pdn_CNtCN,
                              pdn_CNtF,
                              pdn_FtCK,
                              pdn_FtCN,
                              pdn_FtF,
                              pup_CKtCK,
                              pup_CKtCN,
                              pup_CKtF,
                              pup_CNtCK,
                              pup_CNtCN,
                              pup_CNtF,
                              pup_FtCK,
                              pup_FtCN,
                              pup_FtF,
                              bb_sol,
                              bb_the,
                              ar_theComp1_ck,
                              ar_theComp1_cn,
                              ar_theComp1_f,
                              ar_theComp2_ck,
                              ar_theComp2_cn,
                              ar_theComp2_f,
                              S_ck,
                              S_cn,
                              S_f,
                              Edir_ck,
                              Edir_cn,
                              Edir_f,
                              Eup_ck,
                              Eup_cn,
                              Eup_f,
                              Edn_ck,
                              Edn_cn,
                              Edn_f,
                              bb,
                              xx,
                              AA);
      return -8;
    }  //e-if
  }    //e-if

  // Build vector bb_the:
  if (flagThermal) {
    iStatus = buildVectorBthe_twomaxrnd3C (nlev,
                                           Ag,
                                           Bg,
                                           ar_theComp1_ck,
                                           ar_theComp1_cn,
                                           ar_theComp1_f,
                                           ar_theComp2_ck,
                                           ar_theComp2_cn,
                                           ar_theComp2_f,
                                           cf,
                                           cf_ck,
                                           cf_cn,
                                           bb_the);
    if (iStatus != 0) {
      fprintf (stderr, "buildVectorBthe ERROR=%d \n", iStatus);
      freeMemory_twomaxrnd3C (nlev,
                              ar_a11_ck,
                              ar_a11_cn,
                              ar_a11_f,
                              ar_a12_ck,
                              ar_a12_cn,
                              ar_a12_f,
                              ar_a13_ck,
                              ar_a13_cn,
                              ar_a13_f,
                              ar_a23_ck,
                              ar_a23_cn,
                              ar_a23_f,
                              ar_a33_ck,
                              ar_a33_cn,
                              ar_a33_f,
                              pdn_CKtCK,
                              pdn_CKtCN,
                              pdn_CKtF,
                              pdn_CNtCK,
                              pdn_CNtCN,
                              pdn_CNtF,
                              pdn_FtCK,
                              pdn_FtCN,
                              pdn_FtF,
                              pup_CKtCK,
                              pup_CKtCN,
                              pup_CKtF,
                              pup_CNtCK,
                              pup_CNtCN,
                              pup_CNtF,
                              pup_FtCK,
                              pup_FtCN,
                              pup_FtF,
                              bb_sol,
                              bb_the,
                              ar_theComp1_ck,
                              ar_theComp1_cn,
                              ar_theComp1_f,
                              ar_theComp2_ck,
                              ar_theComp2_cn,
                              ar_theComp2_f,
                              S_ck,
                              S_cn,
                              S_f,
                              Edir_ck,
                              Edir_cn,
                              Edir_f,
                              Eup_ck,
                              Eup_cn,
                              Eup_f,
                              Edn_ck,
                              Edn_cn,
                              Edn_f,
                              bb,
                              xx,
                              AA);
      return -9;
    }  //e-if
  }    //e-if

  // Assign bb_sol or bb_the to the final vector bb:
  for (ilev = 0; ilev < NFLUX * nlev; ilev++) {
    bb[ilev] = 0;

    if (flagSolar) {
      bb[ilev] += bb_sol[ilev];
    }  //e-if

    if (flagThermal) {
      bb[ilev] += bb_the[ilev];
    }  //e-if

  }  //e-for

  // Make vector -bb and save it to the same vector bb:
  iStatus = makeVectorMinusB (nlev, NFLUX, bb);

  if (iStatus != 0) {
    fprintf (stderr, "makeVectorMinusB ERROR=%d \n", iStatus);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -10;
  }  //e-if

  // Solve AA*xx=bb to obtain xx:
  double **matr, **mlu;
  long*    indx;
  int      i = 0, nn = NFLUX * nlev, m1 = 8, m2 = 8;  // n=2*input.nlyr+2;
  matr = calloc (nn, sizeof (double*));
  mlu  = calloc (nn, sizeof (double*));
  for (i = 0; i < nn; i++) {
    matr[i] = calloc (m1 + m2 + 1, sizeof (double));
    mlu[i]  = calloc (m2, sizeof (double));
  }
  indx = calloc (nn, sizeof (long));

  matrix_to_band_form (AA, nn, m1, m2, matr);  // Transform matrix AA to band form matr
  bandec_M (matr, nn, m1, m2, mlu, indx);
  banbks_M (matr, nn, m1, m2, mlu, indx, bb);

  for (i = 0; i < nn; i++) {
    xx[i] = bb[i];
  }

  for (i = 0; i < nn; i++) {
    free (mlu[i]);
    free (matr[i]);
  }
  free (indx);
  free (mlu);
  free (matr);

  if (iStatus != 0) {
    fprintf (stderr, "Error %d solving equation system\n", iStatus);
    freeMemory_twomaxrnd3C (nlev,
                            ar_a11_ck,
                            ar_a11_cn,
                            ar_a11_f,
                            ar_a12_ck,
                            ar_a12_cn,
                            ar_a12_f,
                            ar_a13_ck,
                            ar_a13_cn,
                            ar_a13_f,
                            ar_a23_ck,
                            ar_a23_cn,
                            ar_a23_f,
                            ar_a33_ck,
                            ar_a33_cn,
                            ar_a33_f,
                            pdn_CKtCK,
                            pdn_CKtCN,
                            pdn_CKtF,
                            pdn_CNtCK,
                            pdn_CNtCN,
                            pdn_CNtF,
                            pdn_FtCK,
                            pdn_FtCN,
                            pdn_FtF,
                            pup_CKtCK,
                            pup_CKtCN,
                            pup_CKtF,
                            pup_CNtCK,
                            pup_CNtCN,
                            pup_CNtF,
                            pup_FtCK,
                            pup_FtCN,
                            pup_FtF,
                            bb_sol,
                            bb_the,
                            ar_theComp1_ck,
                            ar_theComp1_cn,
                            ar_theComp1_f,
                            ar_theComp2_ck,
                            ar_theComp2_cn,
                            ar_theComp2_f,
                            S_ck,
                            S_cn,
                            S_f,
                            Edir_ck,
                            Edir_cn,
                            Edir_f,
                            Eup_ck,
                            Eup_cn,
                            Eup_f,
                            Edn_ck,
                            Edn_cn,
                            Edn_f,
                            bb,
                            xx,
                            AA);
    return -11;
  }  //e-if

  for (ilev = 0; ilev < nlev; ilev++) {
    Eup_f[ilev]  = xx[NFLUX * ilev];
    Eup_ck[ilev] = xx[NFLUX * ilev + 1];
    Eup_cn[ilev] = xx[NFLUX * ilev + 2];

    Edn_f[ilev]  = xx[NFLUX * ilev + 3];
    Edn_ck[ilev] = xx[NFLUX * ilev + 4];
    Edn_cn[ilev] = xx[NFLUX * ilev + 5];
  }  //e-for

  for (ilev = 0; ilev < nlev; ilev++) {
    Edir_ck[ilev] = S_ck[ilev] * mu0;
    Edir_cn[ilev] = S_cn[ilev] * mu0;
    Edir_f[ilev]  = S_f[ilev] * mu0;
  }  //e-for

  // Sum up the irradiances of the three regions to obtain the final result:
  for (ilev = 0; ilev < nlev; ilev++) {
    (*Edir)[ilev] = Edir_ck[ilev] + Edir_cn[ilev] + Edir_f[ilev];
    (*Eup)[ilev]  = Eup_ck[ilev] + Eup_cn[ilev] + Eup_f[ilev];
    (*Edn)[ilev]  = Edn_ck[ilev] + Edn_cn[ilev] + Edn_f[ilev];
    (*Lavg)[ilev] = NAN;
  }  //e-for

  freeMemory_twomaxrnd3C (nlev,
                          ar_a11_ck,
                          ar_a11_cn,
                          ar_a11_f,
                          ar_a12_ck,
                          ar_a12_cn,
                          ar_a12_f,
                          ar_a13_ck,
                          ar_a13_cn,
                          ar_a13_f,
                          ar_a23_ck,
                          ar_a23_cn,
                          ar_a23_f,
                          ar_a33_ck,
                          ar_a33_cn,
                          ar_a33_f,
                          pdn_CKtCK,
                          pdn_CKtCN,
                          pdn_CKtF,
                          pdn_CNtCK,
                          pdn_CNtCN,
                          pdn_CNtF,
                          pdn_FtCK,
                          pdn_FtCN,
                          pdn_FtF,
                          pup_CKtCK,
                          pup_CKtCN,
                          pup_CKtF,
                          pup_CNtCK,
                          pup_CNtCN,
                          pup_CNtF,
                          pup_FtCK,
                          pup_FtCN,
                          pup_FtF,
                          bb_sol,
                          bb_the,
                          ar_theComp1_ck,
                          ar_theComp1_cn,
                          ar_theComp1_f,
                          ar_theComp2_ck,
                          ar_theComp2_cn,
                          ar_theComp2_f,
                          S_ck,
                          S_cn,
                          S_f,
                          Edir_ck,
                          Edir_cn,
                          Edir_f,
                          Eup_ck,
                          Eup_cn,
                          Eup_f,
                          Edn_ck,
                          Edn_cn,
                          Edn_f,
                          bb,
                          xx,
                          AA);

  return 0;
}  //e-twostream_maxrand3C

//=========================================
//FUNCTION: calc_pdn_transport_coefficients
//=========================================
//INPUT:  nlev, cf, cf_ck, cf_cn;
//OUTPUT: pdn_CKtCK, pdn_CKtCN, pdn_CKtF, pdn_CNtCK, pdn_CNtCN, pdn_CNtF, pdn_FtCK, pdn_FtCN, pdn_FtF;

static int calc_pdn_transport_coefficients_twomaxrnd3C (int     nlev,
                                                        double* cf,
                                                        double* cf_ck,
                                                        double* cf_cn,
                                                        double* pdn_CKtCK,
                                                        double* pdn_CKtCN,
                                                        double* pdn_CKtF,
                                                        double* pdn_CNtCK,
                                                        double* pdn_CNtCN,
                                                        double* pdn_CNtF,
                                                        double* pdn_FtCK,
                                                        double* pdn_FtCN,
                                                        double* pdn_FtF) {
  int ilyr;
  int nlyr;

  double a;  //upper-layer cloud thick;
  double b;  //upper-layer cloud thin;
  double c;  //upper-layer clear-sky;
  double x;  //this-layer cloud thick;
  double y;  //this-layer cloud thin;
  double z;  //this-layer clear-sky;

  nlyr = nlev - 1;

  ilyr = 0;  //SPECIAL CASE

  pdn_CKtCK[ilyr] = 1.0;
  pdn_CKtCN[ilyr] = 0.0;
  pdn_CKtF[ilyr]  = 0.0;

  pdn_CNtCK[ilyr] = 0.0;
  pdn_CNtCN[ilyr] = 1.0;
  pdn_CNtF[ilyr]  = 0.0;

  pdn_FtCK[ilyr] = 0.0;
  pdn_FtCN[ilyr] = 0.0;
  pdn_FtF[ilyr]  = 1.0;

  for (ilyr = 1; ilyr < nlyr; ilyr++) {

    a = cf_ck[ilyr - 1];
    b = cf_cn[ilyr - 1];
    c = 1.0 - cf[ilyr - 1];
    x = cf_ck[ilyr];
    y = cf_cn[ilyr];
    z = 1.0 - cf[ilyr];

    if (cf[ilyr - 1] < 0.0000000001) {  //SPECIAL CASE: In the layer above entirely clear-sky (i.e., a=0, b=0, c=1)
      pdn_CKtCK[ilyr] = 1.0;
      pdn_CKtCN[ilyr] = 0.0;
      pdn_CKtF[ilyr]  = 0.0;

      pdn_CNtCK[ilyr] = 0.0;
      pdn_CNtCN[ilyr] = 1.0;
      pdn_CNtF[ilyr]  = 0.0;

      pdn_FtCK[ilyr] = x;
      pdn_FtCN[ilyr] = y;
      pdn_FtF[ilyr]  = z;

      continue;
    }  //-eif

    if (cf[ilyr] < 0.0000000001) {  // SPECIAL CASE: In this layer entirely clear-sky (i.e., x=0, y=0, z=1);
      pdn_CKtCK[ilyr] = 0.0;
      pdn_CKtCN[ilyr] = 0.0;
      pdn_CKtF[ilyr]  = 1.0;

      pdn_CNtCK[ilyr] = 0.0;
      pdn_CNtCN[ilyr] = 0.0;
      pdn_CNtF[ilyr]  = 1.0;

      pdn_FtCK[ilyr] = 0.0;  //x
      pdn_FtCN[ilyr] = 0.0;  //y
      pdn_FtF[ilyr]  = 1.0;  //z

      continue;
    }  //e-if

    if (cf[ilyr] > cf[ilyr - 1]) {  // Case 1: Cloud in this layer larger than cloud above;
                                    //       = Thick cloud in this layer larger than thick cloud above;
                                    // (implies also: c not equal to 0.0, but: c > 0.0);

      if (x < (a + b)) {  // Case 1a: Thick cloud in this layer smaller than the entire cloud above;
        pdn_CKtCK[ilyr] = 1.0;
        pdn_CKtCN[ilyr] = 0.0;
        pdn_CKtF[ilyr]  = 0.0;

        if (b > 0.0) {
          pdn_CNtCK[ilyr] = (x - a) / b;
          pdn_CNtCN[ilyr] = (a + b - x) / b;
        } else {
          fprintf (stderr, "pdn, Case 1a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -101;
        }  //e-if
        pdn_CNtF[ilyr] = 0.0;

        pdn_FtCK[ilyr] = 0.0;
        if (c > 0.0) {
          pdn_FtCN[ilyr] = (x + y - a - b) / c;
          pdn_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pdn, Case 1a: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -102;
        }  //e-if

      } else {  // Case 1b: Thick cloud in this layer larger than or equal to the entire cloud above;
        pdn_CKtCK[ilyr] = 1.0;
        pdn_CKtCN[ilyr] = 0.0;
        pdn_CKtF[ilyr]  = 0.0;

        pdn_CNtCK[ilyr] = 1.0;
        pdn_CNtCN[ilyr] = 0.0;
        pdn_CNtF[ilyr]  = 0.0;

        if (c > 0.0) {
          pdn_FtCK[ilyr] = (x - a - b) / c;
          pdn_FtCN[ilyr] = y / c;
          pdn_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pdn, Case 1b: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -103;
        }  //e-if

      }  //e-if-Case1a1b

    } else {  // Case 2: Cloud in this layer smaller than or equal to the cloud above;
              //       = Thick cloud in this layer smaller than or equal to the thick cloud above;

      if (a < (x + y)) {  // Case 2a: Thick cloud above smaller than the entire cloud in this layer;
        if (a > 0.0) {
          pdn_CKtCK[ilyr] = x / a;
          pdn_CKtCN[ilyr] = (a - x) / a;
        } else {
          fprintf (stderr, "pdn, Case 2a: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -104;
        }  //e-if
        pdn_CKtF[ilyr] = 0.0;

        pdn_CNtCK[ilyr] = 0.0;
        if (b > 0.0) {
          pdn_CNtCN[ilyr] = (x + y - a) / b;
          pdn_CNtF[ilyr]  = (a + b - x - y) / b;
        } else {
          fprintf (stderr, "pdn, Case 2a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -105;
        }  //e-if

        pdn_FtCK[ilyr] = 0.0;
        pdn_FtCN[ilyr] = 0.0;
        pdn_FtF[ilyr]  = 1.0;

      } else {  // Case 2b: Thick cloud above larger than or equal to the entire cloud in this layer;
        if (a > 0.0) {
          pdn_CKtCK[ilyr] = x / a;
          pdn_CKtCN[ilyr] = y / a;
          pdn_CKtF[ilyr]  = (a - x - y) / a;
        } else {
          fprintf (stderr, "pdn, Case 2b: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -106;
        }  //e-if

        pdn_CNtCK[ilyr] = 0.0;
        pdn_CNtCN[ilyr] = 0.0;
        pdn_CNtF[ilyr]  = 1.0;

        pdn_FtCK[ilyr] = 0.0;
        pdn_FtCN[ilyr] = 0.0;
        pdn_FtF[ilyr]  = 1.0;

      }  //e-if-Case2a2b
    }    //e-if-Case12

  }  //e-for over ilyr;

  return 0;
}  //e-calc_pdn_transport_coefficients

//=========================================
//FUNCTION: calc_pup_transport_coefficients
//=========================================
//INPUT:  nlev, cf, cf_ck, cf_cn;
//OUTPUT: pup_CKtCK, pup_CKtCN, pup_CKtF, pup_CNtCK, pup_CNtCN, pup_CNtF, pup_FtCK, pup_FtCN, pup_FtF;

static int calc_pup_transport_coefficients_twomaxrnd3C (int     nlev,
                                                        double* cf,
                                                        double* cf_ck,
                                                        double* cf_cn,
                                                        double* pup_CKtCK,
                                                        double* pup_CKtCN,
                                                        double* pup_CKtF,
                                                        double* pup_CNtCK,
                                                        double* pup_CNtCN,
                                                        double* pup_CNtF,
                                                        double* pup_FtCK,
                                                        double* pup_FtCN,
                                                        double* pup_FtF) {
  int ilyr;
  int nlyr;

  double x;  //this-layer cloud thick;
  double y;  //this-layer cloud thin;
  double z;  //this-layer clear-sky;
  double a;  //bottom-layer cloud thick;
  double b;  //bottom-layer cloud thin;
  double c;  //bottom-layer clear-sky;

  nlyr = nlev - 1;

  ilyr = (nlyr - 1);  // SPECIAL CASE

  pup_CKtCK[ilyr] = 1.0;
  pup_CKtCN[ilyr] = 0.0;
  pup_CKtF[ilyr]  = 0.0;

  pup_CNtCK[ilyr] = 0.0;
  pup_CNtCN[ilyr] = 1.0;
  pup_CNtF[ilyr]  = 0.0;

  pup_FtCK[ilyr] = 0.0;
  pup_FtCN[ilyr] = 0.0;
  pup_FtF[ilyr]  = 1.0;

  for (ilyr = 0; ilyr < (nlyr - 1); ilyr++) {

    x = cf_ck[ilyr];
    y = cf_cn[ilyr];
    z = 1.0 - cf[ilyr];
    a = cf_ck[ilyr + 1];
    b = cf_cn[ilyr + 1];
    c = 1.0 - cf[ilyr + 1];

    if (cf[ilyr + 1] < 0.0000000001) {  // SPECIAL CASE: In the layer below entirely clear-sky (i.e., a=0, b=0, c=1);
      pup_CKtCK[ilyr] = 1.0;
      pup_CKtCN[ilyr] = 0.0;
      pup_CKtF[ilyr]  = 0.0;

      pup_CNtCK[ilyr] = 0.0;
      pup_CNtCN[ilyr] = 1.0;
      pup_CNtF[ilyr]  = 0.0;

      pup_FtCK[ilyr] = x;
      pup_FtCN[ilyr] = y;
      pup_FtF[ilyr]  = z;

      continue;
    }  //e-if

    if (cf[ilyr] < 0.0000000001) {  // SPECIAL CASE: In this layer entirely clear-sky (i.e., x=0, y=0, z=1);
      pup_CKtCK[ilyr] = 0.0;
      pup_CKtCN[ilyr] = 0.0;
      pup_CKtF[ilyr]  = 1.0;

      pup_CNtCK[ilyr] = 0.0;
      pup_CNtCN[ilyr] = 0.0;
      pup_CNtF[ilyr]  = 1.0;

      pup_FtCK[ilyr] = 0.0;  //x
      pup_FtCN[ilyr] = 0.0;  //y
      pup_FtF[ilyr]  = 1.0;  //z

      continue;
    }  //e-if

    if (cf[ilyr] > cf[ilyr + 1]) {  // Case 1: Cloud in this layer larger than cloud below;
                                    //       = thick cloud in this layer larger than thick cloud below;
                                    // (implies also: c not equal to 0.0, but: c > 0.0);

      if (x < (a + b)) {  // Case 1a: Thick cloud in this layer smaller than the entire cloud below;
        pup_CKtCK[ilyr] = 1.0;
        pup_CKtCN[ilyr] = 0.0;
        pup_CKtF[ilyr]  = 0.0;

        if (b > 0.0) {
          pup_CNtCK[ilyr] = (x - a) / b;
          pup_CNtCN[ilyr] = (a + b - x) / b;
        } else {
          fprintf (stderr, "pup, Case 1a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -201;
        }  //e-if
        pup_CNtF[ilyr] = 0.0;

        pup_FtCK[ilyr] = 0.0;
        if (c > 0.0) {
          pup_FtCN[ilyr] = (x + y - a - b) / c;
          pup_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pup, Case 1a: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -202;
        }  //e-if

      } else {  // Case 1b: Thick cloud in this layer larger than or equal to the entire cloud below;
        pup_CKtCK[ilyr] = 1.0;
        pup_CKtCN[ilyr] = 0.0;
        pup_CKtF[ilyr]  = 0.0;

        pup_CNtCK[ilyr] = 1.0;
        pup_CNtCN[ilyr] = 0.0;
        pup_CNtF[ilyr]  = 0.0;

        if (c > 0.0) {
          pup_FtCK[ilyr] = (x - a - b) / c;
          pup_FtCN[ilyr] = y / c;
          pup_FtF[ilyr]  = z / c;
        } else {
          fprintf (stderr, "pup, Case 1b: This should NOT happen: c = 0.0 !!!, ilyr = %d \n", ilyr);
          return -203;
        }  //e-if

      }  //e-if-Case1a1b

    } else {  // Case 2: Cloud in this layer smaller than or equal to the cloud below;
              //       = thick cloud in this layer smaller than or equal to the thick cloud below;

      if (a < (x + y)) {  // Case 2a: Thick cloud below smaller than the entire cloud in this layer;
        if (a > 0.0) {
          pup_CKtCK[ilyr] = x / a;
          pup_CKtCN[ilyr] = (a - x) / a;
        } else {
          fprintf (stderr, "pup, Case 2a: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -204;
        }  //e-if
        pup_CKtF[ilyr] = 0.0;

        pup_CNtCK[ilyr] = 0.0;
        if (b > 0.0) {
          pup_CNtCN[ilyr] = (x + y - a) / b;
          pup_CNtF[ilyr]  = (a + b - x - y) / b;
        } else {
          fprintf (stderr, "pup, Case 2a: This should NOT happen: b = 0.0 !!!, ilyr = %d \n", ilyr);
          return -205;
        }  //e-if

        pup_FtCK[ilyr] = 0.0;
        pup_FtCN[ilyr] = 0.0;
        pup_FtF[ilyr]  = 1.0;

      } else {  // Case 2b: Thick cloud below larger than or equal to the entire cloud in this layer;
        if (a > 0.0) {
          pup_CKtCK[ilyr] = x / a;
          pup_CKtCN[ilyr] = y / a;
          pup_CKtF[ilyr]  = (a - x - y) / a;
        } else {
          fprintf (stderr, "pup, Case 2b: This should NOT happen: a = 0.0 !!!, ilyr = %d \n", ilyr);
          return -206;
        }  //e-if

        pup_CNtCK[ilyr] = 0.0;
        pup_CNtCN[ilyr] = 0.0;
        pup_CNtF[ilyr]  = 1.0;

        pup_FtCK[ilyr] = 0.0;
        pup_FtCN[ilyr] = 0.0;
        pup_FtF[ilyr]  = 1.0;

      }  //e-if-Case2a2b
    }    //e-if-Case12

  }  //e-for over ilyr;

  return 0;
}  //e-calc_pup_transport_coefficient

//=======================
// FUNCTION: buildMatrixA
//=======================
// INPUT:  nlev, Ag, ar_a11_ck, ar_a11_cn, ar_a11_f, ar_a12_ck, ar_a12_cn, ar_a12_f,
//         pdn_CKtCK, pdn_CKtCN, pdn_CKtF, pdn_CNtCK, pdn_CNtCN, pdn_CNtF, pdn_FtCK, pdn_FtCN, pdn_FtF,
//         pup_CKtCK, pup_CKtCN, pup_CKtF, pup_CNtCK, pup_CNtCN, pup_CNtF, pup_FtCK, pup_FtCN, pup_FtF;
// OUTPUT: matrixA;

static int buildMatrixA_twomaxrnd3C (int      nlev,
                                     double   Ag,
                                     double*  ar_a11_ck,
                                     double*  ar_a11_cn,
                                     double*  ar_a11_f,
                                     double*  ar_a12_ck,
                                     double*  ar_a12_cn,
                                     double*  ar_a12_f,
                                     double*  pdn_CKtCK,
                                     double*  pdn_CKtCN,
                                     double*  pdn_CKtF,
                                     double*  pdn_CNtCK,
                                     double*  pdn_CNtCN,
                                     double*  pdn_CNtF,
                                     double*  pdn_FtCK,
                                     double*  pdn_FtCN,
                                     double*  pdn_FtF,
                                     double*  pup_CKtCK,
                                     double*  pup_CKtCN,
                                     double*  pup_CKtF,
                                     double*  pup_CNtCK,
                                     double*  pup_CNtCN,
                                     double*  pup_CNtF,
                                     double*  pup_FtCK,
                                     double*  pup_FtCN,
                                     double*  pup_FtF,
                                     double** matrixA) {
  int i;     // position in levels (e.g. for nlev=21, i in range: 0 - 20)
  int iRow;  // row position in matrix A (e.g. for nlev=21, iRow in range: 0 - 125)

  int nextColEup;  // index of column for next Eup data
  int nextColEdw;  // index of column for next Edn data; nextColEdn = nextColEup - 6

  // Set values for initial six rows of matrix A:
  // First row:
  matrixA[0][3] = ar_a12_f[0] * pdn_FtF[0];
  matrixA[0][4] = ar_a12_f[0] * pdn_CKtF[0];
  matrixA[0][5] = ar_a12_f[0] * pdn_CNtF[0];

  matrixA[0][6] = ar_a11_f[0] * pup_FtF[0];
  matrixA[0][7] = ar_a11_f[0] * pup_CKtF[0];
  matrixA[0][8] = ar_a11_f[0] * pup_CNtF[0];

  // Second row:
  matrixA[1][3] = ar_a12_ck[0] * pdn_FtCK[0];
  matrixA[1][4] = ar_a12_ck[0] * pdn_CKtCK[0];
  matrixA[1][5] = ar_a12_ck[0] * pdn_CNtCK[0];

  matrixA[1][6] = ar_a11_ck[0] * pup_FtCK[0];
  matrixA[1][7] = ar_a11_ck[0] * pup_CKtCK[0];
  matrixA[1][8] = ar_a11_ck[0] * pup_CNtCK[0];

  // Third row:
  matrixA[2][3] = ar_a12_cn[0] * pdn_FtCN[0];
  matrixA[2][4] = ar_a12_cn[0] * pdn_CKtCN[0];
  matrixA[2][5] = ar_a12_cn[0] * pdn_CNtCN[0];

  matrixA[2][6] = ar_a11_cn[0] * pup_FtCN[0];
  matrixA[2][7] = ar_a11_cn[0] * pup_CKtCN[0];
  matrixA[2][8] = ar_a11_cn[0] * pup_CNtCN[0];

  // 4th row is already zero; (needs to be zero due to upper boundary condition);
  // 5th row is already zero; (needs to be zero due to upper boundary condition);
  // 6th row is already zero; (needs to be zero due to upper boundary condition);

  nextColEup = 9;  // index of column for Eup(level1)
  nextColEdw = 3;  // index of column for Edn(level1)

  for (i = 1; i < nlev; i++) {
    iRow = NFLUX * i;

    if (i == (nlev - 1)) {  // lower boundary condition for the forth and third row from bottom up of matrix A;
      matrixA[iRow][nextColEup]         = Ag;
      matrixA[iRow + 1][nextColEup + 1] = Ag;
      matrixA[iRow + 2][nextColEup + 2] = Ag;
    } else {
      matrixA[iRow][nextColEup]     = ar_a12_f[i] * pdn_FtF[i];
      matrixA[iRow][nextColEup + 1] = ar_a12_f[i] * pdn_CKtF[i];
      matrixA[iRow][nextColEup + 2] = ar_a12_f[i] * pdn_CNtF[i];

      matrixA[iRow][nextColEup + 3] = ar_a11_f[i] * pup_FtF[i];
      matrixA[iRow][nextColEup + 4] = ar_a11_f[i] * pup_CKtF[i];
      matrixA[iRow][nextColEup + 5] = ar_a11_f[i] * pup_CNtF[i];

      matrixA[iRow + 1][nextColEup]     = ar_a12_ck[i] * pdn_FtCK[i];
      matrixA[iRow + 1][nextColEup + 1] = ar_a12_ck[i] * pdn_CKtCK[i];
      matrixA[iRow + 1][nextColEup + 2] = ar_a12_ck[i] * pdn_CNtCK[i];

      matrixA[iRow + 1][nextColEup + 3] = ar_a11_ck[i] * pup_FtCK[i];
      matrixA[iRow + 1][nextColEup + 4] = ar_a11_ck[i] * pup_CKtCK[i];
      matrixA[iRow + 1][nextColEup + 5] = ar_a11_ck[i] * pup_CNtCK[i];

      matrixA[iRow + 2][nextColEup]     = ar_a12_cn[i] * pdn_FtCN[i];
      matrixA[iRow + 2][nextColEup + 1] = ar_a12_cn[i] * pdn_CKtCN[i];
      matrixA[iRow + 2][nextColEup + 2] = ar_a12_cn[i] * pdn_CNtCN[i];

      matrixA[iRow + 2][nextColEup + 3] = ar_a11_cn[i] * pup_FtCN[i];
      matrixA[iRow + 2][nextColEup + 4] = ar_a11_cn[i] * pup_CKtCN[i];
      matrixA[iRow + 2][nextColEup + 5] = ar_a11_cn[i] * pup_CNtCN[i];
    }  //e-if

    nextColEup += NFLUX;

    matrixA[iRow + 3][nextColEdw]     = ar_a11_f[i - 1] * pdn_FtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 1] = ar_a11_f[i - 1] * pdn_CKtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 2] = ar_a11_f[i - 1] * pdn_CNtF[i - 1];

    matrixA[iRow + 3][nextColEdw + 3] = ar_a12_f[i - 1] * pup_FtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 4] = ar_a12_f[i - 1] * pup_CKtF[i - 1];
    matrixA[iRow + 3][nextColEdw + 5] = ar_a12_f[i - 1] * pup_CNtF[i - 1];

    matrixA[iRow + 4][nextColEdw]     = ar_a11_ck[i - 1] * pdn_FtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 1] = ar_a11_ck[i - 1] * pdn_CKtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 2] = ar_a11_ck[i - 1] * pdn_CNtCK[i - 1];

    matrixA[iRow + 4][nextColEdw + 3] = ar_a12_ck[i - 1] * pup_FtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 4] = ar_a12_ck[i - 1] * pup_CKtCK[i - 1];
    matrixA[iRow + 4][nextColEdw + 5] = ar_a12_ck[i - 1] * pup_CNtCK[i - 1];

    matrixA[iRow + 5][nextColEdw]     = ar_a11_cn[i - 1] * pdn_FtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 1] = ar_a11_cn[i - 1] * pdn_CKtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 2] = ar_a11_cn[i - 1] * pdn_CNtCN[i - 1];

    matrixA[iRow + 5][nextColEdw + 3] = ar_a12_cn[i - 1] * pup_FtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 4] = ar_a12_cn[i - 1] * pup_CKtCN[i - 1];
    matrixA[iRow + 5][nextColEdw + 5] = ar_a12_cn[i - 1] * pup_CNtCN[i - 1];

    nextColEdw += NFLUX;
  }  //e-for

  return 0;
}  //e-buildMatrixA

//=========================
// FUNCTION buildVectorBsol
//=========================
// INPUT:  nlev, Ag, mu0, ar_a13_ck, ar_a13_cn, ar_a13_f, ar_a23_ck, ar_a23_cn, ar_a23_f, ar_S_ck, ar_S_cn, ar_S_f,
//         pdn_CKtCK, pdn_CKtCN, pdn_CKtF, pdn_CNtCK, pdn_CNtCN, pdn_CNtF, pdn_FtCK,  pdn_FtCN,  pdn_FtF;
// OUTPUT: vectB;

static int buildVectorBsol_twomaxrnd3C (int     nlev,
                                        double  Ag,
                                        double  mu0,
                                        double* ar_a13_ck,
                                        double* ar_a13_cn,
                                        double* ar_a13_f,
                                        double* ar_a23_ck,
                                        double* ar_a23_cn,
                                        double* ar_a23_f,
                                        double* ar_S_ck,
                                        double* ar_S_cn,
                                        double* ar_S_f,
                                        double* pdn_CKtCK,
                                        double* pdn_CKtCN,
                                        double* pdn_CKtF,
                                        double* pdn_CNtCK,
                                        double* pdn_CNtCN,
                                        double* pdn_CNtF,
                                        double* pdn_FtCK,
                                        double* pdn_FtCN,
                                        double* pdn_FtF,
                                        double* vectB) {
  int i;  // position in levels
  int j;  // position in vector B

  // Set initial six values:
  vectB[0] = ar_a13_f[0] * ar_S_f[0];
  vectB[1] = ar_a13_ck[0] * ar_S_ck[0];
  vectB[2] = ar_a13_cn[0] * ar_S_cn[0];
  vectB[3] = 0.0;  // upper boundary condition
  vectB[4] = 0.0;  // upper boundary condition
  vectB[5] = 0.0;  // upper boundary condition

  for (i = 1; i < (nlev - 1); i++) {
    j            = NFLUX * i;
    vectB[j]     = ar_a13_f[i] * (pdn_FtF[i] * ar_S_f[i] + pdn_CKtF[i] * ar_S_ck[i] + pdn_CNtF[i] * ar_S_cn[i]);
    vectB[j + 1] = ar_a13_ck[i] * (pdn_FtCK[i] * ar_S_f[i] + pdn_CKtCK[i] * ar_S_ck[i] + pdn_CNtCK[i] * ar_S_cn[i]);
    vectB[j + 2] = ar_a13_cn[i] * (pdn_FtCN[i] * ar_S_f[i] + pdn_CKtCN[i] * ar_S_ck[i] + pdn_CNtCN[i] * ar_S_cn[i]);
    vectB[j + 3] =
      ar_a23_f[i - 1] * (pdn_FtF[i - 1] * ar_S_f[i - 1] + pdn_CKtF[i - 1] * ar_S_ck[i - 1] + pdn_CNtF[i - 1] * ar_S_cn[i - 1]);
    vectB[j + 4] =
      ar_a23_ck[i - 1] * (pdn_FtCK[i - 1] * ar_S_f[i - 1] + pdn_CKtCK[i - 1] * ar_S_ck[i - 1] + pdn_CNtCK[i - 1] * ar_S_cn[i - 1]);
    vectB[j + 5] =
      ar_a23_cn[i - 1] * (pdn_FtCN[i - 1] * ar_S_f[i - 1] + pdn_CKtCN[i - 1] * ar_S_ck[i - 1] + pdn_CNtCN[i - 1] * ar_S_cn[i - 1]);
  }  //e-for

  // Treat last six values separately:
  j            = NFLUX * (nlev - 1);
  vectB[j]     = Ag * mu0 * ar_S_f[nlev - 1];   // lower boundary condition
  vectB[j + 1] = Ag * mu0 * ar_S_ck[nlev - 1];  // lower boundary condition
  vectB[j + 2] = Ag * mu0 * ar_S_cn[nlev - 1];  // lower boundary condition
  vectB[j + 3] = ar_a23_f[nlev - 2] * (pdn_FtF[nlev - 2] * ar_S_f[nlev - 2] + pdn_CKtF[nlev - 2] * ar_S_ck[nlev - 2] +
                                       pdn_CNtF[nlev - 2] * ar_S_cn[nlev - 2]);
  vectB[j + 4] = ar_a23_ck[nlev - 2] * (pdn_FtCK[nlev - 2] * ar_S_f[nlev - 2] + pdn_CKtCK[nlev - 2] * ar_S_ck[nlev - 2] +
                                        pdn_CNtCK[nlev - 2] * ar_S_cn[nlev - 2]);
  vectB[j + 5] = ar_a23_cn[nlev - 2] * (pdn_FtCN[nlev - 2] * ar_S_f[nlev - 2] + pdn_CKtCN[nlev - 2] * ar_S_ck[nlev - 2] +
                                        pdn_CNtCN[nlev - 2] * ar_S_cn[nlev - 2]);

  return 0;
}  //e-buildVectorBsol

//=========================
// FUNCTION buildVectorBthe
//=========================
// INPUT:  nlev, Ag, Bg, ar_theComp1_ck, ar_theComp1_cn, ar_theComp1_f,
//                       ar_theComp2_ck, ar_theComp2_cn, ar_theComp2_f, cf, cf_ck, cf_cn;
// OUTPUT: bb_the;

// thermal component 1 = upward emission
// thermal component 2 = downward emission

static int buildVectorBthe_twomaxrnd3C (int     nlev,
                                        double  Ag,
                                        double  Bg,
                                        double* ar_theComp1_ck,
                                        double* ar_theComp1_cn,
                                        double* ar_theComp1_f,
                                        double* ar_theComp2_ck,
                                        double* ar_theComp2_cn,
                                        double* ar_theComp2_f,
                                        double* cf,
                                        double* cf_ck,
                                        double* cf_cn,
                                        double* vectB) {
  int    i;  // position in levels
  int    j;  // position in vector B
  double constPi = 3.141593;

  // Set initial six values:
  vectB[0] = (1.0 - cf[0]) * ar_theComp1_f[0];
  vectB[1] = cf_ck[0] * ar_theComp1_ck[0];
  vectB[2] = cf_cn[0] * ar_theComp1_cn[0];
  vectB[3] = 0.0;  // upper boundary condition
  vectB[4] = 0.0;  // upper boundary condition
  vectB[5] = 0.0;  // upper boundary condition

  for (i = 1; i < (nlev - 1); i++) {
    j            = NFLUX * i;
    vectB[j]     = (1.0 - cf[i]) * ar_theComp1_f[i];
    vectB[j + 1] = cf_ck[i] * ar_theComp1_ck[i];
    vectB[j + 2] = cf_cn[i] * ar_theComp1_cn[i];
    vectB[j + 3] = (1.0 - cf[i - 1]) * ar_theComp2_f[i - 1];
    vectB[j + 4] = cf_ck[i - 1] * ar_theComp2_ck[i - 1];
    vectB[j + 5] = cf_cn[i - 1] * ar_theComp2_cn[i - 1];
  }  //e-for

  // Treat last six values seperately:
  // i=nlev-1; // bottom (ground) level;
  j            = NFLUX * (nlev - 1);
  vectB[j]     = (1.0 - cf[nlev - 2]) * (1.0 - Ag) * constPi * Bg;  // lower boundary condition
  vectB[j + 1] = cf_ck[nlev - 2] * (1.0 - Ag) * constPi * Bg;       // lower boundary condition
  vectB[j + 2] = cf_cn[nlev - 2] * (1.0 - Ag) * constPi * Bg;       // lower boundary condition
  vectB[j + 3] = (1.0 - cf[nlev - 2]) * ar_theComp2_f[nlev - 2];
  vectB[j + 4] = cf_ck[nlev - 2] * ar_theComp2_ck[nlev - 2];
  vectB[j + 5] = cf_cn[nlev - 2] * ar_theComp2_cn[nlev - 2];

  return 0;
}  //e-buildVectorBthe

//=====================
// FUNCTION: freeMemory
//=====================
// Function to free all allocated memory

static void freeMemory_twomaxrnd3C (int      nlev,
                                    double*  ar_a11_ck,
                                    double*  ar_a11_cn,
                                    double*  ar_a11_f,
                                    double*  ar_a12_ck,
                                    double*  ar_a12_cn,
                                    double*  ar_a12_f,
                                    double*  ar_a13_ck,
                                    double*  ar_a13_cn,
                                    double*  ar_a13_f,
                                    double*  ar_a23_ck,
                                    double*  ar_a23_cn,
                                    double*  ar_a23_f,
                                    double*  ar_a33_ck,
                                    double*  ar_a33_cn,
                                    double*  ar_a33_f,
                                    double*  pdn_CKtCK,
                                    double*  pdn_CKtCN,
                                    double*  pdn_CKtF,
                                    double*  pdn_CNtCK,
                                    double*  pdn_CNtCN,
                                    double*  pdn_CNtF,
                                    double*  pdn_FtCK,
                                    double*  pdn_FtCN,
                                    double*  pdn_FtF,
                                    double*  pup_CKtCK,
                                    double*  pup_CKtCN,
                                    double*  pup_CKtF,
                                    double*  pup_CNtCK,
                                    double*  pup_CNtCN,
                                    double*  pup_CNtF,
                                    double*  pup_FtCK,
                                    double*  pup_FtCN,
                                    double*  pup_FtF,
                                    double*  bb_sol,
                                    double*  bb_the,
                                    double*  ar_theComp1_ck,
                                    double*  ar_theComp1_cn,
                                    double*  ar_theComp1_f,
                                    double*  ar_theComp2_ck,
                                    double*  ar_theComp2_cn,
                                    double*  ar_theComp2_f,
                                    double*  S_ck,
                                    double*  S_cn,
                                    double*  S_f,
                                    double*  Edir_ck,
                                    double*  Edir_cn,
                                    double*  Edir_f,
                                    double*  Eup_ck,
                                    double*  Eup_cn,
                                    double*  Eup_f,
                                    double*  Edn_ck,
                                    double*  Edn_cn,
                                    double*  Edn_f,
                                    double*  bb,
                                    double*  xx,
                                    double** AA) {
  int i;

  free (ar_a11_ck);
  free (ar_a11_cn);
  free (ar_a11_f);
  free (ar_a12_ck);
  free (ar_a12_cn);
  free (ar_a12_f);
  free (ar_a13_ck);
  free (ar_a13_cn);
  free (ar_a13_f);
  free (ar_a23_ck);
  free (ar_a23_cn);
  free (ar_a23_f);
  free (ar_a33_ck);
  free (ar_a33_cn);
  free (ar_a33_f);

  free (pdn_CKtCK);
  free (pdn_CKtCN);
  free (pdn_CKtF);
  free (pdn_CNtCK);
  free (pdn_CNtCN);
  free (pdn_CNtF);
  free (pdn_FtCK);
  free (pdn_FtCN);
  free (pdn_FtF);

  free (pup_CKtCK);
  free (pup_CKtCN);
  free (pup_CKtF);
  free (pup_CNtCK);
  free (pup_CNtCN);
  free (pup_CNtF);
  free (pup_FtCK);
  free (pup_FtCN);
  free (pup_FtF);

  if (ar_theComp1_ck != 0)
    free (ar_theComp1_ck);
  if (ar_theComp1_cn != 0)
    free (ar_theComp1_cn);
  if (ar_theComp1_f != 0)
    free (ar_theComp1_f);
  if (ar_theComp2_ck != 0)
    free (ar_theComp2_ck);
  if (ar_theComp2_cn != 0)
    free (ar_theComp2_cn);
  if (ar_theComp2_f != 0)
    free (ar_theComp2_f);
  if (bb_sol != 0)
    free (bb_sol);
  if (bb_the != 0)
    free (bb_the);

  free (S_ck);
  free (S_cn);
  free (S_f);
  free (Edir_ck);
  free (Edir_cn);
  free (Edir_f);
  free (Eup_ck);
  free (Eup_cn);
  free (Eup_f);
  free (Edn_ck);
  free (Edn_cn);
  free (Edn_f);
  free (bb);
  free (xx);

  for (i = 0; i < NFLUX * nlev; i++)
    free (AA[i]);

}  //e-freeMemory
