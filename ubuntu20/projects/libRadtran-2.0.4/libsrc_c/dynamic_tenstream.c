/*-------------------------------------------------------------------
 * This file is part of libRadtran.
 * Copyright (c) 1997-2018 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras,
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
#include "dynamic_tenstream.h"
#include "solver.h"

int dynamic_tenstream (int                 n_iterations,
                       caoth3d_out_struct* caoth3d,
                       int                 n_caoth,
                       float*              dtau_clr_org,
                       float*              omega0_clr_org,
                       float*              g_clr_org,
                       int                 nlev,
                       double              S0,
                       double              mu0,
                       double              Ag,
                       int                 planck,
                       int                 delta,
                       int                 nzout,
                       float*              zd,
                       float*              temper,
                       float               btemp,
                       float               wvnmlo,
                       float               wvnmhi,
                       float*              zout,
                       float*              fldn,
                       float*              flup,
                       float*              fldir,
                       float*              uavg) {

  int                found3D = 0;
  caoth3d_out_struct cld3D;
  int                nlyr = nlev - 1;

  fprintf (stderr, " ... *** Richard Maier's dynamic tenstream solver with %d iterations\n", n_iterations);

  /* output quantities at user levels */
  for (int lu = 0; lu < nzout; lu++) {
    fldn[lu]  = 0.0;
    flup[lu]  = 0.0;
    fldir[lu] = 0.0;
    uavg[lu]  = 0.0;
  }

  for (int isp = MCCAOTH_FIR; isp <= n_caoth; isp++) {
    int ispo = isp - MCCAOTH_FIR;

    /* transfer indices of water and ice from coath3d to atmos structure 
     for subsequent lookups since caoth arrangement can be arbitrary */
    if (strncasecmp (caoth3d[ispo].name, "wc", 2) == 0) {

      found3D = 1;

      /* only 3D clouds with more than zero 3D layers are considered */
      if (caoth3d[ispo].nlyr == 0 || cld3D.nthreed == 0)
        found3D = 0;
      else /* assign cloud caoth to cld3D */
        cld3D = caoth3d[ispo];
    }
  }

  /* if 3D cloud was found, found3D is set to 1 and the 3D cloud is accessible in caoth3d_out_struct cld3D */

  if (!found3D)
    fprintf (stderr, " ... no 3D cloud found\n");
  else
    fprintf (stderr, " ... found 3D water cloud with %d 3D layers\n", cld3D.nthreed);

  if (found3D) {

    /* check if number of levels identical for atmosphere and 3D cloud */
    if (cld3D.nlyr + 1 != nlev) {
      fprintf (stderr,
               "Fatal error, different number of layers in atmosphere (%d) and 3D cloud profile (%d)\n",
               nlev - 1,
               cld3D.nlyr);
      return -1;
    }

    /* !!! atmospheric profile top to bottom, 3D caoth bottom to top !!! */
    //    for (int ilev=0; ilev<nlev; ilev++)
    //      fprintf (stderr, "%d %f %f\n", ilev, cld3D.zd[nlev-1-ilev], zd[ilev]);
    fprintf (stderr, "# z[km]    tau_clr    LWC reff  tau_cld   omega0     asym\n");

    for (int ilyr = 0; ilyr < nlyr; ilyr++) {
      fprintf (stderr, "%7.2f %10.6f  ", zd[ilyr + 1], dtau_clr_org[ilyr]);
      int icld = nlyr - 1 - ilyr;
      if (cld3D.threed[icld]) {
        double tau = cld3D.ext[icld][0][0] * (zd[ilyr] - zd[ilyr + 1]) * 1000.0;
        fprintf (stderr,
                 "%5.3f %4.1f  %7.2f %.6f %.6f\n",
                 cld3D.lwc[icld][0][0],
                 cld3D.reff[icld][0][0],
                 tau,
                 cld3D.ssa[icld][0][0],
                 cld3D.g1[icld][0][0]);
      } else
        fprintf (stderr, "%5.3f %4.1f  %7.2f %.6f %.6f\n", 0.0, 0.0, 0.0, 0.0, 0.0);
    }
  }

  return 0;
}
