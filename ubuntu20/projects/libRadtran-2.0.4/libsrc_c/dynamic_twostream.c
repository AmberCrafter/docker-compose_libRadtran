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
#include "dynamic_twostream.h"
#include "solver.h"

int dynamic_twostream (int    n_iterations,
                       float* dtau_org,
                       float* omega0_org,
                       float* g_org,
                       float* dtau_clr_org,
                       float* omega0_clr_org,
                       float* g_clr_org,
                       float* cf,
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

  /* output quantities at user levels */
  for (int lu = 0; lu < nzout; lu++) {
    fldn[lu]  = 0.0;
    flup[lu]  = 0.0;
    fldir[lu] = 0.0;
    uavg[lu]  = 0.0;
  }

  fprintf (stderr, " ... *** Richard Maier's dynamic twostream solver with %d iterations\n", n_iterations);

  return 0;
}
