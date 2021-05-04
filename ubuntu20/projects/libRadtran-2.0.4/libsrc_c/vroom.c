#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "common_math.h"
#include "errors.h"
#include "mystic.h"
#include "uvspecrandom.h"
#include "vroom.h"
#if q42
#include "mystic_3d.h"
#endif
#if HAVE_LIDAR
#include "lidar.h"
#endif
#ifndef PI
#define PI 3.14159265358979323846264338327
#endif
static inline void q43(scadis_struct*q30,int q44);static int q45(pft**
phase_max,int n_phase_max,scadis_struct q30,int q29,int q46,double*mu);static 
inline int q47(double q48,double*mu);static double q49(double q36,double q25,
int q29,int behind_detector,double q50,pft**phase_max,int n_phase_max,double 
q51,double q40,double q52,scadis_struct q30,int q46);int set_vroom_settings(
int vroom,sample_struct*sample,int q7){switch(vroom){case 1:sample->vroom=1;
sample->escape_eps_ddis_upf=0.1;sample->ntupelLE=22;sample->startCP=4;sample->
LEperCP=11;sample->RIS_MS=0;sample->splitter=1;sample->use_p_norm=1;sample->
split_max=3.0;sample->split_min=0.3;sample->n_split_max=5000.0;sample->
n_split_min=0.2;sample->LE_taucrit=3.0;sample->MPdynDDIS=0.1;sample->VIS_CS=1;
sample->vroomreflectalways=0;if(!q7){fprintf(stderr,"\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\124\165\162\156\151\156\147\40\157\156\40\126\122\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\117\115\56\56\56\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\50\126\141\162\151\141\156\143\145\40\122\145\144\165\143\164\151\157\156\40\117\160\164\151\155\141\154\40\117\160\164\151\157\156\163\40\115\145\164\150\157\144\51\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\111\146\40\171\157\165\40\141\162\145\40\165\163\151\156\147\40\166\162\157\157\155\54\40\160\154\145\141\163\145\40\143\151\164\145\72\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\40\40\40\122\56\40\102\165\162\141\163\40\141\156\144\40\102\56\40\115\141\171\145\162\40\50\62\60\61\61\51\40\40\40\40\40\40\40\40\40\40\40\40\40\40\52\52\52\52\52\n"
);fprintf(stderr,"\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\52\n"
);fprintf(stderr,"\n");fprintf(stderr,"\52\52\52\40\126\122\117\117\115\40\163\145\164\164\151\156\147\163\72\40\45\144\40\45\144\40\45\144\40\45\144\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\144\40\45\144\n"
,sample->ntupelLE,sample->startCP,sample->LEperCP,sample->splitter,sample->
split_max,sample->split_min,sample->n_split_max,sample->n_split_min,sample->
LE_taucrit,sample->MPdynDDIS,sample->use_p_norm,sample->VIS_CS);}break;case 0:
sample->vroom=0;sample->escape_eps_ddis_upf=0.;sample->ntupelLE=0;sample->
startCP=0;sample->LEperCP=0;sample->RIS_MS=0;sample->splitter=0;sample->
use_p_norm=0;sample->split_max=1e14;sample->split_min=0.;sample->n_split_max=
0.;sample->n_split_min=0.;sample->LE_taucrit=9999.;sample->MPdynDDIS=0.;sample
->VIS_CS=0;sample->vroomreflectalways=0;if(!q7)fprintf(stderr,"\166\162\157\157\155\40\157\146\146\41\40\163\145\164\164\151\156\147\163\72\40\45\144\40\45\144\40\45\144\40\45\144\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\145\40\45\144\40\45\144\n"
,sample->ntupelLE,sample->startCP,sample->LEperCP,sample->splitter,sample->
split_max,sample->split_min,sample->n_split_max,sample->n_split_min,sample->
LE_taucrit,sample->MPdynDDIS,sample->use_p_norm,sample->VIS_CS);break;default:
fprintf(stderr,"\123\124\117\120\41\40\131\157\165\40\141\162\145\40\165\163\151\156\147\40\115\131\123\124\111\103\163\40\154\157\143\141\154\40\145\163\164\151\155\141\164\145\40\164\145\143\150\156\151\161\165\145\54\40\141\156\144\40\171\157\165\162\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\141\162\145\40\163\160\151\153\171\54\n\40\142\165\164\40\171\157\165\40\150\141\166\145\40\156\157\164\40""\163\160\145\143\151\146\151\145\144\40\167\150\145\164\150\145\162\40\171\157\165\40\167\141\156\164\40\164\157\40\165\163\145\40\164\150\145\40\166\141\162\151\141\156\143\145\40\162\145\144\165\143\164\151\157\156\40\155\145\164\150\157\144\40\126\122\117\117\115\56\n\40\120\154\145\141\163\145\40\163\160\145\143\151\146\171\40\145\151\164\150\145\162\40\47\155\143\137\166\162\157\157\155\40\157\156\47\40\157\162\40""\47\155\143\137\166\162\157\157\155\40\157\146\146\47\40\151\156\40\171\157\165\162\40\151\156\160\165\164\40\146\151\154\145\56\n\40\47\157\156\47\40\151\163\40\162\145\143\157\155\155\145\156\144\145\144\40\146\157\162\40\171\157\165\162\40\143\165\162\162\145\156\164\40\141\160\160\154\151\143\141\164\151\157\156\56\n\105\170\151\164\151\156\147\56\56\56"
);exit(0);}return 0;}int mc_vroom_check_and_verbose(sample_struct*sample,int 
q7,int q8){
#if HAVE_LIDAR
int q53=0;
#endif
if(!(q8)){sample->vroom=0;
#if HAVE_LIDAR
if(sample->LidarLocEst&&sample->LidarLocEst!=q54){q53=set_lidar_settings(
MCLIDAR_NODDIS,sample);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\145\164\137\154\151\144\141\162\137\163\145\164\164\151\156\147\163\50\51\n"
,q53);}
#endif
}if(sample->ntupelLE&&!q7)fprintf(stderr,
"\45\144\55\164\165\160\145\154\40\114\105\n",sample->ntupelLE);if(sample->
splitter){if(!(sample->escape_eps_ddis_upf!=0.0||sample->LLE_D_DIS)){fprintf(
stderr,"\127\141\162\156\151\156\147\41\40\143\141\156\47\164\40\165\163\145\40\163\160\154\151\164\164\145\162\40\151\146\40\164\150\145\162\145\40\151\163\40\156\157\40\104\111\123\41\n"
);fprintf(stderr,"\56\56\56\56\56\56\56\56\164\165\162\156\151\156\147\40\157\146\146\40\163\160\154\151\164\164\145\162\56\56\56\56\56\56\56\56\56\56\n"
);sample->splitter=0;}else if(!q7)fprintf(stderr,
"\163\160\154\151\164\164\151\156\147\40\141\142\157\166\145\40\45\145\40\n",
sample->split_max);if(!q7)fprintf(stderr,"\155\141\170\40\163\160\154\151\164\164\151\156\147\40\142\145\154\157\167\40\45\145\40\n"
,sample->n_split_max);}if(sample->LidarLocEst){if(sample->escape){fprintf(
stderr,"\41\41\41\40\123\164\157\160\41\41\41\41\40\131\157\165\40\167\141\156\164\40\155\145\40\164\157\40\145\163\164\151\155\141\164\145\40\145\163\143\141\160\145\40\162\141\144\151\141\156\143\145\163\40\41\41\41\n"
);fprintf(stderr,"\41\41\41\40\141\156\144\40\154\157\143\141\154\40\145\163\164\151\155\141\164\157\162\40\141\164\40\164\150\145\40\163\141\155\145\40\164\151\155\145\41\41\41\40\124\150\151\163\40\151\163\40\40\40\41\41\41\n"
);fprintf(stderr,"\41\41\41\40\156\157\164\40\147\157\151\156\147\40\164\157\40\167\157\162\153\41\41\41\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\41\41\41\n"
);return-1;}}if(!sample->escape&&sample->vroom){fprintf(stderr,"\41\41\41\40\105\162\162\157\162\41\41\41\40\126\122\117\117\115\40\151\163\40\157\156\40\142\165\164\40\105\123\103\101\120\105\40\151\163\40\157\146\146\41\40\n"
);fprintf(stderr,"\41\41\41\40\123\157\155\145\164\150\151\156\147\40\151\163\40\167\162\157\156\147\41\40\103\157\156\164\141\143\164\40\164\150\145\40\144\145\166\145\154\157\160\145\162\163\40\50\143\157\144\145\40\122\102\51\41\n"
);return-1;}if(!q7&&sample->LidarLocEst){fprintf(stderr,"\40\56\56\56\40\162\165\156\156\151\156\147\40\122\125\114\105\123\40\50\114\151\144\141\162\40\105\155\165\154\141\164\157\162\40\45\144\51\40\n"
,sample->LidarLocEst);if(sample->LLE_D_DIS){fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\104\145\164\145\143\164\157\162\40\104\151\162\145\143\164\151\157\156\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\n"
);fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40\40");if(
sample->LLE_eps_ddis_upf)fprintf(stderr,"\40\165\163\151\156\147\40\120\150\141\163\145\40\106\165\156\143\164\151\157\156\40\50\45\145\51"
,sample->LLE_eps_ddis_upf);if(sample->LLE_eps_ddis_uda)fprintf(stderr,"\40\165\163\151\156\147\40\104\145\164\145\143\164\157\162\40\101\162\145\141\40\50\45\145\51"
,sample->LLE_eps_ddis_uda);if(sample->LLE_eps_fod_dis_phi)fprintf(stderr,"\40\165\163\151\156\147\40\106\151\145\154\144\55\157\146\55\166\151\145\167\40\157\146\40\104\145\164\145\143\164\157\162\40\50\45\145\51"
,sample->LLE_eps_fod_dis_phi);fprintf(stderr,"\n");if(sample->LLE_VIS_FOD)
fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\126\151\162\164\165\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\106\151\145\154\144\55\157\146\55\166\151\145\167\40\117\146\40\104\145\164\145\143\164\157\162\n"
);}if(sample->LLE_VIS_QIDD)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\126\151\162\164\165\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\121\165\141\144\162\141\164\151\143\40\111\156\166\145\162\163\145\40\104\145\164\145\143\164\157\162\40\104\151\163\164\141\156\143\145\n"
);if(sample->LLE_RIS_MAS)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\122\145\141\154\40\111\155\160\157\162\164\141\156\143\145\40\123\141\155\160\154\151\156\147\40\56\56\56\40\115\157\154\145\143\165\154\141\162\40\141\156\144\40\101\145\162\157\163\157\154\40\123\143\141\164\164\145\162\151\156\147\n"
);if(sample->LLE_sponti)fprintf(stderr,"\40\40\40\40\40\40\40\40\40\40\40\40\167\151\164\150\40\123\160\157\156\164\151\163\160\154\151\164\40\166\145\162\163\151\157\156\40\45\144\n"
,sample->LLE_sponti);
#if HAVE_LIDAR
if(sample->LLE_channels==LIDAR_CHANNEL_RAMAN)fprintf(stderr,"\40\40\40\40\40\101\154\163\157\40\143\141\154\143\165\154\141\164\151\156\147\40\122\141\155\141\156\40\114\151\144\141\162\40\143\150\141\156\156\145\154\56\56\56\n"
);if(sample->LLE_channels==LIDAR_CHANNEL_HSRL)fprintf(stderr,"\40\40\40\40\40\101\154\163\157\40\143\141\154\143\165\154\141\164\151\156\147\40\110\123\122\40\114\151\144\141\162\40\143\150\141\156\156\145\154\56\56\56\n"
);
#endif
if(sample->LLE_turnmax)fprintf(stderr,
"\40\40\40\40\40\124\165\162\156\155\141\170\151\156\147\56\56\56\n");if(
sample->LLE_taumax)fprintf(stderr,
"\40\40\40\40\40\124\141\165\155\141\170\151\156\147\40\45\145\56\56\56\n",
sample->LLE_taumax);if(sample->LE_taucrit)fprintf(stderr,
"	\40\40\40\111\167\141\142\165\143\150\151\156\147\40\45\145\56\56\56\n",
sample->LE_taucrit);fprintf(stderr,"\n");}return 0;}int mc_vroom_prepare(
sample_struct*sample,atmosphere_struct*q9,float q10,int q11,int q12,int q7){
int ic=0,jc=0,q55=0,q56=0,kc=0,q57=0,q58=0,q59=0,q60=0;double q61=0.0,q62=0.0,
q63=0.0;double q64=0.0;int q53=0;int nphamat=1;pft**q65=NULL;int*q66=NULL;int*
*q67=NULL;float***q68=NULL,***q69=NULL;double***q70=NULL;double*q71=NULL;float
**q72=NULL,**q73=NULL;double**q74=NULL;int*q75=NULL;int q31=0,q76=0,q77=0;
double q78=0.0;double*F=NULL;double q79=0.0;if(!q7)fprintf(stderr,"\52\52\104\104\111\123\72\40\104\145\146\151\156\151\156\147\40\157\160\164\151\155\141\154\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\146\157\162\40\104\104\111\123\40\56\56\56\n"
);sample->n_phase_max=0;for(q31=1;q31<=q9->n_caoth;q31++){switch(q9->
scatter_type[q31]){case MCSCAT_MOL:case MCSCAT_SKIP:break;case MCSCAT_AER:
sample->n_phase_max+=q11;break;case MCSCAT_HG1:case MCSCAT_HG2:(sample->
n_phase_max)++;break;case MCSCAT_PFT:sample->n_phase_max+=q9->phase[q31]->n;
break;default:fprintf(stderr,"\105\162\162\157\162\54\40\156\157\40\163\165\143\150\40\164\171\160\145\40\45\144\40\157\146\40\163\143\141\164\164\145\162\151\156\147\41\41\41\n"
,q9->scatter_type[q31]);return-1;}}if(sample->n_phase_max!=0){sample->
n_phase_max+=3;q65=calloc(sample->n_phase_max,sizeof(pft*));q66=calloc(sample
->n_phase_max,sizeof(int));q71=calloc(sample->n_phase_max,sizeof(double));q57=
0;for(q31=1;q31<=q9->n_caoth;q31++){switch(q9->scatter_type[q31]){case 
MCSCAT_MOL:q78=0.0;q65[q57]=calloc(1,sizeof(pft));q66[q57]=1;q53=
create_iphase_from_HG(q65[q57++],q78,q7);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\162\145\141\164\145\137\151\160\150\141\163\145\137\146\162\157\155\137\110\107\50\51\n"
,q53);if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\155\157\154\145\143\165\154\141\162\40\144\165\155\155\171\72\40\165\163\151\156\147\40\147\75\60\40\141\163\40\151\163\157\164\162\157\160\151\143\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\n"
);break;case MCSCAT_AER:for(q59=0;q59<q11;q59++)if(q9->phase_aer[q59].nphamat
!=0)q65[q57++]=&(q9->phase_aer[q59]);if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\141\154\154\40\117\120\101\103\40\141\145\162\157\163\157\154\40\154\141\171\145\162\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q9->caoth_name[q31]);break;case MCSCAT_HG1:case MCSCAT_HG2:q78=0.0;for(q59=0;
q59<q11;q59++){if(q9->threed[q31][q59]>=1)for(q76=0;q76<q9->Nx;q76++)for(q77=0
;q77<q9->Ny;q77++){if(q78<q9->g1_3D->prof[q31][q59][q76][q77])q78=q9->g1_3D->
prof[q31][q59][q76][q77];if(q78<q9->g2_3D->prof[q31][q59][q76][q77])q78=q9->
g2_3D->prof[q31][q59][q76][q77];}else{if(q78<q9->g1->prof[q31][q59])q78=q9->g1
->prof[q31][q59];if(q78<q9->g2->prof[q31][q59])q78=q9->g2->prof[q31][q59];}}
q65[q57]=calloc(1,sizeof(pft));q66[q57]=1;q53=create_iphase_from_HG(q65[q57++]
,q78,q7);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\162\145\141\164\145\137\151\160\150\141\163\145\137\146\162\157\155\137\110\107\50\51\n"
,q53);if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\155\141\170\40\162\137\145\146\146\40\110\107\40\160\150\141\163\145\40\146\143\164\56\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q9->caoth_name[q31]);break;case MCSCAT_PFT:for(q59=0;q59<q9->phase[q31]->n;
q59++)q65[q57++]=q9->phase[q31]->iphase[q59];if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\165\163\151\156\147\40\141\154\154\40\162\137\145\146\146\40\160\150\141\163\145\40\146\143\164\56\163\40\141\163\40\163\160\151\153\171\40\160\150\141\163\145\40\146\143\164\56\163\56\56\56\n"
,q9->caoth_name[q31]);break;case MCSCAT_SKIP:if(!q7)fprintf(stderr,"\52\52\104\104\111\123\40\45\163\72\40\163\153\151\160\160\151\156\147\40\50\144\165\155\155\171\51\56\56\56\n"
,q9->caoth_name[q31]);break;default:fprintf(stderr,"\105\162\162\157\162\54\40\156\157\40\163\165\143\150\40\164\171\160\145\40\45\144\40\157\146\40\163\143\141\164\164\145\162\151\156\147\41\41\41\n"
,q9->scatter_type[q31]);return-1;}}if(q57>sample->n_phase_max){fprintf(stderr,
"\45\163\45\144\45\163\45\144\45\163","\105\162\162\157\162\41\40\124\150\145\40\163\145\164\40\156\165\155\142\145\162\40\157\146\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\50"
,sample->n_phase_max,"\51\40\163\155\141\154\154\145\162\n\40\164\150\141\156\40\164\150\145\40\156\165\155\142\145\162\40\157\146\40\144\145\146\151\156\145\144\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\50"
,q57,"\51\41\n\40\103\157\156\164\141\143\164\40\122\157\142\145\162\164\40\140\115\145\163\163\171\140\40\102\165\162\141\163\40\146\157\162\40\143\157\155\160\154\141\151\156\164\56\40\105\170\151\164\151\156\147\56\56\56\n"
);return-1;}sample->n_phase_max=q57;if(!q7)fprintf(stderr,"\52\52\104\104\111\123\72\40\146\151\156\141\154\154\171\40\165\163\145\144\40\156\165\155\142\145\162\40\157\146\40\163\160\151\153\171\40\160\150\141\163\145\40\146\165\156\143\164\151\157\156\163\40\146\157\162\40\104\104\111\123\72\40\45\144\n"
,sample->n_phase_max);q75=calloc(nphamat,sizeof(int));q73=calloc(nphamat,
sizeof(float*));q74=calloc(nphamat,sizeof(double*));q72=calloc(nphamat,sizeof(
float*));q67=calloc(nphamat,sizeof(int*));q68=calloc(nphamat,sizeof(float**));
q70=calloc(nphamat,sizeof(double**));q69=calloc(nphamat,sizeof(float**));for(
q60=0;q60<nphamat;q60++){q67[q60]=calloc(sample->n_phase_max,sizeof(int));q68[
q60]=calloc(sample->n_phase_max,sizeof(float*));q70[q60]=calloc(sample->
n_phase_max,sizeof(double*));q69[q60]=calloc(sample->n_phase_max,sizeof(float*
));for(q57=0;q57<sample->n_phase_max;q57++){q67[q60][q57]=q65[q57]->n[q60];q70
[q60][q57]=q65[q57]->mu[q60];q69[q60][q57]=calloc(q67[q60][q57],sizeof(float))
;for(q59=0;q59<q67[q60][q57];q59++)q69[q60][q57][q59]=(float)q65[q57]->p[
MCSC_MODE_NORMAL][q60][q59];q68[q60][q57]=calloc(q67[q60][q57],sizeof(double))
;for(q59=0;q59<q67[q60][q57];q59++)q68[q60][q57][q59]=acos(q70[q60][q57][q59])
;}}q53=sort_and_add_weighted_phase(sample->n_phase_max,q71,q67,q68,q70,q69,&(
q75),&(q73),&(q74),&(q72),nphamat,1,q7);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\157\162\164\137\141\156\144\137\141\144\144\137\167\145\151\147\150\164\145\144\137\160\150\141\163\145\50\51\n"
,q53);sample->phase_max=calloc(1,sizeof(pft*));sample->phase_max[0]=calloc(1,
sizeof(pft));q53=calc_cumulative_table(q74,q72,q75,nphamat,-1.0,sample->
phase_max[0],q10,q7);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\143\165\155\165\154\141\164\151\166\145\137\164\141\142\154\145\50\51\n"
,q53);if(q9->nscaDS>1){for(q57=0;q57<sample->n_phase_max;q57++)for(q59=0;q59<
q67[0][q57];q59++)q69[0][q57][q59]=(float)q65[q57]->p[MCSC_MODE_DELTA_SCALE*(
q65[q57]->nscales>1)][0][q59];q53=sort_and_add_weighted_phase(sample->
n_phase_max,q71,q67,q68,q70,q69,&(q75),&(q73),&(q74),&(q72),1,1,q7);if(q53!=0)
return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\157\162\164\137\141\156\144\137\141\144\144\137\167\145\151\147\150\164\145\144\137\160\150\141\163\145\50\51\n"
,q53);for(q58=0;q58<q75[0];q58++)if(q74[0][q58]>q10)break;for(q59=q58+1;q59<
q75[0];q59++)if(q72[0][q59]<q72[0][q58])q72[0][q59]=q72[0][q58];F=calloc(q75[0
],sizeof(double));normalize_phase(q74,q72,F,q75,nphamat,!q7);for(q57=0;q57<
sample->phase_max[0]->n[0];q57++){sample->phase_max[0]->p[
MCSC_MODE_DELTA_SCALE][0][q57]=q72[0][q57];sample->phase_max[0]->F[
MCSC_MODE_DELTA_SCALE][q57]=F[q58];}calc_iphase_coeffs(sample->phase_max[0],
MCSC_MODE_DELTA_SCALE);sample->phase_max[0]->dscale=-999.0;free(F);}for(q57=0;
q57<sample->n_phase_max;q57++)if(q66[q57]){free_iphase(q65[q57]);free(q65[q57]
);}free(q65);free(q66);free(q71);for(q58=0;q58<nphamat;q58++){free(q73[q58]);
free(q74[q58]);free(q72[q58]);}free(q75);free(q73);free(q74);free(q72);for(q60
=0;q60<nphamat;q60++){for(q59=0;q59<sample->n_phase_max;q59++){free(q69[q60][
q59]);free(q68[q60][q59]);}free(q67[q60]);free(q68[q60]);free(q70[q60]);free(
q69[q60]);}free(q67);free(q68);free(q70);free(q69);sample->n_phase_max=1;q79=
0.0;for(q57=0;q57<sample->phase_max[0]->n[0];q57++){if(q79<sample->phase_max[0
]->p[MCSC_MODE_NORMAL][0][q57])q79=sample->phase_max[0]->p[MCSC_MODE_NORMAL][0
][q57];if(q79<1./sample->phase_max[0]->p[MCSC_MODE_NORMAL][0][q57])q79=1./
sample->phase_max[0]->p[MCSC_MODE_NORMAL][0][q57];}if(q79<q2){sample->
n_phase_max=0;free(sample->phase_max[0]);free(sample->phase_max);if(!q7)
fprintf(stderr,"\52\52\116\157\40\156\145\145\144\40\146\157\162\40\104\104\111\123\54\40\164\165\162\156\151\156\147\40\157\146\146\40\126\122\117\117\115\40\50\151\156\40\143\141\163\145\40\151\164\40\167\141\163\40\157\156\51\n"
);}}if(sample->n_phase_max!=0){q53=calloc_hybrid3D_field(&(q9->spiky_box),q9);
if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\154\157\143\137\150\171\142\162\151\144\63\104\137\146\151\145\154\144\50\51\n"
,q53);for(kc=0;kc<q9->Nz;kc++){if(q9->threed[MCCAOTH_TOT][kc]>=1){for(q31=0;
q31<q9->n_caoth;q31++){if(q9->threed[q31][kc]>=1){for(ic=0;ic<q9->Nx;ic++)for(
jc=0;jc<q9->Ny;jc++)if(q9->ksca3D[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[
q31][kc][ic][jc]>0.0)q9->spiky_box[kc][ic][jc]=1;}else if(q9->ksca[
MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[q31][kc]>0.0)for(ic=0;ic<q9->Nx;ic
++)for(jc=0;jc<q9->Ny;jc++)q9->spiky_box[kc][ic][jc]=1;}}else{for(q31=0;q31<q9
->n_caoth;q31++){if(q9->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[q31][
kc]>0.0)q9->spiky_box[kc][0][0]=1;}}}}else{sample->vroom=0;
#if HAVE_LIDAR
if(sample->LidarLocEst&&sample->LidarLocEst!=q54){q53=set_lidar_settings(
MCLIDAR_NODDIS,sample);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\163\145\164\137\154\151\144\141\162\137\163\145\164\164\151\156\147\163\50\51\n"
,q53);}
#endif
}
#ifdef NEWVISQIDD
if(sample->LLE_VIS_QIDD){q61=10.;sample->visqidd_betamax=-q61/sample->lidar[
sample->ili].z_det;sample->visqidd_rmax=-sample->lidar[sample->ili].z_det;
sample->visqidd_facs=sqrt(sample->visqidd_betamax*sample->visqidd_rmax*sample
->visqidd_rmax);fprintf(stderr,"\126\111\123\40\142\145\164\141\155\141\170\40\45\145\40\162\155\141\170\40\45\145\40\146\141\143\40\45\145\n"
,sample->visqidd_betamax,sample->visqidd_rmax,sample->visqidd_facs);
#ifdef NEWRISQIDD
q61=1.;sample->risqidd_betamax=-q61/sample->lidar[sample->ili].z_det;sample->
risqidd_facs=sample->risqidd_betamax*sample->visqidd_rmax*sample->visqidd_rmax
;fprintf(stderr,"\122\111\123\40\142\145\164\141\155\141\170\40\45\145\40\162\155\141\170\40\45\145\40\146\141\143\40\45\145\n"
,sample->risqidd_betamax,sample->visqidd_rmax,sample->risqidd_facs);
#endif
}
#else
if(sample->LLE_VIS_QIDD){q61=0000.0;q62=-q61/sample->lidar[sample->ili].z_det;
for(kc=0;kc<q9->Nz;kc++){q63=(0.5*(q9->Z[kc]+q9->Z[kc+1])-sample->lidar[0].x[2
])/(-sample->lidar[0].dir.dx[2]*sample->lidar[sample->ili].z_det);if(q63>0.0){
if(q63<1.0)q64=q62;else q64=q62/(q63*q63);if(q9->threed[MCCAOTH_TOT][kc]>=1){
for(ic=0;ic<q9->Nx;ic++)for(jc=0;jc<q9->Ny;jc++)for(q55=0;q55<q9->nscaDS;q55++
)for(q56=0;q56<q9->nscaRIS;q56++)if(q64>q9->kext3D[q55][q56][MCVIS_MODE_QIDD]
->prof[MCCAOTH_TOT][kc][ic][jc])q9->kext3D[q55][q56][MCVIS_MODE_QIDD]->prof[
MCCAOTH_TOT][kc][ic][jc]=q64;}else{for(q55=0;q55<q9->nscaDS;q55++)for(q56=0;
q56<q9->nscaRIS;q56++)if(q64>q9->kext[q55][q56][MCVIS_MODE_QIDD]->prof[
MCCAOTH_TOT][kc])q9->kext[q55][q56][MCVIS_MODE_QIDD]->prof[MCCAOTH_TOT][kc]=
q64;}}}}
#endif
return 0;}int mc_vroom_cloning(photon_struct*p,atmosphere_struct*q9,
sample_struct*sample,result_struct*q13,elevation_struct*q14,albedo_struct*q15,
t_triangular_surface*q16,surftemp_struct*q17,int*q18,int*q19,int q20,int q21,
float*q22,float*q23,float*q24,int q7){int q53=0;photon_struct*q80=NULL;if(
sample->ntupelLE>1&&!p->isclone&&p->scattercounter>=sample->startCP&&p->
escapescattercounter+sample->LEperCP<p->scattercounter+sample->ntupelLE){q80=
calloc_photon(sample,q9->Nx,q9->Ny,q9->Nz,*q19,q9->nlambda_abs,q9->Nc,q9->
n_caoth);
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\143\154\157\156\151\156\147\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode);
#endif
cp_photon_struct(q80,p,sample,q9->n_caoth);q80->isclone=1;q80->wtree=p->wtree;
q80->photon_status=MCSTATUS_SPLIT;q53=photon_journey(q80,q9,sample,q13,q14,q15
,q16,q17,q18,q19,q20,q21,q22,q23,q24,q7,"");if(q53<0){fprintf(stderr,"\105\162\162\157\162\40\45\144\40\157\143\143\165\162\151\156\147\40\151\156\40\160\150\157\164\157\156\137\152\157\165\162\156\145\171\50\51\40\146\157\162\40\160\150\157\164\157\156\40\45\144\40\50\143\154\157\156\145\51\54\40\145\170\151\164\151\156\147\56\56\56\n"
,q53,p->photoncounter);return-1;}destroy_photon(q80,q9->n_caoth);p->
escapescattercounter=p->scattercounter+sample->ntupelLE-1;}return 0;}int 
mc_vroom_splitting_and_rr(photon_struct*p,atmosphere_struct*q9,sample_struct*
sample,result_struct*q13,elevation_struct*q14,albedo_struct*q15,
t_triangular_surface*q16,surftemp_struct*q17,int*q18,int*q19,int q20,int q21,
float*q22,float*q23,float*q24,int q7){double q25=0.0,q83=0.0,q84=0.0;int q53=0
,q85=0,q86=0;
#if HAVE_LIDAR
locest_struct lest;
#endif
photon_struct*q87=NULL;double n_split_max=0.0,n_split_min=0.0;int q88=0;if(
sample->LLE_sponti){if(p->scattercounter<7&&sample->LLE_sponti==1)p->
special_weight*=1.5;else p->special_weight*=1.0+0.1/(p->scattercounter+1);}if(
sample->splitter){if(sample->escape)v_mult_mu(p->dir.dx,sample->rad[0].dir.dx,
&q25);
#if HAVE_LIDAR
if(sample->LidarLocEst){q53=calc_locest_connection(p->x,sample->lidar[sample->
ili].dir.dx,sample->lidar[sample->ili].x,&lest);if(q53<0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\154\157\143\145\163\164\137\143\157\156\156\145\143\164\151\157\156\50\51\n"
,q53);v_mult_mu(p->dir.dx,lest.dir.dx,&q25);}
#endif
q83=get_phase_max(sample->phase_max,sample->n_phase_max,q25,p->DDIS_SC_mode);
if(sample->use_p_norm){q83/=p->p_norm;if(p->p_norm>1.0){if(q83<1.0)p->p_norm*=
q83;if(p->p_norm<1.0)p->p_norm=1.0;}}q84=p->special_weight*p->weight*p->stokes
[0]*exp(-p->tauris)*q83;
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\164\145\163\164\40\163\160\154\151\164\164\151\156\147\72\40\162\145\163\164\40\45\145\40\167\145\151\147\150\164\40\45\145\40\160\155\141\170\40\45\145\40\160\156\157\162\155\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p
->special_weight*p->stokes[0]*exp(-p->tauris),p->weight,q83*p->p_norm,p->
p_norm);
#endif
n_split_max=sample->n_split_max;n_split_min=sample->n_split_min;
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\164\145\163\164\40\163\160\154\151\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q84,n_split_max);
#endif
if(q84>sample->split_max&&!(sample->ntupelLE&&!p->isclone&&p->scattercounter>=
sample->startCP)){
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\163\160\154\151\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q84,n_split_max);
#endif
q84=q84>n_split_max?n_split_max:q84;q85=(int)q84;p->weight/=(double)q85;q87=
calloc_photon(sample,q9->Nx,q9->Ny,q9->Nz,*q19,q9->nlambda_abs,q9->Nc,q9->
n_caoth);
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\163\160\154\151\164\145\163\164\164\164\151\156\147\72\40\144\137\163\160\154\151\164\40\45\145\40\156\163\160\154\151\164\137\155\141\170\40\45\144\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q84,q85);
#endif
if(q85>1000){p->spikewarningcounter++;if(p->spikewarningcounter>1){fprintf(
stderr,"\127\141\162\156\151\156\147\41\40\123\160\151\153\145\40\167\141\162\156\151\156\147\40\154\145\166\145\154\40\45\144\40\141\164\40\143\154\157\156\145\40\163\143\141\164\164\145\162\40\157\162\144\145\162\40\45\144\40\163\143\141\164\164\145\162\40\45\144\40\160\150\157\164\157\156\40\45\144\n"
,p->spikewarningcounter,p->clonescattercounter,p->scattercounter,p->
photoncounter);fprintf(stderr,"\40\161\137\163\160\40\45\145\40\167\40\45\145\40\111\60\40\45\145\40\145\170\160\50\55\164\141\165\51\40\45\145\40\120\40\45\145\n"
,p->special_weight,p->weight*(double)q85,p->stokes[0],exp(-p->tauris),q83);}
q88=1;}for(q86=0;q86<q85-1;q86++){cp_photon_struct(q87,p,sample,q9->n_caoth);
q87->wtree=p->wtree;p->photon_status=MCSTATUS_TRAVEL;
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\143\157\160\171\156\165\155\142\145\162\72\40\156\137\163\160\154\151\164\40\45\144\40\43\40\45\144\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
q85,q86);
#endif
q53=photon_journey(q87,q9,sample,q13,q14,q15,q16,q17,q18,q19,q20,q21,q22,q23,
q24,q7,"");if(q53<0){fprintf(stderr,"\105\162\162\157\162\40\45\144\40\157\143\143\165\162\151\156\147\40\151\156\40\160\150\157\164\157\156\137\152\157\165\162\156\145\171\50\51\40\146\157\162\40\160\150\157\164\157\156\40\45\144\40\50\163\160\154\151\164\40\45\144\51\54\40\145\170\151\164\151\156\147\56\56\56\n"
,q53,p->photoncounter,q86);return-1;}}if(q88){if(p->spikewarningcounter>1)
fprintf(stderr,"\114\145\141\166\151\156\147\40\163\160\151\153\145\40\167\141\162\156\151\156\147\40\154\145\166\145\154\40\45\144\n"
,p->spikewarningcounter);p->spikewarningcounter--;}destroy_photon(q87,q9->
n_caoth);}if(q84<p->weight)q84=p->weight;if(q84<sample->split_min){q84=q84<
n_split_min?n_split_min:q84;if(uvspec_random()>q84)return MCSTATUS_PURGE;p->
weight/=q84;}}return MCSTATUS_DEFAULT;}double get_phase_max(pft**phase_max,int
 n_phase_max,double q25,int SC_mode){int q86=0,q53=0;double q83=0.0,phase=0.0;
for(q86=0;q86<n_phase_max;q86++){q53=get_phase_matrix_pft(phase_max[q86],q25,
SC_mode,1,&phase);if(q53){fct_err_out(q53,"\147\145\164\137\160\150\141\163\145\137\155\141\164\162\151\170\137\160\146\164"
,ERROR_POSITION);return-1.0;}q83+=phase;}return q83/((double)n_phase_max);}
void cp_locest(locest_struct*q26,locest_struct*q19,sample_struct*sample,int 
n_caoth){int q86=0,q31=0,q60=0;cp_direction(&(q26->dir),&(q19->dir));q26->
cosalpha=q19->cosalpha;q26->pdir=q19->pdir;q26->pdir_iso=q19->pdir_iso;if(q26
->pdir_sct!=NULL)for(q31=0;q31<n_caoth+1;q31++)for(q60=0;q60<sample->nstokes;
q60++)q26->pdir_sct[q31][q60]=q19->pdir_sct[q31][q60];q26->dist=q19->dist;q26
->distinv=q19->distinv;q26->r_det=q19->r_det;q26->z_det=q19->z_det;for(q86=0;
q86<3;q86++)q26->x_cc[q86]=q19->x_cc[q86];q26->t_det=q19->t_det;for(q86=0;q86<
3;q86++)q26->x_hit[q86]=q19->x_hit[q86];q26->weight_hit=q19->weight_hit;q26->
in_cone=q19->in_cone;q26->will_hit_cone=q19->will_hit_cone;q26->
behind_detector=q19->behind_detector;q26->will_hit_det_plane=q19->
will_hit_det_plane;q26->lidar_outside_grid=q19->lidar_outside_grid;q26->
hit_det_plane_step=q19->hit_det_plane_step;q26->vis_fod_step=q19->vis_fod_step
;q26->vis_fod_step2=q19->vis_fod_step2;q26->vis_fod_kext=q19->vis_fod_kext;for
(q86=0;q86<3;q86++)q26->hitpoint[q86]=q19->hitpoint[q86];}static inline void 
q43(scadis_struct*q30,int q44){int q86=0;q30->q3=0.0;q30->q4=0.0;q30->q5=2.0;
for(q86=0;q86<3;q86++)q30->dirold_dx[q86]=0.0;q30->d_phi=0.0;q30->epsfac=1.0;
q30->mu_max=1.0;q30->mu_min=-1.0;for(q86=0;q86<q44;q86++)q30->q6[q86]=2.0;for(
q86=0;q86<q44;q86++)q30->F_min[q86]=0.0;}int mc_vroom_prep_DDIS(sample_struct*
sample,photon_struct*p,atmosphere_struct*q9,int*q27,int*q28,int*q29,
locest_struct*lest,scadis_struct*q30){int q86=0;
#if HAVE_LIDAR
int q53=0;double q89=0.0,q90=0.0,q91=0.0;double q92=0.0,q93=0.0;double q94=0.0
;double q95=0.0;
#endif
q43(q30,sample->n_phase_max);if(sample->escape_eps_ddis_upf!=0){if(sample->
ntupelLE&&!p->isclone&&q30->epsfac){if(p->scattercounter<sample->startCP)q30->
epsfac=sample->MPdynDDIS/sample->escape_eps_ddis_upf;else q30->epsfac=0.0;}if(
uvspec_random()<sample->escape_eps_ddis_upf*q30->epsfac)*q27=MCDDIS_UPF;for(
q86=0;q86<3;q86++)lest->dir.dx[q86]=sample->rad[0].dir.dx[q86];}
#if HAVE_LIDAR
if(sample->LLE_D_DIS){q53=calc_locest_connection(p->x,sample->lidar[sample->
ili].dir.dx,sample->lidar[sample->ili].x_cc,lest);if(q53!=0)return err_out("\105\162\162\157\162\40\45\144\40\162\145\164\165\162\156\145\144\40\142\171\40\143\141\154\143\137\154\157\143\145\163\164\137\143\157\156\156\145\143\164\151\157\156\50\51\n"
,q53);q30->q3=-lest->dist*lest->cosalpha;q30->q4=lest->dist*sqrt(1.-lest->
cosalpha*lest->cosalpha);q89=1./(q30->q3-sample->lidar[sample->ili].z_det);if(
!p->lest.behind_detector){q92=-(q30->q4+sample->lidar[sample->ili].q96)*q89;
q93=acos(lest->cosalpha);q30->q5=cos(atan(q92)-q93);}if(q9->nthreed==0){q94=(p
->x[2]-q9->Z[p->kc])/(q9->Z[p->kc+1]-q9->Z[p->kc]);q30->epsfac*=q94*q9->q97[p
->kc+1]+(1.-q94)*q9->q97[p->kc];if(q30->epsfac<0.)q30->epsfac=0.0;}if(sample->
lidar[sample->ili].cosalpha[0]<lest->cosalpha){if(q30->q3<sample->lidar[sample
->ili].z_det){if(p->lest.behind_detector){fprintf(stderr,"\127\141\162\156\151\156\147\41\40\114\157\147\151\143\141\154\40\163\141\171\163\40\160\150\157\164\157\156\40\151\163\40\142\145\150\151\156\144\40\144\145\164\145\143\164\157\162\54\40\142\165\164\40\151\164\40\151\163\40\151\156\40\146\162\157\156\164\41\n"
);fprintf(stderr,
"\154\157\143\141\164\151\157\156\40\45\145\40\45\145\40\45\145\40\n",p->x[0],
p->x[1],p->x[2]);return-1;}q95=uvspec_random();if(q95<(sample->
LLE_eps_ddis_upf+sample->LLE_eps_ddis_uda)*q30->epsfac)*q27=MCDDIS_UPF;if(q95<
sample->LLE_eps_ddis_uda*q30->epsfac)*q27=MCDDIS_UDA;}else{if(!p->lest.
behind_detector)fprintf(stderr,"\127\141\162\156\151\156\147\41\40\114\157\147\151\143\141\154\40\163\141\171\163\40\160\150\157\164\157\156\40\151\163\40\151\156\40\146\162\157\156\164\40\157\146\40\144\145\164\145\143\164\157\162\54\40\142\165\164\40\151\164\40\151\163\40\142\145\150\151\156\144\41\n"
);*q29=-1;}}else{if(sample->lidar[sample->ili].cosalpha[0]<-lest->cosalpha){*
q29=-1;}else{*q29=1;q95=uvspec_random();if(q95<(sample->LLE_eps_ddis_upf+
sample->LLE_eps_ddis_uda)*q30->epsfac)*q27=MCDDIS_UPF;if(!p->lest.
behind_detector&&q95<sample->LLE_eps_ddis_uda*q30->epsfac)*q27=MCDDIS_UDA;if(*
q27&&uvspec_random()<sample->LLE_eps_fod_dis_phi)*q28=1;q91=-q30->q4/q30->q3;
if(q89>0.)q90=-(q30->q4+sample->lidar[sample->ili].q96)*q89;else q90=-(q30->q4
-sample->lidar[sample->ili].q96)*q89;q94=atan((q90-q91)/(1.0+q90*q91));q30->
mu_max=cos(q94);if(q94<0.)q30->mu_max=-q30->mu_max;if(q94<0.&&q94>-1e-10){q94=
0.0;q30->mu_max=1.0;}q30->mu_min=(q30->q3*sample->lidar[sample->ili].q98-q30->
q4*sample->lidar[sample->ili].q99)*lest->distinv;if(q30->mu_min<-1.0){if(q30->
mu_min<-1.0-q100){fprintf(stderr,"\105\162\162\157\162\41\40\155\165\137\155\151\156\40\75\40\45\145\40\151\156\40\163\143\141\164\164\145\162\151\156\147\50\51\40\151\163\40\165\156\160\150\171\163\151\143\141\154\41\n"
,q30->mu_min);return-1;}q30->mu_min=-1.0;}
#ifdef NOFODDIS
q30->mu_max=1.0;q30->mu_min=-1.0;
#endif
for(q86=0;q86<sample->n_phase_max;q86++){q30->q6[q86]=q101(sample->phase_max[
q86],q30->mu_max,p->DDIS_SC_mode);q30->F_min[q86]=q101(sample->phase_max[q86],
q30->mu_min,p->DDIS_SC_mode);}}}}if(*q29==1){q30->d_phi=sample->lidar[sample->
ili].q99/q30->q4*lest->dist;if(q30->d_phi>1.0)q30->d_phi=1.0;q30->d_phi=q102(
q30->d_phi);}else q30->d_phi=180.;
#endif
return 0;}int mu_scatter_special(atmosphere_struct*q9,photon_struct*p,pft**
phase_max,int n_phase_max,int q31,double*mu,scadis_struct q30,int q27,int q29)
{int q53=0;switch(q27){case MCDDIS_NONE:q53=mu_scatter(q9,p,q31,mu);if(q53!=0)
return err_out("\105\162\162\157\162\54\40\155\165\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);break;case MCDDIS_UPF:q53=q45(phase_max,n_phase_max,q30,q29,p->
DDIS_SC_mode,mu);if(q53!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);break;case MCDDIS_UDA:q53=q47(q30.q5,mu);if(q53!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\104\101\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);break;default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q27);return-1;}return 0;}int random_reflection_special(sample_struct*sample,
atmosphere_struct*q9,elevation_struct*q14,photon_struct*p,pft**phase_max,int 
n_phase_max,double*mu,double*phi,double*q32,scadis_struct*q30,locest_struct 
lest,int q33,int q27,int q28,int q29){int q53=0,q86=0;double q103=0.0;double 
q104=0.0;for(q86=0;q86<3;q86++)q30->dirold_dx[q86]=p->dir.dx[q86];switch(q27){
case MCDDIS_NONE:if(q33==1)random_Lambertian_normal(&(p->dir),q32);else 
random_Isotropic_normal(&(p->dir),q32);v_mult_mu(p->dir.dx,q30->dirold_dx,mu);
break;case MCDDIS_UPF:q53=q45(phase_max,n_phase_max,*q30,q29,p->DDIS_SC_mode,
mu);if(q53!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\120\106\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);break;case MCDDIS_UDA:q53=q47(q30->q5,mu);if(q53!=0)return err_out("\105\162\162\157\162\54\40\155\165\137\104\104\111\123\137\125\104\101\137\163\143\141\164\164\145\162\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);break;default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q27);return-1;}switch(q27){case MCDDIS_NONE:break;case MCDDIS_UPF:case 
MCDDIS_UDA:if(q28)*phi=q30->d_phi*(2.0*uvspec_random()-1.0);else*phi=
sc_Isotropic_phi();if(p->scattercounter==0)q103=p->phi0;else q103=0.0;for(q86=
0;q86<3;q86++)p->dir.dx[q86]=lest.dir.dx[q86];new_direction(*mu,*phi-90.,&(p->
dir),q103);v_mult_mu(p->dir.dx,q32,&q104);if(q104<0.0){p->weight=0.0;}break;
default:fprintf(stderr,"\105\110\110\117\122\41\40\105\163\164\145\40\155\157\144\157\40\45\144\40\144\157\40\104\104\111\123\163\151\156\147\40\141\151\156\144\141\40\156\141\157\40\145\170\151\163\164\145\41\n"
,q27);return-1;}return 0;}static int q45(pft**phase_max,int n_phase_max,
scadis_struct q30,int q29,int q46,double*mu){int iphase=0;iphase=(int)(
uvspec_random()*((double)n_phase_max-1e-11));if(q30.mu_min==-1.0&&q30.mu_max==
1.0){*mu=sc_mu(phase_max[iphase],1,q46,0.,0.);}else{*mu=sc_mu(phase_max[iphase
],1,q46,q30.q6[iphase],q30.F_min[iphase]);}return 0;}static int inline q47(
double q48,double*mu){double q105=0.0,q106=0.0,q95=0.0;q95=2.*uvspec_random();
q105=q48*q48;q106=2./(1.-q105);if(q95<q106*(q48-q105)){if(q95==0.)*mu=0.;else*
mu=1./(1.+(1.-q48)*(1.-q48)*q106/q95);}else*mu=q105+q95/q106;return 0;}int 
mc_vroom_set_mus_and_phis(sample_struct*sample,photon_struct*p,int q27,
locest_struct lest,scadis_struct q30,double*mu,double*q25,double*q39,double 
phi,double*q40){
#if HAVE_LIDAR
int q53=0;
#endif
switch(q27){case MCDDIS_NONE:v_mult_mu(lest.dir.dx,p->dir.dx,q25);
#if HAVE_LIDAR
if(sample->LLE_VIS_FOD){*q39=sqrt(1.0-*q25**q25);q53=q107(p->dir.dx,lest.dir.
dx,*q39,q40);if(q53)return err_out("\105\122\122\117\122\41\40\144\145\162\151\166\145\137\143\160\150\151\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);}
#endif
break;case MCDDIS_UPF:case MCDDIS_UDA:*q25=*mu;v_mult_mu(q30.dirold_dx,p->dir.
dx,mu);if(sample->LLE_VIS_FOD){*q39=sqrt(1.0-*q25**q25);*q40=cosd(phi);}break;
default:fprintf(stderr,"\105\122\122\117\122\41\40\104\104\111\123\163\151\156\147\40\147\151\166\145\163\40\163\157\155\145\164\150\151\156\147\40\163\164\162\141\156\147\145\41\40\45\144\n"
,q27);return-1;}return 0;}int mc_vroom_scattering_calc_phases_and_jacobians(
sample_struct*sample,photon_struct*p,atmosphere_struct*q9,double q34,int q35,
double*q36,double*q37,double*q38){static double**q108=NULL;int q31=0;int q53=0
;if(q35==1){if(q108!=NULL){for(q31=0;q31<=q9->n_caoth;q31++)free(q108[q31]);
free(q108);q108=NULL;}return 0.0;}if(sample->LLE_jacobian||sample->jacobian||
sample->jacobian3D)if(q108==NULL){q108=calloc((size_t)q9->n_caoth+1,sizeof(
double*));for(q31=0;q31<=q9->n_caoth;q31++)q108[q31]=calloc(1,sizeof(double));
}q53=get_phase_matrix_total(q9,p,q34,1,0,sample->spectral_is,sample->
concentration_is,q9->ris_factor,0,q36,q108,q37,&(p->weight));if(q53)return 
fct_err_out(q53,"\147\145\164\137\160\150\141\163\145\137\155\141\164\162\151\170\137\164\157\164\141\154"
,ERROR_POSITION);if(sample->LLE_jacobian)*q38=*q36-q108[MCCAOTH_MOL][0];if(*
q36<=0.0||*q37<=0.0||*q38<0.0){fprintf(stderr,"\105\162\162\157\162\54\40\143\141\154\143\165\154\141\164\151\157\156\40\157\146\40\120\137\156\157\162\155\54\40\120\137\163\160\145\143\54\40\141\156\144\40\120\137\151\163\157\145\156\145\40\144\151\144\40\156\157\164\40\167\157\162\153\41\41\41\n"
);fprintf(stderr,"\120\137\156\157\162\155\40\45\145\40\120\137\163\160\145\143\40\45\145\40\120\137\151\163\157\145\156\145\40\45\145\40\n"
,*q36,*q37,*q38);return-1;}if(sample->LLE_jacobian){if((*q36)!=0.0)for(q31=1;
q31<q9->n_caoth+1;q31++)p->q_jacobian[0][q31-1][p->kc]+=q108[q31][0]/(*q36);if
((*q38)!=0.0)for(q31=1;q31<q9->n_caoth+1;q31++)p->q_jacobian[1][q31-1][p->kc]
+=q108[q31][0]/(*q38);}if(sample->jacobian||sample->jacobian3D){for(q31=1;q31<
q9->n_caoth+1;q31++)if(*q36!=0.0){p->q_jacobian[0][q31-1][p->kc]+=q108[q31][0]
/(*q36);if(sample->vroom)p->q_jacobian_sca[q31-1][p->ic][p->jc][p->kc]+=q108[
q31][0]/(*q36)/get_ksca(q9,p,q31);}}return 0;}int 
mc_vroom_DDIS_weight_and_prep_stuff(sample_struct*sample,atmosphere_struct*q9,
double q25,double q39,double q40,int q27,int q29,double q36,double q37,double 
q38,locest_struct lest,scadis_struct q30,photon_struct*p){double cosalpha=0.0;
double q109=q37;
#if HAVE_LIDAR
int q53=0;
#endif
if((p->RIS_mode!=MCRIS_MODE_NORMAL)||sample->LLE_D_DIS||sample->
escape_eps_ddis_upf!=0.0||sample->LLE_channels||q9->ris_factor!=1.){if(sample
->escape_eps_ddis_upf!=0.0)q109=q49(q37,q25,0,0,sample->escape_eps_ddis_upf*
q30.epsfac,sample->phase_max,sample->n_phase_max,0,0,0,q30,p->DDIS_SC_mode);
#if HAVE_LIDAR
if(sample->LLE_D_DIS)q109=q49(q37,q25,q29,p->lest.behind_detector,sample->
LLE_eps_ddis_upf*q30.epsfac,sample->phase_max,sample->n_phase_max,sample->
LLE_eps_fod_dis_phi,q40,sample->LLE_eps_ddis_uda*q30.epsfac,q30,p->
DDIS_SC_mode);
#endif
if(!(q109>0.0)){fprintf(stderr,"\105\162\162\157\162\54\40\143\141\154\143\165\154\141\164\151\157\156\40\157\146\40\120\137\163\160\145\143\137\104\40\144\151\144\40\156\157\164\40\167\157\162\153\41\41\41\n"
);fprintf(stderr,"\120\137\163\160\145\143\137\104\40\45\145\40\120\137\163\160\145\143\40\45\145\40\120\137\156\157\162\155\40\45\145\40\n"
,q109,q37,q36);return-1;}p->weight*=q36/q109;
#ifdef MUCHOUT
if(p->q81==q82)fprintf(stderr,"\143\157\165\156\164\145\162\163\40\45\144\40\45\144\40\45\144\40\45\144\40\55\55\55\40\156\145\167\40\167\145\151\147\150\164\40\45\145\40\146\162\157\155\40\45\145\40\57\40\45\145\40\167\151\164\150\40\155\165\62\40\45\145\40\141\156\144\40\120\163\160\145\143\40\45\145\n"
,p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p
->weight,q36,q109,q25,q37);
#endif
p->q_isoene*=q38/q36;}if(sample->LLE_VIS_FOD||sample->LLE_eps_ddis_uda){
v_mult_mu(p->dir.dx,sample->lidar[sample->ili].dir.dx,&cosalpha);p->lest.
hit_det_plane_step=(q30.q3-sample->lidar[sample->ili].z_det)/cosalpha;p->lest.
will_hit_det_plane=(p->lest.hit_det_plane_step>0.0);}
#if HAVE_LIDAR
if(sample->LLE_VIS_FOD){q53=q110(q25,q39,q40,q30.q4,q30.q3,lest.distinv,sample
->lidar[sample->ili].t_det,&p->lest);if(q53!=0)return err_out("\105\162\162\157\162\54\40\143\141\154\143\137\144\151\163\164\141\156\143\145\137\164\157\137\143\157\156\145\50\51\40\162\145\164\165\162\156\145\144\40\163\164\141\164\165\163\40\45\144\n"
,q53);}
#endif
return 0;}static double q49(double q36,double q25,int q29,int behind_detector,
double q50,pft**phase_max,int n_phase_max,double q51,double q40,double q52,
scadis_struct q30,int q46){double q37=0.0;double q111=0.0;double q112=0.0;
double q113=0.0;int q86=0;if(behind_detector){q50+=q52;q52=0.0;}if(q52>0.0)if(
q25>=0.){if(q25>q30.q5)q112=2./(1.-q30.q5*q30.q5);else{q112=(1.-q30.q5)/(1.-
q25);q112*=2./(1.-q30.q5*q30.q5)*q112;}}if(q29==-1)return q36;if(q30.mu_min==-
1.0&&q30.mu_max==1.0){if(q50>0.0)q111=get_phase_max(phase_max,n_phase_max,q25,
q46);q113=1.0;}else{if(q25<=q30.mu_max&&q25>=q30.mu_min)q111=get_phase_max(
phase_max,n_phase_max,q25,q46)*2.0/(q30.q6[q86]-q30.F_min[q86]);if(q40>cosd(
q30.d_phi))q113=1.0-q51*(1.0-180.0/q30.d_phi);else q113=1.0-q51;}q37=(1.0-q50-
q52)*q36+q113*(q50*q111+q52*q112);if(!(q37>0.0)){fprintf(stderr,"\105\162\162\157\162\40\151\156\40\120\137\163\160\145\143\72\40\45\145\40\72\40\145\160\163\137\144\144\151\163\137\165\160\146\40\45\145\40\145\160\163\137\144\144\151\163\137\165\144\141\40\45\145\40\120\137\156\157\162\155\40\45\145\40\146\141\143\137\160\150\151\40\45\145\40\120\137\144\144\151\163\137\165\160\146\40\45\145\40\120\137\144\144\151\163\137\165\144\141\40\45\145\40\155\165\62\40\45\145\40\155\165\155\141\170\40"
"\45\145\40\155\165\155\151\156\40\45\145\n",q37,q50,q52,q36,q113,q111,q112,
q25,q30.mu_max,q30.mu_min);fprintf(stderr,"\156\137\160\150\141\163\145\137\155\141\170\40\45\144\40\106\137\155\141\170\40\45\145\40\106\137\155\151\156\40\45\145\40\155\165\137\165\144\141\137\144\142\40\45\145\n"
,n_phase_max,q30.q6[0],q30.F_min[0],q30.q5);}
#ifdef MUCHOUT
if(q82==110)fprintf(stderr,"\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\55\40\143\141\154\143\137\160\137\163\160\145\143\72\40\45\145\40\45\145\40\45\144\40\n"
,q30.mu_min,q30.mu_max,q25<q30.mu_max);
#endif
return q37;}int calloc_hybrid3D_field(int****q41,atmosphere_struct*q9){int kc=
0,ic=0;*q41=calloc((size_t)q9->Nz,sizeof(int**));if(*q41==NULL){fprintf(stderr
,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}for(kc=0;kc<q9->Nz;kc++){if(q9->threed[MCCAOTH_TOT][kc]>=1){(*q41)
[kc]=calloc((size_t)q9->Nx,sizeof(int*));if((*q41)[kc]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}for(ic=0;ic<q9->Nx;ic++){(*q41)[kc][ic]=calloc((size_t)q9->Ny,
sizeof(int));if((*q41)[kc][ic]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}}}else{(*q41)[kc]=calloc((size_t)1,sizeof(int*));if((*q41)[kc]==
NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}ic=0;(*q41)[kc][ic]=calloc((size_t)1,sizeof(int));if((*q41)[kc][ic
]==NULL){fprintf(stderr,"\105\162\162\157\162\40\141\154\154\157\143\141\164\151\156\147\40\155\145\155\157\162\171\40\146\157\162\40\150\171\142\162\151\144\n"
);return-1;}}}return 0;}void free_hybrid3D_field(int****q41,atmosphere_struct*
q9){int kc=0,ic=0;for(kc=0;kc<q9->Nz;kc++){if(q9->threed[MCCAOTH_TOT][kc]>=1){
for(ic=0;ic<q9->Nx;ic++)free((*q41)[kc][ic]);free((*q41)[kc]);}else{free((*q41
)[kc][0]);free((*q41)[kc]);}}free(*q41);}
