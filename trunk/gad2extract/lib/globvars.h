
#define  GAMMA (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define  PI  3.1415926

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS 9.10953e-28
#define  THOMPSON     6.65245e-25
#define  ELECTRONCHARGE  4.8032e-10


#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  HYDROGEN_MASSFRAC  0.76

#define  HUBBLE      1.0     /* Hubble constant in 100km/sec/Mpc */  


/***********************************/
extern char COM_FILE[500],OutputDir[500], OutputFile[500];
extern int SNAP_START, SNAP_END,SNAP_STEP;
extern char SnapDir[500], SnapMask[500], AFOF_MASK[512];
extern int TYPE;
extern double RMIN, RMAX;
///////////// PROFILERS PARAMS
extern double profRMAX;
extern double xStep;
extern int fm;
extern int NCOEF;
extern double dTime;
/***********************************/
#include <string>
void read_parameterfile(const char *fname);
void GenName(std::string &fname,const unsigned short isnap, std::string ext="_");
void GetFileName(std::string &fname,const unsigned short isnap);
///////////////////////////////////
#define VERBOSE_flag 0
