
char COM_FILE[500],OutputDir[500], OutputFile[500];
int SNAP_START, SNAP_END,SNAP_STEP;
char SnapDir[500], SnapMask[500], AFOF_MASK[512];

int TYPE;
double RMIN, RMAX;
/////////////////////
double profRMAX;
double xStep;
int fm;
int NCOEF;
double dTime=0.005;

#include <string>
using std::string;
void GenName(string &fname,const unsigned short isnap, string ext="_")
{
  char buf[16000], snapname[512];
  string extmask;
  string str(SnapMask);
  size_t found=str.find("%");
  extmask.assign(str,found, 4);
  sprintf(snapname, SnapMask, isnap);
  sprintf(buf, "%s../OTHER/%s%s", OutputDir, ext.c_str(),snapname);
  fname.assign(buf);
}

void GetFileName(string &fname,const unsigned short isnap)
{
  char buf[16000], snapname[512];
  string extmask;
  string str(SnapMask);
  size_t found=str.find("%");
  extmask.assign(str,found, 4);
  sprintf(snapname, SnapMask, isnap);
  sprintf(buf, "%s/%s", SnapDir, snapname);
  fname.assign(buf);
}

/*	std::sort(idu.begin(),idu.end());// stl vector
	std::sort(snapu.begin(),snapu.end());// stl vector
	idu.erase(unique(idu.begin(), idu.end()),  idu.end());// stl vector
	snapu.erase(unique(snapu.begin(), snapu.end()),  snapu.end());// stl vector
*/	
