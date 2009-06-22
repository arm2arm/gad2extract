#ifndef MY_EXTRACTOR_
#define MY_EXTRACTOR_
#include <string>
#include <map>
#include <set>
#include <vector>
#include <iostream>    
#include <cmath>
#include "utils.h"
typedef struct tagCOM{
double  pos[3];
  //double  vel[3];
 friend std::ostream& operator<<(std::ostream& os, const tagCOM  e)
    {
     
      // os<<e.pos[0];//<<" "<<e.pos[1]<<" "<<e.pos[2]<<std::endl;
    }
} typeCom;
typedef struct tagidInReg{
  std::vector<int> type_idx[6];
  // mapping Type of particle and its position in the file;
} id_maps;

typedef struct tagIDXinFile
{
  std::vector<int> mergedIDX, IDXWithMass,IDXgas, IDXsfr, IDXage, IDXz;
} TIDXinFile;

using std::cout;
using std::endl;

class CExtractor
	{
	public:
		CExtractor(std::string path, std::map<std::string, std::string> ini);
		~CExtractor(void);
		void GetIdToTrace(void);
		void ParseIniFile(void);
		void ReportHaloStats()
		  {
		    for(int ih=0;ih<m_nhalo;ih++)
		      {
			std::cout<<"#Halo "<<ih<<"\n#---------------\n";
			for(int it=0;it<6;it++)
			  std::cout<<"#\t T"<<it<<" \t Np = "<<
			    id_loc[ih].type_idx[it].size()<<"\n";
			std::cout<<"\n#---------------\n"<<std::endl;			
		      }
		  }
		bool isInReg(float x, float y,float  z, int hcur=0, int type=4)
		  {
		    //		    double x=pos[0];
		    //		    double y=pos[1];
		    //		    double z=pos[2];
		
		    double xl=m_com[hcur].pos[0]-m_Rad;
		    double yl=m_com[hcur].pos[1]-m_Rad;
		    double zl=m_com[hcur].pos[2]-m_Rad;

		    double xr=m_com[hcur].pos[0]+m_Rad;
		    double yr=m_com[hcur].pos[1]+m_Rad;
		    double zr=m_com[hcur].pos[2]+m_Rad;
		    if(  
		       xl<x && x<xr &&
		       yl<y && y<yr &&
		       zl<z && z<zr 
			)
		      {
		
			return true;
		      }
		    return false;

		  }
	private:
		std::map<std::string, std::string> m_inifile; 
		std::string m_path;

		std::vector<typeCom> m_com;
		int m_nhalo;
		double m_Rad;
		double m_Box, m_Box2;
		std::string idfile;
		std::string m_file;
		std::string m_dumpfile;
		std::vector<std::vector<int> > m_ID;//IDs to be traced for COM
		std::vector<id_maps> id_loc; // will hold particle location per halo
		std::vector<std::string> which_blocks;//which blocks are present??
		std::vector<TIDXinFile> floc;//halos Ids location per HALO;
	};

#endif

