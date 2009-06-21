#include "Extractor.h"
#include "globvars.h"
#include "utils.h"
#include "data_readers.h"
#include <stdio.h>

#include <fstream>
#include <iostream>

typedef unsigned char UCHAR;

//////////////////////////////
class CGetData{


public:
  CGetData(string file, int type=4):m_infile(file),
				    m_type(type),
				    pID(NULL),   
				    pPOS(NULL), pTYPE(NULL)
  {        
    if(! ReadData()){
      cout<<"Fail to read data...exiting"<<endl;EXIT_FAIL;}
  };
  ~CGetData(){
    safeFree(pID);
    safeFree(pPOS);
    safeFree(pTYPE);
    delete m_gin;
  };
  void GetCom(std::vector<int> &trID, double *COM)
  {
    int i, ip;
    std::map<int, int>::iterator it;	  
    int nprod=trID.size();
    COM[0]=0.0 ; COM[1]=0.0 ; COM[2]=0.0;
    cout<<"# Get Com"<<endl;
    for(i=0;i<nprod;i++)
      {
	
	it=ID.find(trID[i]);
	if(it==ID.end()){std::cerr<<ID.size()<<" ERROR: "<<
			   trID[i]<<" not in set"<<
			   std::endl;exit(0);};
	ip=(*it).second*3;
	COM[0]+=(double)pPOS[ip];
	COM[1]+=(double)pPOS[ip+1];
	COM[2]+=(double)pPOS[ip+2];
      }	
   
    COM[0]/=double(nprod);
    COM[1]/=double(nprod);
    COM[2]/=double(nprod);
  };
  
  unsigned int GetNp(){ return nelem; };
  inline double GetBoxSize(){return m_gin->myhead.BoxSize;}
 
protected:  
  bool ReadData()
  {
    m_gin=new CGadget(m_infile, false);
    bool ans;
   
    ans = GetIDBlock(m_gin);
    ans = ans && GetPosBlock(m_gin);
    ans = ans && GetTypeBlock(m_gin);
    cout<<"#Fill Map"<<endl;
    for(int i=0;i<nelem;i++)
      {
	ID.insert(std::make_pair(pID[i],i));

      }
    
    
    return ans;
  };  

  bool GetIDBlock(CGadget *g)
  {
    bool retval=true;
    
    nelem=0;
    safeFree(pID);
    nelem=g->read_block(pID, 
			"ID  ",
			m_type);
    
    if(nelem==0)
      {  cout<<"Error reading block: "<<"ID  "<<
	   " In file: "<<g->m_filename<<endl;
      retval=false;
      }
#ifdef VERBOSE
    else      cout<<"Getting: "<<"ID  "<<" with "<<nelem<<" elems.."<<endl;
#endif
    return retval;
    
  };
 
  bool GetTypeBlock(CGadget *g)
  {
    ntotal=0;
    for(int itype=0;itype<6;itype++)
      {
	ntotal+=g->myhead.npart[itype];
      }
    pTYPE=new UCHAR[ntotal];
    int ip=0;
    for(int itype=0;itype<6;itype++)
      {
	for(int i=0;i<g->myhead.npart[itype];i++)
	  {	  
	    pTYPE[ip]=itype;
	    ip++;
	  }
      }
    return true;
  }
  bool GetPosBlock(CGadget *g)
  {
    bool retval=true;
    
    nelem=0;
    safeFree(pPOS);
    string blname = "POS ";
    nelem=g->read_blockv3(pPOS, 
			"POS ",
			m_type);
    
    if(nelem==0)
      {  cout<<"Error reading block:"<<"POS"<<
	   " In file: "<<g->m_filename<<endl;
      retval=false;
      }
#ifdef VERBOSE
    else cout<<"Getting: "<<"POS"<<" with "<<nelem<<" elems.."<<endl;
#endif
    return retval;
    
  };
public:
  bool GetAllPos()
  {
    bool retval=true;
    
    nelem=0;
    safeFree(pPOS);
    string blname = "POS ";
    nelem=m_gin->read_whole_block((char*&)pPOS, 
			"POS ")/3;
    
    if(nelem==0)
      {  cout<<"Error reading block:"<<"POS"<<
	   " In file: "<<m_gin->m_filename<<endl;
      retval=false;
      }
#ifdef VERBOSE
    else cout<<"Getting: "<<"POS"<<" with "<<nelem<<" elems.."<<endl;
#endif
    return retval;
    
  };
  void blocks(std::vector<std::string> &which_blocks)
  {
    const char* allblocks[]={"POS ", "VEL ","ID  ", 
			     "MASS","U   ","RHO ",
			     "NE  ", "NH  ", "HSML",
			     "SFR ", "AGE ", "Z   ", 
			     "POT ", "ACCE", "TSTP", 
			     "COOR", "ENDT", "BFLD", 
			     "DBDT", "DIVB", "ABVC", 
			     "CONR", "BFSM", "EGYP",
			     "EGYC", "CRC0", "CRQ0", 
			     "BHMA","BHMD", "MACH", 
			     "DISSE"
    };
    int nb =32;
    for(int i=0;i<nb;i++)
      if(m_gin->find_block(&m_gin->m_file,allblocks[i])>1 )
	{
	  cout<<allblocks[i]<<endl;
	  which_blocks.push_back(allblocks[i]);
	}
    
  }
protected:
  int m_type;
  unsigned int nelem;
  string m_infile;
  CGadget *m_gin;
  int m_isnap;
  std::map<int, int> ID;
  int ntotal;
public:
  int *pID;
  UCHAR *pTYPE;
  float *pPOS;

};

//////////////////////////////

CExtractor::CExtractor(std::string path, std::map<std::string, std::string> ini):
m_inifile(ini), m_path(path)
	{
	ParseIniFile();
	GetIdToTrace();
	/////////////////////////////////////
	UCHAR myType;
	CGetData *data= new CGetData(m_file);
	
	m_Box=data->GetBoxSize();
	m_Box2=data->GetBoxSize()*0.5;
	cout<<"#BoxSize/2="<<m_Box2<<endl;
	m_com.resize(m_nhalo);
	id_loc.resize(m_nhalo);
	for(int hcur=0;hcur<m_nhalo;hcur++)
	  {
	    cout<<"# Anal:"<<hcur<<endl;
	    data->GetCom(m_ID[hcur],(double *) &(m_com[hcur]) );
	    cout<<"# COM is :"<<m_com[hcur].pos[0]<<
	      " "<<m_com[hcur].pos[1]<<
	      " "<<m_com[hcur].pos[2]<<endl;
	   
	  }
	cout<<"#get Whole data"<<endl;
	data->GetAllPos();
	cout<<"# Fill ID location info with Np="<<data->GetNp()<<endl;
	for(int hcur=0;hcur<m_nhalo;hcur++)
	  {
	    for(int i=0;i<data->GetNp();i++)
	      {
		if( isInReg(&data->pPOS[i*3],hcur,data->pTYPE[i] ) )
		  {
		    myType=data->pTYPE[i];
		    id_loc[hcur].type_idx[myType].push_back(i);
		    
		  }

	      }	   
	  }
	      	    
	ReportHaloStats();
	//////////////////////////////////////////////////////
	//which blocks are present??
	data->blocks(which_blocks);	
	
	//////////////////////////////////////////////////////
	cout<<"# start to write dumps "<< m_dumpfile<<
	  " found "<<which_blocks.size()<<" blocks."
	    <<endl;
		std::string snap=m_file.substr(m_file.find_last_of("_")+1,  
			 m_file.find_last_of("_")+3);
	char dumpname[256];
	for(int hcur=0;hcur<m_nhalo;hcur++)
	  {	
	    sprintf(&dumpname[0], m_dumpfile.c_str(),hcur,int(m_Rad), snap.c_str());
	    cout<<"# "<<dumpname<<endl;
	  }
	
	/////////////////////////////////////
	delete data;
	/////////////////////////////////////
	cout<<"# DONE CExtractor"<<endl;
	}

CExtractor::~CExtractor(void)
	{
	}

void CExtractor::GetIdToTrace()
	{
	std::ifstream inidfile(idfile.c_str());
	if (!inidfile.bad())
		{
		  ////////////////////////
		  string oneline;
		  getline(inidfile,oneline);
		  cout<<oneline<<endl;
		  inidfile>>m_nhalo;
		  cout<<"# number of haloes to be tracked = "<<m_nhalo<<endl;
		  ////////////////////////
		  m_ID.resize(m_nhalo);
		  ////////////////////////
		  for(int currHalo=0;currHalo<m_nhalo;currHalo++)
		    {
		      int np, tID;
		      inidfile>>np;
		      cout<<"\t#reding H"<<currHalo<<" with np = "<<np<<endl;
		      for(int ip=0;ip<np;ip++)
			{
			  inidfile>>tID;
			  m_ID[currHalo].push_back(tID);
			
			}
		      
		    }
		  /////////////////////////
		inidfile.close();
		}
	else 
		{
		std::cerr<<"Error cannot open ID file for input:\n"<<idfile<<std::endl;exit(0);
		}
	}


void CExtractor::ParseIniFile()
	{
	char strParam[256]="";
	string blkName="EXTRACTREGION.FILE";
	if(GetParam(m_inifile,blkName,(void*)strParam))
		{
		m_file=strParam;
		}
	else
		{
		std::cerr<<"Error cannot find "<< blkName <<" block in the ini file:\n"<<std::endl;
		exit(0);
		}
	blkName="EXTRACTREGION.DUMPFILE";
	if(GetParam(m_inifile,blkName,(void*)strParam))
		{
		m_dumpfile=strParam;
		}
	else
		{
		std::cerr<<"Error cannot find "<< blkName <<" block in the ini file:\n"<<std::endl;
		exit(0);
		}
	blkName="EXTRACTREGION.DUMPFILE";
	if(GetParam(m_inifile,blkName,(void*)strParam))
		{
		m_dumpfile=strParam;
		}
	else
		{
		std::cerr<<"Error cannot find "<< blkName <<" block in the ini file:\n"<<std::endl;
		exit(0);
		}
	blkName="EXTRACTREGION.IDFILE";
	if(GetParam(m_inifile,blkName,(void*)strParam))
		{
		idfile=strParam;
		}
	else
		{
		std::cerr<<"Error cannot find "<< blkName <<" block in the ini file:\n"<<std::endl;
		exit(0);
		}
	blkName="EXTRACTREGION.R";
	if(!GetParam(m_inifile,blkName,(void*)&m_Rad,2))
		{
		std::cerr<<"Error cannot find "<< blkName <<" block in the ini file:\n"<<std::endl;
		exit(0);
		}
		

	}
