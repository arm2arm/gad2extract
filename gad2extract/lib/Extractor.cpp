#include "Extractor.h"
#include "globvars.h"
#include "utils.h"
#include "data_readers.h"
#include <stdio.h>

#include <fstream>
#include <iostream>

typedef int UCHAR;

//////////////////////////////
class CGetData{


public:
  CGetData(string file, int type=4):m_infile(file),
				    m_type(type),
				    pID(NULL),   
				    pPOS(NULL), pTYPE(NULL), 
				    pMASS(NULL), pVEL(NULL)
  {        
    if(! ReadData()){
      cout<<"Fail to read data...exiting"<<endl;EXIT_FAIL;}
  };
  ~CGetData(){
    safeFree(pID);
    safeFree(pPOS);
    safeFree(pTYPE);
    safeFree(pVEL);
    safeFree(pMASS);
    safeFreeOne( m_gin);
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
    safeFree(pVEL);
    safeFree(pID);
    safeFree(pMASS);
    string blname = "POS ";
    nelem=m_gin->read_whole_block((char*&)pPOS, 
			"POS ")/3;
    nelem=m_gin->read_whole_block((char*&)pVEL, 
			"VEL ")/3;
    m_gin->read_whole_block((char*&)pID, 
			"ID  ");
    nwithmass=m_gin->read_whole_block((char*&)pMASS, 
			   "MASS");

    
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
    int nb =31;
    for(int i=0;i<nb;i++)
      if(m_gin->find_block(&m_gin->m_file,allblocks[i])>1 )
	{
	  cout<<"#PRESENT =>  "<<allblocks[i]<<endl;
	  which_blocks.push_back(allblocks[i]);
	}
    
  }

 CGadget *m_gin;
protected:
  int m_type;
  unsigned int nelem;
  string m_infile;
  int m_isnap;
  std::map<int, int> ID;
  int ntotal;
public:
  int *pID;
  UCHAR *pTYPE;
  float *pPOS;
  float *pVEL;
  float *pMASS;
  int nwithmass;
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
	floc.resize(m_nhalo);
	for(int hcur=0;hcur<m_nhalo;hcur++)
	  {
	    int im=0;// mass block counter
	    int ig=0;// ig block counter	
	    int iz=0;// iZ block counter
	    int iage=0;// iage block counter

	    for(int i=0;i<data->GetNp();i++)
	      {	
		myType=data->pTYPE[i];
		//first setup counters
		if(data->m_gin->myhead.mass[myType]==0 && 
		       data->m_gin->myhead.npart[myType]>0)
		  {
		    im++;
		  }
		if(myType==0)ig++;
		if(myType==0||myType==4)iz++; 
		if(myType==4)iage++;
		//check region
		if( isInReg(data->pPOS[i*3], data->pPOS[i*3+1], data->pPOS[i*3+2],hcur) )
		  {		    
		    id_loc[hcur].type_idx[myType].push_back(i);
		    floc[hcur].mergedIDX.push_back(i);
		    
		   if(data->m_gin->myhead.mass[myType]==0 && 
		       data->m_gin->myhead.npart[myType]>0)
		     {
		        floc[hcur].IDXWithMass.push_back(im);
		     }
		   if(myType==0)
		      floc[hcur].IDXgas.push_back(ig);
		   if(myType==0 || myType==4)
		      floc[hcur].IDXz.push_back(iz); 
		   if(myType==4)
		     floc[hcur].IDXage.push_back(iage);  
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
	char dumpname[1024];
	for(int hcur=0;hcur<m_nhalo;hcur++)
	  {	
	    sprintf(&dumpname[0], m_dumpfile.c_str(),hcur,int(m_Rad), snap.c_str());
	    cout<<"# "<<dumpname<<endl;
	    
	    cout<<"# got total in region:"<<floc[hcur].mergedIDX.size()<<endl;
	    cout<<"# got Nwithmass:"<<floc[hcur].IDXWithMass.size()<<endl;
	    cout<<"# got NwithZ:"<<floc[hcur].IDXz.size()<<endl;
	
	    ofstream file;	    
	    file.open(dumpname,ios::binary);
	    CGadget::io_header head;
	    int totpart=0;
	   
	    CGadget *tg=new CGadget();
	    memcpy(&head,&data->m_gin->myhead,sizeof(head));
	    for(int it=0;it<6;it++)
	      {
		head.npart[it]=id_loc[hcur].type_idx[it].size();
		head.npartTotal[it]=head.npart[it];
		totpart+=head.npart[it];
		if(head.npart[it]==0)
		  head.mass[it]=0.0;
	      }
	    head.num_files=1;
	    /////////////////////////////////////////////////////
	    int blsize;
	    blsize=sizeof(head);
	    tg->WriteOneBlock(file,"HEAD", (char*)&head, blsize);
	    
	    //////////////////////////////////////////////////////
	  
	      {
		{
		    int idx=0, k=0;
		    float *pPos=new float[totpart*3];
		    float *pVel=new float[totpart*3];
		    int *pID=new int[totpart];
		   
		
		    for(int ip=0;ip<totpart;ip++)
		      {
			idx=floc[hcur].mergedIDX[ip];
			pPos[k]=data->pPOS[idx*3];
			pPos[k+1]=data->pPOS[idx*3+1];
			pPos[k+2]=data->pPOS[idx*3+2];

			pVel[k]=data->pVEL[idx*3];
			pVel[k+1]=data->pVEL[idx*3+1];
			pVel[k+2]=data->pVEL[idx*3+2];

			pID[ip]=data->pID[idx];
			
			k+=3;
		      }
		    ///////////////////////////////
		    cout<<"try to free some data"<<endl;
		    cout<<"done"<<endl;
		    //////////////////////
		    blsize=sizeof(float)*totpart*3;		   
		    tg->WriteOneBlock(file,string("POS "), (char*)pPos, blsize);	          
		    tg->WriteOneBlock(file,string("VEL "), (char*)pVel, blsize);	 
		    blsize=sizeof(int)*totpart;
		    tg->WriteOneBlock(file,string("ID  "), (char*)pID, blsize);
		    safeFree(pPos);		
		    safeFree(pVel);
		    safeFree(pID);
		    for(int ib=1;ib<which_blocks.size();ib++)
		    {
		    ///////////////////////////////
		    if(which_blocks[ib]=="MASS" && floc[hcur].IDXWithMass.size()>0)
		      {
			float *pMASS=
			  pMASS=new float[floc[hcur].IDXWithMass.size()];
			for(int ip=0;ip<floc[hcur].IDXWithMass.size();ip++)
			  {
			    pMASS[ip]=data->pMASS[floc[hcur].IDXWithMass[ip]];
			  }
			  blsize=sizeof(int)*floc[hcur].IDXWithMass.size();
			  tg->WriteOneBlock(file,string("MASS"), (char*)pMASS, blsize);
			  safeFree( pMASS);		    
		      }
		    if(which_blocks[ib]=="U   "||
		       which_blocks[ib]=="RHO "||
		       which_blocks[ib]=="NE  "||
		       which_blocks[ib]=="NH  "||
		       which_blocks[ib]=="HSML"||
		       which_blocks[ib]=="SFR "
		       )
		      {
			float *pBUF=NULL;//new float[data->m_gin->myhead.npart[0]];
			data->m_gin->read_whole_block((char*&)pBUF, 
			   which_blocks[ib].c_str());			
			for(int ip=0;ip<floc[hcur].IDXgas.size();ip++)
			  {
			    pBUF[ip]=pBUF[floc[hcur].IDXgas[ip]];
			  }
			blsize=sizeof(float)*floc[hcur].IDXgas.size();
			tg->WriteOneBlock(file,which_blocks[ib], (char*)pBUF, blsize);
			safeFree(pBUF);
		      }
		    if(which_blocks[ib]=="AGE ")
		      {
			float *pBUF=NULL;//=new float[data->m_gin->myhead.npart[4]];
			data->m_gin->read_whole_block((char*&)pBUF, 
			   which_blocks[ib].c_str());			
			for(int ip=0;ip<floc[hcur].IDXage.size();ip++)
			  {
			    pBUF[ip]=pBUF[floc[hcur].IDXage[ip]];
			  }
			blsize=sizeof(float)*floc[hcur].IDXage.size();
			tg->WriteOneBlock(file,which_blocks[ib], (char*)pBUF, blsize);
			safeFree(pBUF);
		      }
		    if(which_blocks[ib]=="Z   ")
		      {
			float *pBUF=NULL;//new float[data->m_gin->myhead.npart[0]+data->m_gin->myhead.npart[4]];
			data->m_gin->read_whole_block((char*&)pBUF, 
			   which_blocks[ib].c_str());			
			for(int ip=0;ip<floc[hcur].IDXz.size();ip++)
			  {
			    pBUF[ip]=pBUF[floc[hcur].IDXz[ip]];
			  }
			blsize=sizeof(float)*floc[hcur].IDXz.size();
			tg->WriteOneBlock(file,which_blocks[ib], (char*)pBUF, blsize);
			safeFree(pBUF);
		      }
		    }
		  
		    file.close();		    
		   		    
		  }
	      }	    
	    //////////////////////////////////////////////////////
	    

	    safeFreeOne( tg);
	    
	  }
	
	/////////////////////////////////////
	safeFreeOne( data);
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
