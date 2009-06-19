#include "Extractor.h"
#include "globvars.h"
#include "utils.h"
#include <fstream>
#include <iostream>
CExtractor::CExtractor(std::string path, std::map<std::string, std::string> ini):
m_inifile(ini), m_path(path)
	{
	ParseIniFile();
	GetIdToTrace();
	}

CExtractor::~CExtractor(void)
	{
	}

void CExtractor::GetIdToTrace()
	{
	std::ifstream inidfile(idfile.c_str());
	if (!inidfile.bad())
		{
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