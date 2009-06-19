#ifndef MY_EXTRACTOR_
#define MY_EXTRACTOR_
#include <string>
#include <map>
#include <vector>

typedef struct tagCOM{
double  pos[3];
double  vel[3];
} typeCom;

class CExtractor
	{
	public:
		CExtractor(std::string path, std::map<std::string, std::string> ini);
		~CExtractor(void);
		void GetIdToTrace(void);
		void ParseIniFile(void);
	private:
		std::map<std::string, std::string> m_inifile; 
		std::string m_path;
		std::vector<std::map<int, int> > m_ID;
		std::vector<typeCom> m_com;
		double m_Rad;
		std::string idfile;
		std::string m_file;
		std::string m_dumpfile;
	};

#endif

