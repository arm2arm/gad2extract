#ifndef _CRange_
#define _CRange_
#include <iostream>
#include <algorithm>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include "cycle.h"

using std::string;
using std::min;
using std::max;
using std::cout;
using std::endl;

class CRange{
public:
	float Min;
	float Max;
	float m_sum;
	float m_mean;
	unsigned int m_np;
	CRange(){Reset();};
	~CRange(){};
	void sum(float x){m_sum+=x;m_np++;};

	void getbound(float x)
	{
		Min=min(x, Min);
		Max=max(x, Max);
		sum(x);
	};
void Reset(void ){Min=1e10; Max=-1e10;m_sum=0;m_np=0;};
void print(const char *str){
	cout<<str<<" Min="<<Min<<"\tMax="<<Max<<"\t mean="<<m_sum/m_np<<endl;

};
};

#define printOpenGLError() printOglError(__FILE__, __LINE__)

int printOglError(char *file, int line);
float   stoptimer_( int *flag);
void  starttimer_( void);
float mydrand48(void);
void process_mem_usage(double& vm_usage, double& resident_set);
void my_mem_usage();

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795f
#endif
#define EXIT_FAIL exit(1);

template <class T> void safeFree(T v)
{
  if(v!=NULL)
    {
      delete [] v;
      v=NULL;
    }
}

template <class T> void safeFreeOne(T v)
{
  if(v!=NULL)
    {
      delete  v;
      v=NULL;
    }
}

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip> 
#include <map>
using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using std::map;
template < class T >
string ToString(const T &arg)
{
	ostringstream	out;

	out <<std::setprecision(1)<<std::fixed<< arg;

	return(out.str());
}
//////////////////////////////
//keeper XYZ used for GTS analysis
template < class T > 
struct tagCoord
{
    T x;
    T y;
    T z;
};

//////////////////////////////
//////////////////////
float mydrand48(void);
bool FileExists(string filename);
bool StringToInt(const string &s, int &i);
bool GetSnapPath( string snap, string &path);
bool GetSnapName( string &snap, int &isnap);
bool GetParam(std::map<std::string, std::string> m_inifile, std::string bl, void *val, unsigned short type=0);

//////////////////////
#endif

