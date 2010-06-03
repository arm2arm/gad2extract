//last change Arman Khalatyan akhalatyan@aip.de
//2004-06-18

#include "utils.h"
#include <stdio.h>
#include <string.h>

//int __gxx_personality_v0;

//    void    starttimer_(void);
//    float  stoptimer_(int  *flag);
//    float  gettimer_(void);
//    double drand_(int);


clock_t startm, stopm;
void  starttimer_( void)
    //  Start Timer
{
    
    if ( (startm = clock()) == -1) {printf("Error calling clock\n\n");exit(1);}
}
float   stoptimer_( int *flag)
    //Stop and print
{
    //printf("flag = %d", *flag);
    if ( (stopm = clock()) == -1) {printf("Error calling clock\n\n");exit(1);}
    
    float run_time = ((float)(stopm-startm))/CLOCKS_PER_SEC/60.0f;
    if(*flag!=0)
      printf( "\n%6.4f min used by the processor.\n", run_time);
   
    return run_time;
}

float   gettimer_( void)
    //Stop and print
{
       
    float timer;
    if ( (timer = float(clock())) == -1) 
    {
	printf("Error calling clock\n\n");
	exit(1);
    }

    // printf( "\nclock = %6.4f\n", timer);
   
    return timer/CLOCKS_PER_SEC;
}
double drand_(int seed)
{
 srand(seed);
 return double(rand())/float(RAND_MAX);
}

float mydrand48(void)
{
float x = rand()/float(RAND_MAX);
return x;

}

//////////////// MEM UTILS ///////////////////
#ifndef WIN32
#include <unistd.h>
#endif
#include <ios>
#include <iostream>
#include <fstream>
#include <string>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set)
{
#ifndef WIN32
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / (1024.0*1024.0);
   resident_set = rss * page_size_kb;
#endif
}

void my_mem_usage()
{

  using std::cout;
  using std::endl;
  double vm, rss;

  process_mem_usage(vm, rss);
  std::cout.precision(2);
#ifndef WIN32

  std::cout <<std::fixed << " vMEM USAGE: " << vm << " Mb ;" << std::endl;
#else
    cout << " MEMORY USAGE NOT IMPLEMENTED ON WINDOWS!!!"<< endl;
#endif
  
}

// File Utils ////////
/////////////////////////////
// Checks if a file exists //
/////////////////////////////
bool FileExists(string filename) {

	std::ifstream file;
	file.open(filename.c_str());

	// Check if the file exists
	if (file.is_open() == true) {
		file.close();
		return true;
		}
	file.close();
	return false;

	}
#include <map>
#include <string>
///Just for parameter file parser
bool GetParam(std::map<std::string, std::string> m_inifile, std::string bl, void *val, unsigned short type)
{
	if(m_inifile.count(bl))
	{
	  if(type==0)//STRING
	    {
	      strcpy((char*)val,m_inifile[bl.c_str()].c_str());		
	    }
	 else
	   if(type==1)//INT
	     {
	       int ival=atoi(m_inifile[bl.c_str()].c_str());
	       memcpy(val, &ival, sizeof(int));
	     }else//DOUBLE
	       {
		 double dval=atof(m_inifile[bl.c_str()].c_str());
		 memcpy(val, &dval, sizeof(double));
	       }
	}else
	  {
	    cout<<"# Cannot find block in INI: "<<bl<<endl;
	    return false;
	  }
	return true;	
}

