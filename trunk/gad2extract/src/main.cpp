//#include "UseThreads.h"
#include <iostream>
#include <string>
#include <set>
#include <sstream>
#include <exception>
#include <fstream>

#include <boost/config.hpp>
#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "utils.h"
#include "Extractor.h"

using namespace boost::program_options;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

namespace pod = boost::program_options::detail;

int RunAnalysis(string path,string inifile)
{
	std::ifstream config(inifile.c_str());
	//parameters
	std::set<std::string> options;
	std::map<std::string, std::string> parameters;
	options.insert("*");

	try
	{      
		for (pod::config_file_iterator i(config, options), e ; i != e; ++i)
		{
		  std::cout<<"# " << i->string_key <<" "<<i->value[0] << std::endl;
			parameters[i->string_key] = i->value[0];
		}
		//	std::cout <<"# testing maping(sould be same as previos line): "<< parameters["OTHER.LOGFILE"] << std::endl;
	}
	catch(std::exception& e)    
	{
		std::cerr<<"Exception: "<<e.what()<<std::endl;
	}
	//////////////////////////////////
	CExtractor *extract=new CExtractor(path, parameters);
	delete extract;
	//////////////////////////////////
	return 0;
}

int RunBenchmark(unsigned int initmem)
{
	unsigned short  NUM_CPU=2;
	//CUseThreads * analysis=new CUseThreads(NUM_CPU, initmem);
	
	//delete analysis;
	return 0;
}


int  main(int ac, char* av[])
{     
	try
	{
		std::string inifile="extractor.ini";
		std::string path;
		std::string inifname;
		unsigned int vecsize=128; 
		options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce a help screen")
			("version,v", "print the version number")
			("bench",value<unsigned int>(&vecsize),
			"Run scalability Benchmark\n arg size in the Mb")
			("simulation-path,p", value<std::string>(),
			"simulation root path to be analysed")
			("inifile", value<string>(&inifile)->default_value("extractor.ini"),
			" Setup custom ini file name ")	       
			;

		variables_map vm;
		store(parse_command_line(ac, av, desc), vm);
		////////////////////////////////////////////////
		if (vm.count("help")||ac<2) {
			cout << "Usage: analyzer [options]\n";
			cout << desc;
			return 0;
		}
		if (vm.count("version")) {
			cout << "Version 1.\n";
			return 0;
		}

       
		if (vm.count("simulation-path")) {
			path=vm["simulation-path"].as<string>();
			cout << "# The analyzer will use simulation:  \""
				<< path << "\"\n";
			inifname=path+"/extractor.ini";
			if (vm.count("inifile"))
				inifname=path+vm["inifile"].as<string>();
			if(!FileExists(inifname))
			{
				inifname=path+"extractor.ini";//lets try default ini file
				cout << "# INI file : "<<inifname<<" ..... not found "<<endl;
					if(!FileExists(inifname))
					{
						throw "Cannot find any INI file for the model.";
					}
			}
			cout<<"# Found INI file: \n# "<<inifname<<endl;
		}
		////////////////////////////////////////////////
		RunAnalysis(path, inifname);
		////////////////////////////////////////////////
	}
	catch(std::exception& e)    
	{
		std::cerr<<"Exception: "<<e.what()<<std::endl;
		return 0;
	}
	return 0;
}
