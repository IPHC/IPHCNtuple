/* BASH COLORS */
#define RST   "[0m"
#define KRED  "[31m"
#define KGRN  "[32m"
#define KYEL  "[33m"
#define KBLU  "[34m"
#define KMAG  "[35m"
#define KCYN  "[36m"
#define KWHT  "[37m"
#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST
#define BOLD(x) "[1m" x RST
#define UNDL(x) "[4m" x RST

#include "TString.h"

#include <iostream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <sstream>
#include <fstream>

#include <cassert> 	//Can be used to terminate program if argument is not true.
//Ex : assert(test > 0 && "Error message");
#include <sys/stat.h> // to be able to use mkdir

using namespace std;

//Convert a double into a TString
// precision --> can choose if TString how many digits the TString should display
TString Convert_Number_To_TString(double number, int precision=10)
{
	stringstream ss;
	ss << std::setprecision(precision) << number;
	TString ts = ss.str();
	return ts;
}

//Convert a TString into a double
double Convert_TString_To_Number(TString ts)
{
	double number = 0;
	string s = ts.Data();
	stringstream ss(s);
	ss >> number;
	return number;
}



/**
 * Executes bash command and returns the output as a Tstring
 * !! PROBLEMS WHEN PIPES INVOLVED (errno=ENOMEM, memory-related issue) --
 */
TString GetStdoutFromCommand(TString cmd_ts)
{
	string output = "";
	FILE* stream=0;
	const int max_buffer = 500; //Lack of buffer memory causes errors (errno=12)
	char buffer[max_buffer]; //Create buffer
	cmd_ts+= " 2>&1"; //Get both stdout and stderr outputs
	string cmd = cmd_ts.Data();

	stream = popen(cmd.c_str(), "r"); //Open read-only stream, run command

	if(stream)
	{
		while(!feof(stream))
		{
			if(fgets(buffer, max_buffer, stream) != NULL) output.append(buffer); //Get the output
		}

		pclose(stream);
	}
	else //If stream was not opened <-> insufficient buffer memory. Retry with more !
	{
		pclose(stream);

		return "";
	}

	output.erase(std::remove(output.begin(), output.end(), '\n'), output.end()); //Remove the "end of line" character

	TString ts_output(output); //Convert string to TString

	return ts_output;
}

TString  Rename(TString name)
{
	if(name.Contains("THQ") ) {return "tHq.root";}
	else if(name.Contains("THW") ) {return "tHW.root";}
	else if(name.Contains("TTTT") ) {return "TTTT.root";}
	else if(name.Contains("WZZ") ) {return "WZZ.root";}
	else if(name.Contains("WWZ") ) {return "WWZ.root";}
	else if(name.Contains("TTW") ) {return "ttW.root";}
	else if(name.Contains("TTZ") ) {return "ttZ.root";}
	else if(name.Contains("WZTo3") ) {return "WZ.root";}
	else if(name.Contains("ZZTo") ) {return "ZZ.root";}
	else if(name.Contains("ZZZ") ) {return "ZZZ.root";}
	else if(name.Contains("tZq") ) {return "tZq.root";}
	else if(name.Contains("ttH") ) {return "ttH.root";}
	else if(name.Contains("DY") ) {return "DY.root";}
	else if(name.Contains("tWll") ) {return "tWZ.root";}
	else if(name.Contains("TTG") ) {return "TTG.root";}
}

void Rename_Samples(TString path)
{
	mkdir("samples_renamed", 0777);

	TString command = "ls " + path + " | wc -l";
	int nfiles = Convert_TString_To_Number( GetStdoutFromCommand(command) );

	for(int ifile=1; ifile<=nfiles; ifile++)
	{
		command = "ls " + path + " | sed -n '" + Convert_Number_To_TString(ifile) + "p'"; //Look to each filename, one at a time

		TString name = GetStdoutFromCommand(command);

		command = "cp " + path + "/" + name + " " + "samples_renamed/" + Rename(name);
		cout<<"--- Will : "<<command<<endl;

		system(command.Data() );
	}

}

int main(int argc, char **argv)
{
    if(argc == 1 || argc > 2) {cout<<FRED("Error : specify directory")<<endl; return 1;}

	TString path = argv[1];

  	Rename_Samples(path);
}
