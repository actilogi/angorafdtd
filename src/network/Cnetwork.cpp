//Definition of the network object "Cnetwork"
//An network object may comprise multiple port objects, represented by Cport.

#include "headers.h"

#include "Cnetwork.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

extern bool check_mode;

extern int create_path(const string& path);


Cnetwork::Cnetwork()
	: NumOfPorts(0), Ports(0)		//initialize with no ports
{
	NetworkOutputDir = "./network";
	if (!check_mode)
	{
		create_path(NetworkOutputDir);
	}
}

int Cnetwork::AddPort(const PortDataType& MyData,
					   const string& InputCoaxMovieFileName, const string& InputCoaxVoltageFileName)
{
	//if not specified, name the output files according to the index in the port collection
	string CoaxMovieFileName, CoaxVoltageFileName;
	if (InputCoaxMovieFileName == "")
	{
		ostringstream stream;
		//build and assign filename
		stream << NetworkOutputDir << "/CoaxMovie_" << NumOfPorts+1 << ".movie" ;
		CoaxMovieFileName = stream.str();
	}
	else
	{
		ostringstream stream;
		//just append the output directory
		stream << NetworkOutputDir << "/" << InputCoaxMovieFileName;
		CoaxMovieFileName = stream.str();
	}
	if (InputCoaxVoltageFileName == "")
	{
		ostringstream stream;
		//build and assign filename
		stream << NetworkOutputDir << "/CoaxVoltage_" << NumOfPorts+1 << ".voltage" ;
		CoaxVoltageFileName = stream.str();
	}
	else
	{
		ostringstream stream;
		//just append the output directory
		stream << NetworkOutputDir << "/" << InputCoaxVoltageFileName;
		CoaxVoltageFileName = stream.str();
	}
	NumOfPorts++;
	Ports.resizeAndPreserve(NumOfPorts);
	Ports(NumOfPorts-1) = new Cport(MyData,CoaxMovieFileName,CoaxVoltageFileName,NumOfPorts-1);

	return NumOfPorts-1;	//indices of ports begin with 0
}

void Cnetwork::UpdateV(const int& n)
{
	for (int i=0; i<NumOfPorts; i++)
	{
		Ports(i)->UpdateV(n);		//update the port voltage
	}
}

void Cnetwork::UpdateI(const int& n)
{
	for (int i=0; i<NumOfPorts; i++)
	{
		Ports(i)->UpdateI(n);		//update the port current
	}
}

void Cnetwork::RecordV()
{
	for (int i=0; i<NumOfPorts; i++)
	{
		Ports(i)->RecordV();		//record the port variables (movie and voltage at a point)
	}
}

Cnetwork::~Cnetwork()
{//releases dynamic memory allocated for each port
	for (int i=0; i<NumOfPorts; i++)
	{
		delete Ports(i);
	}
}
