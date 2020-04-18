#ifndef CNETWORK_H
#define CNETWORK_H

//Declaration of the network object "Cnetwork"

#include "Cport.h"


class Cnetwork
{
 public:
	 Cnetwork();		//constructor
	 ~Cnetwork();		//destructor

	 int AddPort(const PortDataType& MyData, 
		 const string& InputCoaxMovieFileName = "", const string& InputCoaxVoltageFileName = "");
	 void UpdateV(const int& n);
	 void UpdateI(const int& n);
	 void RecordV();	//record the voltage (on the whole line and/or at one point)

	 int NumberOfPorts()
	 {
		 return NumOfPorts;
	 }
	 bool PortIsInNode(const int& PortIndex)
	 {
		 return Ports(PortIndex)->PortBelongsToNode;
	 }

 private:	 
	 string NetworkOutputDir;	//network output directory

	 int NumOfPorts;		//number of ports in the network
	 Array<Cport*,1> Ports;	//array of pointers to the ports in the network

};

#endif
