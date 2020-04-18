//Includes the routine that reads the network definitions

#include "headers.h"

#include "read_network.h"

//definition of Cnffft needed
#include "Cnetwork.h"


void read_network(Cnetwork &Network)
{
	// 	PortDataType PortData1;
	// 	ostringstream stream;
	// 	stream << "CoaxMovie_" << GridIndex << ".movie";
	// 	string CoaxMovieName = stream.str();
	// 	stream.str("");
	// 	stream << "CoaxVoltage_" << GridIndex << ".voltage" ;
	// 	string CoaxVoltageFileName = stream.str();
	// 	PortData1.Gap_back = MyPortData.Gap_back;
	// 	PortData1.Gap_front = MyPortData.Gap_front;
	// 	PortData1.Gap_left = MyPortData.Gap_left;	
	// 	PortData1.Gap_right = MyPortData.Gap_right;	
	// 	PortData1.Gap_lower = MyPortData.Gap_lower;
	// 	PortData1.Gap_upper = MyPortData.Gap_upper;
	// 	PortData1.vp = MyPortData.vp;
	// 	PortData1.Zc = MyPortData.Zc;	
	// 	//PortData1.WaveShape = "modulatedsinegaussian";
	// 	//PortData1.tau = 64.4e-12;
	// 	PortData1.WaveShape = MyPortData.WaveShape;
	// 	PortData1.tau = MyPortData.tau;
	// 	PortData1.delay = MyPortData.delay;
	// 	//PortData1.w_0 = 2*M_PI*6.85e9;
	// 	PortData1.V0 = MyPortData.V0;
	// 	
	// 	//if (GridIndex==NumberOfGrids-1)
	// 	//{
	// 		//PortData1.mode = "transmit";	
	// 	//}
	// 	//else
	// 	//{
	// 	//	PortData1.mode = "receive";
	// 	//}
	// 	//SlowestGrid = NumberOfGrids-1;	//the NFFFT node is probably the slowest
	// 	PortData1.mode = MyPortData.mode;
	
	// 	int Port1 = Network.AddPort(PortData1,CoaxMovieName,CoaxVoltageFileName);
		//if (GridIndex==0)
		//{
		//	if (Network.PortIsInNode(Port1))
		//	{
		//		cout << "At the port node, iback is " << iback << " and ifront is " << ifront << endl 
		//			<< "The gap is between " << PortData1.Gap_back << " and " << PortData1.Gap_front << endl;	
		//		cout << "At the port node, jleft is " << jleft << " and jright is " << jright << endl 
		//			<< "The gap is between " << PortData1.Gap_left << " and " << PortData1.Gap_right << endl;
		//		cout << "At the port node, klower is " << klower << " and kupper is " << kupper << endl 
		//			<< "The gap is between " << PortData1.Gap_lower << " and " << PortData1.Gap_upper << endl;
		//	}
		//}
}
