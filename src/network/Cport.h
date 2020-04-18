#ifndef CPORT_H
#define CPORT_H

//Declaration of the port object "Cport"

#include <fstream>

extern int OriginX,OriginY,OriginZ;


struct PortDataType
{
	PortDataType():
		mode("transmit"),
		vp(c),
		Zc(50.0),
		Gap_lower(OriginZ+1),	//lowermost z cell in the gap
		Gap_upper(OriginZ),		//uppermost z cell in the gap
		Gap_back(OriginX-1),	//rearmost x cell in the gap
		Gap_front(OriginX),		//foremost x cell in the gap
		Gap_left(OriginY),		//leftmost y cell in the gap
		Gap_right(OriginY),		//rightmost y cell in the gap
		V0(1.0),
		WaveShape("gaussian"),
		tau(40e-12),
		w_0(2*M_PI*6.32e9),
		delay(6)
	{};

	string mode;	//transmitting or receiving mode
	double vp;	//velocity of propagation in the coaxial line
	double Zc;	//characteristic impedance of the coaxial line
	int Gap_lower,Gap_upper,Gap_back,Gap_front,Gap_left,Gap_right;
	double V0;	//maximum voltage
	string WaveShape;	//shape of the waveform
	double tau;	//time constant for the Gaussian waveform
	double w_0;	//angular modulation frequency (if modulated)
	double delay;	//delay in the waveform (in tau)
};

class Cport
{
 private:
	 int PortIndex;	//index of the port

	 friend class Cnetwork;

	 Cport(const PortDataType& MyData, const string& CoaxMovieFileName, const string& CoaxVoltageFileName, 
		 const int& Index);		//constructor
	 void UpdateV(const int& n);
	 void UpdateI(const int& n);
	 void RecordV();	//record the voltage (on the whole line and/or at one point)

	 //basic waveform shape
	 double waveform(double t, double tau);

	 const PortDataType Data;

	 bool PortBelongsToNode;		//is port area included in the node?

	 int length;	//length of the coaxial line
	 Array<double,1> voltage,current;	//voltage and current in the coaxial line
	 int feed_point;	//location of the one-way injector on the coaxial line
	 int obs_point;		//location of the voltage observation point on the coaxial line
	 double observed_voltage; //observed voltage value at that point
	 Array<double,1> waveformV,waveformI;	//incident voltage and current waveforms on the coaxial line

	 double S_coax;	//courant number in the coaxial line

	 int Loop_back,Loop_front,Loop_left,Loop_right,Loop_lower,Loop_upper;

	 bool clockwiseloop;	//is the loop clockwise? (y axis pointing away from the clock)

	 double PreviousV;	//storage variable used in the 1st order ABC

	 int GapWidth;	//width of the gap (in dx)
	 double GapVoltage;			//voltage applied at the gap
	 double TerminalCurrent;	//terminal current

	 ofstream CoaxMovieFile,CoaxVoltageFile;

};

#endif
