//Defines the port object "Cport"

#include "headers.h"

#include "Cport.h"

// //Uses Gaussian-type waveforms
// #include "waveforms/waveforms.h"
// 
// extern double courant,dt,dx;
// extern int NSTEPS;
// extern int OriginX,OriginY,OriginZ;
// 
// extern Array<double,3> Ex,Ey,Ez;
// extern Array<double,3> Hx,Hy,Hz;
// 
// extern int GridIndex;
// extern int rank;
// extern int iback,ifront;
// extern int jleft,jright;
// extern int klower,kupper;
// 
// extern double MaxFieldValue;
// 
// extern int i,j,k;
// 

Cport::Cport(const PortDataType& MyData, const string& CoaxMovieFileName, const string& CoaxVoltageFileName,
			 const int& Index): Data(MyData), PortIndex(Index), PortBelongsToNode(false)
{
// 	length = 100;
// 	feed_point = length/2;
// 	obs_point = feed_point;
// 
// 	voltage.resize(Range(1,length));
// 	current.resize(Range(1,length));
// 	voltage = 0;
// 	current = 0;
// 
// 	S_coax = courant/1.73205*(Data.vp/c);	//courant number in the coaxial line = vp*dt/dx = (S/sqrt(3))*(vp/c)
// 
// 	//width of the gap (in dx)
// 	GapWidth = (Data.Gap_right-Data.Gap_left+1);
// 	//position of the loop around the current sheet 
// 	Loop_back = Data.Gap_back-1;			//rearmost x cell in the loop around the current sheet
// 	Loop_front = Data.Gap_front+1;			//foremost x cell in the loop around the current sheet
// 	Loop_left = Data.Gap_left;				//y position of the leftmost loop
// 	Loop_right = Data.Gap_right;			//y position of the rightmost loop
// 	Loop_lower = Data.Gap_lower-1;			//lowermost z cell in the loop around the current sheet
// 	Loop_upper = Data.Gap_upper+1;			//uppermost z cell in the loop around the current sheet
// 
// 	clockwiseloop = false;
// 
// 	if ((Data.Gap_front>=iback)&&(Data.Gap_front<=ifront)
// 		&&(Data.Gap_right>=jleft)&&(Data.Gap_right<=jright)
// 		&&(Data.Gap_upper>=klower)&&(Data.Gap_lower<=kupper))
// 	{
// 		PortBelongsToNode = true;
// 	}
// 
// 	if ((Data.mode!="transmit")&&(Data.mode!="receive")&&(Data.mode!="probe"))
// 	{
// 		if (rank==0)
// 		{
// 			cout << "Invalid port type (" << Data.mode << ") for port " 
// 				<< PortIndex << " in grid " << GridIndex << endl;
// 			exit(-1);
// 		}
// 	}
// 
// 	if ((Data.mode=="transmit")||(Data.mode=="probe"))
// 	{
// 		waveformV.resize(NSTEPS);
// 		waveformI.resize(NSTEPS);
// 		waveformV = 0;
// 		waveformI = 0;
// 
// 		double E0 = Data.V0/dx;
// 		if (E0>MaxFieldValue)	//update MaxFieldValue for better appearance in the movie, not in the far field
// 		{
// 			MaxFieldValue=E0;
// 		}
// 		double distance = dx/2;
// 		double time_offset = dt; //ensures that the corrected voltage element in the line, and not waveformV, 
// // has the desired Gaussian shape. Since V is updated first in the coaxial line, Vinc is dt ahead of the voltage. Therefore,
// // the first element of the waveformV array effectively corrects the voltage that is dt ahead in time.
// 		for (int n=0; n<NSTEPS; n++)
// 		{
// 			waveformV(n) = Data.V0*waveform(n*dt-Data.delay*Data.tau+time_offset,Data.tau);
// 			waveformI(n) = Data.V0/Data.Zc*waveform((n-0.5)*dt-Data.delay*Data.tau-distance/Data.vp+time_offset,Data.tau);
// 			//Since the voltage is updated first in the coaxial line, the incident V is 0.5dt ahead of the incident I. 
// 			//I is in the total-field region, dx/2 above V. Therefore, it lags V by dx/2/vp.
// 		}
// 	}
// 
// 	//ofstream tempfile("tempfile.file",ios::binary);
// 	//for (int n=0; n<NSTEPS; n++)
// 	//{
// 	//	double temp = waveformV(n);
// 	//	tempfile.write((char*)&temp,sizeof(temp));
// 	//}
// 	//tempfile.close();
// 
// 	if (PortBelongsToNode&&((Data.mode=="receive")||(Data.mode=="probe")))
// 	{
// 		//open the file for recording the voltage or current waveform on the entire coaxial line
// 		CoaxMovieFile.open(CoaxMovieFileName.c_str(),ios::binary);
// 		if (!CoaxMovieFile)
// 		{
// 			cout << "Error opening coaxial line movie file!" << endl << endl;
// 			exit(-1);
// 		}
// 		CoaxMovieFile.write((char*)&(Data.V0),sizeof(Data.V0));
// 		CoaxMovieFile.write((char*)&length,sizeof(length));
// 		CoaxMovieFile.write((char*)&NSTEPS,sizeof(NSTEPS));
// 
// 		//open the file for recording the voltage waveform at "z=obs_point" on the coaxial line
// 		CoaxVoltageFile.open(CoaxVoltageFileName.c_str(),ios::binary);
// 		if (!CoaxVoltageFile)
// 		{
// 			cout << "Error opening voltage value recorder file!" << endl << endl;
// 			exit(-1);
// 		}
// 		CoaxVoltageFile.write((char*)&(Data.V0),sizeof(Data.V0));
// 		CoaxVoltageFile.write((char*)&NSTEPS,sizeof(NSTEPS));
// 		CoaxVoltageFile.write((char*)&dt,sizeof(dt));
// 	}
// 
// 	PreviousV = 0;
}

void Cport::UpdateV(const int& n)
{
// 	if (PortBelongsToNode)
// 	{
// 		//update the voltage in the coaxial line
// 		for (k=2; k<=length; k++)
// 		{
// 			voltage(k) = voltage(k) - dt*Data.vp*Data.Zc/dx*(current(k)-current(k-1));
// 		}
// 
// 		//apply 1st order ABC
// 		voltage(1) = PreviousV + (S_coax-1)/(S_coax+1)*(voltage(2)-voltage(1));
// 		PreviousV = voltage(2);
// 
// 		//if the port is in transmitting mode, inject the incident voltage wave
// 		if ((Data.mode=="transmit")||(Data.mode=="probe"))
// 		{
// 			//voltage(feed_point) is in the scattered field region
// 			voltage(feed_point) = voltage(feed_point) + dt*Data.vp*Data.Zc/dx*waveformI(n);
// 		}
// 
// 		//finally, send terminal voltage to the 3D grid
// 		GapVoltage = voltage(length);
// 
// 		for (i=Data.Gap_back; i<=Data.Gap_front+1; i++)
// 		{
// 			for (j=Data.Gap_left; j<=Data.Gap_right; j++)
// 			{
// 				for (k=Data.Gap_lower; k<=Data.Gap_upper+1; k++)
// 				{
// 					Ey(i,j,k) = GapVoltage/(GapWidth*dx);
// 				}
// 			}
// 		}
// 
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 		//typedef unsigned short ElectricMaterialIndexType;
// 
// 		//extern Array<ElectricMaterialIndexType,3> Media_Ey;
// 		//extern Array<double,1> Cb_Y;
// 		//double DiffGaussian(double t, double tau);
// 		//k=Data.Gap_z;
// 		//for (i=Data.Gap_back; i<=Data.Gap_front+1; i++)
// 		//{
// 		//	for (j=Data.Gap_left; j<=Data.Gap_right; j++)
// 		//	{
// 		//		Ey(i,j,k) += -Cb_Y(Media_Ey(i,j,k))*dx*
// 		//			(1e-4/(3*dx*dx))*DiffGaussian((n+0.5)*dt-6*40e-12,40e-12); //volume current waveform
// 		//	}
// 		//}
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 	}
}

void Cport::UpdateI(const int& n)
{
// 	if (PortBelongsToNode)
// 	{
// 		//update the current in the coaxial line
// 		for (k=1; k<=length-1; k++)
// 		{
// 			current(k) = current(k) - dt*Data.vp/Data.Zc/dx*(voltage(k+1)-voltage(k));
// 		}
// 
// 		//if the port is in transmitting mode, inject the incident current wave
// 		if ((Data.mode=="transmit")||(Data.mode=="probe"))
// 		{
// 			//current(feed_point) is in the total field region
// 			current(feed_point) = current(feed_point) + dt*Data.vp/Data.Zc/dx*waveformV(n);
// 		}
// 
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 		//Loop_y = OriginY-1;
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 		//!!!!!temporary!!!!!
// 
// 		//receive terminal current from the 3D grid
// 		//calculate the line integrals along the H-loops, then take the average
// 		//(loops are counterclockwise, since the feed draws current from the negative terminal)
// 		TerminalCurrent = 0;
// 		
// 		for (j=Loop_left; j<=Loop_right; j++)
// 		{
// 			//upper arm
// 			for (i=Loop_front; i>=Loop_back+1; i--)
// 			{
// 				TerminalCurrent += Hx(i,j,Loop_upper)*(-dx);	//(-dx) represents	 ->   ->   ->    ->
// 								//					 x  . dx = x  . -(x) dx = -dx
// 			}
// 			//lower arm
// 			for (i=Loop_back+1; i<=Loop_front; i++)
// 			{
// 				TerminalCurrent += Hx(i,j,Loop_lower)*(+dx);
// 			}
// 			//back arm
// 			for (k=Loop_upper; k>=Loop_lower+1; k--)
// 			{
// 				TerminalCurrent += Hz(Loop_back,j,k)*(-dx);
// 			}
// 			//front arm
// 			for (k=Loop_lower+1; k<=Loop_upper; k++)
// 			{
// 				TerminalCurrent += Hz(Loop_front,j,k)*(+dx);
// 			}
// 		}
// 
// 		if (clockwiseloop)
// 		{//if loop is clockwise, flip the sign of the current
// 			TerminalCurrent = -TerminalCurrent;
// 		}
// 
// 		current(length) = TerminalCurrent/GapWidth;	//assign current at the end of the coaxial line
// 	}
}

void Cport::RecordV()
{//writes the voltage information into files
// 	
// 	if (PortBelongsToNode&&((Data.mode=="receive")||(Data.mode=="probe")))
// 	{
// 		CoaxMovieFile.write((char*)voltage.data(),voltage.size()*sizeof(voltage(1)));
// 		//write the voltage value at "z=obs_point" into file
// 		observed_voltage = voltage(obs_point);	//record the normalized voltage
// 		CoaxVoltageFile.write((char*)&observed_voltage,sizeof(observed_voltage));
// 	}
}

// double Cport::waveform(double t, double tau)
// {//returns the basic waveform of the incident voltage wave
// 	if (Data.WaveShape=="gaussian")
// 	{
// 		return Gaussian(t,tau);	//normalized to 1
// 	}
// 	else if (Data.WaveShape=="diffgaussian")
// 	{
// 		return DiffGaussian(t,tau);		//normalized to 1
// 	}
// 	else if (Data.WaveShape=="doublediffgaussian")
// 	{
// 		return DoubleDiffGaussian(t,tau);		//normalized to 1
// 	}
// 	else if (Data.WaveShape=="triplediffgaussian")
// 	{
// 		return TripleDiffGaussian(t,tau);		//normalized to 1
// 	}
// 	else if (Data.WaveShape=="modulatedsinegaussian")
// 	{
// 		return ModulatedSineGaussian(t,tau,Data.w_0);
// 	}
// 	else if (Data.WaveShape=="diffmodulatedsinegaussian")
// 	{
// 		return DiffModulatedSineGaussian(t,tau,Data.w_0);
// 	}
// 	else if (Data.WaveShape=="modulatedcosinegaussian")
// 	{
// 		return ModulatedCosineGaussian(t,tau,Data.w_0);
// 	}
// 	else if (Data.WaveShape=="diffmodulatedcosinegaussian")
// 	{
// 		return DiffModulatedCosineGaussian(t,tau,Data.w_0);
// 	}
// 	else
// 	{
// 		if (rank==0)
// 		{
// 			cout << "Invalid wave shape (" << Data.WaveShape << ") for port " << PortIndex << " in node " 
// 				<< GridIndex << endl;
// 			exit(-1);
// 		}
// 	}
// }
