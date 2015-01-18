// pupillometry-Sequential.cpp : main project file.

#include "stdafx.h"
#include <ctime>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <opencv\cv.h>
#include <opencv\highgui.h>
#include <msclr/marshal_cppstd.h>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "dataStructs.h"
#include <windows.h>


using namespace System;
using namespace cv;
using namespace std;

#define BLINK_THRESHOLD 98
#define DEBUG_ON 1
#define CROOP_REPETITION 0
#define GRAPHIC_MODE 1
#define LOCK_CONDITION 40
#define LONG_BLINK_THRESHOLD 40
#define PIXEL_MM 0.05
#define VERBOSE_MODE 0
#define FRAME_WIDTH 320
#define FRAME_HEIGHT 240
#define M_PI 3.14159265358979323846


//Global Variables
//bool
bool _G_blinkCounterEnable = false;
bool _G_blinkFinished = false;
bool _G_blingStarted = false;
//int
int REMOVE_BLINK = 0;
int _G_blinkDelay = 0;
int _G_blinkDuration = 0;
int _G_lastBlinkDuration = 0;
int _longBlinkCounter = 0;
int _minutes = 0;
int _milisecond = 0;
int _seconds = 0;
//string Variables
string _finalFileName;

struct blinkStorage
{
	//int variables
	int _duration;
	int _startFrame;
};

struct fftSerie
{
	vector <double> _realBlinkArray;
	vector <double> _complexBlinkArray;
	int _totalElements;
};

struct ellipticFilterSerie
{
	vector <double> _inputData;
	vector <double> _outputData;
	int _totalElements;
};

//vector Variables
vector <blinkStorage> _blinkArray;

class signalFiltering
{
 
 ofstream _outputSS;
 public:
	signalFiltering::signalFiltering(vector<double> _inputArray,int __size, string _fileName)
	{
	  vector<double> __outputArray;
	  _outputSS.open(_fileName + ".ss", ofstream::out | ofstream::trunc); 
	  float __dbuffer[7];
	  int __k;
	  int __j;
	  static const float __dv0[7] = { 0.14751796914954551, -0.041634098790647008,
		-0.3454446156844046, -6.5511138358307581E-17, 0.34544461568440449,
		0.04163409879064732, -0.14751796914954565 };

	  static const float __dv1[7] = { 1.0, -1.9363880106582814, 1.8480252780314839,
		-1.5693896876190729, 1.2968629637769831, -0.62059713142908479,
		0.1607169988431125 };

	  for (__k = 0; __k < 6; __k++) 
	  {
		__dbuffer[__k + 1] = 0.0;
	  }

	  for (__j = 0; __j < __size; __j++) 
	  {
		for (__k = 0; __k < 6; __k++) 
		{
		  __dbuffer[__k] = __dbuffer[__k + 1];
		}

		__dbuffer[6] = 0.0;
		for (__k = 0; __k < 7; __k++) 
		{
		  __dbuffer[__k] += _inputArray[__j] * __dv0[__k];
		}

		for (__k = 0; __k < 6; __k++) 
		{
		  __dbuffer[__k + 1] -= __dbuffer[0] * __dv1[__k + 1];
		}

		__outputArray.push_back(__dbuffer[0]);
	    }
	  	//Perform SS
		long float __ss;
		for(int __i=0; __i < __size; __i++)
		{
			long float __d = __outputArray[__i]*__outputArray[__i];
			__ss = __ss+__d;
		}
		_outputSS<<__ss<<endl;
		_outputSS.close();
	}

};

class frameTick
{
public:
	void getTick()
	{
		_milisecond++;

		if(_milisecond > 120)
		{
			_milisecond = 0;
			_seconds++;
		}

		if(_seconds > 60)
		{
			_seconds = 0;
			_minutes++;
		}

		if(_minutes > 60)
			_minutes = 0;

	}
};

struct valueFrequency
{
	int _value;
	int _frequency;

	valueFrequency::valueFrequency()
	{
		_frequency = 0;
		_value = 0;
	}
};

//Perform FFT 
class computeFFT
{
public:
		computeFFT::computeFFT(vector<double> &_real, vector<double> &_imag)
		{
			transform(_real,_imag);

		}
	
private:
	void transform(vector<double> &_real, vector<double> &_imag) 
	{
		if (_real.size() != _imag.size())
		{
			cout<<"Mismatched Lengths..."<<endl;
			cout<<"ERROR 001"<<endl;
			exit(0);
		}
	
		size_t __n = _real.size();
		if (__n == 0)
			return;
		else if ((__n & (__n - 1)) == 0) 
			transformRadix2(_real, _imag);
		else 
			transformBluestein(_real, _imag);
	}


	void inverseTransform(vector<double> &_real, vector<double> &_imag) 
	{
		transform(_imag, _real);
	}


	void transformRadix2(vector<double> &_real, vector<double> &_imag) 
	{

		if (_real.size() != _imag.size())
		{
			cout<<"Mismatched Lengths..."<<endl;
			cout<<"ERROR 001"<<endl;
			exit(0);
		}

		size_t __n = _real.size();
		unsigned int __levels;
		{
			size_t __temp = __n;
			__levels = 0;
			while (__temp > 1) 
			{
				__levels++;
				__temp >>= 1;
			}
			if (1u << __levels != __n)
			{
			cout<<"The array size isn't power of 2..."<<endl;
			cout<<"ERROR 002"<<endl;
			exit(0);
			}
		}
	
		vector<double> __cosTable(__n / 2);
		vector<double> __sinTable(__n / 2);
		for (size_t __i = 0; __i < __n / 2; __i++) 
		{
			__cosTable[__i] = cos(2 * M_PI * __i / __n);
			__sinTable[__i] = sin(2 * M_PI * __i / __n);
		}
	
		for (size_t __i = 0; __i < __n; __i++) 
		{
			size_t __j = reverseBits(__i, __levels);
			if (__j > __i) 
			{
				double __temp = _real[__i];
				_real[__i] = _real[__j];
				_real[__j] = __temp;
				__temp = _imag[__i];
				_imag[__i] = _imag[__j];
				_imag[__j] = __temp;
			}
		}
	
		for (size_t __size = 2; __size <= __n; __size *= 2) 
		{
			size_t __halfsize = __size / 2;
			size_t __tablestep = __n / __size;
			
			for (size_t __i = 0; __i < __n; __i += __size) 
			{
				for (size_t __j = __i, __k = 0; __j < __i + __halfsize; __j++, __k += __tablestep) 
				{
					double __tpre =  _real[__j+__halfsize] * __cosTable[__k] + _imag[__j+__halfsize] * __sinTable[__k];
					double __tpim = -_real[__j+__halfsize] * __sinTable[__k] + _imag[__j+__halfsize] * __cosTable[__k];
					_real[__j + __halfsize] = _real[__j] - __tpre;
					_imag[__j + __halfsize] = _imag[__j] - __tpim;
					_real[__j] += __tpre;
					_imag[__j] += __tpim;
				}
			}
			if (__size == __n)
				break;
		}
	}


	void transformBluestein(vector<double> &_real, vector<double> &_imag) 
	{

		if (_real.size() != _imag.size())
		{
			cout<<"Mismatched Lengths..."<<endl;
			cout<<"ERROR 001"<<endl;
			exit(0);
		}

		size_t __n = _real.size();
		size_t __m;
		{
			size_t __target;
			if (__n > (SIZE_MAX - 1) / 2)
			{
			cout<<"Maximum size exceeded..."<<endl;
			cout<<"ERROR 003"<<endl;
			exit(0);
			}
			
			__target = __n * 2 + 1;
			for (__m = 1; __m < __target; __m *= 2) 
			{
				if (SIZE_MAX / 2 < __m)
				{
				cout<<"Maximum size exceeded..."<<endl;
				cout<<"ERROR 003"<<endl;
				exit(0);
				}
			}
		}
	
		vector<double> __cosTable(__n), __sinTable(__n);
		for (size_t __i = 0; __i < __n; __i++) 
		{
			double __temp = M_PI * (size_t)((unsigned long long)__i * __i % ((unsigned long long)__n * 2)) / __n;
			__cosTable[__i] = cos(__temp);
			__sinTable[__i] = sin(__temp);
		}
	
		vector<double> __areal(__m), __aimag(__m);
		for (size_t __i = 0; __i < __n; __i++) 
		{
			__areal[__i] =  _real[__i] * __cosTable[__i] + _imag[__i] * __sinTable[__i];
			__aimag[__i] = -_real[__i] * __sinTable[__i] + _imag[__i] * __cosTable[__i];
		}

		vector<double> __breal(__m), __bimag(__m);
		__breal[0] = __cosTable[0];
		__bimag[0] = __sinTable[0];
		for (size_t __i = 1; __i < __n; __i++) 
		{
			__breal[__i] = __breal[__m - __i] = __cosTable[__i];
			__bimag[__i] = __bimag[__m - __i] = __sinTable[__i];
		}
	
		vector<double> __creal(__m), __cimag(__m);
		convolve(__areal, __aimag, __breal, __bimag, __creal, __cimag);
	
		for (size_t __i = 0; __i < __n; __i++) 
		{
			_real[__i] =  __creal[__i] * __cosTable[__i] + __cimag[__i] * __sinTable[__i];
			_imag[__i] = -__creal[__i] * __sinTable[__i] + __cimag[__i] * __cosTable[__i];
		}
	}


	void convolve(const vector<double> &__x, const vector<double> &__y, vector<double> &__outArray) 
	{
		if (__x.size() != __y.size() || __x.size() != __outArray.size())
		{
			cout<<"Mismatched Lengths..."<<endl;
			cout<<"ERROR 001"<<endl;
			exit(0);
		}

		size_t __n = __x.size();
		vector<double> __ximag(__n), __yimag(__n), __zimag(__n);
		convolve(__x, __ximag, __y, __yimag, __outArray, __zimag);
	}


	void convolve(const vector<double> &__xreal, const vector<double> &__ximag, const vector<double> &__yreal, const vector<double> &__yimag, vector<double> &__outreal, vector<double> &__outimag) 
	{
		if (__xreal.size() != __ximag.size() || __xreal.size() != __yreal.size() || __yreal.size() != __yimag.size() || __xreal.size() != __outreal.size() || __outreal.size() != __outimag.size())
		{
			cout<<"Mismatched Lengths..."<<endl;
			cout<<"ERROR 001"<<endl;
			exit(0);
		}
	
		size_t __n = __xreal.size();
		vector<double> __xr(__xreal);
		vector<double> __xi(__ximag);
		vector<double> __yr(__yreal);
		vector<double> __yi(__yimag);
	
		transform(__xr, __xi);
		transform(__yr, __yi);
		for (size_t __i = 0; __i < __n; __i++) 
		{
			double __temp = __xr[__i] * __yr[__i] - __xi[__i] * __yi[__i];
			__xi[__i] = __xi[__i] * __yr[__i] + __xr[__i] * __yi[__i];
			__xr[__i] = __temp;
		}

		inverseTransform(__xr, __xi);
		for (size_t __i = 0; __i < __n; __i++) 
		{
			__outreal[__i] = __xr[__i] / __n;
			__outimag[__i] = __xi[__i] / __n;
		}
	}


	static size_t reverseBits(size_t __x, unsigned int __n) 
	{
		size_t __result = 0;
		unsigned int __i;
		for (__i = 0; __i < __n; __i++, __x >>= 1)
			__result = (__result << 1) | (__x & 1);
		return __result;
	}
};

class parametersCalculator
{
	//declaringFFT variables
	fftSerie _fourierSeire;
	//declaring filter series
	ellipticFilterSerie _filteredSerie;
	//String variables
	string _fileLine;
	string _outputFftFileName;
	//vector
	vector <double> _series;
	vector <double> _puiAverageSerie;
	//Float
	float _diameter;
	float _standardAverage;
	//ofstream Variables
	ofstream _outputFFT;
	ofstream _params;

	ofstream _powerSpectrum;
	ofstream _puiFile;
	ofstream _pdrFile;
	ofstream _pvrFile;
	ofstream _pfrFile;
	//int variables
	long float _puiCounter;
	long float	_puiAverage;
	long float	_pui;
	long float	_puiAverageSummation;
	int _totalNumberofFrames;

public:
	parametersCalculator::parametersCalculator(string _fileName)
	{
		ifstream _fileOpen(_fileName + ".plb", ofstream::in);

		_diameter = 0;
		_puiCounter = 0;
		_puiAverage = 0;
		_pui = 0;
		_totalNumberofFrames = 0;
		_puiAverageSummation = 0;
		_outputFftFileName = _fileName + ".pws";
		_outputFFT.open(_outputFftFileName, ofstream::out | ofstream::trunc); 
		//_params.open(_fileName + ".par", ofstream::out | ofstream::trunc);

		//_powerSpectrum.open(_fileName + ".pws", ofstream::out | ofstream::trunc);//power spectrum file
		_puiFile.open(_fileName + ".pui", ofstream::out | ofstream::trunc);//pui file
		_pdrFile.open(_fileName + ".pdr", ofstream::out | ofstream::trunc);//pdr file
		_pvrFile.open(_fileName + ".pvr", ofstream::out | ofstream::trunc);//pvr of squares file
		_pfrFile.open(_fileName + ".pfr", ofstream::out | ofstream::trunc);//pfr of squares file
		_standardAverage = 0;

		while(getline(_fileOpen, _fileLine))
		{
			_diameter = atof(_fileLine.c_str());
			_series.push_back(_diameter);
			_fourierSeire._realBlinkArray.push_back((double)(_diameter));
			_filteredSerie._inputData.push_back((double)(_diameter));
			_fourierSeire._complexBlinkArray.push_back(0);
			_standardAverage += _diameter;

			if(_diameter == -1)
				break;

			_totalNumberofFrames++;

			//PUI calculation
			if (_puiCounter <= 16)
			{
				_puiAverage += _diameter*PIXEL_MM;
				_puiCounter++;
			}
			else
			{
				_puiAverageSerie.push_back(_puiAverage);
				_puiAverage = 0;
				_puiCounter = 0;
			}
		}
		//serie average
		_standardAverage = _standardAverage/_totalNumberofFrames;
		_filteredSerie._totalElements = _totalNumberofFrames;
		//Performing FFT Power Spectrum
		computeFFT __fftCalculator(_fourierSeire._realBlinkArray,_fourierSeire._complexBlinkArray);
		//elliptic Filter Variable
		signalFiltering  elliptic(_filteredSerie._inputData,_filteredSerie._totalElements,_fileName);

		for(int __i=0; __i < _fourierSeire._realBlinkArray.size(); __i++)
		{
			_outputFFT<<_fourierSeire._realBlinkArray[__i]<<","<<_fourierSeire._complexBlinkArray[__i]<<endl;
		}

		_outputFFT.close();

		//PUI absolute difference of the average series
		for(int __i=1; __i < _puiAverageSerie.size(); __i++)
		{
			_puiAverageSummation += abs(_puiAverageSerie[__i] - _puiAverageSerie[__i-1]);
		}

		//compute the average from the PUI series
		_pui = _puiAverageSummation/(1+(_totalNumberofFrames-16)*((1/120)));

		_puiFile<<_pui<<endl;

		//Pupilare Fatique Ratio calculation (PFR)
		float __pfr = 0;
		float __pvr = 0;
		float __sumDeviation=0.0;
		float __standardDeviation = 0;
		float __seriesVariance = 0;
		float __segmentSize = 100;
		float __segmentCounter = 0;
		float __firstSegmentDiameterMean = 0;
		float __firstSegmentDiameterSD = 0;
		float __firstSegmentDiameterVariance = 0;
		float __firstSegmentDiameterSummation = 0;
		float __segmentDiameterSD = 0;
		float __segmentDiameterVariance = 0;
		float __segmentDiameterSummation = 0;

		float __segmentMean = 0;
		//vector <float> __prvArray;

		if (__standardDeviation == 0)
			__standardDeviation = 1;

		if (__firstSegmentDiameterMean == 0)
			__firstSegmentDiameterMean = 1;

		//PFR
		//_params<<"pfr"<<endl;
		for(int __i=0; __i<_totalNumberofFrames;++__i)
			__sumDeviation+=(_series[__i] -_standardAverage)*(_series[__i] -_standardAverage);

		__standardDeviation =sqrt(__sumDeviation/(_totalNumberofFrames-1));   
		__pfr = _standardAverage/__standardDeviation;
		_pfrFile<<__pfr/__firstSegmentDiameterMean<<endl;
		__seriesVariance = __sumDeviation/(_totalNumberofFrames-1);

		//PDR
		//_params<<"pdr,pvr"<<endl;
		for(int __i=0; __i<(int)_totalNumberofFrames/__segmentSize;++__i)
		{
			__firstSegmentDiameterMean += _series[__i];
			__firstSegmentDiameterSummation+=(_series[__i] -_standardAverage)*(_series[__i] -_standardAverage);
		}
			__firstSegmentDiameterSD =sqrt(__firstSegmentDiameterSummation/(_totalNumberofFrames-1));  
			__firstSegmentDiameterVariance = __firstSegmentDiameterSummation/(_totalNumberofFrames-1);
			__firstSegmentDiameterMean = __firstSegmentDiameterMean/(_totalNumberofFrames/__segmentSize);

		for(int __i=(int)_totalNumberofFrames/__segmentSize; __i<(int)_totalNumberofFrames;++__i)
		{
			if (__segmentCounter<(_totalNumberofFrames/__segmentSize))
			{
				__segmentCounter++;
				__segmentMean += _series[__i];
				__segmentDiameterSummation+=(_series[__i] -_standardAverage)*(_series[__i] -_standardAverage);
			}
			else
			{
				__segmentMean = __segmentMean/(_totalNumberofFrames/__segmentSize);
				__segmentDiameterVariance =  __segmentDiameterSummation/(_totalNumberofFrames-1);
				if(__firstSegmentDiameterMean == 0)
				__firstSegmentDiameterMean = 1;
				if(__firstSegmentDiameterVariance == 0)
				__firstSegmentDiameterVariance = 1;
				_pdrFile<<__segmentMean/__firstSegmentDiameterMean<<endl;
				_pvrFile<<__segmentDiameterVariance/__firstSegmentDiameterVariance<<endl;
				__segmentMean = 0;
				__segmentDiameterSummation = 0;
				__segmentDiameterVariance = 0;
				__segmentCounter = 0;
			}
		}
	}
	
};

class fileHandler
{
private:
	//Vector Variables
	vector <valueFrequency> _serieValues;
	vector <double> _series;
	//String Variables
	string _fileLine;
	string _outputFileName;
	//Int Variables
	int _diameter;
	int _previousDiameter;
	int _puiCounter;
	double _puiAverage;
	double _puiAverageSummation;
	int _repetitionCounter;
	int _totalElements;
	int _totalNumberofFrames;
	//floar variables
	long float _pui;
	//DebugVariables
	#if DEBUG_ON == 1
	//ofStream variables
	ofstream _outputFile2;


	#endif
	
	bool isOnVector(int _value)
	{
		for(int __i=0;__i<_totalElements;__i++)
		{
			if(_serieValues[__i]._value == _value)
			{
				_serieValues[__i]._frequency++;
				return true;
			}
		}
		return false;
	}

	bool isRemovable(int _value, int _threshold, int _minimumCount)
	{
		for(int __i=0;__i<_totalElements;__i++)
		{
			if(_value == _serieValues[__i]._value && _value <= _threshold && _serieValues[__i]._frequency < _minimumCount)
				return true;
		}

		return false;
	}

	bool isOnBlink(int _pointPosition)
	{
		for (int __i = 0; __i < _longBlinkCounter; __i++)
		{
			if(_pointPosition >= _blinkArray[__i]._startFrame)
			{
				if( _pointPosition <= _blinkArray[__i]._duration+_blinkArray[__i]._startFrame)
				{
					#if DEBUG_ON == 1
					_outputFile2<<"Current Point: "<<_pointPosition<<" start Point: "<<_blinkArray[__i]._startFrame<<" Final Point:"<< _blinkArray[__i]._duration+_blinkArray[__i]._startFrame<<endl;
					#endif
					return true;
				}
			}	
			else
				return false;
		}
	}

public:
	fileHandler::fileHandler(string _fileName, string _inputFile, string _fftName)
	{
		//debugVariable
		int _removalCounter = 0;
		//declaring ofStream
		ofstream _outputFile;
		//strings
		string _finalFileName = _fileName;
		string _fftFinalName = _fftName;
		//initializing variables
		_totalElements = 0;
		_totalNumberofFrames = 0;
		_outputFileName = _fftFinalName + ".plb";
		_outputFile.open(_outputFileName, ofstream::out | ofstream::trunc); 
		ifstream _fileOpen(_fileName, ofstream::in);
		valueFrequency _tempValue;
		//Analyzing Data
		while(getline(_fileOpen, _fileLine))
		{
			_diameter = atoi(_fileLine.c_str());
			_series.push_back(_diameter);
			_tempValue._value = _diameter;
			if(!isOnVector(_diameter))
			{
				_serieValues.push_back(_tempValue);
				_totalElements++;
			}
		}
		//sorting values using Bubble Sort
		//int variables
		valueFrequency _aux;
		int _k = _totalElements - 1 ;
		//Sorting
		for(int _i = 0; _i < _totalElements; _i++)
		{
			for(int _j = 0; _j < _k; _j++)
			{
				if(_serieValues[_j]._value > _serieValues[_j+1]._value)
				{
					_aux = _serieValues[_j];
					_serieValues[_j] = _serieValues[_j+1];
					_serieValues[_j+1] = _aux;
				}
			}
			_k--;
		}
		//New Data Analysis
		//int Variables;
		int _mean = 0;
		for(int _j=0; _j<_totalElements;_j++)
		{
			#if VERBOSE_MODE == 1
			cout<<_serieValues[_j]._value<<","<<_serieValues[_j]._frequency<<endl;
			#endif
			_mean += _serieValues[_j]._value;

		}
		_mean = _mean/_totalElements;

		#if VERBOSE_MODE == 1
		cout<<endl<<endl<<"Mean: "<<_mean<<" Total: "<<_totalElements<<endl<<endl;
		#endif
		//Remove Data
		int _newFileCounter = 0;
		ifstream _fileOpen2(_fileName, ofstream::in);
		while(getline(_fileOpen2, _fileLine))
		{
			_totalNumberofFrames++;	
			_diameter = atoi(_fileLine.c_str());
			_tempValue._value = _diameter;
			
			if(isRemovable(_diameter,0,0))
			{
				_removalCounter++;
				continue;
			}
			else
			{
				if( REMOVE_BLINK == 1)
				{
					if(_previousDiameter == _diameter)
						_repetitionCounter++;
					else
						_repetitionCounter=0;

					if(_repetitionCounter < 10)
					{
					_outputFile<<_diameter*PIXEL_MM<<endl;
					}
					_previousDiameter = _diameter ;
				}
				else
				{
					_outputFile<<_diameter*PIXEL_MM<<endl;
				}
			}

			_newFileCounter++;
		}

		_outputFile<<"-1"<<endl;
		std::string _separator = ":";
		std::string _string1 = std::to_string(_minutes) + _separator + std::to_string(_seconds) + _separator + std::to_string(_milisecond);
		_outputFile<<_string1<<endl;
		std::string _convertedString =msclr::interop::marshal_as< std::string >(DateTime::Today.ToString());
		_outputFile<<_convertedString<<endl;
		_outputFile<<_fileName<<endl;
		_outputFile<<_inputFile<<endl;
		_outputFile<<"-1"<<endl;

		//adding blink data to the file.
		for(int __i=0; __i < _longBlinkCounter; __i++)
		{
			_outputFile<<_blinkArray[__i]._duration<<endl;
			_outputFile<<_blinkArray[__i]._startFrame<<endl;
		}
				parametersCalculator parameters(_fftFinalName);

	}
	
};

class lockDetector
{
private:
	int _lockArray;
	int _repeatTime;
	int _previousPoint;

public:

	bool _detectLock(int _point)
	{
		if(_previousPoint == _point)
		{
			_repeatTime++;
		}
		else
		{
			_repeatTime = 0;
			_previousPoint = _point;
		}

		if(_repeatTime > LOCK_CONDITION)
			return true;
		else
			return false;
	}

};

class standardDeviator
{
private:
	int _size;
	int _totalPoints;
	float _data[3];
	float _dataPoint;
	float _previousStandardDeviation;

public:
	void _classConstructor()
	{
		_size = 0;
		_totalPoints = 0;
		_dataPoint = 0;
		_previousStandardDeviation = 0;
	}

	float _standardDeviation()
	{
		float __mean=0.0;
		float __sumDeviation=0.0;
		
		for(int __i=0; __i<_totalPoints;++__i)
		{
			__mean+=_data[__i];
		}

		__mean=__mean/_totalPoints;

		for(int __i=0; __i<_totalPoints;++__i)
			__sumDeviation+=(_data[__i]-__mean)*(_data[__i]-__mean);

		_previousStandardDeviation = sqrt(__sumDeviation/_totalPoints);  
		return 10*(__sumDeviation/_totalPoints);//sqrt(__sumDeviation/_totalPoints);   
	}


	void _addPoint(float _dataPoint)
	{
		_data[_size] = _dataPoint;
		_size++;

		if(_totalPoints < 3)
			_totalPoints++;

		if(_size == 3)
			_size = 0;
	}
};

//sdandardDeviator Variables
standardDeviator _G_sd;

class blinkDetector
{
private:
	//Binary Mat
	cv::Mat _blinkMat;
	//int variable
	int _height;
	int _whiteCounter;
	int _whitePercentage;
	int _width;

public:
	blinkDetector::blinkDetector(cv::Mat _frame,int _frameWidth, int _frameHeight)
	{
		_blinkMat = _frame;
		_height = _frameHeight;
		_width = _frameWidth;
	}

	bool _blinkDetector()
	{
		IplImage *__roiImage=  new IplImage(_blinkMat);
		cvSetImageROI(__roiImage,cvRect(0,0,320,240));
		_whiteCounter = cvCountNonZero(__roiImage);
		_whitePercentage = ((_whiteCounter*100) / (_width*_height));

		if (_whitePercentage >= BLINK_THRESHOLD)
		{
			_G_blinkFinished = false;
			_G_blingStarted = true;
			_G_blinkDuration++;
			_G_blinkDelay = 0;
			return true;
		}
		else
		{
			_G_blinkFinished=true;
			_G_blingStarted = false;
			_G_lastBlinkDuration = _G_blinkDuration;
			_G_blinkDuration = 0;
			return false;
		}
	}
};

class pupilEstimationDrawer
{
private:
	//diameter variables
	diameter _diameterDrawer;
	//Mat Variables
	cv::Mat _drawedMat;
		
public:
	
	void drawFrame(cv::Mat _frame, diameter _diameter, int _minutes, int _seconds, int _miliseconds, int _frameNumber, int _totalNumberOfFrames)
	{
		_diameterDrawer = _diameter;
		line(_frame, Point(_diameterDrawer.__smallerY+(_diameterDrawer.__diameter/2)-10,  _diameterDrawer.__smallerX+(_diameterDrawer.__diameter/2)), Point(_diameterDrawer.__smallerY+(_diameterDrawer.__diameter/2)+10, _diameterDrawer.__smallerX+(_diameterDrawer.__diameter/2)),cvScalar(0,255,255));
		line(_frame, Point(_diameterDrawer.__smallerY+(_diameterDrawer.__diameter/2),  _diameterDrawer.__smallerX+(_diameterDrawer.__diameter/2)-10), Point(_diameterDrawer.__smallerY+(_diameterDrawer.__diameter/2),  _diameterDrawer.__smallerX+(_diameterDrawer.__diameter/2)+10),cvScalar(0,255,255));
		line(_frame, Point(_diameterDrawer.__smallerY, _diameterDrawer.__smallerX), Point(_diameterDrawer.__smallerY, _diameterDrawer.__smallerX+(_diameterDrawer.__diameter)),cvScalar(0,0,255));
		line(_frame, Point(_diameterDrawer.__smallerY+_diameterDrawer.__diameter, _diameterDrawer.__smallerX), Point(_diameterDrawer.__smallerY+_diameterDrawer.__diameter, _diameterDrawer.__smallerX+(_diameterDrawer.__diameter)),cvScalar(0,0,255));
		line(_frame, Point(_diameterDrawer.__smallerY, _diameterDrawer.__smallerX+(_diameterDrawer.__diameter)), Point(_diameterDrawer.__smallerY+(_diameterDrawer.__diameter), _diameterDrawer.__smallerX+(_diameterDrawer.__diameter)),cvScalar(0,0,255));
		line(_frame, Point(_diameterDrawer.__smallerY, _diameterDrawer.__smallerX), Point(_diameterDrawer.__smallerY+_diameterDrawer.__diameter, _diameterDrawer.__smallerX),cvScalar(0,0,255));
		circle(_frame, Point(_diameterDrawer.__smallerY+(_diameterDrawer.__diameter/2),_diameterDrawer.__smallerX+(_diameterDrawer.__diameter/2)), (_diameterDrawer.__diameter/2) +1,cvScalar(0, 255, 0),1, 8, 0);

		std::string separator = ":";
		std::string _minuteString = std::to_string(_minutes);
		std::string _secondString = std::to_string(_seconds);
		std::string _frameNumberString = std::to_string(_totalNumberOfFrames);
		std::string _milisecondString = std::to_string(_miliseconds);
		std::string _string1 = "Time: " + _minuteString + separator + _secondString + separator + _milisecondString;
		std::string _string2 = "Frame Number: ";
		std::string _diameterString = std::to_string(_diameter.__diameter);
		std::string _numberString = std::to_string(_frameNumber);
		std::string _completeFrameNumberString = _string2 + _numberString + " of " +_frameNumberString;
		std::string _completeDiameterString =  "Diameter: " + _diameterString;

		//cv::putText(_frame, _string1, cvPoint(6,230), cv::FONT_HERSHEY_SIMPLEX, 0.3, CV_RGB(255,255,0), 1, 8);
		//cv::putText(_frame, _completeFrameNumberString, cvPoint(6,220), FONT_HERSHEY_SIMPLEX, 0.3, CV_RGB(255,255,0), 1, 8);
		//cv::putText(_frame, _completeDiameterString, cvPoint(6,210), FONT_HERSHEY_SIMPLEX, 0.3, CV_RGB(255,255,0), 1, 8);

		imshow("Estimated Diameter",_frame);
		cvWaitKey(1000/120);			
	}
};

class contourExtractor
{
private:
	//Mat Variables
	Mat _contourMat;
	//int variables
	int _currentPointX;
	int _currentPointY;
	int _diameter;
	int _findCountourTimeout;
	int _firstPointX;
	int _firstPointY;
	int _frameWidth;
	int _frameHeight;
	int _greaterY;
	int _previousPointX;
	int _previousPointY;
	int _smallerX;
	int _smallerY;
	//bool variables
	bool _firstPointContour;
	bool _pointFound;
	//diameter variables
	diameter _outPutDiameter;

public:
   contourExtractor(Mat _frame, int _width, int _height)
   {
		_contourMat = _frame;
		_frameHeight = _height;
		_frameWidth = _width;
		_firstPointContour = false;
		_pointFound = false;
		_greaterY = 0;
		_smallerX = 1000;
		_smallerY = 1000;
		_previousPointX = 0;
		_previousPointY = 0;
		_currentPointX = 0;
		_currentPointY = 0;
		_diameter = 0;
		_findCountourTimeout = 0;
   }

   void findFirstPoint()
	{
		for(int __j=0;__j<_contourMat.rows-1;__j++)
		{
			for(int __i=_contourMat.cols-1;__i>80;__i--)
			{
				if(_contourMat.at<unsigned char>(__j, __i) != 255)
				{
					continue;
				}
				else
				{
					_firstPointX = __j;
					_firstPointY = __i;
					_previousPointX = 0;
					_previousPointY = 0;
					_currentPointX = __j;
					_currentPointY = __i;
					_pointFound = true;
					break;
				}
			}

			if(_pointFound == true)
				break;
		}
	}

	diameter diameterExtractor()
	{
		findFirstPoint();
		if(_pointFound == true)
		{
			for(;;)
			{
				if(_currentPointX == _firstPointX && _currentPointY == _firstPointY && _firstPointContour == true)
					break;

				if(_findCountourTimeout > 500)
					break;
			
				if(_currentPointX-1 > 0 && _contourMat.at<unsigned char>(_currentPointX-1, _currentPointY) == 255  &&  _previousPointX != _currentPointX-1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointX = _currentPointX-1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}

				if( _currentPointX-1 > 0 && _currentPointY+1 < _contourMat.cols && _contourMat.at<unsigned char>(_currentPointX-1, _currentPointY+1) == 255  &&  _previousPointX != _currentPointX-1  &&  _previousPointY != _currentPointY+1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointX = _currentPointX-1;
					_currentPointY = _currentPointY+1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}

				if(_currentPointY+1 < _contourMat.cols && _contourMat.at<unsigned char>(_currentPointX, _currentPointY+1) == 255   &&  _previousPointY != _currentPointY+1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointY = _currentPointY+1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}
			
				if(_currentPointX+1 < _contourMat.rows && _currentPointY+1 < _contourMat.cols && _contourMat.at<unsigned char>(_currentPointX+1, _currentPointY+1) == 255  &&  _previousPointX != _currentPointX+1  &&  _previousPointY != _currentPointY+1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointX = _currentPointX+1;
					_currentPointY = _currentPointY+1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}

				if(_currentPointX+1 < _contourMat.rows && _contourMat.at<unsigned char>(_currentPointX+1, _currentPointY) == 255  &&  _previousPointX != _currentPointX+1 )
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointX = _currentPointX+1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}
	
				if(_currentPointX+1 < _contourMat.rows &&  _currentPointY-1 > 0  && _contourMat.at<unsigned char>(_currentPointX+1, _currentPointY-1) == 255  &&  _previousPointX != _currentPointX+1  &&  _previousPointY != _currentPointY-1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointX = _currentPointX+1;
					_currentPointY = _currentPointY-1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}

				if( _currentPointY-1 > 0  && _contourMat.at<unsigned char>(_currentPointX, _currentPointY-1) == 255   &&  _previousPointY != _currentPointY-1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointY = _currentPointY-1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}

				if(_currentPointX-1  > 0 &&  _currentPointY-1 > 0  && _contourMat.at<unsigned char>(_currentPointX-1, _currentPointY-1) == 255  &&  _previousPointX != _currentPointX-1  &&  _previousPointY != _currentPointY-1)
				{
					_contourMat.at<unsigned char>(_currentPointX, _currentPointY) =0;

					_previousPointX = _currentPointX;
					_previousPointY = _currentPointY;
					_currentPointX = _currentPointX-1;
					_currentPointY = _currentPointY-1;
					_firstPointContour = true;
					if(_currentPointY  > _greaterY )
						_greaterY  = _currentPointY;
					if(_currentPointY  < _smallerY )
						_smallerY  = _currentPointY;
					if(_currentPointX  < _smallerX )
						_smallerX  = _currentPointX;
					_findCountourTimeout++;
					continue;
				}
			break;
			}
			_outPutDiameter.__diameter =  (_greaterY -_smallerY );
			_outPutDiameter.__smallerX = _smallerX;
			_outPutDiameter.__smallerY = _smallerY;
			_outPutDiameter.__diameterDetectionFail = false;

			if(_outPutDiameter.__diameter < 41)
			{
				_outPutDiameter.__diameter = -2;
				_outPutDiameter.__diameterDetectionFail = true;
			}
		}
		else
		{
			_outPutDiameter.__diameter = -2;
			_outPutDiameter.__diameterDetectionFail = true;
		}

	return _outPutDiameter;
	}
};


int main(array<System::String ^> ^args)
{
	//diameter variables
	diameter _diameter;
	diameter _previousDiameter;
	//int variables
	int _blinkCounter = 0;
	int _frameHeight = 0;
	int _frameRate = 0;
	int _frameWidth = 0;
	int _totalNoFrames = 0;
	int __k = 0;
	//string variables
	string _fileName;
	//ofStream variables
	ofstream _outputFile;
	//timer variables
	clock_t startTime;
	//cvMat variables
	Mat _frame;
	Mat _postProcessingFrame;
	//bool variables
	bool _blinking = false;
	bool _diameterDetectionFail = false;
	bool _longBlinkFinished = false;
	//frame clock
	frameTick _clockTick;
	//lock detector
	lockDetector _lockDetector;
	//blink storage variables
	blinkStorage _longBlink;
	//initializing variables
	_previousDiameter.__diameter = 0;
	_previousDiameter.__smallerX = 0;
	_previousDiameter.__smallerY = 0;
	_previousDiameter.__diameterDetectionFail = false;
	_G_sd._classConstructor();
	//starting program
	cout<<"================================================================================"<<endl;
	cout<<"                     PUPILAB V 1.O - SEQUENTIAL ALGORITHM                       "<<endl;
	cout<<"================================================================================"<<endl;
	string _fftName;
	VideoCapture _capture(0);
	_capture.set(CV_CAP_PROP_FRAME_WIDTH, FRAME_WIDTH);
	_capture.set(CV_CAP_PROP_FRAME_HEIGHT, FRAME_HEIGHT);

	if(!_capture.isOpened())
	{
		_capture.open(0);
		if (!_capture.isOpened()) 
		{
		cout<<"Invalid Video File. Please check the input file"<<endl;
		return 0;
		}
	}
	else
	{
		
		cout<<"Please enter the output file name: ";
		cin>>_fileName;
		_finalFileName = _fileName;
		_fftName = _fileName;
		_fileName += ".raw";
		_outputFile.open(_fileName, ofstream::out | ofstream::trunc); 
		//starting timmer
		string _answer = "";
		bool _validOption = false;
		while(true)
		{
			cout<<"Remove blink sequences from the output file (remove repetitions) (y/n)? ";
			cin>>_answer;
			cout<<endl;
			if(_answer == "n")
			{
				REMOVE_BLINK = 0;
				break;
			}

			if(_answer == "y")
			{
				REMOVE_BLINK = 1;
				break;
			}

		}
		startTime = clock();
		//starting detection;
		cout<<"Starting Video Analysis. Please Wait..."<<endl;
		while(1)
		{
			__k++;
			_clockTick.getTick();
			_capture>>_frame;
			cvtColor(_frame, _postProcessingFrame, CV_BGR2GRAY);
			GaussianBlur(_postProcessingFrame,  _postProcessingFrame,  Size(9, 9), 2, 2); 
			threshold(_postProcessingFrame,  _postProcessingFrame, 21, 250, CV_THRESH_BINARY);
			blinkDetector _detector(_postProcessingFrame,FRAME_WIDTH,FRAME_HEIGHT);
			_blinking = _detector._blinkDetector();
			//cheking long blink
			if(_blinking == true)
			{
				_longBlink._startFrame = __k;
				_blinkCounter++;
				_longBlinkFinished = false;
			}
			else
			{
				_longBlinkFinished = true;
			}

			if(_blinkCounter > LONG_BLINK_THRESHOLD)
			{
				if(_longBlinkFinished == true)
				{
				_longBlinkFinished = false;
				_longBlink._duration =_blinkCounter - LONG_BLINK_THRESHOLD;
				_blinkArray.push_back(_longBlink);
				_longBlinkCounter++;
				_blinkCounter = 0;
				#if VERBOSE_MODE == 1
				cout<<"Long Blink Detected - Duration: "<<_blinkArray[_longBlinkCounter-1]._duration<<" Started At: "<<_blinkArray[_longBlinkCounter-1]._startFrame<<endl;
				#endif
				}
			}

			dilate(_postProcessingFrame,  _postProcessingFrame, Mat(), Point(-1, -1), 2, 1, 1);
			erode(_postProcessingFrame,  _postProcessingFrame, Mat(), Point(-1, -1), 2, 1, 1);
			Canny(_postProcessingFrame,  _postProcessingFrame, 20, 20*2, 3 );
			imshow("Threshold",_postProcessingFrame);
			contourExtractor _extractor(_postProcessingFrame,_frameWidth,_frameHeight);
			_diameter = _extractor.diameterExtractor();

			if(__k == 0)
				_previousDiameter.__diameter =  _diameter.__diameter;
			
			if(_diameter.__diameterDetectionFail == true)
			{
				 _G_sd._addPoint(_diameter.__diameter);
				_diameter = _previousDiameter;
				
			}
			else
			{
				_G_sd._addPoint(_diameter.__diameter);
				if(_G_sd._standardDeviation() > 10)
				{
					if(_lockDetector._detectLock(_previousDiameter.__diameter) == true)
					{
						if(_blinking ==false)
						{
							_previousDiameter = _diameter;
							#if CROOP_REPETITION == 1
							_outputFile<<_diameter.__diameter*PIXEL_MM<<endl;
							#endif
						}
						else
						{
							_diameter = _previousDiameter;
						}
					}
					else
					{
					_diameter = _previousDiameter;
					}
				}
				else
				{
					if(abs(_diameter.__diameter - _previousDiameter.__diameter) < 3)
					{
						if(_blinking ==false)
						{
							_previousDiameter = _diameter;
							#if CROOP_REPETITION == 1
							_outputFile<<_diameter.__diameter<<endl;
							#endif
						}
						else
						{
							_diameter = _previousDiameter;
						}
					}
					else
					{
						if(_lockDetector._detectLock(_previousDiameter.__diameter) == true)
						{
							_previousDiameter = _diameter;
							#if CROOP_REPETITION == 1
							_outputFile<<_diameter.__diameter<<endl;
							#endif
						}
						else
						{
							_diameter = _previousDiameter;
						}
					}
				}
			}

		#if CROOP_REPETITION == 0
			_outputFile<<_diameter.__diameter<<endl;
		#endif

		#if GRAPHIC_MODE == 1
			pupilEstimationDrawer _drawer;
			_drawer.drawFrame(_frame,_diameter,_minutes,_seconds,_milisecond,__k,_totalNoFrames);
		#endif

		#if VERBOSE_MODE == 1
			cout<<_diameter.__diameter*PIXEL_MM<<" "<<_blinking<<" "<<__k<<" "<<_G_sd._standardDeviation()<<" "<<endl;
		#endif

		char _key = cvWaitKey(10);

        if (char(_key) == 27)
		{
            break;  
        }


		}
	}
	
	cout<<endl<<"Saving Output File, Please Wait..."<<endl<<endl;
	//fileHandler Variable;
	fileHandler _handler(_fileName,"Camera Data",_fftName);
	remove(_fileName.c_str());

	system("pause");
    return 0;
}
