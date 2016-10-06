#ifndef TIME_H
#define TIME_H

/*
 * @author Vasileios Zois
 * @email vzois@usc.edu
 *
 * Utility functions for measuring execution time.
 */

#include<ratio>
#include<ctime>
#include<chrono>
#include <iostream>
#include<string>

#include<ratio>
#include<climits>
#include<chrono>
#include<ctime>

typedef std::pair<uint64_t, uint64_t> arr2D;

typedef std::chrono::duration<double, std::ratio<1, 1000000000>> nanos;
typedef std::chrono::duration<double, std::ratio<1,1000>> millis;
typedef std::chrono::duration<double, std::ratio<1,1> > secs;
typedef std::chrono::duration<double, std::ratio<60, 1> > mins;
typedef std::chrono::duration<double, std::ratio<3600, 1> > hours;
typedef std::chrono::duration<double, std::ratio<60 * 60 * 24, 1> > days;
typedef std::chrono::duration<double, std::ratio<60 * 60 * 24 * 30, 1> > months;

template<class M>
class Time{
public:
	Time(){};
	~Time(){};
	void start();
	void reset();
	double lap(std::string);
	double lap();
private:
	std::chrono::high_resolution_clock::time_point tp1;
	std::chrono::high_resolution_clock::time_point tp2;
};

template<class M>
void Time<M>::start(){
	this->tp1 = std::chrono::high_resolution_clock::now();
}

template<class M>
void Time<M>::reset(){
	this->tp1 = std::chrono::high_resolution_clock::now();
}

template<class M>
double Time<M>::lap(){
	return this->lap("");
}

template<class M>
double Time<M>::lap(std::string comment){
	this->tp2 = std::chrono::high_resolution_clock::now();
	M time_span = std::chrono::duration_cast<M>(this->tp2 - this->tp1);

	double tt= time_span.count();
	std::cout << "Elapsed Time ( " << comment << " ): " << tt << std::endl;
	this->tp1 = std::chrono::high_resolution_clock::now();

	return tt;
}

#endif
