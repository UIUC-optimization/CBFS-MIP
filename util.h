// util.h: Wenda Zhang, modified from David R. Morrison
// Utility functions

#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <set>
#include <algorithm>
#include <math.h>
#include <exception>
#include <string>
#include <sstream>
#include <iomanip>
#include <limits>
#include <time.h>
using std::vector;
using std::string;
using std::stringstream;
using std::setprecision;

enum Mode { Disable = -1, CplexOnly, Weighted, LBContour, RandCont, DiveComp, ContSel };
enum TieBreak { FIFO, LIFO, OG, ARB };
enum ContSelMode { Subtree, WBranch };
enum ContScoreMode { UCB1, SVM };

// options
struct opts
{
	Mode m;
	ContSelMode cm;
	double posW, nullW;
	const char* json_filename;
	bool disableAdvStart;
	int timelimit;
	int mob;
	int cPara;
	int maxDepth;
	int probInterval;
	int randSeed;
	int earlyTerm;
	bool jsDetail;
};

/*** Constants and Types ***/
const double Tolerance = 0.000000001;
const double Infinity = std::numeric_limits<double>::max();
const int MaxInt = std::numeric_limits<int>::max();

/*** Error handling ***/
class Error : public std::exception
{
public:
	Error(const char* file, const int line) 
	{ _msg << setprecision(10) << "ERROR (" << file << ":" << line << "): "; }
	Error(const Error& e) { _msg << e._msg.str(); }
	virtual const char* what() const throw() { return (_msg.str() + "\n").c_str(); }
	~Error() throw () {}

	Error& operator<<(const char* s) { _msg << s; return *this; }
	Error& operator<<(const string s) { _msg << s; return *this; }
	Error& operator<<(const int i) { _msg << i; return *this; }
	Error& operator<<(const double f) { _msg << f; return *this; }
	Error& operator<<(const bool b) { _msg << (b ? "true" : "false"); return *this; }
	Error& operator<<(const void* p) { _msg << p; return *this; }

private:
	stringstream _msg;
};

#define ERROR Error(__FILE__, __LINE__)

/*** Miscellaneous utility functions ***/
template <typename ContainerT, typename U>
bool contains(const ContainerT& con, const U& el) 
{ 
	typedef typename ContainerT::value_type ElementT;
	return find(con.begin(), con.end(), static_cast<ElementT>(el)) != con.end(); 
}

template <typename T>
void sort(T& con) { sort(con.begin(), con.end()); }

template <typename T>
void reverse(T& con) { reverse(con.begin(), con.end()); }

template <typename T>
void append(T& c1, const T& c2) { c1.insert(c1.end(), c2.begin(), c2.end()); }

template <typename U>
int vec_index(const vector<U>& vec, const U& el) 
{ return find(vec.begin(), vec.end(), el) - vec.begin(); }


#endif // UTIL_H



