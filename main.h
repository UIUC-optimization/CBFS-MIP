// Main.h: David R. Morrison
// Include file for the dual network graph visualizer

#ifndef MAIN_H
#define MAIN_H

#define BOOL "b"
#define ARG "a"
#define COMMENT "c"
#define BREAK "x"
#define SEP "h"

#include "util.h"

// Check for gcc version < 4.6, in which case nullptr needs to be implemented
#if __GNUC__ == 4 && __GNUC_MINOR__ < 6
    const                        // this is a const object...
    class {
    public:
      template<class T>          // convertible to any type
        operator T*() const      // of null non-member
        { return 0; }            // pointer...
      template<class C, class T> // or any type of null
        operator T C::*() const  // member pointer...
        { return 0; }
    private:
      void operator&() const;    // whose address can't be taken
    } nullptr = {};              // and whose name is nullptr
#endif

struct opts
{
	Mode m;
	double posW, nullW;
	const char* json_filename;
	bool disableAdvStart;
	int timelimit;
	int mob;
	int cPara;
};

static const char* optStrings[][3] =
{
	{BOOL, "C", "Cplex only (no CBFS)"},
	{BOOL, "L", "LBContour (default 50)"},
	{BOOL, "D", "Disable callbacks entirely (default)"},
	{BOOL, "b", "Use no contours (BFS)"},
	{BOOL, "w", "Use weighted contour rule (default pw = 1, nw = 1)"},
	{BOOL, "W", "Disable warm starts and cut generation"},
	{ARG, "t", "Time limit"},
	{SEP, "", ""},
	{ARG, "M", "Measure of best: 1 - estimate, 2 - lower bound, "
		"3 - worst estimate, 4 - worst lower bound"},
	{ARG, "l", "LB contour number"},
	{ARG, "P", "Weight for the positive assignments"},
	{ARG, "N", "Weight for the null assignments"},
	{ARG, "j", "JSON output file"},
	{SEP, "", ""},
	{BOOL, "h", "Help"}
};

int main(int argc, char* argv[]);
const char* parseOpts(int argc, char* argv[], opts& options);
void usage(const char* name);

#endif // MAIN_H


