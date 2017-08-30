// main.cpp: David R. Morrison
// Entry point 

#include "main.h"
#include "cplex.h"

#include <cstdio>
#include <cstdlib>
// Removed for Windows compatibility (also not currently being used)
//#include <unistd.h>

using namespace std;

int main(int argc, char* argv[])
{
	try
	{
		opts options;
		const char* infile = parseOpts(argc, argv, options);

		FILE* jsonFile = nullptr;
		if (options.json_filename)
		{
			jsonFile = fopen(options.json_filename, "w");
			if (!jsonFile) throw ERROR << "Could not open " << options.json_filename << 
				" for writing";
		}

		CbfsData* cbfs = new CbfsData(options.m, options.posW, options.nullW, options.mob, options.cPara, options.maxDepth, options.probInterval);
		Cplex cplex(infile, jsonFile, cbfs, options.timelimit, options.disableAdvStart, options.randSeed);
		cplex.solve();

		if (jsonFile)
			fclose(jsonFile);
	}
	catch (IloCplex::Exception& e) { printf("ERROR: %s\n", e.getMessage()); }

	return 0;
}

const char* parseOpts(int argc, char* argv[], opts& options)
{
	options.m = Disable;
	options.posW = 1;
	options.nullW = 1;
	options.json_filename = nullptr;
	options.disableAdvStart = false;
	options.timelimit = 3600;
	options.mob = 3;
	options.cPara = 50;
	options.maxDepth = -1;
	options.probInterval = -1;
	options.randSeed = 0;

	int len = sizeof(optStrings)/(3 * sizeof(char*));
	int numOpts = 0;
	for (int i = 0; i < len; ++i)
		if (optStrings[i][0][0] == BOOL[0] || optStrings[i][0][0] == ARG[0]) ++numOpts;

	char flags[numOpts * 2 + 1];
	int pos = 0;
	for (int i = 0; i < len; ++i)
	{
		if (optStrings[i][0][0] != BOOL[0] && optStrings[i][0][0] != ARG[0]) continue;
		flags[pos++] = optStrings[i][1][0];
		if (optStrings[i][0][0] == ARG[0]) flags[pos++] = ':';
	}
	flags[pos] = '\0';

	int opt;
	stringstream str; 
	while ((opt = getopt(argc, argv, flags)) != -1)
	{
		switch (opt)
		{
		case 'D':
			options.m = Disable;
			break;
		case 'C':
			options.m = CplexOnly;
			break;
		case 'w':
			options.m = Weighted;
			break;
		case 'W':
			options.disableAdvStart = true;
			break;
		case 'L':
			options.m = LBContour;
			break;
		case 'R':
			options.m = RandCont;
			break;
		case 'l':
			options.cPara = atoi(optarg);
			break;
		case 'd':
			options.maxDepth = atoi(optarg);
			break;
		case 'p':
			options.probInterval = atoi(optarg);
			break;
		case 'r':
			options.randSeed = atoi(optarg);
			break;
		case 't':
			options.timelimit = atoi(optarg);
			break;
		case 'M':
			options.mob = atoi(optarg);
			break;
		case 'P':
			options.posW = atof(optarg);
			break;
		case 'N':
			options.nullW = atof(optarg);
			break;
		case 'j':
			options.json_filename = optarg;
			break;

		case 'h':
			usage(argv[0]);
			exit(0);
		default:
			usage(argv[0]);
			exit(-1);
		}
	}

	if (optind == argc - 1)
		return argv[optind];
	else
	{
		fprintf(stderr, "No file specified\n");
		usage(argv[0]);
		exit(-1);
	}
}

void usage(const char* name)
{
	int len = sizeof(optStrings)/(3 * sizeof(char*));
	printf("%s <options> [filename]\n", name);
	for (int i = 0; i < len; ++i)
	{
		if (optStrings[i][0][0] == BREAK[0]) 
			printf("\n");
		else if (optStrings[i][0][0] == SEP[0])
			printf("---------------------------------------------\n");
		else if (optStrings[i][0][0] == COMMENT[0])
			printf("\t\t%s\n", optStrings[i][2]);
		else
			printf("\t-%s: %s\n", optStrings[i][1], optStrings[i][2]);
	}
	printf("\n");
}




