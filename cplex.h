#ifndef CPLEX_H
#define CPLEX_H 1

/* 
 * cplex.h: David R. Morrison, Aug. 2013
 *
 * A wrapper around CPLEX to use CBFS for MIPs
 */

#include "util.h"

#include <map>
#include <unordered_set>
#include <list>
#include <utility>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>

using namespace std;

class CbfsNodeData;
typedef IloCplex::MIPCallbackI::NodeId NID;
typedef multimap<double, NID> CbfsHeap;
typedef map<int, CbfsHeap> ContourMap;
typedef map<int, pair<int,int>> RangeMap;
typedef list<CbfsNodeData*> CbfsDive;
struct opts;

class CbfsData
{
public:
	CbfsData(Mode m, double posW, double nullW, int mob, int cPara, int maxDepth, int probInterval) : 
		mMode(m), mPosW(posW), mNullW(nullW), bestLB(-INFINITY), bestUB(INFINITY), nIters(0),
		mCurrContour(mContours.begin()), mMob(mob), mContPara(cPara), mDiveCount(0), mDiveStart(0),
		mMaxDepth(maxDepth), mProbStep(0), mProbInterval(probInterval)
	{
		mDiveStatus = (maxDepth > 0) ? true : false;
		mProbStatus = (probInterval > 0) ? true : false;
		mNInfeasibleCont = 20;
		mContScores.resize(cPara, 0);
		mTieBreak = OG;
		srand(time(0)); // Randomness is used in random contour
	}
	~CbfsData();

	void addNode(CbfsNodeData* nodeData);
	void delNode(CbfsNodeData* nodeData);
	void updateContScores();
	void updateContScores(int contID, int score);
	void updateBounds(double lb, double ub, double gap);
	void probStep();
	int calContour(CbfsNodeData* nodeData);

	NID getNextNode();
	ContourMap::iterator getNextCont();

	Mode getMode() { return mMode; }
	void setSense(IloObjective::Sense s) { mSense = s; }
	long getNumIterations() { return nIters; }
	double getBestUB() { return bestUB; }
	void setBestUB(double ub) { bestUB = ub; }
	void setNumIterations(long i) { nIters = i; }
	void updateContourBegin() { mCurrContour = mContours.begin(); }

private:
	Mode mMode;
	TieBreak mTieBreak;
	IloObjective::Sense mSense;
	double mPosW, mNullW, bestLB, bestUB, mReOptGap;
	long nIters;

	ContourMap mContours;
	ContourMap::iterator mCurrContour;
	CbfsDive mDiveCand, mProbLeft, mProbRight, mProbPre;
	int mMob, mContPara;
	int mDiveCount, mDiveStart, mMaxDepth, mProbStep, mProbInterval; // Diving and Probing parameters
	int mNIntVars, mNInfeasibleCont;								 // Infeasible variable contour parameters
	bool mDiveStatus, mProbStatus;									 // Diving and Probing triggers

	vector<int> mContScores;
};

class Cplex
{
public:
	Cplex(const char* filename, FILE* jsonFile, CbfsData* cbfs, int timelimit, 
			bool disableAdvStart, int rand);

	void solve();
	IloObjective::Sense getSense() { return mObj.getSense(); }

private:
	IloEnv mEnv;
	IloModel mModel;
	IloCplex mCplex;
	IloObjective mObj;

	FILE* mJsonFile;

};

class CbfsBranchCallback : public IloCplex::BranchCallbackI
{
public:
	CbfsBranchCallback(IloEnv env, CbfsData* cbfs, FILE* jsonFile) : 
		IloCplex::BranchCallbackI(env), mCbfs(cbfs), mJsonFile(jsonFile) { }
	IloCplex::CallbackI* duplicateCallback() const 
		{ return (new (getEnv()) CbfsBranchCallback(*this)); }
	void main();

private:
	CbfsData* mCbfs;
	FILE* mJsonFile;
};

class CbfsNodeCallback : public IloCplex::NodeCallbackI
{
public:
	CbfsNodeCallback(IloEnv env, CbfsData* cbfs) : 
		IloCplex::NodeCallbackI(env), mCbfs(cbfs) { }
	IloCplex::CallbackI* duplicateCallback() const 
		{ return (new (getEnv()) CbfsNodeCallback(*this)); }
	void main();

private:
	CbfsData* mCbfs;
};

class CbfsNodeData : public IloCplex::MIPCallbackI::NodeData
{
public:
	int depth, numPos, numNull, numInfeasibles;

	int contour;
	NID id;
	double estimate, lpval, incumbent;
	RangeMap rmap;
	int parent;
	bool explored;

	CbfsNodeData(CbfsData* cbfs, FILE* jsonFile) : 
		explored(false), mCbfs(cbfs), mJsonFile(jsonFile) {}
	~CbfsNodeData();

private:
	CbfsData* mCbfs;
	FILE* mJsonFile;
};
#endif // CPLEX_H
