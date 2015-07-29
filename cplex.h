#ifndef CPLEX_H
#define CPLEX_H 1

/* 
 * cplex.h: David R. Morrison, Aug. 2013
 *
 * A wrapper around CPLEX to use CBFS for MIPs
 */

#include "util.h"

#include <map>
#include <list>
#include <utility>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>

using namespace std;

class CbfsNodeData;
typedef IloCplex::MIPCallbackI::NodeId NID;
typedef multimap<double, CbfsNodeData*> CbfsHeap;
typedef map<int, CbfsHeap> ContourMap;
typedef map<int, pair<int,int>> RangeMap;
typedef list<CbfsNodeData*> CbfsDive;
struct opts;

class CbfsData
{
public:
	CbfsData(Mode m, double posW, double nullW, int mob, int cPara) : 
		mMode(m), mPosW(posW), mNullW(nullW), nIters(0), mCurrContour(mContours.begin()), mMob(mob),
		bestLB(-INFINITY), bestUB(INFINITY), mcontPara(cPara), mDiveCount(0), mDiveStatus(true), mProbStep(0), mDiveStart(0)
	{ }
	~CbfsData();

	void addNode(CbfsNodeData* nodeData);
	void delNode(CbfsNodeData* nodeData);
	void updateContour();
	void updateBounds(double lb, double ub, double gap);
	void updateProfile(int pb, int nb);
	void estimateTree();
	void checkTermi(double time);
	void probStep();
	int calContour(double lb);

	NID getNextNode();

	Mode getMode() { return mMode; }
	void setSense(IloObjective::Sense s) { mSense = s; }
	long getNumIterations() { return nIters; }
	void setNumIterations(long i) { nIters = i; }
	void setTime(double time) { startTime = time; }
	void updateContourBegin() { mCurrContour = mContours.begin(); }

private:
	Mode mMode;
	IloObjective::Sense mSense;
	double mPosW, mNullW, bestLB, bestUB, startTime, mReOptGap;
	long nIters;
	vector<int> treeProfile;

	ContourMap mContours;
	ContourMap::iterator mCurrContour;
	CbfsDive mDiveCand, mProbLeft, mProbRight, mProbPre;
	int mMob, mcontPara, count, mDiveCount, mProbStep, mTreeDepth, mDiveStart;
	bool mDiveStatus;
};

class Cplex
{
public:
	Cplex(const char* filename, FILE* jsonFile, CbfsData* cbfs, int timelimit, 
			bool disableAdvStart);

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
	int depth, numPos, numNull;

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
