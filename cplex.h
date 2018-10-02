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
typedef multimap<double, CbfsNodeData*> CbfsHeap;
typedef map<int, CbfsHeap> ContourMap;
typedef map<int, pair<int,int>> RangeMap;
typedef list<CbfsNodeData*> CbfsDive;
struct opts;

class CbfsData
{
public:
	CbfsData(Mode m, double posW, double nullW, int mob, int lbPara, int maxDepth, int probInterval, int term) : 
		mMode(m), mPosW(posW), mNullW(nullW), bestLB(-INFINITY), bestUB(INFINITY), nIters(0),
		mCurrContour(mContours.begin()), mMob(mob), mLBContPara(lbPara),
		mDiveMaxDepth(maxDepth), mProbInterval(probInterval),
		earlyTermIter(term)
	{
		// WRA Contour Score Matrix Initialization
		mContScores.resize(lbPara, 0);
		// Tie-breaking Rule Initialization
		mTieBreak = OG;
		// Dive Initialization
		mDiveStatus = (maxDepth > 0) ? true : false;
		// Check if use CPLEX Dive
		mUseCplexDive = (maxDepth == -1) ? true : false;
		preId = prepreId = 0;

		mNumDive = 0; mBFSInterval = 1000; mDiveMinDepth = 0;
		mCurDiveCount = mDiveStartCont = 0;
		mDiveTol = 0.25; mExpdMaxDepth = 0;
		// Prob Initialization
		mProbStatus = (probInterval > 0) ? true : false;
		mProbStep = mPreProbEnd = 0;
		// Set random seed for CBFS
		srand(time(0));
		// EXPERIMENT: Hybrid Best Parameter
		mHybridBestPara = 0.1;
		// EXPERIMENT: Dive Contour
		mDiveContPara = 200;
		prtIDLstInstedNode = -1;
		penaltyPara = 1; penaltyOn = false;
		// EXPERIMENT: Early Termination
		mNumIntInfsblUnexpl = mSumIntInfsblUnexpl = 0;
		mEarlyTermOn = (earlyTermIter > 0) ? true : false;
		//geoMeanShift = 1;
		
	}
	~CbfsData();

	void addNode(CbfsNodeData* nodeData);
	void delNode(CbfsNodeData* nodeData);
	void updateContScores();
	void updateContScores(int contID, int score);
	void updateBounds(double lb, double ub, double gap);
	void updateExpdMaxDepth(int depth);
	void probStep();
	void terminateProb();
	void printScores();
	void printCriteria();
	int calContour(CbfsNodeData* nodeData);
	int getNumNodesInCont(int contID);

	NID getNextNode();
	NID getNextNodePenalty();
	NID getNextNodeBFS();
	ContourMap::iterator getNextCont();

	Mode getMode() { return mMode; }
	void setSense(IloObjective::Sense s) { mSense = s; }
	long getNumIterations() { return nIters; }
	double getBestUB() { return bestUB; }
	void setBestUB(double ub) { bestUB = ub; }
	void setNumIterations(long i) { nIters = i; }
	void updateContourBegin() { mCurrContour = mContours.begin(); }

	// EXPERIMENT: Early termination and criteria for best contour strategy
	void calCriteria();
	bool isEarlyTermOn() { return mEarlyTermOn; }
	bool isEarlyTermToStart() { return (nIters == earlyTermIter); }
	void setVarArray(IloNumVarArray input) { mAllVars = input; }
	void setIntVarArray(IloNumVarArray input) { mAllIntVars = input; }
	void setVarCount(int input) { mNVars = input; }
	void setIntVarCount(int input) { mNIntVars = input; }
	void setDiveMaxDepth(int input) { mDiveMaxDepth = input; }
	void setRootLB(double lb) { rootLB = lb; }
	void setSplxIterInfo(int numSplxIter) { mSumIterSplx = numSplxIter; }
	void setJsonFile(FILE* json) { mJsonFile = json; }
	void setRootTime() { clock_gettime(CLOCK_MONOTONIC, &mRoot); }
	timespec getRootTime() { return mRoot; }
	IloNumVarArray getIntVarArray() { return mAllIntVars; }
	CbfsDive rcntAddedNodes;

	// EXPERIMENT: Use CPLEX Dive
	bool isCplexDiveOn() { return mUseCplexDive; }
	bool isInSameDive(int id);
	void updateDivePre(int id);

private:
	Mode mMode;
	TieBreak mTieBreak;
	IloObjective::Sense mSense;
	double mPosW, mNullW, bestLB, bestUB, mReOptGap;
	long nIters;

	ContourMap mContours;
	ContourMap::iterator mCurrContour;
	CbfsDive mDiveCand, mProbOneSide, mProbCand;
	int mMob, mLBContPara;
	int mNumDive, mBFSInterval;
	int mCurDiveCount, mDiveStartCont;										// Diving parameters
	int mDiveMinDepth, mDiveMaxDepth;										
	double mDiveTol;
	int mProbStep, mProbInterval, mPreProbEnd;								// Probing parameters
	bool mDiveStatus, mProbStatus;											// Diving and Probing flags
	vector<int> mContScores;

	// EXPERIMENT: Dive contour, deviating path penalty
	int mDiveContPara, prtIDLstInstedNode;
	int preNodeID, preContID;												// ID of the node explored just prior to current one
	double preNodeLB, penaltyPara;
	bool penaltyOn;

	// EXPERIMENT: Hybrid estimate and value
	double mHybridBestPara;

	// Record the maximum depth in current run
	int mExpdMaxDepth;

	// EXPERIMENT: Early termination and criteria for best contour strategy
	bool mEarlyTermOn;
	int earlyTermIter;
	double rootLB;
	//double geoMeanShift;
	int mNumUnexplNodes;													// Criteria 1: Number of unexplored nodes
	double mLbImprv;														// Criteria 2: LB percentage improvement from root LB
	int mNumIntInfsblUnexpl;												// Criteria 3: Sum of number of integer-infeasible variables of unexplored ndoes
	int mSumDpthUnexpl;														// Criteria 4: Sum of depth of unexplored nodes
	double mSumIntInfsblUnexpl;												// Criteria 5: Sum of integer infeasibilites of unexplored ndoes
	int mSumIterSplx;														// Criteria 8: Sum of simplex iterations
	int mNVars, mNIntVars;													// Num of Decision Variables
	timespec mRoot;

	// Use CPLEX Dive
	bool mUseCplexDive;
	int preId, prepreId;

	IloNumVarArray mAllVars, mAllIntVars;

	FILE* mJsonFile;
};

class Cplex
{
public:
	Cplex(const char* filename, FILE* jsonFile, CbfsData* cbfs, int timelimit, 
			bool disableAdvStart, int rand, bool jsonDetail);

	void solve();
	IloObjective::Sense getSense() { return mObj.getSense(); }
	IloNumVarArray getVars();
	IloNumVarArray getIntVars();

private:
	IloEnv mEnv;
	IloModel mModel;
	IloCplex mCplex;
	IloObjective mObj;

	FILE* mJsonFile;
	bool mJsonDetail;

	int varsCount;
	int intVarsCount;
	
	CbfsData* mCbfs;
};

class CbfsBranchCallback : public IloCplex::BranchCallbackI
{
public:
	CbfsBranchCallback(IloEnv env, CbfsData* cbfs, FILE* jsonFile, bool jsonDetail) : 
		IloCplex::BranchCallbackI(env), mCbfs(cbfs), mJsonFile(jsonFile), mJsonDetail(jsonDetail) { }
	IloCplex::CallbackI* duplicateCallback() const 
		{ return (new (getEnv()) CbfsBranchCallback(*this)); }
	void main();

private:
	CbfsData* mCbfs;
	FILE* mJsonFile;
	bool mJsonDetail;
	IloArray<IloCplex::ControlCallbackI::IntegerFeasibility> varsFeasibility;
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

	int infNum;
	double infSum;

	CbfsNodeData(CbfsData* cbfs, FILE* jsonFile, bool jsonDetail) : 
		explored(false), mCbfs(cbfs), mJsonFile(jsonFile), mJsonDetail(jsonDetail) {}
	~CbfsNodeData();

private:
	CbfsData* mCbfs;
	FILE* mJsonFile;
	bool mJsonDetail;
};
#endif // CPLEX_H
