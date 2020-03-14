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
#include <queue>
#include <utility>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <gmp.h>
#include <gmpxx.h>

using namespace std;

class CbfsNodeData;
typedef IloCplex::MIPCallbackI::NodeId NID;
typedef multimap<double, CbfsNodeData*> CbfsHeap;
typedef map<int, CbfsHeap> ContourMap;
typedef map<int, pair<int,int>> RangeMap;
typedef list<CbfsNodeData*> CbfsDive;
//struct opts;

class CbfsData
{
public:
	CbfsData(Mode m, ContSelMode cm, double posW, double nullW, int mob, int lbPara, int maxDepth, int probInterval, 
			int term, int numCont, double ucb) : 
		mMode(m), mPosW(posW), mNullW(nullW), bestLB(-INFINITY), bestUB(INFINITY), nIters(0),
		mCurrContour(mContours.begin()), mMob(mob), mLBContPara(lbPara),
		mDiveMaxDepth(maxDepth), mProbInterval(probInterval),
		earlyTermIter(term), mContSelMode(cm), mMaxNumConts(numCont), mBiasPara(ucb)
	{
		// WRA Contour Score Matrix Initialization
		//mContScores.resize(lbPara, 0);
		// Tie-breaking Rule Initialization
		mTieBreak = OG;
		// Dive Initialization
		mIsCustomDiveOn = (maxDepth > 0) ? true : false;
		mDiveStatus = (maxDepth > 0) ? true : false;
		// Check if use CPLEX Dive
		mIsCplexDiveOn = (maxDepth == -1) ? true : false;
		preId = prepreId = 0;

		mNumUnExpd = 0;
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
		mNumIntInfsblUnexplr = mSumIntInfsblUnexplr = 0;
		mEarlyTermOn = (earlyTermIter > 0) ? true : false;
		//geoMeanShift = 1;
		// Cont selection
		mPreDepth = 0;
		mMinNumContSel = 100;
		//mMaxNumConts = 8;
		//mBiasPara = 1.4;
		mIsContInitd = false;
		mIsRepopulate = false;
		switch (mContSelMode)
		{
		case Subtree:
			mContSelMulti = 1;
			break;
		case WBranch:
			mContSelMulti = 1;
			break;
		default:
			mContSelMulti = 1;
		}
		if (mMode == ContSel)
			turnOffDive();

		mRexSolCountOn = true;
		checkAgainst = 31.87;
	}
	CbfsData(opts* o);
	~CbfsData();

	void addNode(CbfsNodeData* nodeData);
	void delNode(CbfsNodeData* nodeData);
	void updateBounds(double lb, double ub, double gap);
	void updateExpdMaxDepth(int depth);
	void probStep();
	void terminateProb();
	void printScores();
	void printCriteria();
	int calContour(CbfsNodeData* nodeData);
	int getNumNodesInCont(int contID);

	NID getNextNode();
	CbfsNodeData* getNextNodePenalty();
	CbfsNodeData* getNextNodeBFS();
	ContourMap::iterator getNextCont();

	Mode getMode() { return mMode; }
	void setSense(IloObjective::Sense s) { mSense = s; }
	long getNumIterations() { return nIters; }
	double getBestUB() { return bestUB; }
	void setBestUB(double ub) { bestUB = ub; }
	void setNumIterations(long i) { nIters = i; }
	void updateContourBegin() { mCurrContour = mContours.begin(); }

	// FUNCTIONS: Early termination and criteria for best contour strategy
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

	// FUNCTIONS: Use CPLEX Dive
	bool isCplexDiveOn() { return mIsCplexDiveOn; }
	void setCplexDive(bool flag) { mIsCplexDiveOn = flag; }
	bool isInSameDive(int id);
	void updateDivePre(int id);

	// FUNCTIONS: Contour selection
	void updateCurContPayoff();
	void updateOneStep(double value);
	void updateContScores();
	void initContScores();
	void populateConts();
	void resetSubtrees();
	void setContScoresMP(int contID, mpf_class score) { mContScoresMP[contID] = score; }
	void setContScores(int contID, double score) { mContScores[contID] = score; }
	void resetCurPayoff() { mCurPayoff = 0; mCurNumNodesExplrd = 0; }
	void resetCurPayoffMP() { mCurPayoffMP = 0; mCurNumNodesExplrd = 0; }
	void turnOnDive();
	void turnOffDive() { mIsCustomDiveOn = false; mDiveStatus = false; mIsCplexDiveOn = false; }
	int getContNum() { return mContours.size(); }

	void updateNodeCount();
	int getNumLNode() { return mNumLNode; }
	int getNumENode() { return mNumENode; }
	int getNumGNode() { return mNumGNode; }

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
	int mNumUnExpd;
	// VARIABLES: Diving and probing
	bool mIsCustomDiveOn;
	int mNumDive, mBFSInterval;
	int mCurDiveCount, mDiveStartCont;										// Diving parameters
	int mDiveMinDepth, mDiveMaxDepth;										
	double mDiveTol;
	int mProbStep, mProbInterval, mPreProbEnd;								// Probing parameters
	bool mDiveStatus, mProbStatus;											// Diving and Probing flags
	
	// Dive contour, deviating path penalty
	int mDiveContPara, prtIDLstInstedNode;
	int preNodeID, preContID;												// ID of the node explored just prior to current one
	double preNodeLB, penaltyPara;
	bool penaltyOn;
	// VARIABLES: Hybrid estimate and value
	double mHybridBestPara;
	// VARIABLES: Early termination and criteria for best contour strategy
	bool mEarlyTermOn;
	int earlyTermIter;
	double rootLB;
	//double geoMeanShift;
	int mExpdMaxDepth;														// Record the maximum depth in current run
	int mNumUnexplNodes;													// Criteria 1: Number of unexplored nodes
	double mLbImprv;														// Criteria 2: LB percentage improvement from root LB
	int mNumIntInfsblUnexplr;												// Criteria 3: Sum of number of integer-infeasible variables of unexplored ndoes
	int mSumDpthUnexplr;													// Criteria 4: Sum of depth of unexplored nodes
	double mSumIntInfsblUnexplr;											// Criteria 5: Sum of integer infeasibilites of unexplored ndoes
	int mSumIterSplx;														// Criteria 8: Sum of simplex iterations
	int mNVars, mNIntVars;													// Num of Decision Variables
	timespec mRoot;
	// VARIABLES: Use CPLEX Dive
	bool mIsCplexDiveOn;
	int preId, prepreId;
	// VARIABLES: ContourSelection
	ContSelMode mContSelMode;
	ContScoreMode mContScoreMode;
	int mContSelMulti;
	bool mIsContInitd;
	bool mIsRepopulate;
	int mPreDepth;
	int mMaxNumConts, mNumContVisits, mCurNumNodesExplrd;
	int mMinNumContSel;
	double mBiasPara, mCurPayoff;
	map<int, int> mContIndexMap;
	vector<double> mContScores;
	vector<double> mContPayoffs;
	vector<int> mContVisits;
	vector<mpf_class> mContScoresMP;
	vector<mpf_class> mContPayoffsMP;
	mpf_class mCurPayoffMP;

	bool mRexSolCountOn;
	map<double, int> mRexSolCount;
	int mNumLNode, mNumENode, mNumGNode;
	double checkAgainst;

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
	int lastBranch;			// 0 if branching down (left), 1 if branching up (right)
	int contour;
	NID id;
	double estimate, lpval, incumbent;
	RangeMap rmap;
	int parent, parentCont;
	int nodeContID, parentNodeContID;
	int subtree, parentSubtree;
	int weiContour;
	bool explored;

	int infNum;
	double infSum;

	CbfsNodeData(CbfsData* cbfs, FILE* jsonFile, bool jsonDetail) : 
		explored(false), mCbfs(cbfs), mJsonFile(jsonFile), mJsonDetail(jsonDetail) {}
	~CbfsNodeData();
	void calSubtree();

private:
	CbfsData* mCbfs;
	FILE* mJsonFile;
	bool mJsonDetail;
};

#endif // CPLEX_H
