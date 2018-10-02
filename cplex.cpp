// cplex.cpp: David R. Morrison, Aug. 2013
//            Wenda Zhang, 2015 - 
// Implementation of the CPLEX wrapper for CBFS on MIPs

#include "cplex.h"

// Import the model and set up the CPLEX solver
Cplex::Cplex(const char* filename, FILE* jsonFile, CbfsData* cbfs, int timelimit, 
		bool disableAdvStart, int rand, bool jsonDetail) :
	mModel(mEnv),    //Initialize model
	mCplex(mEnv),    //Initialize algorithm(These two require enverioment instance mEnv to initialize)
	mJsonFile(jsonFile),
	mJsonDetail(jsonDetail)
{
	// Import the model
	mCplex.importModel(mModel, filename);   //Import model from file "filename"
	mCplex.extract(mModel);                 //Extract
	mObj = mCplex.getObjective();   
	cbfs->setSense(mObj.getSense()); //To specify whether the invoking objective is a maximization (Maximize) or minimization

	// Set CPLEX parameters
	mCplex.setParam(IloCplex::MIPSearch, 1);		// Disable dynamic search for non-callback version
	mCplex.setParam(IloCplex::Threads, 1);			// Limit to 1 thread for all versions
	mCplex.setParam(IloCplex::TiLim, timelimit);	// 1 hour time limit
	mCplex.setParam(IloCplex::NodeFileInd, 0);		// Callbacks won't work with compressed node files
//	mCplex.setParam(IloCplex::NodeSel, 2);          // A best estimate node selection strategy
//	mCplex.setParam(IloCplex::BBInterval, 0);       // Never select best bound when using best estimate strategy
//	mCplex.setParam(IloCplex::VarSel, -1);			// Use the ''maximum infeasibility'' rule for variable selection
	if (rand > 1000)
	    mCplex.setParam(IloCplex::RandomSeed, rand);// Set random seed
	if (disableAdvStart)
	{
		printf("Disabling advanced start methods and cut generation\n");
		mCplex.setParam(IloCplex::AdvInd, 0);		// Disable use of warm-start info
		mCplex.setParam(IloCplex::EachCutLim, 0); 	// Disable cut generation
		mCplex.setParam(IloCplex::FracCuts, -1); 	// Do not generate Gomory cuts
	}

	// Initialize callbacks for branching and node selection
	if (cbfs->getMode() != Disable)
	{
		mCplex.use(new (mEnv) CbfsBranchCallback(mEnv, cbfs, mJsonFile, mJsonDetail));  //IloCplex::BranchCallbackI::CbfsBranchCallback
		mCplex.use(new (mEnv) CbfsNodeCallback(mEnv, cbfs));               //IloCplex::NodeCallbackI::CbfsNodeCallback
	}

	// Storing decision variables info in cbfs
	//cbfs->setVarArray(getVars());
	//cbfs->setVarCount(varsCount);
	// Storing integer variables info in cbfs
	//cbfs->setIntVarArray(getIntVars());
	//cbfs->setIntVarCount(intVarsCount);
	//cbfs->setDiveMaxDepth(countVars);
	cbfs->setJsonFile(mJsonFile);
	mCbfs = cbfs;
}

// Use CPLEX to solve the problem
void Cplex::solve() 
{ 
	// Write output about solution procedure to JSON file
	if (mJsonFile && mJsonDetail) 
	{
		fprintf(mJsonFile, "{\"start_time\": %ld, \"clocks_per_sec\": %ld, "
			"\"num_vars\": %ld, \"num_int_vars\": %ld, \"num_bin_vars\": %ld}\n", 
			clock(), CLOCKS_PER_SEC, mCplex.getNcols(), mCplex.getNintVars(), mCplex.getNbinVars());
	}

	printf("Starting solve!\n\n");

	timespec start, end;
	clock_gettime(CLOCK_MONOTONIC, &start);
	mCplex.solve(); 
	clock_gettime(CLOCK_MONOTONIC, &end);
	double elapsed = (end.tv_sec - start.tv_sec);
	elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
	double elpToRoot = (mCbfs->getRootTime()).tv_sec - start.tv_sec;
	elpToRoot += ((mCbfs->getRootTime()).tv_nsec - start.tv_nsec) / 1000000000.0;

	printf("DATA_START\n{\n");

	// Determine the solution found by CPLEX
	auto status = mCplex.getCplexStatus();
	printf("\t\"cplex_status\" : %d,\n", status);
	printf("\t\"lower_bound\" : "); try { printf("%0.2f", mCplex.getBestObjValue()); } 
	catch (const IloCplex::Exception& e) { printf("\"inf\""); }
	printf(",\n");
	printf("\t\"incumbent\" : "); try { printf("%0.2f", mCplex.getObjValue()); } 
	catch (const IloCplex::Exception& e) { printf("\"inf\""); }
	printf(",\n");
	printf("\t\"gap\" : "); try { printf("%0.2f", mCplex.getMIPRelativeGap()); } 
	catch (const IloCplex::Exception& e) { printf("\"inf\""); }
	printf(",\n");
	printf("\t\"nodes_to_best\" : "); try { printf("%d", mCplex.getIncumbentNode()); } 
	catch (const IloCplex::Exception& e) { printf("-1"); }
	printf(",\n");
	printf("\t\"total_nodes\" : %d,\n", mCplex.getNnodes());
	printf("\t\"total_time\" : %0.2f\n", elapsed);
	printf("}\nDATA_END\n");
	switch (status)
	{
		case IloCplex::Optimal:
			printf("Solution was optimal\n");
			break;
		case IloCplex::OptimalTol:
			printf("Solution was optimal within tolerance\n");
			break;
		case IloCplex::Infeasible:
			printf("Problem was infeasible\n");
			break;
		case IloCplex::AbortTimeLim:
			printf("Reached time limit\n");
			break;
		default:
			printf("Unknown CPLEX status: %d\n", mCplex.getCplexStatus());
	}
	printf("Best lower bound: %0.2f\n", mCplex.getBestObjValue());
	printf("Value of incumbent: "); try { printf("%0.2f", mCplex.getObjValue()); }
	catch (const IloCplex::Exception& e) { printf("\"inf\""); }
	printf(",\n");
	printf("Solution gap: "); try { printf("%0.2f", mCplex.getMIPRelativeGap()); }
	catch (const IloCplex::Exception& e) { printf("\"inf\""); }
	printf(",\n");
	printf("Processed nodes: %ld\n", mCplex.getNnodes());

	// Write output about solution status to JSON file
	if (mJsonFile) 
	{
		if (status != IloCplex::Infeasible)
		{
			fprintf(mJsonFile, "{\"cpl_status\": %d, ", status);
			fprintf(mJsonFile, "\"obj_value\": "); try { fprintf(mJsonFile, "%0.2f, \"opt_gap\": %0.2f, \"nodes_to_best\": %d, ", 
																mCplex.getObjValue(), mCplex.getMIPRelativeGap(), mCplex.getIncumbentNode()); }
			catch (const IloCplex::Exception& e) { fprintf(mJsonFile, "inf, "); }
			fprintf(mJsonFile, "\"total_nodes\": %ld, \"end_time\": %0.2f, \"num_vars\": %d, ", mCplex.getNnodes(), elapsed, mCplex.getNintVars());
			fprintf(mJsonFile, "\"root_time\": %0.2f}\n", elpToRoot);
		}
		else
			fprintf(mJsonFile, "{\"cpl_status\": %d, \"total_nodes\": %ld, \"end_time\": %0.2f, \"num_vars\": %d, \"root_time\": %0.2f}\n",
			status, mCplex.getNnodes(), elapsed, mCplex.getNintVars(), elpToRoot);
	}

}

// Attempt to extract all variables from CPLEX model
IloNumVarArray Cplex::getVars() 
{
	IloModel::Iterator iter(mModel);
	unordered_set<int> vars;
	IloNumVarArray allVars(mEnv);
	int count = 0;
	IloExpr expr;
	// Extract all variables and store all integer variables
	while (iter.ok()) 
	{
		IloExtractable extr = *iter;
		if (extr.isVariable()) 
		{
			IloNumVar temp = extr.asVariable();
			vars.insert(temp.getId());
			if (vars.size() > count) 
			{
				allVars.add(temp);
				++count;
			}
			//printf("Variable detected. %d\n", count);
		}
		++iter;
	}
	//printf("Count %d variables in total.\n", count);
	varsCount = count;
	return allVars;
}

// Attempt to extract integer variables from CPLEX model
IloNumVarArray Cplex::getIntVars()
{
	IloModel::Iterator iter(mModel);
	unordered_set<int> vars;
	IloNumVarArray allIntVars(mEnv);
	int count = 0;
	IloExpr expr;
	// Extract all variables and store all integer variables
	while (iter.ok())
	{
		IloExtractable extr = *iter;
		if (extr.isVariable())
		{
			IloNumVar temp = extr.asVariable();
			if (temp.getType() != 2)
			{
				vars.insert(temp.getId());
				if (vars.size() > count)
				{
					allIntVars.add(temp);
					++count;
				}
			}
			//printf("Variable detected. %d\n", count);
		}
		++iter;
	}
	//printf("Count %d variables in total.\n", count);
	intVarsCount = count;
	return allIntVars;
}

// Compute the contour that a new node belongs to, and add it to the data structure
void CbfsData::addNode(CbfsNodeData* nodeData)
{

	// Compute the contour for the node
	nodeData->contour = calContour(nodeData);

	// Store diving candidate if diving is on
	if (mDiveStatus) 
		mDiveCand.push_back(nodeData);

	// Store the recently added nodes
	if (mEarlyTermOn)
		rcntAddedNodes.push_back(nodeData);

	// 1,2: Best Estimate/LB; 3,4: Worst Estimate/LB; 5: hybrid between Best Estimate and LB
	double best;
	switch (mMob)
	{
		case 1: 
			best = (mSense == IloObjective::Maximize) ? nodeData->estimate : -nodeData->estimate;
			break;
		case 2:
			best = (mSense == IloObjective::Maximize) ? nodeData->lpval : -nodeData->lpval;
			break;
		case 3:
			best = (mSense == IloObjective::Maximize) ? -nodeData->estimate : nodeData->estimate;
			break;
		case 4:
			best = (mSense == IloObjective::Maximize) ? -nodeData->lpval : nodeData->lpval;
			break;
		case 5:
			best = nodeData->estimate * mHybridBestPara + nodeData->lpval * (1 - mHybridBestPara);
			best = (mSense == IloObjective::Maximize) ? -best : best;
			break;
		default:
			best = (mSense == IloObjective::Maximize) ? -nodeData->estimate : nodeData->estimate;
			break;
	}
	// The contour heap now store pointers to the nodeData instances, instead of NID.
	mContours[nodeData->contour].insert({best, nodeData});
}

void CbfsData::delNode(CbfsNodeData* nodeData)
{
	// If mDiveCand has elements, see if node to be deleted is there
	if (!mDiveCand.empty())
	{
		if (mDiveCand.size() > 2) throw ERROR << "More than 2 elements in mDiveCand.";
		if (nodeData->id == mDiveCand.front()->id) mDiveCand.pop_front();
		else if (nodeData->id == mDiveCand.back()->id) mDiveCand.pop_back();
	}
	// If during probing, then check stored probing nodes
	if (mProbStep == 1 && !mProbCand.empty())
	{
		if (mProbCand.front()->id == nodeData->id) mProbCand.clear();
	}
	else if (mProbStep == 2 && !mProbOneSide.empty())
	{
		if (mProbOneSide.front()->id == nodeData->id) mProbOneSide.pop_front();
		else if (mProbOneSide.back()->id == nodeData->id) mProbOneSide.pop_back();
	}
	// If early termination info recording is on, then update infeasibility info
	if (mEarlyTermOn)
	{
		// If rcntAddedNodes has elements, see if node to be deleted is there
		if (!rcntAddedNodes.empty())
		{
			if (rcntAddedNodes.size() > 2) throw ERROR << "More than 2 elements in rcntAddedNodes.";
			if (nodeData->id == rcntAddedNodes.front()->id) rcntAddedNodes.pop_front();
			else if (nodeData->id == rcntAddedNodes.back()->id) rcntAddedNodes.pop_back();
		}
	}
	// Look up the node by the measure of best key
	double search;
	// 1,2: Best Estimate/LB; 3,4: Worst Estimate/LB; 5: hybrid between Best Estimate and LB
	switch (mMob)
	{
		case 1: 
			search = (mSense == IloObjective::Maximize) ? nodeData->estimate : -nodeData->estimate;
			break;
		case 2:
			search = (mSense == IloObjective::Maximize) ? nodeData->lpval : -nodeData->lpval;
			break;
		case 3:
			search = (mSense == IloObjective::Maximize) ? -nodeData->estimate : nodeData->estimate;
			break;
		case 4:
			search = (mSense == IloObjective::Maximize) ? -nodeData->lpval : nodeData->lpval;
			break;
		case 5:
			search = nodeData->estimate * mHybridBestPara + nodeData->lpval * (1 - mHybridBestPara);
			search = (mSense == IloObjective::Maximize) ? -search : search;
			break;
		default:
			search = (mSense == IloObjective::Maximize) ? -nodeData->estimate : nodeData->estimate;
			break;
	}
	auto range = mContours[nodeData->contour].equal_range(search);

	// Look at all nodes in the current contour with this measure of best, and delete
	// the first one we find with the same id.  We shouldn't have multiple nodes with 
	// the same id, so this should be fine.  Maybe at some point it's worth making this
	// a little cleaner.
	for (auto i = range.first; i != range.second; ++i)
	{
		// pointers to nodeData are stored in contour heap
		if (i->second->id == nodeData->id) 
		{ 
			mContours[nodeData->contour].erase(i); 
			break;
		}
	}

	// If the contour is empty, remove it from the map
	if (mContours[nodeData->contour].empty())
	{
		/*************************************/
		// When contour becomes empty, clear corresponding score
		//TODO: Use a separate function for contour score updates
		updateContScores(nodeData->contour, 0);
		/*************************************/
		// If the current iterator points to the contour we're deleting, increment it
		if (mCurrContour == mContours.find(nodeData->contour))
		{
			if (mCurrContour == mContours.begin())
				mCurrContour = mContours.end();
			--mCurrContour;
		}
		mContours.erase(nodeData->contour);
	}

	// Update the max depth of processed nodes
	updateExpdMaxDepth(nodeData->depth);
}

// Select the next node to explore from the CBFS data structure
NID CbfsData::getNextNode()
{
	NID id; id._id = -1;

	if (!mDiveCand.empty())
	{
		// Perform probing either when specific condition satisfied or probing already started
		if (mProbStep != 0 || (mProbStatus && (mCurDiveCount == 0 
						   || (mCurDiveCount - mPreProbEnd) == mProbInterval)))
			probStep();
		if (mProbStep != 2) mCurDiveCount++;							// Avoid an extra depth update when probing

		CbfsNodeData* diveNode = mDiveCand.front();
		id = diveNode->id;
		mDiveCand.clear();

		if (mCurDiveCount >= mDiveMaxDepth) mDiveStatus = false;		// If maximal depth of a dive is reached, terminate current dive
		else if (mCurDiveCount >= mDiveMinDepth && bestUB != INFINITY)
		{
			// If current LB too close to global UB, terminate current dive
			if ((diveNode->lpval - bestLB) / (bestUB - bestLB) > mDiveTol) mDiveStatus = false;
			// If no change in relaxed solution, terminate diving with probability
			//int rnum = rand() % mDiveMaxDepth;
			//if (mDiveCand.front()->lpval == preNodeLB && mDiveCount > rnum)
			//	mDiveStatus = false;
			//preNodeLB = mDiveCand.front()->lpval;
		}
		return id;
	}
	// Hhandle empty successor list mid-probing
	else if (mProbStatus && mProbStep != 0)
	{
		if (mProbStep == 1) mDiveCand = mProbCand;		// If child nodes of first candidate are pruned
		else mDiveCand = mProbOneSide;
		if (mDiveCand.empty())
		{
			terminateProb();
			mCurDiveCount = 0; mPreProbEnd = 0;
		}
		else
		{
			id = mDiveCand.front()->id;
			mDiveCand.clear();
			terminateProb();
			mPreProbEnd = mCurDiveCount;
			mCurDiveCount++;
			return id;
		}
	}
	// Two possibilities for mDiveCand to be empty: 
	// 1. maximal depth reached in last iteration
	// 2. node in previous iteration pruned somehow
	// mMaxDepth > 0 means diving is on, then this block enables diving for the next contour.
	else if (mDiveMaxDepth > 0)
	{
		mCurDiveCount = 0; mPreProbEnd = 0;
		mDiveStatus = true;
	}

	// Increment the iterator into the heap structure
	//++mCurrContour;
	mCurrContour = getNextCont();
	if (mCurrContour == mContours.end())
		mCurrContour = mContours.begin();

	// Get one node from contour mCurrContour points to
	if (mCurrContour != mContours.end()) 
	{
		// Tie breaking rules
		double search = mCurrContour->second.begin()->first;
		auto range = mCurrContour->second.equal_range(search);
		switch (mTieBreak)
		{
		case FIFO:
			id = range.first->second->id;
			break;
		case LIFO:
			id = (--range.second)->second->id;
			break;
		case OG:
			id = mCurrContour->second.begin()->second->id;
			break;
		default:
			id = mCurrContour->second.begin()->second->id;
		}

		// EXPERIMENT: Penalized deviating from path
		if (penaltyOn && mCurrContour != mContours.begin())
			id = getNextNodePenalty();
	}

	// There should always be a next node when we call this
	if (id._id == -1) 
	{
		if (mContours.empty())
			printf("No contours exist");
		for (auto i = mContours.begin(); i != mContours.end(); ++i)
			printf("Contour %d has size %ld\n", i->first, i->second.size());
		throw ERROR << "Node ID should not be -1!";
	}

	if (mDiveStatus) 
	{
		mDiveStartCont = mCurrContour->first;		// record contour number of the dive start node
		preNodeLB = -INFINITY;
		mNumDive++;
		if (mNumDive % mBFSInterval == 0)
			id = getNextNodeBFS();
	}
	preNodeID = id._id;
	preContID = mCurrContour->first;

	return id;
}

NID CbfsData::getNextNodeBFS()
{
	double minLB = mCurrContour->second.begin()->second->lpval;
	NID minID = mCurrContour->second.begin()->second->id;
	for (auto iter = mCurrContour->second.begin(); iter != mCurrContour->second.end(); iter++)
	{
		if (iter->second->lpval < minLB) 
		{
			minLB = iter->second->lpval;
			minID = iter->second->id;
		}
	}
	return minID;
}

// EXPERIMENT: Penalized deviating from path
// Currently only with minimizing and best-estimate search
NID CbfsData::getNextNodePenalty()
{
	double minPSearch = INFINITY;
	double curPSearch;
	NID minID;
	// Define penalty for deviating
	double penalty = penaltyPara;
	if (bestUB != INFINITY)
		penalty = mReOptGap * penaltyPara / 2;

	minID = mCurrContour->second.begin()->second->id;
	for (auto iter = mCurrContour->second.begin(); iter != mCurrContour->second.end(); iter++)
	{
		if (preNodeID == iter->second->parent)
			curPSearch = iter->second->estimate;
		else
			curPSearch = iter->second->estimate + penalty * fabs(iter->second->estimate);
		if (curPSearch < minPSearch)
		{
			minPSearch = curPSearch;
			minID = iter->second->id;
		}
	}
	return minID;
}

// Select the next contour to be explored
ContourMap::iterator CbfsData::getNextCont()
{
	ContourMap::iterator iterToBe;
	int sumScore = 0;
	int randScore;
	switch (mMode) 
	{
	case LBContour:
		iterToBe = mContours.begin();
		for (int i = 0; i < mContScores.size(); i++)
		{
			sumScore += mContScores[i];
		}
		randScore = rand() % sumScore;
		sumScore = 0;
		for (int i = 0; i < mContScores.size(); i++, iterToBe++)
		{
			sumScore += mContScores[i];
			if (mContScores[i] != 0 && randScore < sumScore)
				break;
		}
		break;
	default:
		iterToBe = mCurrContour;
		iterToBe++;
	}
	return iterToBe;
}

void CbfsData::updateContScores()
{
	switch (mMode)
	{
	case LBContour:
		mContScores[mDiveStartCont]++;
		break;
	default:
		break;
	}
}

void CbfsData::updateContScores(int contID, int score)
{
	switch (mMode) 
	{
	case LBContour:
		mContScores[contID] = score;
		break;
	default:
		break;
	}
}

// Update lower and upper bound, optimality gap
void CbfsData::updateBounds(double lb, double ub, double gap)
{
	bestLB = (bestLB < lb) ? lb : bestLB;
	bestUB = (bestUB > ub) ? ub : bestUB;
	mReOptGap = gap;
}

void CbfsData::updateExpdMaxDepth(int depth)
{
	mExpdMaxDepth = (depth > mExpdMaxDepth) ? depth : mExpdMaxDepth;
}

int CbfsData::calContour(CbfsNodeData* nodeData)
{
	int contour, numInBaseCont;
	int lb = nodeData->lpval;
	switch (mMode) {
	/*----------------------------------------------------*/
	// Assign nodes by the weight of branch on their path from root
	case Weighted:
		contour = mPosW * nodeData->numPos + mNullW * nodeData->numNull;
		break;
	/*----------------------------------------------------*/
	// Assign nodes by their lower bound
	case LBContour:
		if (bestUB == INFINITY)
		{
			// The labeling function here could lead to extremely large contour numbers that might be a problem.
			contour = 0;
		}
		else
		{
			contour = int(floor(fabs((lb - bestLB) / (bestUB - bestLB) * mLBContPara)));
		}
		if (contour > (mLBContPara - 1))
			contour = mLBContPara - 1;
		/*************************************/
		// When new non-empty contour appears, update corresponding score
		//TODO: Use a separate function for contour score updates
		if (mContScores[contour] == 0)
			updateContScores(contour, 1);
		/*************************************/
		break;
	/*----------------------------------------------------*/
	// Assign randomly a contour to each node
	case RandCont:
		contour = rand() % 5;
		break;
	/*----------------------------------------------------*/
	// Assign nodes to the next contour or the first one
	case DiveComp:
		numInBaseCont = getNumNodesInCont(0);
		if (preContID < mDiveContPara && prtIDLstInstedNode != nodeData->parent)
		{
			// First of the two child nodes generated will go to the next contour
			// as long as it does not exceed the contour limit
			contour = preContID + 1;
			prtIDLstInstedNode = nodeData->parent;
		}
		else
		{
			// Second of the two child nodes may go to the first contour or the next
			if (preContID >= mDiveContPara / 2)
				contour = 0;
			else if (numInBaseCont <= mDiveContPara / 10)
				contour = 0;
			else
				contour = preContID + 1;
		}
		break;
	/*----------------------------------------------------*/
	// CPLEX default node selection
	case CplexOnly:
		contour = -1;
		break;
	/*----------------------------------------------------*/
	default:
		contour = 0;
	}
	return contour;
}

// Main procedure of probing step: store the two parent node, explore both and compare results
void CbfsData::probStep()
{
	if (mProbStep == 0 && mDiveCand.size() < 2) return;
	switch (mProbStep)
	{
	case 0:
		// Store prob candidate (the remaining branch)
		mProbCand = mDiveCand;
		mProbCand.pop_front();
		mProbStep = 1;
		break;
	case 1:
		if (mProbCand.empty())
		{
			terminateProb();
			mPreProbEnd = mCurDiveCount;
			break;
		}
		// Store child nodes from one probing candidate and direct search to the other candidate (binary branching)
		mProbOneSide = mDiveCand;
		mDiveCand = mProbCand;
		mProbStep = 2;
		break;
	case 2:
		if (mProbOneSide.empty()) {}
		// Check estimate of child nodes on both candidates, choose the better one
		else if (mProbOneSide.front()->estimate < mDiveCand.front()->estimate) mDiveCand = mProbOneSide;
		terminateProb();
		mPreProbEnd = mCurDiveCount;
		break;
	default:
		throw ERROR << "Invalid probing parameter.";
	}
}

void CbfsData::terminateProb()
{
	if (mProbStep != 0) 
	{
		mProbCand.clear();
		mProbOneSide.clear();
		mProbStep = 0;
	}
}

void CbfsData::printScores()
{
	if (mJsonFile)
	{
		if (mMode == LBContour)
		{
			fprintf(mJsonFile, "{ ");
			for (int i = 0; i < mContScores.size(); i++)
				fprintf(mJsonFile, "%d ", mContScores[i]);
			fprintf(mJsonFile, "}\n");
		}
	}
}

int CbfsData::getNumNodesInCont(int contID)
{
	auto iter = mContours.find(contID);
	if (iter == mContours.end())
		return 0;
	else
		return iter->second.size();
}

void CbfsData::calCriteria()
{
	ContourMap::iterator iterCont;
	double lbGap;
	mNumUnexplNodes = 0;
	mLbImprv = 1;
	mSumDpthUnexpl = 0;
	mNumIntInfsblUnexpl = mSumIntInfsblUnexpl = 0;
	for (iterCont = mContours.begin(); iterCont != mContours.end(); iterCont++)
	{
		for (auto iterNode = iterCont->second.begin(); iterNode != iterCont->second.end(); iterNode++)
		{
			mNumUnexplNodes++;
		}
	}
	for (iterCont = mContours.begin(); iterCont != mContours.end(); iterCont++)
	{
		for (auto iterNode = iterCont->second.begin(); iterNode != iterCont->second.end(); iterNode++)
		{
			// LB Percentage Improvement from Root LB
			//mLbImprv *= pow((iterNode->second->lpval - rootLB) / (fabs(rootLB) + 0.0001) * 100, 1 / mNumUnexplNodes);
			//mLbImprv *= (iterNode->second->lpval - rootLB) / (fabs(rootLB) + 0.0001) * 100;
			lbGap = (iterNode->second->lpval - rootLB) / (fabs(rootLB) + 0.0001) * 100;
			fprintf(mJsonFile, "%0.2f\n", lbGap);
			mLbImprv *= lbGap;
			//mLbImprv += (iterNode->second->lpval - rootLB) / (fabs(rootLB) + 0.0001) * 100;
			//printf("Current mLbImprv: %0.2f, LB: %0.2f, Root: %0.2f\n", mLbImprv, iterNode->second->lpval, rootLB);

			// Sum of number of integer-infeasible variables of unexplored ndoes
			mNumIntInfsblUnexpl += iterNode->second->infNum;
			mSumIntInfsblUnexpl += iterNode->second->infSum;

			// Unexplored Nodes Depth
			mSumDpthUnexpl += iterNode->second->depth;
		}
	}
	//mLbImprv = mLbImprv / mNumUnexplNodes;
	mLbImprv = pow(mLbImprv, 1 / mNumUnexplNodes);
}

void CbfsData::printCriteria()
{
	if (mJsonFile)
	{
		fprintf(mJsonFile, "{\"numUnexpNode\": %d, \"impvMean\": %0.2f, \"numIntInf\": %d, \"sumDepth\": %d, ", mNumUnexplNodes, mLbImprv, mNumIntInfsblUnexpl, mSumDpthUnexpl);
		fprintf(mJsonFile, "\"sumIntInf\": %0.2f, \"bestLB\": %0.2f, \"bestUB\": %0.2f, \"numIterSplx\": %d}\n", mSumIntInfsblUnexpl, bestLB, bestUB, mSumIterSplx);
	}
}

bool CbfsData::isInSameDive(int id)
{
	return ((id == preId) || (id == prepreId));
}

void CbfsData::updateDivePre(int id)
{
	prepreId = preId;
	preId = id;
}

// Do branching and store the generated nodes
void CbfsBranchCallback::main()
{
	// Output

	// If there aren't any branches, just prune and return
	if (getNbranches() == 0)
	{
		if (mJsonFile && mJsonDetail)
		{
			fprintf(mJsonFile, "{\"node_id\": %lld, \"is_int_feas\": ", getNodeId()._id);
			if (isIntegerFeasible()) fprintf(mJsonFile, "true, ");
			else fprintf(mJsonFile, "false, ");
			fprintf(mJsonFile, "\"upper_bound\": %0.2f}\n", getObjValue());
		}
		// If the solution is feasible, then check if this is the new incumbent solution
		if (isIntegerFeasible()) 
		{
			int ub = mCbfs->getBestUB();
			//int curSol = getObjValue();
			//if (ub != INFINITY && ub > curSol)
			if (ub != INFINITY)
				mCbfs->updateContScores();
			mCbfs->terminateProb();
			mCbfs->printScores();
		}
		prune();
		return;
	}

	// Update lower bound and upper bound, and update contour if necessary
	double incumVal = hasIncumbent() ? getIncumbentObjValue() : INFINITY;
	double bestObj = getBestObjValue();
	double optGap = getMIPRelativeGap();
	mCbfs->updateBounds(bestObj, incumVal, optGap);

	// myNodeData will be valid except at the root of the search
	CbfsNodeData* myNodeData = (CbfsNodeData*)getNodeData();
	if (!myNodeData)
	{
		mCbfs->updateContourBegin();
		mCbfs->setRootLB(getObjValue());
		mCbfs->setRootTime();
	}
	if (mJsonFile && mJsonDetail)
	{
		double incumbentVal = INFINITY;
		if (hasIncumbent()) incumbentVal = getIncumbentObjValue();
		fprintf(mJsonFile, "{\"node_id\": %lld, \"explored_at\": %ld, \"lower_bound\": %0.2f, "
			"\"upper_bound\": %0.2f, ", getNodeId()._id, getNnodes(), getObjValue(), incumbentVal);
		if (!myNodeData)
			fprintf(mJsonFile, "\"contour\": 0, ");
	}

	// Loop through all of the branches produced by CPLEX
	for (int i = 0; i < getNbranches(); ++i)
	{
		// Get information about the current branch
		IloNumVarArray vars(getEnv());
		IloNumArray bounds(getEnv());
		IloCplex::BranchDirectionArray dirs(getEnv());
		double estimate = getBranch(vars, bounds, dirs, i);

		if (vars.getSize() > 1) throw ERROR << "More than one variable branched on";
		else if (vars.getSize() == 0) throw ERROR << "No variables available to branch on";
		int varId = vars[0].getId();

		// Output
		if (i == 0 && mJsonFile && mJsonDetail)
			fprintf(mJsonFile, "\"branching_var\": \"%s\"}\n", vars[0].getName());

		// Populate the CbfsNodeData structure corresponding to the new branch
		CbfsNodeData* newNodeData;
		if (myNodeData) newNodeData = new CbfsNodeData(*myNodeData);
		else
		{ 
			newNodeData = new CbfsNodeData(mCbfs, mJsonFile, mJsonDetail); 
			newNodeData->depth = 0; newNodeData->numPos = 0; newNodeData->numNull = 0;
			newNodeData->contour = 0; newNodeData->parent = 0;
		}
		newNodeData->infNum = 0; newNodeData->infSum = 0;
		newNodeData->estimate = estimate;
		newNodeData->lpval = getObjValue();
		newNodeData->incumbent = hasIncumbent() ? getIncumbentObjValue() : INFINITY;
		newNodeData->parent = getNodeId()._id;
		//newNodeData->numInfeasibles = countInfeasibles;
		newNodeData->depth++;

		// Compute the values of the new bounds for the variable we're branching on
		auto& rmap = newNodeData->rmap;
		int oldLB = rmap.find(varId) != rmap.end() ? rmap[varId].first : vars[0].getLB();
		int oldUB = rmap.find(varId) != rmap.end() ? rmap[varId].second : vars[0].getUB();
		if (dirs[0] == IloCplex::BranchUp) 
		{
			// Branching up means that we're imposing a new lower bound on the variable
			newNodeData->numPos++;
			if (rmap.find(varId) != rmap.end()) rmap[varId].first = bounds[0];
			else rmap[varId] = make_pair(bounds[0], oldUB);
		}
		if (dirs[0] == IloCplex::BranchDown) 
		{
			// Branching down means that we're imposing a new upper bound on the variable
			newNodeData->numNull++;
			if (rmap.find(varId) != rmap.end()) rmap[varId].second = bounds[0];
			else rmap[varId] = make_pair(oldLB, bounds[0]);
		}

		// Tell CPLEX to create the branch and return the id of the new node
		newNodeData->id = makeBranch(vars, bounds, dirs, estimate, newNodeData);

		// Insert the newly-created node into our CBFS data structure
		//if (mCbfs->getMode() != CplexOnly)
		//	mCbfs->addNode(newNodeData);
		//else newNodeData->contour = -1;

		// For the sake of collecting info, addNode is used with CplexOnly as well
		mCbfs->addNode(newNodeData);

		// Output
		if (mJsonFile && mJsonDetail)
		{
			fprintf(mJsonFile, "{\"node_id\": %lld, \"child\": %lld, \"branch_dir\": %d, "
				"\"estimate\": %0.2f, \"contour\": %d}\n", getNodeId()._id, newNodeData->id._id,
				dirs[0], newNodeData->estimate, newNodeData->contour);
		}
	}

	// Mark this node as explored so we don't print output when it gets deleted
	if (myNodeData) myNodeData->explored = true;
}

// Select a new node to branch on according to the CBFS rule
void CbfsNodeCallback::main()
{
	NID nodeId;
	if (mCbfs->getMode() != CplexOnly)
	{
		if (mCbfs->isCplexDiveOn())
		{
			nodeId = getNodeId(0);
			CbfsNodeData* nodeData = (CbfsNodeData*)getNodeData(nodeId);
			if (!mCbfs->isInSameDive(nodeData->parent))
				nodeId = mCbfs->getNextNode();
			mCbfs->updateDivePre((nodeData->id)._id);
		} 
		else
			nodeId = mCbfs->getNextNode();
	}
	else nodeId = getNodeId(0);

	mCbfs->setNumIterations(getNnodes());
	if (mCbfs->isEarlyTermOn())
	{
		while (!mCbfs->rcntAddedNodes.empty())
		{
			CbfsNodeData* node = mCbfs->rcntAddedNodes.front();
			node->infNum = getNinfeasibilities(node->id);
			node->infSum = getInfeasibilitySum(node->id);
			mCbfs->rcntAddedNodes.pop_front();
		}
		// Remove the node to be explored from info
		//mCbfs->updateInfsblInfo(-getNinfeasibilities(nodeId), -getInfeasibilitySum(nodeId));
		if (mCbfs->isEarlyTermToStart())
		{
			mCbfs->setSplxIterInfo(getNiterations());
			mCbfs->calCriteria();
			mCbfs->printCriteria();
			abort();
			return;
		}
	}

	selectNode(nodeId);
}

// Clean up a node from the CBFS data structure once it's deleted from the tree
CbfsNodeData::~CbfsNodeData()
{
	//if (mCbfs->getMode() != CplexOnly)
		mCbfs->delNode(this);

	// If the node hasn't been explored yet, mark it as pruned in the JSON file
	if (mJsonFile && !explored && mJsonDetail)
	{
		fprintf(mJsonFile, "{\"node_id\": %lld, \"pruned_at\": %ld, \"lower_bound\": %0.2f, "
			"\"upper_bound\": %0.2f}\n", id._id, mCbfs->getNumIterations(), lpval, incumbent);
	}
}