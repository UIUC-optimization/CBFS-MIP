// cplex.cpp: David R. Morrison, Aug. 2013
// Implementation of the CPLEX wrapper for CBFS on MIPs

#include "cplex.h"

// Import the model and set up the CPLEX solver
Cplex::Cplex(const char* filename, FILE* jsonFile, CbfsData* cbfs, int timelimit, 
		bool disableAdvStart) :
	mModel(mEnv),    //Initialize model
	mCplex(mEnv),    //Initialize algorithm(These two require enverioment instance mEnv to initialize)
	mJsonFile(jsonFile)
{
	// Import the model
	mCplex.importModel(mModel, filename);   //Import model from file "filename"
	mCplex.extract(mModel);                 //Extract
	mObj = mCplex.getObjective();   
	cbfs->setSense(mObj.getSense()); //To specify whether the invoking objective is a maximization (Maximize) or minimization

	// Set CPLEX parameters
//	mCplex.setParam(IloCplex::MIPSearch, 1);		// Disable dynamic search for non-callback version
	mCplex.setParam(IloCplex::Threads, 1);			// Limit to 1 thread for all versions
	mCplex.setParam(IloCplex::TiLim, timelimit);	// 1 hour time limit
	mCplex.setParam(IloCplex::NodeFileInd, 0);		// Callbacks won't work with compressed node files
//	mCplex.setParam(IloCplex::NodeSel, 2);          // A best estimate node selection strategy
//	mCplex.setParam(IloCplex::BBInterval, 0);       // Never select best bound when using best estimate strategy
//	mCplex.setParam(IloCplex::VarSel, -1);			// Use the ''maximum infeasibility'' rule for variable selection
	mCplex.setParam(IloCplex::RandomSeed, 20150624);// Set random seed
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
		mCplex.use(new (mEnv) CbfsBranchCallback(mEnv, cbfs, mJsonFile));  //IloCplex::BranchCallbackI::CbfsBranchCallback
		mCplex.use(new (mEnv) CbfsNodeCallback(mEnv, cbfs));               //IloCplex::NodeCallbackI::CbfsNodeCallback
	}
}

// Use CPLEX to solve the problem
void Cplex::solve() 
{ 
	// Write output about solution procedure to JSON file
	if (mJsonFile) 
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
		    fprintf(mJsonFile, "{\"cpl_status\": %d, \"obj_value\": %0.2f, \"opt_gap\": %0.2f, \"nodes_to_best\": %d, \"total_nodes\": %ld, \"end_time\": %0.2f}\n", 
			status, mCplex.getObjValue(), mCplex.getMIPRelativeGap(), mCplex.getIncumbentNode(), mCplex.getNnodes(), elapsed);
		else
			fprintf(mJsonFile, "{\"cpl_status\": %d, \"total_nodes\": %ld, \"end_time\": %0.2f}\n",
			status, mCplex.getNnodes(), elapsed);
	}

}

// Compute the contour that a new node belongs to, and add it to the data structure
void CbfsData::addNode(CbfsNodeData* nodeData)
{
	//Store diving candidate if diving is on
	if (mDiveStatus)
		mDiveCand.push_back(nodeData);

	// Compute the contour for the node
	switch (mMode)
	{
		case Weighted:
			nodeData->contour = mPosW * nodeData->numPos+ mNullW * nodeData->numNull;
			break;
		case LBContour:
			//Keep in mind, this is not the final version, yet. One improvement is to prune all contours with keys greater than 1
			//if we have an incumbent solution, and not let new subproblem like that be added to contours.
			//But this procedure can be done during branching callback and not here.
			nodeData->contour = calContour(nodeData->lpval);
			break;
	}

	// TODO - should we use estimate or lower bound here?  I think estimate's the right thing, but
	// I'm not 100% sure of this.  We should try the other as well
	// Here we could try to use lower bound as a measure of best 
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
	}
	// The contour heap now store pointers to the nodeData instances, instead of NID. Below, changes are made accordingly.
	mContours[nodeData->contour].insert({best, nodeData});
}

void CbfsData::delNode(CbfsNodeData* nodeData)
{
	//If mDiveCand has elements, then first see if node to be deleted is in here.
	if (!mDiveCand.empty())
	{
		if (nodeData->id == mDiveCand.front()->id) mDiveCand.pop_front();
		else if (nodeData->id == mDiveCand.back()->id) mDiveCand.pop_back();
		if (mDiveCand.size() > 2) throw ERROR << "More than 2 elements in mDiveCand.";
	}
	// Look up the node by the measure of best key
	double search;
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
	}
	auto range = mContours[nodeData->contour].equal_range(search);

	// Look at all nodes in the current contour with this measure of best, and delete
	// the first one we find with the same id.  We shouldn't have multiple nodes with 
	// the same id, so this should be fine.  Maybe at some point it's worth making this
	// a little cleaner.
	for (auto i = range.first; i != range.second; ++i)
	{
		// pointers to nodeData are stored in contour heap
		if ((i->second)->id == nodeData->id) 
		{ 
			mContours[nodeData->contour].erase(i); 
			break; 
		}
	}

	// If the contour is empty, remove it from the map
	if (mContours[nodeData->contour].empty())
	{
		// If the current contour pointer points to the thing we're deleting, increment it
		if (mCurrContour == mContours.find(nodeData->contour))
		{
			if (mCurrContour == mContours.begin())
				mCurrContour = mContours.end();
			--mCurrContour;
		}
		mContours.erase(nodeData->contour);
	}
}

// Select the next node to explore from the CBFS data structure
NID CbfsData::getNextNode()
{
	NID id; id._id = -1;
	if (!mDiveCand.empty())
	{
		id = mDiveCand.front()->id;
		mDiveCand.clear();
		mDiveCount++;
		if (mDiveCount >= 500) mDiveStatus = false;     //If maximal depth of a dive is reached, then stop.
		return id;
	}
	else
	{
		mDiveCount = 0;
		mDiveStatus = true;
	}
	// Increment the iterator into the heap structure
	++mCurrContour;
	if (mCurrContour == mContours.end())
		mCurrContour = mContours.begin();

	// If the heap has a non-empty contour, return the best thing in that contour
	// pointers to nodeData are stored in contour heap
	if (mCurrContour != mContours.end())
		id = (mCurrContour->second.begin()->second)->id;

	// There should always be a next node when we call this
	if (id._id == -1) 
	{
		if (mContours.empty())
			printf("No contours exist");
		for (auto i = mContours.begin(); i != mContours.end(); ++i)
			printf("Contour %d has size %ld\n", i->first, i->second.size());
		throw ERROR << "Node ID should not be -1!";
	}

	return id;
}

// Do branching and store the generated nodes
void CbfsBranchCallback::main()
{
	mCbfs->setNumIterations(getNnodes());
	// Output

	// If there aren't any branches, just prune and return
	if (getNbranches() == 0)
	{
		prune();
		return;
	}

	// myNodeData will be valid except at the root of the search
	CbfsNodeData* myNodeData = (CbfsNodeData*)getNodeData();

	//Set starting time of estimation period
	//if (!myNodeData)
	//	mCbfs->setTime(getCplexTime());

	// Update lower bound and upper bound, and update contour if necessary
	double incumVal = hasIncumbent() ? getIncumbentObjValue() : INFINITY;
	double bestObj = getBestObjValue();
	double optGap = getMIPRelativeGap();
	mCbfs->updateBounds(bestObj, incumVal, optGap);

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

		// Populate the CbfsNodeData structure corresponding to the new branch
		CbfsNodeData* newNodeData;
		if (myNodeData) newNodeData = new CbfsNodeData(*myNodeData);
		else
		{ 
			newNodeData = new CbfsNodeData(mCbfs, mJsonFile); 
			newNodeData->depth = 0; newNodeData->numPos = 0; newNodeData->numNull = 0;
			newNodeData->contour = 0; newNodeData->parent = 0;
		}
		newNodeData->estimate = estimate;
		newNodeData->lpval = getObjValue();
		newNodeData->incumbent = hasIncumbent() ? getIncumbentObjValue() : INFINITY;
		newNodeData->parent = getNodeId()._id;
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

		//Update the tree profile
		//mCbfs->updateProfile(newNodeData->numPos, newNodeData->numNull);

		// Tell CPLEX to create the branch and return the id of the new node
		newNodeData->id = makeBranch(vars, bounds, dirs, estimate, newNodeData);

		// Insert the newly-created node into our CBFS data structure
		if (mCbfs->getMode() != CplexOnly)
			mCbfs->addNode(newNodeData);
		else newNodeData->contour = -1;

		// Output
	}

	//Check if ready for estimation
	//mCbfs->checkTermi(getCplexTime());

	// Mark this node as explored so we don't print output when it gets deleted
	if (myNodeData) myNodeData->explored = true;
}

// Select a new node to branch on according to the CBFS rule
void CbfsNodeCallback::main()
{
	NID nodeId;
	if (mCbfs->getMode() != CplexOnly)
		nodeId = mCbfs->getNextNode();
	else nodeId = getNodeId(0);

	selectNode(nodeId);
}

// Clean up a node from the CBFS data structure once it's deleted from the tree
CbfsNodeData::~CbfsNodeData()
{
	if (mCbfs->getMode() != CplexOnly)
		mCbfs->delNode(this);

	// If the node hasn't been explored yet, mark it as pruned in the JSON file
}

// Update lower and upper bound, as well as contPara. The first time contPara changes,
// all contour will be updated to accomodate new contPara.
void CbfsData::updateBounds(double lb, double ub, double gap)
{
	bestLB = (bestLB < lb) ? lb : bestLB;
	bestUB = (bestUB > ub) ? ub : bestUB;
	//if (ub != INFINITY)
	//{
	//	if (gap <= 0.01 && contPara != 2)
	//	{
	//		contPara = 2;
	//		updateContour();
	//	}
	//	if (gap <= 0.05 && gap > 0.01 && contPara != 5)
	//	{
	//		contPara = 5;
	//		updateContour();
	//	}
	//	if (gap <= 0.1 && gap > 0.05 && contPara != 10)
	//	{
	//		contPara = 10;
	//		updateContour();
	//	}
	//	if (gap <= 0.2 && gap > 0.1 && contPara != 25)
	//	{
	//		contPara = 25;
	//		updateContour();
	//	}
	//}
}

// Contour update function. *This is not very efficient, but we should be fine for the moment because
// it will only be called at most 4 times.
void CbfsData::updateContour()
{
	ContourMap tempContours;
	int contour;

	//Go through all nodes in existing contours and determine their new contour
	for (ContourMap::iterator i = mContours.begin(); i != mContours.end(); i++)
	{
		for (CbfsHeap::iterator j = ((*i).second).begin(); j != ((*i).second).end(); j++)
		{
			contour = calContour((*j).second->lpval);
			((*j).second)->contour = contour;
			tempContours[contour].insert({ (*j).first, (*j).second });
		}
	}
	mContours.clear();
	mContours = tempContours;
	mCurrContour = mContours.begin();
}

int CbfsData::calContour(double lb)
{
	int contour;
	if (bestUB == INFINITY)
	{
		contour = (bestLB != 0) ? int(floor(fabs((lb - bestLB) / bestLB) * mcontPara)) : int(floor(fabs(lb)));
	}
	else
	{
		contour = int(floor(fabs((lb - bestLB) / (bestUB - bestLB) * mcontPara)));
	}
	return contour;
}

//Profile of the partial tree produced by the bnb algorithm, record number of nodes in each level
void CbfsData::updateProfile(int pb, int nb)
{
	int level = pb + nb;
	if (treeProfile.size() < level)
		treeProfile.resize(level);
	treeProfile[level - 1]++;
}

//Produce the estimate tree size with existing partial tree
void CbfsData::estimateTree()
{
	int dT, bT;
	int lT = 0;
	int waistMin, waistMax;
	int waistVal = -1;

	dT = treeProfile.size() + 1;
	for (int ind = 0; ind < dT - 2; ind++)
	{
		if (lT == 0 && (treeProfile[ind + 1] / treeProfile[ind]) < 2)
			lT = ind + 1;
		if (waistVal < treeProfile[ind])
		{
			waistMin = ind;
			waistMax = ind;
			waistVal = treeProfile[ind];
		}
		else if (waistVal == treeProfile[ind])
		{
			waistMax = ind;
		}
	}
	bT = ceil((waistMin + waistMax) / 2.0);

	vector<double> lambdaSeq;
	for (int ind = 0; ind <= dT; ind++)
	{
		if (ind <= lT - 1)
			lambdaSeq.push_back(2);
		else if (ind >= lT && ind <= bT - 1)
			lambdaSeq.push_back(2 - (ind - lT + 1) / (bT - lT + 1));
		else
			lambdaSeq.push_back(1 - (ind - bT + 1) / (dT - bT + 1));
	}

	double estSum = 1;
	for (int i = 1; i <= dT; i++)
	{
		double estMulti = 1;
		for (int j = 0; j < i; j++)
			estMulti = estMulti * lambdaSeq[j];
		estSum += estMulti;
	}

	printf("Estimate tree size: %0.2f\n", estSum);
}

void CbfsData::checkTermi(double time)
{
	if (time - startTime > 5 && count != 1)
	{
		int treeSize = 1;
		for (int ind = 0; ind < treeProfile.size(); ind++)
			treeSize += treeProfile[ind];
		if (treeSize / treeProfile.size() >= 20)
		{
			estimateTree();
			count = 1;
		}
	}
}
