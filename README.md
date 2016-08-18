# CBFS-MIP
Search strategy research for MIP
Current features:
1. CPLEX default mode
2. CPLEX with dummy callbacks for comparison with CBFS
3. CBFS using Weighted Branch contours
4. CBFS using Lower Bound contours
5. CBFS using number of infeasible variable as contours
6. CBFS using random contours (randomly assign nodes to contours)
Experimental features:
1. Probing steps
2. Disable diving when optimality gap is small (disabled)
3. Update the contour of every existing nodes if big change in optimality gap (disabled)
4. Infeasible variable contour can be computed in many ways:
```
contour = double(mNIntVars - nodeData->numInfeasibles) / (double)(mNIntVars) * 10;
```
or
```
contour = mNIntVars - nodeData->numInfeasibles;
double def = double(mNIntVars - nodeData->numInfeasibles) / (double)(mNIntVars);
if (def <= 0.5) 
    contour = 0;
else {
	def = (def - 0.5) / 0.5;
	contour = floor(def * mNInfeasibleCont + 1);
}
```