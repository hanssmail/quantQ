.quantQ.lof.L2distance:{[point;dataVector]
    // point -- reference point
    // dataVector -- vector of data
    :sqrt t wsum t:flip[dataVector]-point;
 };

.quantQ.lof.getKNNProperties:{[k;data;point]
    // k -- the parameter of the kNN
    // data -- the data set with one column x, each observation is a list of features
    // point -- reference point
    // calculate distance to point and sort it
    data: `distanceToPoint xasc update distanceToPoint: .quantQ.lof.L2distance[point;x] from data;
    // distance to k-th
    distanceToK: data[`distanceToPoint][min[(k-1;count[data])]];
    // number of elements within the distance
    numberToK: count select from data where distanceToPoint<=distanceToK;
    // return values
    :(`distanceToK`numberToK)!(distanceToK;numberToK);
 };

.quantQ.lof.reachabilityDistance:{[X;indexY;data;k]
    // X -- point to find distance, can be out of data
    // indexY -- index of the point Y, it has to come from data
    // data -- the data set with one column x, each observation is a row of features
    // k -- the parameter of the kNN
    // extract point Y
    Y: first exec x from data where i=indexY;
    // component 1: distance between X and Y
    distanceXY: first .quantQ.lof.L2distance[X;enlist Y];
    // component 2: k-th distance of Y within dataset
    kDistanceY: .quantQ.lof.getKNNProperties[k;data;Y][`distanceToK];
    // return the reachability distance
    :max[(distanceXY;kDistanceY)];
    };

.quantQ.lof.localReachabilityDensity:{[X;data;k]
    // X -- point to find distance, can be out of data
    // data -- the data set with one column x, each observation is a row of features
    // k -- the parameter of the kNN
    // update distance to point to data
    data: update distanceToPoint: .quantQ.lof.L2distance[X;x] from data;
    // find k-neighbourhood of X
    neighbourhood: .quantQ.lof.getKNNProperties[k;data;X];
    // indices of k-neighbours -- contains neighbourhood[`numberToK] elements
    indicesB: exec i from data where distanceToPoint<=neighbourhood[`distanceToK];
    // return density
    :1.0%(sum .quantQ.lof.reachabilityDistance[X; ;data;k] each indicesB)%neighbourhood[`numberToK];
 };

.quantQ.lof.LOF:{[X;data;k]
    // X -- point to find distance, can be out of data
    // data -- the data set with one column x, each observation is a row of features
    // k -- the parameter of the kNN
    // update distance to point to data
    data: update distanceToPoint: .quantQ.lof.L2distance[X;x] from data;
    // find k-neighbourhood of X
    neighbourhood: .quantQ.lof.getKNNProperties[k;data;X];
    // indices of k-neighbours -- contains neighbourhood[`numberToK] elements
    indicesY: exec i from data where distanceToPoint<=neighbourhood[`distanceToK];
    // return LOF
    :sum[.quantQ.lof.localReachabilityDensity[ ;data;k] each (exec x from data)[indicesY]]
        %neighbourhood[`numberToK]*.quantQ.lof.localReachabilityDensity[X;data;k];
 };

.quantQ.lof.LOFdata:{[data;k]
    // data -- the data set with one column x, each observation is a row of features
    // k -- the parameter of the kNN
    // LOF for every row of data
    LOFarray:{[dataSet;k;j] .quantQ.lof.LOF[(exec x from dataSet)[j];select from dataSet where i<>j;k]}
        [data;k;]  each til count[data];
    // update table and publish
    :update LOF:LOFarray from data;
 };












