.quantQ.rf.sampleTree:{[s;n]
    // s -- sample dictionary with predictor and predicted variables
    // n -- sample size
    z:`x`y!(s[`x][;i];s[`y]i:n?n);
    :z,`oobi`ibi!((til n) except distinct i; i);
 };

.quantQ.rf.bootstrapTree:{[params;m;n;B]
    // params -- same as the ones required by .quantQ.trees.learnTree
    // m -- select m of the features in each split
    // n -- sample size
    // B -- size of bootstrap: create B bootstrap sample trees
    z: .quantQ.rf.sampleTree[`x`y#params;n];
    tree_b:   .quantQ.trees.learnTree @[params;`x`y;:;z`x`y],enlist[`m]!enlist m;
    tree_oob: raze .quantQ.trees.predictOnTree[tree_b]each flip params[`x;;z`oobi];
    tree_oob: update pred_error: abs obs_y-{first x where y}[z`y]each bitmap from
        update obs_y:params[`y]z`oobi from tree_oob;
    :`tree`oob`ibi!
        (`B xcols update B from tree_b;
        `B xcols update B from tree_oob;
        enlist[B]!enlist z`ibi);
 };

.quantQ.rf.randomForest:{[params]
    // params dictionary,
    ensemble: .quantQ.rf.bootstrapTree[params;params`m;params`n] peach til params`B;
    :raze each flip ensemble;
 };

.quantQ.rf.predictOnRF:{[y;ensemble;data]
    // y -- predicted variable vector
    // ensemble -- a random forest: a dictionary of
    //              `tree: list of treetables, a table itself
    //              `oob:  the table of out-of-bag predictions
    //              `ibi:  in-the-bag indices for each sample B
    // data -- a tuple of the features at a data point i, ie X[i]
    // returns a classification of data X based on majority rule
    rf:{[data;tree;ibi;b]
        prediction: .quantQ.trees.predictOnTree[select from tree where B=b] data;
        update ibi: enlist ibi b from prediction
       }[data;tree;ensemble`ibi] each exec distinct B from tree:ensemble`tree;
    prediction: {first where x=max x}count each
        group exec {first x[y] where z}[y]'[ibi; bitmap] from raze rf;
    :`prediction`mean_error`dev_error!
        enlist[prediction],value exec avg pred_error,dev pred_error from ensemble`oob;
 };
