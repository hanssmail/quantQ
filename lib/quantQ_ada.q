.quantQ.ada.stump:{[variable;threshold]
    // variable -- array of variables
    // threshold -- threshold for > operation
    :neg[1]+2*eval (>;variable;threshold);
 };


.quantQ.ada.addNegFeatureTable:{[tab]
    // tab -- input table, 3rd column onwards are features
    colsx: 2_cols tab;
    coly: first 1#cols tab;
    //add negative features 
    :tab,'flip raze{[tab;yy;name] flip
     eval (?;tab;();0b;(enlist`$string[name],"n")!enlist (-:;name))
     }[tab;coly;] each colsx;
 };

.quantQ.ada.weightedError:{[tab;prediction]
    // tab -- table with data and current weights
    // prediction -- array of predictions
    :sum tab[`weights] where tab[`y]<>prediction;
 };

.quantQ.ada.runOneAdaBoost:{[bucket]
    // bucket -- dictionary with tab and ada tables
    tab:bucket[`tab];
    ada:bucket[`ada];
    // get names of features and of the dependent variable
    colsx: 2_cols tab;
    coly: first 1#cols tab;
    //randomly choose feature to split
    feature: first 1?colsx;
    // iterate through possible stumps and find the best one, if more than exists, choose randomly one
    optimalStump: first 1?(flip (t;z)) where (min[t])=t:.quantQ.ada.weightedError[tab;] each
        .quantQ.ada.stump[tab[feature];] each z: distinct asc tab[feature];
    // error of the optimal stump
    err: optimalStump[0];
    //split point
    splitPoint:optimalStump[1];
    // alpha of the adaBoost
    alpha: 0.5*log[(1-err)%err];
    // correction to weights
    corr2w: exp neg alpha*tab[coly]*.quantQ.ada.stump[tab[feature];splitPoint];
    // update data table
    tab:update weights: (weights*corr2w)%sum[weights*corr2w] from tab;
    // update ada table 
    ada: ada,([] i:(1+count ada);alpha: alpha;stump: enlist(feature;splitPoint;err));
    // return bucket
    :((`tab`ada)!(tab;ada));
 };




.quantQ.ada.adaBoostPrediction:{[bucket]
    // bucket -- dictionary with tab and ada tables
    tab:bucket[`tab];
    ada:bucket[`ada];
    // calculate the sum of the individual predictions
    res: sum {[tab;config] config[0]*.quantQ.ada.stump[tab[config[1]];config[2]]}[tab;] each
        (ada[`alpha],'ada[`stump]);
    // add prediction to tab table
    :update isPredicted: prediction*y from update prediction: signum res from ([] res:res; y:tab[`y]);
 };    

.quantQ.ada.adaPredictionStats:{[tabAdaPrediction]
    // tabAdaPrediction -- output of .quantQ.ada.adaBoostPrediction function
    truePositive:sum (tabAdaPrediction[`y]=tabAdaPrediction[`prediction]) and (tabAdaPrediction[`y]=1);
    trueNegative:sum (tabAdaPrediction[`y]=tabAdaPrediction[`prediction]) and (tabAdaPrediction[`y]=-1);
    falsePositive:sum (tabAdaPrediction[`y]=1) and (tabAdaPrediction[`prediction]=-1);
    falseNegative:sum (tabAdaPrediction[`y]=-1) and (tabAdaPrediction[`prediction]=1);
    accuracy:("f"$truePositive+trueNegative)%"f"$count tabAdaPrediction;
    recall:("f"$truePositive)%"f"$(truePositive+falseNegative);
    precision: ("f"$truePositive)%"f"$(truePositive+falsePositive);
    specificity: ("f"$trueNegative)%"f"$falsePositive+trueNegative;    
    F1score:2.0%((1.0%recall)+(1.0%precision));
    // output is dictionary with all measures
    :`truePositive`trueNegative`falsePositive`falseNegative`accuracy`recall`precision`specificity`F1score!(truePositive;trueNegative;falsePositive;falseNegative;accuracy;recall;
        precision;specificity;F1score);
 };



