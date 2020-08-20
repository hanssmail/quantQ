.quantQ.trees.classify:{[breakpoints;y]
    // breakpoints -- points for tree splits
    // y -- variable to classify 
    :asc[breakpoints] binr y;
 };

.quantQ.trees.entropy:{[y;classes]
    // y -- vector of sampled predicted variable
    // classes -- domain we wish to classify y, the distinct classes
    p:{sum[x=y]%count x}[y]each classes;
    :neg p wsum 2 xlog p;
 };

.quantQ.trees.infogain:{[yp;ysplit;classes]
    // yp -- the parent node classification set of attributes
    // ysplit -- the child nodes classifications set of attributes after the split
    // classes -- the set of classes
    uncond: .quantQ.trees.entropy[yp;classes];
    cond:  (wsum). flip ({count[y]%count x}[yp];.quantQ.trees.entropy[;classes])@\:/:ysplit;
    :uncond - cond;
 };

.quantQ.trees.applyRule:{[cbm;bmi;xi;rule;j]
    // cbm -- count of the bitmap: number of sampled points for predicted and predictor variable
    // bmi -- bitmap indicating indices of the sampled variables at this node
    // xi -- the ith predicted variable
    // rule -- rule to split on
    // j -- split xi at point j based on rule
    :(`bitmap`appliedrule!( @[zb; bmi;:;] not b; not);
        `bitmap`appliedrule!( @[zb:cbm#0b; bmi;:;] b:eval (rule;xi;j); (::)));
 };

.quantQ.trees.chooseSplitXi:{[yy;y;cbm;bmi;rule;classes;xi]
    // yy -- predicted variable, all sample points
    // y -- predicted variable at node
    // cbm -- count of the bitmap: number of sampled points for predicted and predictor variable
    // bmi -- bitmap indicating indices of the sampled variables at this node
    // rule -- the logical rule to apply
    // classes -- the k classification classes for y
    // xi -- the ith predicted variable
    j:asc distinct xi;
    info:{[yy;y;cbm;bmi;xi;rule;classes;j]
        split: .quantQ.trees.applyRule[cbm;bmi;xi;rule;j];
        (.quantQ.trees.infogain[y;{x where y}[yy] each split[`bitmap];classes];split)
        }[yy;y;cbm;bmi;xi;rule;classes] each j;
    :(j;info)@\:first idesc flip[info]0;
 };
 
// Wrap the rule -- helper function
.quantQ.trees.runRule:{[rule;i;j;x] 
    :rule[x i;j];
 };

// Choose optimal split for a node
.quantQ.trees.chooseSplit:{[xy;treeinput]
    // treeinput -- dictionary with splitting information
    // xy -- dictionary of
    //         predictor (`x) set
    //         predicted (`y) vector
    //         distinct vector of k classification classes (`classes)
    //         number of predictor vectors to sample (`m)
    bm: treeinput`bitmap;
    x: xy[`x][;bmi: where bm];
    y: xy[`y][bmi];
    classes: xy`classes;
    m: xy`m;
    // the logical rule to apply
    rule: treeinput`rule;rulepath: treeinput`rulepath;
    cx:count x;
    info: .quantQ.trees.chooseSplitXi[xy`y;y;count bm;bmi;rule;classes]peach x@sampled:asc neg[m]?cx;
    is:info[;1;0];
    summary: (`infogains`xi`j`infogain!enlist[sampled!is],i,(-1_r)),/:last r:raze info i:first idesc is;
    cnt: count summary;
    rulepathfunc:{[rp;ar;r;i;j] rp,enlist (ar;(`.quantQ.trees.runRule;r;i;j))};
    :update rule:rule,
        rulepath:rulepathfunc[rulepath]'[appliedrule;rule;xi;j] from summary;
 };

.quantQ.trees.growTree:{[xy;treeinput]
    // treeinput -- dictionary with splitting information
    // xy -- dictionary of
    //       predictor (`x) set
    //       predicted (`y) vector
    //       distinct vector of k classification classes (`classes)
    //       number of predictor vectors to sample (`m)
    //       maximum tree depth (`maxdepth)
    :{[xy;r]
        if[(1>=count distinct xy[`y]where r`bitmap)
            or xy[`maxdepth]<=count r`rulepath;
            :r];
         enlist[r],$[98h<>type rr:.quantQ.trees.growTree[xy;r];raze @[rr;where 99h=type each rr;enlist];rr]
    }[xy]each  r:.quantQ.trees.chooseSplit[xy;treeinput];
 };

.quantQ.trees.learnTree:{[params]
    // params -- dictionary with parameters for learning the tree
    // m -- number of features to randomly select defaults to all features
    // maxdepth -- maximum depth of tree if pure nodes not reached defaults to 100
    params[`m]: $[`m in key params;params`m;count params`x];
    params[`maxdepth]: $[`maxdepth in key params;params`maxdepth;100];
    r0: `infogains`xi`j`infogain`bitmap`x`y`appliedrule`rule`rulepath`classes`m`maxdepth#
        params,
        `infogains`xi`j`infogain`bitmap`appliedrule`rulepath!(()!();0N;0N;0n;count[params`y]#1b;::;());
    tree: enlist[r0],
    $[ 98<>type r:.quantQ.trees.growTree[; r0:dc _r0] (dc:`x`y`classes`m`maxdepth)#r0;
        raze @[r;where 99h=type each r;enlist];
        r];
    tree: update p:{x?-1_'x} rulepath  from tree;
    :`i`p`path xcols update path:{(x scan)each til count x}p,i:i from tree;
 };

.quantQ.trees.leaves:{[tree]
    // tree -- decision tree
    :select from tree where i in til[count p]except p
 };

.quantQ.trees.predictOnTree:{[tree;x]
    // tree -- a previously grown tree
    // x -- a tuple of the features at a data point i, ie X[i]
    :({[x;tree]
        if[1=count tree;:tree];
        @[tree;`rulepath;1_'] where {value y[0],value[y 1]x}[x]each tree[`rulepath][;0]
        }[x]over)[.quantQ.trees.leaves tree];
 };

