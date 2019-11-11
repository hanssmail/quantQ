.quantQ.mcts.MAB:{[choice;mabSpec]
    // mabSpec -- vector of maxBoundK's
    // choice -- choice to be drawn
    :first 1?"f"$mabSpec[choice];
 };

.quantQ.mcts.mabExplorationFirstGA:{[bucket]
    // bucket -- the dictionary with all parameters
    $[bucket[`counter]<=bucket[`epsilon]*bucket[`N];
    // exploration
    [
        // decision randomly
        choice:first 1?bucket[`K];
        bucket[`strategy],:choice;
        payoffs:bucket[`payoffs];
        payoffs[choice],:bucket[`mab] choice;
        bucket[`payoffs]:payoffs
    ];
    // exploitation
    [
        // decision based on the best outcome
        avgReward:bucket[`minBound]^avg each bucket[`payoffs];
        choice: first where avgReward=max avgReward;
        bucket[`strategy],:choice;
        payoffs:bucket[`payoffs];
        payoffs[choice],:bucket[`mab] choice;
        bucket[`payoffs]:payoffs
    ]
    ];    
    // increase counter
    bucket[`counter]+:1;
    // return bucket
    :bucket;
 };    

.quantQ.mcts.mabExplorationRandomGA:{[bucket]
    // bucket -- the dictionary for algorithm
    $[bucket[`epsilon]>first 1?1.0;
    // exploration
    [
        // decision randomly
        choice:first 1?bucket[`K];
        bucket[`strategy],:choice;
        payoffs:bucket[`payoffs];
        payoffs[choice],:bucket[`mab] choice;
        bucket[`payoffs]:payoffs
    ];
    // exploitation
    [
        // decision based on the best outcome
        avgReward:bucket[`minBound]^avg each bucket[`payoffs];
        choice: first where avgReward=max avgReward;
        bucket[`strategy],:choice;
        payoffs:bucket[`payoffs];
        payoffs[choice],:bucket[`mab] choice;
        bucket[`payoffs]:payoffs
    ]];    
    // increase counter
    bucket[`counter]+:1;
    // return bucket
    :bucket;
 };    

.quantQ.mcts.mabExplorationDecreasingGA:{[bucket]
    // bucket -- the dictionary for algorithm
    $[(bucket[`epsilon]%(sqrt 1.0+bucket[`counter]))>first 1?1.0;
    // exploration
    [
        // decision randomly
        choice:first 1?bucket[`K];
        bucket[`strategy],:choice;
        payoffs:bucket[`payoffs];
        payoffs[choice],:bucket[`mab] choice;
        bucket[`payoffs]:payoffs
    ];
    // exploitation
    [
        // decision based on the best outcome
        avgReward:bucket[`minBound]^avg each bucket[`payoffs];
        choice: first where avgReward=max avgReward;
        bucket[`strategy],:choice;
        payoffs:bucket[`payoffs];
        payoffs[choice],:bucket[`mab] choice;
        bucket[`payoffs]:payoffs
    ]];    
    // increase counter
    bucket[`counter]+:1;
    // return bucket
    :bucket;
 };    

.quantQ.mcts.randomGameMove:{[bucket]
    // bucket -- information about the game
    // list of possible moves
    possibleMoves: movesSpace where not movesSpace in bucket[`movesState];
    // update movesState by a random move
    bucket[`movesState]: bucket[`movesState],1?possibleMoves;
    // return updated bucket
    :bucket;    
 };

.quantQ.mcts.checkGameWinner:{[bucket]
    // bucket -- information about the game    
    // extract moves of two players
    player1: bucket[`movesState]{x where 0=mod[x;2]} til count bucket[`movesState];
    player2: bucket[`movesState]{x where 1=mod[x;2]} til count bucket[`movesState];
    // check the winners
    p1win:(sum[player1 in `a1`b2`c3]=3) or (sum[player1 in `a3`b2`c1]=3) or (0<sum 3=sum 
        each player1 in/:{`$string[x] cross \:string[y]}[`a`b`c;`1`2`3]) or (0<sum 3=sum 
        each player1 in/:{`${reverse each x} each string[y] cross \:string[x]}[`a`b`c;`1`2`3]);
    p2win:(sum[player2 in `a1`b2`c3]=3) or (sum[player2 in `a3`b2`c1]=3) or (0<sum 3=sum 
        each player2 in/:{`$string[x] cross \:string[y]}[`a`b`c;`1`2`3]) or (0<sum 3=sum 
        each player2 in/:{`${reverse each x} each string[y] cross \:string[x]}[`a`b`c;`1`2`3]);
    // update winner
    bucket[`isWinner]:$[0=count player1;0;p1win]+($[0=count player2;0;p2win]*2)+$[((p1win=0) 
        or (0=count player1)) and ($[0=count p2win;1;p2win=0] or (0=count player2)) and 
        prd[bucket[`movesState][til 9] in movesSpace];1;0]*3;
    // output bucket
    :bucket;
 };

.quantQ.mcts.oneTurnRun:{[bucket]
    // bucket -- information about the game  
    :.quantQ.mcts.checkGameWinner .quantQ.mcts.randomGameMove bucket;
 };

.quantQ.mcts.visualiseBoard:{[bucket]
    // bucket -- information about the game  
    // player 1 moves
    player1: bucket[`movesState]{x where 0=mod[x;2]} til count bucket[`movesState];
    // player 2 moves
    player2: bucket[`movesState]{x where 1=mod[x;2]} til count bucket[`movesState];
    // visualise the board
    :3 3#{`$x} each ssr[;"0";"."] ssr[;"2";"x"] ssr[;"1";"o"] raze string  (movesSpace in
        player1)+(2*movesSpace in player2);
 };

.quantQ.mcts.uct:{[tab;perspective;c]
    // tab -- table with history of data (gameTree and isWinner columns)
    // perspective -- player making move
    // c -- parameter of the UCT
    // adjust table by perspective -- isWinner
    tab: update w: ?[(isWinner=perspective) and isWinner<3;1;?[(not isWinner=perspective) 
        and isWinner<3;0;0.5] ] from tab;
    // return choice based on UCT
    :first exec gameTree from `uct xdesc select gameTree, uct: (wSum%"f"$nSum)+c*sqrt 
        log["f"$count[tab]]%nSum from select nSum: count w, wSum: sum w by gameTree from tab;
 };

.quantQ.mcts.oneMCTSTurnRun:{[bucket]
    // bucket -- information about the game  
    availableMoves: asc movesSpace where not movesSpace in bucket[`movesState];
    // relevant trees, exception for first iteration 
      relevantTrees: $[0=count bucket[`treeMCTS];bucket[`treeMCTS] ;select index, gameTree, isWinner from (update {[x;y] x~y}[bucket[`movesState]; ] each subTree from update subTree: gameTree[;til count bucket[`movesState]] from bucket[`treeMCTS]) where subTree=1];
    // all moves used in given level of tree, distinct and non-empty
    allExploredMoves:asc t where not null t:distinct {(x,())[y]}[;bucket[`treeDepthCounter]] each exec gameTree from relevantTrees;
    // decide what to do:
    $[(availableMoves~allExploredMoves);
        // use UCT rule, when we have explored all moves we can do
        [
            // decide who is on move -- affects the UCT rule below
            isOnMove:1+mod[count bucket[`movesState];2];
            // use UCT rule to find the best move 
            bestMove: .quantQ.mcts.uct[select {(x,())[y]}[;bucket[`treeDepthCounter]] each gameTree, 
                isWinner from bucket[`treeMCTS];isOnMove;bucket[`c]];
            movesStateTMP:bucket[`movesState];
            bucket[`movesState]:movesStateTMP,bestMove;
            // increase indexInTree when availableMoves>1   
            indexTMP:bucket[`indexInTree];
            bucket[`indexInTree]:indexTMP+"j"$(1<count availableMoves);
        ];
        // use random non-explored node
        [
            movesStateTMP:bucket[`movesState];
            bucket[`movesState]:movesStateTMP,1?availableMoves where not availableMoves in
                allExploredMoves;
        ]
    ];
    // increase the treeDepthCounter by 1
    tmp: bucket[`treeDepthCounter];
    bucket[`treeDepthCounter]:tmp+1;
    // return the bucket with one turn and with updated game check
    :.quantQ.mcts.checkGameWinner bucket;
 };

.quantQ.mcts.runMCTS:{[bucket]
    // decide move first vs further 
    $[0=count bucket[`treeMCTS];
        // randomly finish game
        [
        bucketTMP:(.quantQ.mcts.oneTurnRun/)[{0=x[`isWinner]};(`movesState`isWinner)!
            (bucket[`movesState];0)];
        tabTMP: bucket[`treeMCTS];
        bucket[`treeMCTS]:tabTMP,([] index:1; gameTree:enlist enlist bucketTMP[`movesState]
            [count bucket[`movesState]]; isWinner:bucketTMP[`isWinner])
        ];
        // continue building tree
        [
        // create local bucket and add features    
        bucketTMP:bucket;
        // indexInTree defines what part of the game is stored into the persisitng MCTS table
        bucketTMP[`indexInTree]:0j;
        // treeDepthCounter counts local iterations, to keep track when going through the existing tree
        bucketTMP[`treeDepthCounter]:0j;
        // isWinner determines local winner
        bucketTMP[`isWinner]:0j;
        // run the MCTS game
        bucketTMP:(.quantQ.mcts.oneMCTSTurnRun/)[{0=x[`isWinner]};bucketTMP];
        // extract variables for updated gameTree     
        newGameTree:((count bucket[`movesState])_bucketTMP[`movesState])[til 1+bucketTMP[`indexInTree]];
        newIsWinner:bucketTMP[`isWinner];
        newIndex: count[bucket[`treeMCTS]]+1;
        // update the existing game tree    
        bucketTreeTMP: bucket[`treeMCTS];        
        bucket[`treeMCTS]:bucketTreeTMP, ([] index: newIndex; gameTree:enlist
            newGameTree;isWinner:newIsWinner)
        ]
    ];
    // stop if the maximum number of trees have been reached
    if[(bucket[`numberOfRuns])<=count bucket[`treeMCTS];bucket[`isWinner]:1j];
    // return bucket
    :bucket;
 };

.quantQ.mcts.evaluateMCTS:{[tab]
    // tab -- table with `treeMCTS
    // decide who is on move -- it affects perspective
    isOnMove:1+mod[count bucket[`movesState];2];
    // convert isWinner to numerical value
    tab: update w: ?[(isWinner=1) and isWinner<3;1;?[(not isWinner=1) and isWinner<3;0;0.5] ] from tab;
    // evaluate each possible first move -- take the one with best average value
    tab: 0!select avg w by gameTree from select first each gameTree, w from tab;
    // return one of the best moves
    :first $[isOnMove=1;1?exec gameTree from (0!tab) where w=max[tab `w];1?exec gameTree 
        from (0!tab) where w=min[tab `w]];
 };
