/ Unit tests for lib/quantQ_mcts.q -- multi-armed bandits and Monte-Carlo tree search on tic-tac-toe (.quantQ.mcts)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`mcts;

/ the tic-tac-toe functions read the board from a global called movesSpace (as in the book)
movesSpace:`a1`a2`a3`b1`b2`b3`c1`c2`c3;

/ multi-armed bandit
.qt.true["MAB draws uniformly below the chosen arm's bound";{r:.quantQ.mcts.MAB[;2 5 9] each 200#1; (all r>=0) and (all r<5) and any r>2}];
mctsBandit:{[eps] `counter`epsilon`N`K`strategy`payoffs`mab`minBound!(1;eps;100;3;`long$();(1 1f;5 5f;enlist 2f);{[x] 10f};0f)};
.qt.eq["mabExplorationFirstGA exploits the best average arm";{(.quantQ.mcts.mabExplorationFirstGA mctsBandit 0f)`strategy};enlist 1];
/ counter starts at 1: the first-epsilon rule explores while counter <= epsilon*N, i.e. also at counter 0
.qt.eq["mabExplorationFirstGA records the payoff on the chosen arm";{(.quantQ.mcts.mabExplorationFirstGA mctsBandit 0f)`payoffs};(1 1f;5 5 10f;enlist 2f)];
.qt.eq["mabExplorationFirstGA increments the counter";{(.quantQ.mcts.mabExplorationFirstGA mctsBandit 0f)`counter};2];
.qt.true["mabExplorationFirstGA explores within the first epsilon*N rounds";{s:((.quantQ.mcts.mabExplorationFirstGA/)[50;mctsBandit 1f])`strategy; (50=count s) and all s in til 3}];
.qt.eq["mabExplorationRandomGA with epsilon 0 always exploits";{distinct ((.quantQ.mcts.mabExplorationRandomGA/)[20;mctsBandit 0f])`strategy};enlist 1];
.qt.eq["mabExplorationDecreasingGA with epsilon 0 always exploits";{distinct ((.quantQ.mcts.mabExplorationDecreasingGA/)[20;mctsBandit 0f])`strategy};enlist 1];
.qt.true["mabExplorationRandomGA with epsilon > 1 always explores";{s:((.quantQ.mcts.mabExplorationRandomGA/)[300;mctsBandit 1.1])`strategy; (asc distinct s)~0 1 2}];

/ tic-tac-toe: player 1 plays the even moves, player 2 the odd ones
mctsDraw:`a1`a2`a3`b2`b1`b3`c2`c1`c3;
.qt.approx["checkGameWinner: game in progress";{(.quantQ.mcts.checkGameWinner enlist[`movesState]!enlist `a1`b1)`isWinner};0;0f];
.qt.approx["checkGameWinner: player 1 wins a row";{(.quantQ.mcts.checkGameWinner enlist[`movesState]!enlist `a1`b1`a2`b2`a3)`isWinner};1;0f];
.qt.approx["checkGameWinner: player 1 wins a column";{(.quantQ.mcts.checkGameWinner enlist[`movesState]!enlist `a1`b2`b1`c3`c1)`isWinner};1;0f];
.qt.approx["checkGameWinner: player 1 wins a diagonal";{(.quantQ.mcts.checkGameWinner enlist[`movesState]!enlist `a1`a2`b2`a3`c3)`isWinner};1;0f];
.qt.approx["checkGameWinner: player 2 wins";{(.quantQ.mcts.checkGameWinner enlist[`movesState]!enlist `b1`a1`b2`a2`c3`a3)`isWinner};2;0f];
.qt.approx["checkGameWinner: full board without a line is a draw";{(.quantQ.mcts.checkGameWinner enlist[`movesState]!enlist mctsDraw)`isWinner};3;0f];
.qt.eq["visualiseBoard: player 1 is o, player 2 is x";{.quantQ.mcts.visualiseBoard enlist[`movesState]!enlist mctsDraw};(`o`x`o;`o`x`x;`x`o`o)];

.qt.true["randomGameMove plays one new legal move";{b:.quantQ.mcts.randomGameMove enlist[`movesState]!enlist `a1`b2; m:last b`movesState; (3=count b`movesState) and (m in movesSpace) and not m in `a1`b2}];
.qt.true["oneTurnRun returns a valid game state";{b:.quantQ.mcts.oneTurnRun `movesState`isWinner!(`symbol$();0); b[`isWinner] in 0 1 2 3}];
.qt.true["a random game always terminates with a result";{b:(.quantQ.mcts.oneTurnRun/)[{0=x`isWinner};`movesState`isWinner!(`symbol$();0)]; (b[`isWinner] in 1 2 3) and count[b`movesState] within 5 9}];

/ UCT
.qt.eq["uct with c=0 picks the move with the best win rate";{.quantQ.mcts.uct[([] gameTree:`a1`a1`b2; isWinner:1 2 1);1;0f]};`b2];
.qt.eq["uct with a large c favours the less explored move";{.quantQ.mcts.uct[([] gameTree:`a1`a1`a1`b2; isWinner:1 1 1 2);1;100f]};`b2];

/ MCTS
mctsBucket:`movesState`treeMCTS`numberOfRuns`c`isWinner!(`symbol$();([] index:`long$(); gameTree:(); isWinner:`long$());20;1.4;0);
.qt.eq["runMCTS adds one game to an empty tree";{count (.quantQ.mcts.runMCTS mctsBucket)`treeMCTS};1];
.qt.true["runMCTS builds the tree up to numberOfRuns games";{b:(.quantQ.mcts.runMCTS/)[{0=x`isWinner};mctsBucket]; 20=count b`treeMCTS}];
.qt.true["evaluateMCTS returns a legal first move";{b:(.quantQ.mcts.runMCTS/)[{0=x`isWinner};mctsBucket]; .quantQ.mcts.evaluateMCTS[b`treeMCTS] in movesSpace}];

.qt.finish[];
