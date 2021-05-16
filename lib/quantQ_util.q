// Utility name space with various functions

// Wrapper for functional select
.quantQ.util.selectCol:{[tab;listCols]
    // tab -- table (pass by value or reference)
    // listCols -- array of symbols with columns to select
    :?[tab; (); 0b; ((),listCols)!((),listCols)];
 };
// exa: tab: ([] a1: til 10; a2: til 10; a3: til 10);
// .quantQ.util.selectCol[tab;`a1]
// .quantQ.util.selectCol[`tab;`a1`a2] 

// Wrapper for functional select
.quantQ.util.deleteCol:{[tab;listCols]
    // tab -- table (pass by value or reference)
    // listCols -- array of symbols with columns to select
    :![tab; (); 0b; (),listCols];
 };
// exa: tab: ([] a1: til 10; a2: til 10; a3: til 10);
// .quantQ.util.deleteCol[tab;`a1]

// Generalised prev
.quantQ.util.prev:{[n;x]
    // n -- number lags
    // x -- array
    :(prev/)[n;] x;
 };
// exa .quantQ.util.prev[5] til 10

// Generalised next
.quantQ.util.next:{[n;x]
    // n -- number lags
    // x -- array
    :(next/)[n;] x;
 };
// exa .quantQ.util.next[5] til 10
