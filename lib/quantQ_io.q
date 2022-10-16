.quantQ.io.saveSplayedTab:{[tabPath;pvar;table]
    // tabPath -- path where to store the table
    // pvar -- variable to sort and index on for fast querying
    // table -- name of the table to save
    @[;pvar;`p#] pvar xasc (` sv (tabPath;`;table;`)) set .Q.en[tabPath] get table
 };

.quantQ.io.appendSplayedTab:{[tabPath;pvar;table;table2Add]
    // tabPath -- path where to store the table
    // pvar -- variable to sort and index on for fast querying
    // table -- name of the saved table
    // table2Add -- table to append
    @[;pvar;`p#] pvar xasc (` sv (tabPath;`;table;`)) upsert .Q.en[tabPath] table2Add
 };

k).quantQ.io.dpftsAppend:{[d;p;f;t;s;o]r:+.Q.enxs[$;d;;s]`. . `\:t;{[d;t;o;x]@[d;x;o; t[x]]}[d:.Q.par[d;p;t];r;o]'[!r];
 @[d;`.d;:;f,r@&~f=r:!r];
 @[.q.xasc[f] d;f;`p#];
 t};















