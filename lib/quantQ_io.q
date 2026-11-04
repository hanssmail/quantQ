.quantQ.io.saveSplayedTab:{[tabPath;pvar;table]
    // tabPath -- path where to store the table
    // pvar -- variable to sort and index on for fast querying
    // table -- name of the table to save
    // sorted and given the parted attribute in memory, then written once:
    // the same table on disk as sorting and amending the attribute in place,
    // and it needs only set, which every q implementation supports
    (` sv tabPath,table,`) set @[;pvar;`p#] pvar xasc .Q.en[tabPath] get table
 };

.quantQ.io.appendSplayedTab:{[tabPath;pvar;table;table2Add]
    // tabPath -- path where to store the table
    // pvar -- variable to sort and index on for fast querying
    // table -- name of the saved table
    // table2Add -- table to append
    // the whole table is re-sorted, so it is rewritten: read, append, sort, write
    path:` sv tabPath,table,`;
    path set @[;pvar;`p#] pvar xasc (select from get path),.Q.en[tabPath] table2Add
 };

.quantQ.io.dpftsAppend:{[d;p;f;t;s;o]
    // d -- root of the partitioned database, as a file handle
    // p -- partition, e.g. a date
    // f -- column to sort on and give the parted attribute
    // t -- name of the table to write
    // s -- name of the sym file the symbol columns are enumerated against
    // o -- how the partition is amended: (:) overwrites, (,) appends
    // written in q rather than k, and through set alone, so that it runs on
    // q implementations without k or on-disk amend (PeachQ, m)
    r:$[s~`sym;.Q.en[d] get t;.Q.enxs[$;d;get t;s]];
    path:` sv .Q.par[d;p;t],`;
    // the sort column first, as .Q.dpft lays a partition out
    c:f,cols[r] except f;
    old:$[()~key ` sv .Q.par[d;p;t],`.d;0#r;select from get path];
    path set @[;f;`p#] f xasc o[c#old;c#r];
    t};















