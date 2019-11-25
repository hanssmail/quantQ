.quantQ.upd.trades:{[new]
    // newdata -- the incoming data sent by the upstream feedhandler or tickerplant, in table type format
    // aggregate the new data by time,volume and price
    // then insert each list to each column for each key (record) of the global trades cache
    trades:: trades,'select time,volume,price by sym from new;
    // run a 3-tick moving vwap using each-both on every record
    :select sym, (-3#'volume) wavg'(-3#'price) from trades;
 };

.quantQ.upd.trades:{[new]
    // newdata -- the incoming data sent by the upstream feedhandler or tickerplant, in table type format
    // aggregate the new data by time,volume and price
    // then insert each list to each column for each key (record) of the global trades cache
    trades:: trades,''select time,volume,price by sym from new;
    // run a 3-tick moving vwap
    ret: select sym, (-3#'volume) wavg'(-3#'price) from trades;
    trades:: update -3#'time,-3#'volume,-3#'price from trades;
    :ret;
 };

\

.quantQ.upd.trades newdata