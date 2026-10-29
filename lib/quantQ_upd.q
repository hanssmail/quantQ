.quantQ.upd.trades:{[new]
    // newdata -- the incoming data sent by the upstream feedhandler or tickerplant, in table type format
    // aggregate the new data by time,volume and price
    // then insert each list to each column for each key (record) of the global trades cache
    trades:: trades,'select time,volume,price by sym from new;
    // run a 3-tick moving vwap using each-both on every record
    :select sym, (-3 sublist'volume) wavg'(-3 sublist'price) from trades;
 };

.quantQ.upd.trades:{[new]
    // newdata -- the incoming data sent by the upstream feedhandler or tickerplant, in table type format
    // aggregate the new data by time,volume and price
    // then insert each list to each column for each key (record) of the global trades cache
    trades:: trades,''select time,volume,price by sym from new;
    // run a 3-tick moving vwap
    // -3 sublist (unlike -3#) does not wrap around when fewer than 3 ticks are cached
    ret: select sym, (-3 sublist'volume) wavg'(-3 sublist'price) from trades;
    trades:: update -3 sublist'time,-3 sublist'volume,-3 sublist'price from trades;
    :ret;
 };

\

.quantQ.upd.trades newdata