/ Unit tests for lib/quantQ_upd.q -- real-time update handler with a 3-tick VWAP (.quantQ.upd)
/ The handler keeps its state in a global keyed table called `trades` (as in the book).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`upd;

updReset:{[] `trades set ([sym:`a`b] time:2#enlist `time$(); volume:2#enlist `long$(); price:2#enlist `float$()); };
updNew1:([] time:09:00:00.000 09:00:01.000 09:00:02.000 09:00:03.000; sym:`a`b`a`b; volume:1 2 3 2; price:10 5 20 7f);
updNew2:([] time:09:00:04.000 09:00:05.000 09:00:06.000 09:00:07.000; sym:`a`b`a`b; volume:4 1 2 1; price:30 6 40 8f);
updVwap:{[r] last value flip 0!r};

.qt.approx["trades: VWAP after the first batch";{updReset[]; updVwap .quantQ.upd.trades updNew1};(70%4;24%4);1e-12];
.qt.approx["trades: VWAP uses only the last 3 ticks";{updReset[]; .quantQ.upd.trades updNew1; updVwap .quantQ.upd.trades updNew2};((60+120+80)%9;(14+6+8)%4);1e-12];
.qt.eq["trades: the cache is trimmed to the last 3 ticks";{updReset[]; .quantQ.upd.trades updNew1; .quantQ.upd.trades updNew2; count each exec price from trades};3 3];
.qt.eq["trades: one result row per symbol";{updReset[]; count .quantQ.upd.trades updNew1};2];

.qt.finish[];
