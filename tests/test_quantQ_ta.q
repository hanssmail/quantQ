/ Unit tests for lib/quantQ_ta.q -- technical-analysis indicators (.quantQ.ta)
/ Reference series computed in Python with q's windowing conventions (mavg/mdev/mmax/mmin use the
/ observations available at the start of the series).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`ta;

taO:10 11 12 11 13 14 13 15 16 15 17 18f;
taC:10.5 11.8 11.2 12.4 13.6 13.1 14.2 15.5 15.2 16.4 17.9 17.3;
taH:10.8 12.1 12.3 12.7 13.9 14.3 14.5 15.8 16.3 16.7 18.2 18.3;
taL:9.6 10.6 10.8 10.6 12.6 12.7 12.6 14.6 14.8 14.6 16.6 16.9;
taP:([] price:taO);
taOHLC:([] open:taO; high:taH; low:taL; close:taC);
taM3:enlist[`memory]!enlist 3;

/ moving averages and deviations
.qt.eq["ma adds a column named <col>MA<memory>";{cols .quantQ.ta.ma[`price;taM3;taP]};`price`priceMA3];
.qt.eq["ma default memory is 10";{last cols .quantQ.ta.ma[`price;()!();taP]};`priceMA10];
.qt.approx["ma (mavg 3)";{exec priceMA3 from .quantQ.ta.ma[`price;taM3;taP]};10 10.5 11 11.333333333333334 12 12.666666666666666 13.333333333333334 14 14.666666666666666 15.333333333333334 16 16.666666666666668;1e-12];
.qt.approx["ema with alpha = 2/(memory+1)";{exec priceEMA3 from .quantQ.ta.ema[`price;taM3;taP]};10 10.5 11.25 11.125 12.0625 13.03125 13.015625 14.0078125 15.00390625 15.001953125 16.0009765625 17.00048828125;1e-12];
.qt.approx["msd is the moving standard deviation";{exec priceMSD3 from .quantQ.ta.msd[`price;taM3;taP]};0 0.5 0.816496580927726 0.471404520791032 0.816496580927726 1.247219128924647 0.471404520791032 0.816496580927726 1.247219128924647 0.471404520791032 0.816496580927726 1.247219128924647;1e-12];
.qt.eq["bollingerBands column names";{cols .quantQ.ta.bollingerBands[`price;taM3;taP]};`price`priceMiddleBound3`priceLowerBound3`priceUpperBound3];
.qt.approx["bollingerBands: middle band is the moving average";{exec priceMiddleBound3 from .quantQ.ta.bollingerBands[`price;taM3;taP]};mavg[3;taO];1e-12];
.qt.approx["bollingerBands: bands are 2 standard deviations away";{t:.quantQ.ta.bollingerBands[`price;taM3;taP]; exec priceUpperBound3-priceMiddleBound3 from t};2*mdev[3;taO];1e-12];
.qt.approx["momentum: x_t - x_(t-m), 0 for the first m points";{exec priceMom3 from .quantQ.ta.momentum[`price;taM3;taP]};0 0 0 1 2 2 2 2 2 2 2 2f;1e-12];

/ OHLC indicators
.qt.approx["accBands: upper band";{exec closeUpperBand3 from .quantQ.ta.accBands[`open`high`low`close;taM3;taOHLC]};13.34117647058824 14.319707178025396 14.711406516951996 16.023861284266648 16.466964379049568 17.198572523991047 17.627865821215867 18.183595711980512 18.768617346691325 19.640406837749172 20.724536702641284 21.313426403508476;1e-9];
.qt.approx["accBands: middle band";{exec closeMiddleBand3 from .quantQ.ta.accBands[`open`high`low`close;taM3;taOHLC]};mavg[3;taC];1e-12];
.qt.approx["accBands: lower band";{exec closeLowerBand3 from .quantQ.ta.accBands[`open`high`low`close;taM3;taOHLC]};7.341176470588232 7.569707178025393 7.711406516951993 7.523861284266647 8.300297712382902 8.865239190657713 9.627865821215863 10.350262378647175 11.101950680024656 11.64040683774917 12.057870035974622 12.81342640350848;1e-9];

taIchi:{.quantQ.ta.Ichimoku[`open`high`low`close;`memoryShort`memoryMedium`memoryLong!2 3 4;taOHLC]};
.qt.approx["Ichimoku: turning line";{exec closeTurningLine2x3x4 from taIchi[]};10.2 10.85 11.45 11.65 12.25 13.45 13.55 14.2 15.45 15.65 16.4 17.45;1e-9];
.qt.approx["Ichimoku: standard line";{exec closeStandardLine2x3x4 from taIchi[]};10.2 10.85 10.95 11.65 12.25 12.45 13.55 14.2 14.45 15.65 16.4 16.45;1e-9];
.qt.approx["Ichimoku: leading span 1";{exec closeLeadingSpan1x2x3x4 from taIchi[]};0 0 0 10.2 10.85 11.2 11.65 12.25 12.95 13.55 14.2 14.95;1e-9];
.qt.approx["Ichimoku: leading span 2";{exec closeLeadingSpan2x2x3x4 from taIchi[]};0 0 0 10.2 10.85 10.95 11.15 12.25 12.45 12.55 14.2 14.45;1e-9];

taMacd:{.quantQ.ta.MACD[`price;`memoryFast`memorySlow`memory!2 4 3;taP]};
.qt.approx["MACD line";{exec priceMACD2x4 from taMacd[]};0 0 0.5 0.5 0.25 1 0.75 0.25 1 0.75 0.25 1;1e-12];
.qt.approx["MACD signal line";{exec priceMACDSignalLine2x4x3 from taMacd[]};mavg[3;0 0 0.5 0.5 0.25 1 0.75 0.25 1 0.75 0.25 1];1e-12];
.qt.approx["MACD histogram = MACD - signal";{t:taMacd[]; exec priceMACDHistogram2x4x3-priceMACD2x4-priceMACDSignalLine2x4x3 from t};12#0f;1e-12];

.qt.approx["Aroon up";{exec priceAroonUp3 from .quantQ.ta.Aroon[`price;taM3;taP]};1 1 1 0.5 1 1 0.5 1 1 0.5 1 1;1e-12];
.qt.approx["Aroon down";{exec priceAroonDown3 from .quantQ.ta.Aroon[`price;taM3;taP]};1 0.5 0 1 0.5 0 1 0.5 0 1 0.5 0;1e-12];
.qt.approx["Aroon oscillator = up - down";{exec priceAroonOscillator3 from .quantQ.ta.Aroon[`price;taM3;taP]};0 0.5 1 -0.5 0.5 1 -0.5 0.5 1 -0.5 0.5 1;1e-12];

/ RSI: only gains -> 100, only losses -> 0, equal gains and losses -> 50
.qt.approx["RSI with only gains is 100";{last exec closeRSI3 from .quantQ.ta.RSI[`open`close;taM3;([] open:10#1f; close:10#2f)]};100f;1e-12];
.qt.approx["RSI with only losses is 0";{last exec closeRSI3 from .quantQ.ta.RSI[`open`close;taM3;([] open:10#2f; close:10#1f)]};0f;1e-12];
.qt.approx["RSI with balanced gains and losses is 50";{last exec closeRSI3 from .quantQ.ta.RSI[`open`close;taM3;([] open:10#1 2f; close:10#2 1f)]};50f;1e-12];
.qt.approx["RSI with no movement is 50";{last exec closeRSI3 from .quantQ.ta.RSI[`open`close;taM3;([] open:10#1f; close:10#1f)]};50f;1e-12];

/ ATR (close lies inside [low, high], so the true range is high - low)
.qt.approx["ATR";{exec closeATR3 from .quantQ.ta.ATR[`open`high`low`close;taM3;taOHLC]};mavg[3;taH-taL];1e-12];

.qt.finish[];
