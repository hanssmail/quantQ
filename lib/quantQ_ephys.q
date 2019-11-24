.quantQ.ephys.powerLaw:{[alpha;beta;x]
	// alpha -- the coefficient
	// beta -- the power of the distribution
	// x -- the independent variable 
	:alpha*(x xexp neg beta);
 };

.quantQ.ephys.getSquareAround:{[nGrid;coord]
    // nGrid -- the overall size of the grid for cyclic conditions
    // coord -- 2-dimensional coordinate
    grid:(cross).{((y+-1 0 1)mod x)}[nGrid]each coord;
    // exclude the point itself
    :grid except enlist coord;
 };

.quantQ.ephys.probPositiveSpin:{[beta;localEnergy]
    // beta -- the paratemer proportional to inverse temperature
    // localEnergy -- the local energy of a trader
    :1.0%(1+ exp[-2*beta*localEnergy]);
 };

.quantQ.ephys.oneRunIsing:{[isingModel]
    // isingModel -- object with all data
    // randomly choose one candidate to update
    pivot:2?isingModel[`nGrid];
    // recover price
    price:last isingModel[`priceIsing][`price];
    t:last isingModel[`priceIsing][`t];
    // calculate local energy
    term1: (first exec spin from isingModel[`wallStreet] where coord in enlist pivot)*sum 
        exec spin from isingModel[`wallStreet] where coord in isingModel[`getSquareAround][isingModel[`nGrid];pivot];
    term2:isingModel[`alpha]*abs[price]*first exec spin from isingModel[`wallStreet] 
        where coord in enlist pivot;
    localEnergy:term1-term2;
    // probability for positive spin
    pPos:isingModel[`probPositiveSpin][isingModel[`beta];localEnergy];   
    // set the new spin
    newSpin:$[(first 1?1.0)<=pPos;1;-1];
    isingModel[`wallStreet]:update spin:newSpin from isingModel[`wallStreet] 
        where coord in enlist pivot;
    // update price and time in the priceIsing table
    isingModel[`priceIsing]:isingModel[`priceIsing] upsert (t+1;sum exec spin from 
        isingModel[`wallStreet]);
    :isingModel;
 };
