.quantQ.preml.alpha2Int:{[list]
    // list -- the input list of all non-numeric features, underlying list is alphabetically ordered ordinal variable
    listIndep: asc distinct list;
    // dictionary mapping the underlying list on the integers
    :listIndep!til count listIndep;
 };

.quantQ.preml.distanceCategorical:{[point;list]
    // point --- a reference symbol
    // list -- array of symbols
    :point=list;
 };

.quantQ.preml.zScoreNorm:{[sample]
    // sample -- array of values to normalize
    :(sample - avg[sample])%sdev[sample];
 };

.quantQ.preml.zScoreNormExt:{[sample;params]
    // sample -- array of values to normalize
    // params -- dictionary with `mean and `std
    : (sample - params[`mean])%params[`std];
 };

.quantQ.preml.rangeScaleNorm:{[sample;range]
    // sample -- array of values to normalize
    // range -- two-dimensional array of a and b
    :range[0] + ((sample-min[sample])*(range[1]-range[0]))%(max[sample]-min[sample]);
 };

.quantQ.preml.rangeScaleNormExt:{[sample;range;params]
    // sample -- array of values to normalize
    // range -- two-dimensional array of a and b
    // params -- dictionary with `min and `max
   :range[0] + ((sample-params[`min])*(range[1]-range[0]))%(params[`max]-params[`min]);
 };


