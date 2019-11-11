.quantQ.art.placeCircle: {[diameter; canvas]
    // diameter -- diameter of the circle in units
    // canvas -- table with canvas to be modified   
    // random position
    position: (first 1?1+exec max x from canvas;first 1?1+exec max y from canvas);
    // find a circle around the position and invert colour
    :update colour:colour*neg[1] from canvas where
        diameter>=sqrt[xexp[x-position[0];2]+xexp[y-position[1];2]];
 };

.quantQ.art.placeSquare: {[side; canvas]
    // diameter -- length of side of a square
    // canvas -- table with canvas to be modified   
    // random position
    position: (first 1?1+exec max x from canvas;first 1?1+exec max y from canvas);
    // find a circle around the position and invert colour
    :update colour:colour*neg[1] from canvas where
        (x-position[0])<=0.5*side,
        (x-position[0])>=-0.5*side,
        (y-position[1])<=0.5*side,
        (y-position[1])>=-0.5*side;
 };

.quantQ.art.pivotColours:{[canvas]
    // canvas -- table with canvas and patterns
    :update (-1)^colour1, (-1)^colour2 from
        (select x, colour1:y from canvas where colour=1) uj
        (select x, colour2:y from canvas where colour=-1);
 };



.quantQ.art.placeCircleGradient: {[diameter; canvas]
    // diameter -- diamater of the circle in units
    // canvas -- table with canvas to be modified   
    // random position, position is float
    position: (sqrt first 1? t*t:1+exec max x from canvas;sqrt first 1?t*t:1+exec max y from canvas);
    // find a circle around the position and invert the colour
    :update colour:colour*neg[1] from canvas where
        diameter>=sqrt[xexp[x-position[0];2]+xexp[y-position[1];2]];
 };

