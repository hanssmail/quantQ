.quantQ.dtw.cost:{[a;b]
 dtw:(1+count b)#/:((1+count a)#0w); 
 dtw[0;0]:0f;
 i:1;
 j:1;
 while[i<=count a;
  while[j<=count b;
   cost:abs a[i-1]-b[j-1];
   dtw[i;j]:cost + min(dtw[i-1;j];dtw[i;j-1];dtw[i-1;j-1]);
   j:j+1;
   ];
  i:i+1;
  j:1;
  ];
 dtw}


.quantQ.dtw.optimalWarpingPath:{[a;b]
 dtw:(1+count b)#/:((1+count a)#0w); 
 dtw[0;0]:0f;
 i:1;
 j:1;
 while[i<=count a;
  while[j<=count b;
   cost:abs a[i-1]-b[j-1];
   dtw[i;j]:cost + min(dtw[i-1;j];dtw[i;j-1];dtw[i-1;j-1]);
   j:j+1;
   ];
  i:i+1;
  j:1;
  ];
 dtw[(count a);(count b)]}


.quantQ.dtw.optimalWarpingPathWindow:{[a;b;w]
 dtw:(1+count b)#/:((1+count a)#0w); 
 dtw[0;0]:0f;
 i:1;
 w:max(w, abs (count a)-(count b));
 while[i<=count a;
  j:max(1;i-w);
  while[j< min(count b;i+w);   
   dtw[i;j]:0f;
   j:j+1;
   ];
  i:i+1;
  ]; 
 i:1;
 while[i<=count a;
  j:max(1; i-w);
  while[j<=min(count b;i+w);
   cost:abs a[i-1]-b[j-1];
   dtw[i;j]:cost + min(dtw[i-1;j];dtw[i;j-1];dtw[i-1;j-1]);
   j:j+1;
   ];
  i:i+1;
  ];
 dtw[(count a);(count b)]}


.quantQ.dtw.warpingPath:{[a;b]
 dtw:(1+count b)#/:((1+count a)#0w); 
 dtw[0;0]:0f;
 x:count a; 
 y:count b; 
 i:1;
 j:1;
 while[i<=count a;
  while[j<=count b;
   cost:abs a[i-1]-b[j-1];
   dtw[i;j]:cost + min(dtw[i-1;j];dtw[i;j-1];dtw[i-1;j-1]);
   j:j+1;
   ];
  i:i+1;
  j:1;
  ];
 i:count a;
 j:count b;
 while[(i > 0) and (j>0);
  $[i=0;
   j:j-1;
   $[j=0;
    i:i-1;
    [
     $[dtw[i-1;j]=min(dtw[i-1;j];dtw[i;j-1];dtw[i-1;j-1]);
      i:i-1;
      [
       $[dtw[i;j-1]=min(dtw[i-1;j];dtw[i;j-1];dtw[i-1;j-1]);
	j:j-1;
	[
	 i:i-1;
	 j:j-1;
	 ]
	];
       ]
      ];
     if[(i <> 0) and (j <> 0);
      x:x,i;
      y:y,j;
      ];
     ]
    ]
   ];
  ];
 (reverse x;reverse y)
 }
