 function k = rowwiseLast(A)
 %Finds locations k(i) of first non-zero value in each row A(i,:), with k(i)=NaN if 
 %the entire row is zero.
      m=size(A,2);
      [val,loc] = max(  logical(A),  [],2);
      k=m+1-loc;
       k(val==0)=nan;
 end