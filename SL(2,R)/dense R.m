
// The following file contains code to run the function "IsDense", which
// determines if a two-generator subgroup of SL(2,R) is dense

// input A, B, f	                       where A, B are in SL(2, K) and f is an embedding of
//                                       K into the real numbers R 
//                                       [ i.e. if K=Rationals(), then f=Infinity() ]
// output:
//   true, _                              if G=<A,B> is dense in SL(2, R)
//   false, "reason"                      if G is not dense (and a reason why)

// we first load the Magma package "sl2r" from
// http://www.math.rwth-aachen.de/~Markus.Kirschmer/magma/sl2r.html
Attach("sl2r.m");


// we now define the function
IsDense := function (A, B, f)

    assert Type (A) eq AlgMatElt and Type (B) eq AlgMatElt;
    assert Determinant(A) eq 1;
    assert Determinant(B) eq 1;

    // checking if G fixes a point in hyperbolic plane H^2 or on boundary

    if Trace(A*B*A^-1*B^-2) eq 2 then 
      return false, "global fixed point or fixed point on boundary";
    end if;

    // checking if G preserves a pair of boundary points of H^2

    Id := Matrix(2,[1,0,0,1]);
    if A*(B*A*B^-1)*A^-1*(B*A^-1*B^-1) eq Id or (A*B*A^-1)*B*(A*B^-1*A^-1)*B^-1 eq Id then 
      return false, "invariant pair of boundary points";
    end if;

    // checking if G is discrete

    G:= TwoGeneratorSubgroupSL2R(A, B, f);
    if IsDiscrete(G) then 
      return false, "group is discrete";
    end if;

    return true, _;

end function;
