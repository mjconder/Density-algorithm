// The following file contains code to run the function "IsDense", which
// determines if a two-generator subgroup of SL(2,Qp) is dense

// input A, B, p, (optional: precision)	  where p is prime and A, B in SL(2, Q_p) representing
//					  elts of PSL(2,Q_p), and precision=20 unless o/wise specified
// output:
//   true, _                              if G=<A,B> is dense in SL(2, Q_p)
//   false, "reason"                      if G is not dense (and a reason why)
//   "error", "precision too low".        ord(B);if precision argument needs to be increased

// we first load the Magma package "discrete magma"
load "discrete magma.m";


// we now define the function
IsDense := function (A, B, p : precision:=20)
    K := pAdicField(p, precision);

    assert Type (A) eq AlgMatElt and Type (B) eq AlgMatElt;
    assert Determinant(A) in {1};
    assert Determinant(B) in {1};
    for i,j in [1..2] do
      assert A[i][j] in K and B[i][j] in K;
    end for;

    // checking if G is discrete and if not then setting appropriate generators

    a, b := IsDiscrete(A,B,p);
    if a eq true then
      return false, "group is discrete";
    else 
      X := b[1]; Y := b[2];
    end if;

    // checking if G fixes a point, or preserves a point or boundary point of B-T tree
    if TL(Y,p) eq 0 and TL(X*Y,p) eq 0 then
      return false, "global fixed point";
    elif TL(Y,p) eq 0 and TL(X*Y,p) ne 0 then
      return true, _;
    elif Trace(X*Y*X^-1*Y^-1) eq 2 then
      return false, "fixed point on boundary";
    elif (X*Y*X^-1)*Y*(X*Y^-1*X^-1)*Y^-1 eq Matrix(2, [1 , 0, 0, 1]) then
      return false, "invariant pair of boundary points";
    else 
      return true, _;
    end if;
end function;
