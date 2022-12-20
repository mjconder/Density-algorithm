# Density algorithms

The following contains implementations in Magma of two algorithm described in [1]. It takes as input two elements A and B of either SL(2,R) or SL(2,Q_p) (where Q_p denotes the p-adic numbers) and determines whether or not the subgroup G=<A,B> is dense.

These two algorithms require the discreteness algorithms given in  https://github.com/mjconder/Discrete-algorithm and http://www.math.rwth-aachen.de/~Markus.Kirschmer/magma/sl2r.html respectively

# References

[1] M. J. Conder and J. Schillewaert, Discrete two-generator subgroups of PSL2 over non-archimedean local fields, https://arxiv.org/abs/2208.12404
