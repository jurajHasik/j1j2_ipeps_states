# Infinite PEPS states for J1-J2 model on a square lattice

This repository contains two datasets

* **single-site_pg-C4v-A1**
* **single-site_pg-C4v-A1_internal_U1**, basis for [SciPost Phys. 10, 012 (2021)](https://scipost.org/SciPostPhys.10.1.012)
 
representing translationary invariant infinite PEPS states defined by a single rank-5 tensor *A<sup>s</sup><sub>uldr</sub>*.
The index *s* runs over two states of spin-1/2 degree of freedom and the virtual indices *u,l,d,r*, for up, left, down and right directions on square lattice, run over values 
from 0 to *D*-1 with *D* being the bond dimension.

This tensor transforms as A<sub>1</sub> irreducible representation of C4<sub>v</sub> point group (symmetry group of square lattice) 
under 90° rotations and reflections represented by permutations of virtual indices

*A<sup>s</sup><sub>uldr</sub>*=*A<sup>s</sup><sub>ldru</sub>*=*A<sup>s</sup><sub>ruld</sub>*=*A<sup>s</sup><sub>urdl</sub>*=*A<sup>s</sup><sub>dlur</sub>*.

These single-site translationally invariant PEPS are endowed with antiferromagnetic correlations by
applying a unitary rotation *-iσ<sup>y</sup>* to spins[^1] on every sublattice-B site of a square lattice under bipartite tiling

*B<sup>s</sup><sub>uldr</sub>*=*-iσ<sup>y</sup><sub>sk</sub>A<sup>k</sup><sub>uldr</sub>*.

The tensors in the dataset **single-site_pg-C4v-A1_internal_U1** also possess U(1) symmetry [see [SciPost Phys. 10, 012 (2021)](https://scipost.org/SciPostPhys.10.1.012)]. 
The states are stored in plain text format (JSON).

# Observables

Each state is accompanied by corresponding *.dat file containing the value of selected observables
as evaluated by corner transfer matrix algorithm implemented in [``peps-torch``](https://github.com/jurajHasik/peps-torch) library.
For a set of environment bond dimensions χ, those are
* energy per site of J<sub>1</sub>-J<sub>2</sub> model[^2]
* on-site magnetization m=|⟨S⟩|, with S=(S<sup>z</sup>,S<sup>x</sup>,S<sup>y</sup>) the vector of spin-1/2 operators
* leading eigenvalues λ<sub>0</sub>, λ<sub>1</sub>,... of (width-1) transfer matrix. The spectrum is normalized (λ<sub>0</sub>=1)
  and the leading correlation length can be obtained as ξ = -1/ln(λ<sub>1</sub>)

# Reading and exporting states

To parse the states use the Python script ``ipeps_io.py`` which can export 
the **dense tensor** *A<sup>s</sup><sub>uldr</sub>* to either NumPy's *.npz format or MATLAB's *.mat format (requires SciPy)

```
python ipeps_io.py --instate path/to/json/file --format mat [--out optional/name/for/exported/file]
python ipeps_io.py --instate path/to/json/file --format npz [--out optional/name/for/exported/file]
```

Or access the (NumPy) tensor directly in the interactive mode

```
>>> from ipeps_io import load_from_pepstorch_json_dense
>>> A=load_from_pepstorch_json_dense("single-site_pg-C4v-A1/j20.0/state_1s_A1_j20.0_D3_chi_opt108.json")
>>> type(A)
<class 'numpy.ndarray'>
>>> A.shape
(2, 3, 3, 3, 3)
>>> A[0,0,0,:,:]
array([[-0.49333601,  0.1460243 ,  0.03569442],
       [ 0.1460243 , -0.26334648, -0.01241685],
       [ 0.03569442, -0.01241685,  0.40762404]])
``` 

It is also possible to use ``ipeps_io.py`` script to export U(1)-symmetric states of **single-site_pg-C4v-A1_internal_U1** dataset to commonly adopted
block-sparse form. See ``single-site_pg-C4v-A1_internal_U1/README.md``.

[^1]: In practice it is more convenient to instead rotate either operators or reduced density matrices.
[^2]: both nearest- and next-nearest neighbour spin exchange terms are evaluated in 2x2 patch embedded in CTM environment