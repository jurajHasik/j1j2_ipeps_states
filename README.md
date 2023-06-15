# Infinite PEPS states for J1-J2 model on a square lattice

This repository contains two datasets 

* **single-site_pg-C4v-A1**
* **single-site_pg-C4v-A1_internal_U1**

representing translationary invariant infinite PEPS states defined by a single rank-5 tensor *A<sup>s</sup><sub>uldr</sub>*.
The index *s* runs over two states of spin-1/2 degree of freedom and the virtual indices *u,l,d,r*, for up, left, down and right directions on square lattice, run over values 
from 0 to *D*-1 with *D* being the bond dimension.

This tensor transforms as A<sub>1</sub> irreducible representation of C4<sub>v</sub> point group (symmetry group of square lattice) 
under 90° rotations and reflections represented by permutations of virtual indices

*A<sup>s</sup><sub>uldr</sub>*=*A<sup>s</sup><sub>ldru</sub>*=*A<sup>s</sup><sub>ruld</sub>*=*A<sup>s</sup><sub>urdl</sub>*=*A<sup>s</sup><sub>dlur</sub>*

The tensors in the dataset **single-site_pg-C4v-A1_internal_U1** also possess U(1) symmetry. The states are stored in plain text format (JSON).