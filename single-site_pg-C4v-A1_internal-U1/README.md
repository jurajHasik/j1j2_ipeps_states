## Format

The states in this dataset, **single-site_pg-C4v-A1_internal-U1**, are stored as simple plain text files in JSON format. 
The site tensor *A<sup>s</sup><sub>uldr</sub>* is defined as a linear combination 

*A<sup>s</sup><sub>uldr</sub>*=*Σ<sub>i</sub> c<sub>i</sub>(e<sub>i</sub>)<sup>s</sup><sub>uldr</sub>*,

of elementary tensors *e<sub>i</sub>* stored in ``"elem_tensors"`` array and coefficients *c<sub>i</sub>*
stored in ``"coeffs"`` array. 

The ``"coeffs"`` array holds a single dictionary, with array ``"entries"`` containing
the numerical values of the coefficients.
The ``"elem_tensors"`` list holds different elementary tensors, different representatives of A<sub>1</sub> irreducible representation of C4<sub>v</sub> point group. 
These tensors are given in coordinate-sparse form, whose non-zero entries are listed in array ``"entries"``. 
Each entry is a string composed of 5 integers, the indices of the non-zero element
in order *s,u,l,d,r* corresponding to physical spin-1/2 degree of freedom, up, left, down, and right virtual degrees of freedom, followed by a real number which is the value of the tensor element.

An example of such format is given below

```
"total_u1_charge": 1,
"u1_charges": [1,-1,0,2,0],
"coeffs": [
        {
            "siteId": "A0",
            "numEntries": 12,
            "entries": [
                "0 1.0",
                "1 0.3566495454318669",
                "2 -0.01732843774078798",
                "3 0.4227036534709078",
                "4 0.5107228093236208",
                "5 0.5809497493044713",
                ...
            ]
        }
"elem_tensors": [
        {
            "physDim": 2,
            "auxDim": 3,
            "numEntries": 4,
            "entries": [
                "1 1 2 2 2 0.5",
                "1 2 1 2 2 0.5",
                "1 2 2 1 2 0.5",
                "1 2 2 2 1 0.5"
            ]
        },
        ...
```

### U(1) charges

The site tensor *A<sup>s</sup><sub>uldr</sub>* in this dataset possess U(1) symmetry, see [SciPost Phys. 10, 012 (2021)](https://scipost.org/SciPostPhys.10.1.012) for detailed 
explanation. The first two charges specified in array ``"u1_charges"``
are assigned to spin up (*s*=0) and spin down (*s*=1) state respectively of physical index. The remaining charges, whose number is equal to the dimension of virtual
degrees of freedom, are assigned to their virtual states. The total charge of all non-zero entries adds up to ``"total_u1_charge"``. Hence,
in general for any non-zero element *A<sup>s</sup><sub>uldr</sub>* the following charge conservation rule holds[^1]
```
c_physical[s] + c_virtual[u] + c_virtual[l] + c_virtual[d] + c_virtual[r] = total_u1_charge,
```
where 
```
c_physical, c_virtual = u1_charges[:2], u1_charges[2:].
```

In the example above, the 0<sup>th</sup> and 2<sup>nd</sup> virtual state have charge 0, while the 1<sup>st</sup> virtual
state is charged with charge 2. The total charge of all non-zero entries adds up to ``"total_u1_charge"``, which in the example above is 1.

### Exporting U(1) states in block-sparse form

The site tensors *A<sup>s</sup><sub>uldr</sub>* from this dataset can be also exported in the typical block-sparse format, i.e.
as a dictionary of blocks (dense tensors) indexed by tuples of charges[^2]. For any block, all its elements have the same charges, equal to the block's index.
To export these U(1) symmetric states in block-sparse format, use ``ipeps_io.py`` script
specifying the desired format as either NumPy or MATLAB

```
python ipeps_io.py --instate path/to/state --format npz_blocks --out optional/output/filename
python ipeps_io.py --instate path/to/state --format mat_blocks --out optional/output/filename
```

Or access the block-sparse format in interactive form

```
>>> from ipeps_io import load_from_pepstorch_json_blocksparse
>>> A=load_from_pepstorch_json_blocksparse("single-site_pg-C4v-A1_internal-U1/j20.0/state_1s_A1_U1B_j20.0_D3_chi_opt72.json")
>>> A.keys()
dict_keys([(-1, 2, 0, 0, 0), (-1, 0, 2, 0, 0), (-1, 0, 0, 2, 0), (-1, 0, 0, 0, 2), (1, 0, 0, 0, 0)])
>>> A[(-1, 2, 0, 0, 0)].shape
(1, 1, 2, 2, 2)
>>> A[(-1, 2, 0, 0, 0)]
array([[[[[0.11581714, 0.11530105],
          [0.12885121, 0.19440972]],

         [[0.11530105, 0.16758501],
          [0.19440972, 0.5       ]]]]])
```

[^1]: this convention corresponds to all indices of *A<sup>s</sup><sub>uldr</sub>* having the same signature in sense of formalism developed in works introducing symmetric tensor networks.
[^2]: in the output files the indices are strings made of tuples of charges. This is imposed by the underlying NumPy or SciPy implementations.