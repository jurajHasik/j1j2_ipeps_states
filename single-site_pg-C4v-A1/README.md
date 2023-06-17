## Format

The states in this dataset, **single-site_pg-C4v-A1**, are stored as simple plain text files in JSON format. The field ``"sites"`` is a list holding
a single dictionary, the site tensor *A<sup>s</sup><sub>uldr</sub>*, whose non-zero elements are given in the list ``"entries"``, as shown in the example below

```
"sites": [
        {
            "siteId": "A0",
            "physDim": 2,
            "auxDim": 2,
            "numEntries": 32,
            "entries": [
                "0 0 0 0 0 1.0",
                "0 0 0 0 1 0.017282585382645214",
                "0 0 0 1 0 0.017282585382645214",
                "0 0 0 1 1 -4.472056243694618e-05",
...
```

Each entry is a string composed of 5 integers, the indices of the non-zero element
in order *s,u,l,d,r* corresponding to physical spin-1/2 degree of freedom, up, left, down, and right virtual degrees of freedom, followed by a real number which is the value of the tensor element.