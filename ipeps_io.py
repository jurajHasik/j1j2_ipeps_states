import json
import logging
import numpy as np
import argparse
logger = logging.getLogger("ipeps_io")


parser= argparse.ArgumentParser(description='',allow_abbrev=False)
# additional model-dependent arguments
parser.add_argument("--instate", type=str, default=None, help="state to parse")
parser.add_argument("--out",type=str, default=None, help="output file name")
parser.add_argument("--format", type=str, default="npz", help="desired format",\
    choices=["npz","mat","npz_blocks","mat_blocks"])
args, unknown_args= parser.parse_known_args()


def load_from_pepstorch_json_dense(filename)->np.ndarray:
    r"""
    """
    with open(filename) as f:
        state = json.load(f)

        if "sites" in state:
            logger.info("Processed format detected: peps-torch single dense tensor")
            assert len(state["sites"]) == 1, "Multisite dense not yet implemented"
            site = state["sites"][0]
            dims = site.get("dims", None)
            if dims is None:
                assert (
                    "auxDim" in site and "physDim" in site
                ), "Missing dims or auxDim and physDim fields"
                dims = [site["physDim"]] + [site["auxDim"]] * 4
            A = (
                np.zeros(dims, dtype=np.complex128)
                if len(site["entries"][0].split()) > 6
                else np.zeros(dims, dtype=np.float64)
            )
            for entry in site["entries"]:
                if len(entry.split()) == 7:
                    # Complex
                    A[tuple(int(i) for i in entry.split()[:5])] = float(
                        entry.split()[5]
                    ) + 1j * float(entry.split()[6])
                else:
                    # Real
                    A[tuple(int(i) for i in entry.split()[:5])] = float(
                        entry.split()[5]
                    )
        else:
            for _key_elem_ts in ["su2_tensors", "sym_tensors", "elem_tensors", None]:
                if _key_elem_ts in state.keys():
                    break
            assert not _key_elem_ts is None, "Missing elementary tensors"
            logger.info(
                "Processed format detected: peps-torch linear combination of dense tensor"
            )
            elem_ts = state[_key_elem_ts]
            coeffs = state["coeffs"][0]["entries"]
            # get physical dimension and auxiliary bond dimension from (first) elementary tensor
            dtype = (
                "float64" if not "dtype" in elem_ts[0].keys() else elem_ts[0]["dtype"]
            )
            assert dtype == "float64", "Unexpected dtype"
            pd = elem_ts[0]["physDim"]
            ad = elem_ts[0]["auxDim"]
            A = np.zeros((pd, ad, ad, ad, ad), dtype=dtype)

            for elem_t, coeff in zip(elem_ts, coeffs):
                c = float(coeff.split()[1])
                for entry in elem_t["entries"]:
                    A[tuple(int(i) for i in entry.split()[:5])] = c * float(
                        entry.split()[5]
                    )
        f.close()

    return A


def load_from_pepstorch_json_blocksparse(filename)->dict[tuple[int],np.ndarray]:
    r"""
    """
    with open(filename) as j:
        raw_state = json.load(j)

        assert "total_u1_charge" in raw_state,"Missing total charge data"
        assert "u1_charges" in raw_state,"Missing u1 charges data"
        tot_charge= raw_state["total_u1_charge"]
        charges= raw_state["u1_charges"]
        logger.info(f"charges {charges} total_charge {tot_charge}")

        # read the list of considered U(1)-symmetric tensors
        assert "elem_tensors" in raw_state,"Missing elementary tensors"
        sym_tensor_key= "elem_tensors"

        sym_tensors=[]
        for symt in raw_state[sym_tensor_key]:
            meta=dict({"meta": symt["meta"]})
            dims=[symt["physDim"]]+[symt["auxDim"]]*4
            
            sparse_rep=[]
            sparse_rep_oc=[]
            for elem in symt["entries"]:
                tokens= elem.split(' ')
                inds=tuple([int(i) for i in tokens[0:5]])
                sparse_rep.append(\
                    ( [(charges[:2][inds[0]],inds[0])]+[(charges[2:][inds[i]], inds[i])\
                        for i in range(1,5)] , float(tokens[5]) )
                )
            sym_tensors.append((meta,sparse_rep,sparse_rep_oc))

        # Loop over non-equivalent tensor,coeffs pairs in the unit cell
        coeffs={}
        for ts in raw_state["map"]:
            coord = (ts["x"],ts["y"])

            # find the corresponding tensor of coeffs (and its elements) 
            # identified by "siteId" in the "sites" list
            t = None
            for s in raw_state["coeffs"]:
                if s["siteId"] == ts["siteId"]:
                    t = s
            if t == None:
                raise Exception("Tensor with siteId: "+ts["sideId"]+" NOT FOUND in \"sites\"") 

            X= np.zeros(t["numEntries"])

            # 1) fill the tensor with elements from the list "entries"
            # which list the coefficients in the following
            # notation: Dimensions are indexed starting from 0
            # 
            # index (integer) of coeff, (float) Re, Im  
            for entry in t["entries"]:
                tokens = entry.split()
                X[int(tokens[0])]=float(tokens[1])
            coeffs[coord]=X

    # split charges into physical and auxiliary
    oc_p= charges[:2]
    c_a= charges[2:]
    
    # the initial state creates association index_value->charge
    # 0) we need to sort the charges to create charged sectors with D>1
    #    e.g. D=7 aux-charges (0, 2, -2, 0, 2, -2, 2) -> (-2,-2,0,0,2,2,2) <=> (-2, D=2), (0, D=2), (2, D=3) 
    oc_a= sorted(c_a)
    oc_d= {k: oc_a.count(k) for k in set(oc_a)}

    # 1) we need to map the index_values from unsorted charges to index_values within sectors
    # 1a) while sorting the charges, sort the index_values occordingly
    #    e.g. D=7 (0, 2, -2, 0, 2, -2, 2)->(-2, -2, 0, 0, 2, 2, 2)
    #             [0,1,2,3,4,5,6]        ->[ 2,  5, 0, 3, 1, 4, 6]
    oc_a_i= sorted(range(len(c_a)), key=c_a.__getitem__)

    # 1b) now map the sorted index values according to charge sectors
    i0= 0
    c0= oc_a[0]
    oc_a_si= []
    a_map= dict() # maps original index_value into (charge, sorted_index_value)
    for i in range(len(oc_a)):
        if c0 != oc_a[i]:
            c0= oc_a[i]
            i0= 0
        oc_a_si.append(i0)
        a_map[oc_a_i[i]]= (oc_a[i], i0)
        i0+=1

    # 2) build blocks
    blocks= dict()
    for i,T in enumerate(sym_tensors):
        x= coeffs[(0,0)][i]
        for elem in T[1]: 
            # split list[(charge, index_value)] into list[tuple(charges), tuple(index_values)]
            c, iv= tuple(zip(*elem[0]))
            # check if the charged block (key) exists in blocks
            if c not in blocks:
                blocks[c]= np.zeros([1]+[oc_d[_c] for _c in c[1:]])
                # blocks[c]= torch.zeros([1]+[oc_d[_c] for _c in c[1:]],\
                #     dtype=cfg.global_args.dtype, device=cfg.global_args.device)
                logger.debug(f"Creating block c={c} of D={blocks[c].shape}")
            # map dense index_values to block index_values (physical dimension has always size 1)
            iv_b= tuple([0]+[a_map[v][1] for v in iv[1:]])
            if T[0]["meta"]["pg"]=="A_2":
                assert np.iscomplexobj(blocks[c]), "pg A_2 requires complex dtype"
                blocks[c][iv_b]+= 1.0j * x*elem[1]
            else:
                blocks[c][iv_b]+= x*elem[1] 
            logger.debug(f"elem {c},{iv} -> {iv_b} val {x}*{elem[1]} -> {x*elem[1]}")

    return blocks


if __name__=='__main__':
    if len(unknown_args)>0:
        print("args not recognized: "+str(unknown_args))
        raise Exception("Unknown command line arguments")
    if args.format=="npz":
        outf= args.out if not (args.out is None) else "A.npz"
        np.savez(outf, A=load_from_pepstorch_json_dense(args.instate))
    elif args.format=="mat":
        from scipy.io import savemat
        outf= args.out if not (args.out is None) else "A.mat"
        savemat(outf, {"A": load_from_pepstorch_json_dense(args.instate)})
    elif args.format=="npz_blocks":
        outf= args.out if not (args.out is None) else "A.npz"
        np.savez(outf, **{f"{c}":b for c,b in load_from_pepstorch_json_blocksparse(args.instate).items()})
    elif args.format=="mat_blocks":
        from scipy.io import savemat
        outf= args.out if not (args.out is None) else "A.mat"
        savemat(outf, {f"{c}":b for c,b in load_from_pepstorch_json_blocksparse(args.instate).items()})
    