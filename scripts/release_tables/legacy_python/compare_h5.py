#!/usr/bin/env python
"""compare_h5.py -- bit-level comparison of a rebuilt release table against a reference.

For every dataset present in either file: presence, dtype, shape, chunk shape, filter pipeline,
and the VALUES compared bitwise (raw bytes; NaN payloads and -0.0 count), optionally on a row
slice. Mismatching float columns report the count and the maximum difference in ulps. Every
dataset attribute and root attribute is compared by value AND by HDF5 string flavour (fixed vs
variable length, size, padding, character set). Exit status 0 only if everything matches
(attributes named with --ignore-attr are reported but not counted).

Usage:
  python compare_h5.py REFERENCE CANDIDATE [--rows START:STOP] [--ignore-attr NAME ...]
                       [--only-datasets a,b,...] [--json OUT]
"""
import argparse
import json
import sys

import h5py
import numpy as np

STEP = 4_000_000


def attr_flavour(obj, name):
    aid = h5py.h5a.open(obj.id, name.encode())
    t = aid.get_type()
    d = {"class": int(t.get_class()), "shape": list(aid.get_space().shape)}
    if t.get_class() == h5py.h5t.STRING:
        d.update(vlen=bool(t.is_variable_str()), pad=int(t.get_strpad()), cset=int(t.get_cset()))
        if not t.is_variable_str():
            d["size"] = int(t.get_size())
    else:
        d["size"] = int(t.get_size())
    return d


def attrs_equal(a, b):
    try:
        return bool(np.array_equal(np.asarray(a), np.asarray(b)))
    except Exception:
        return a == b


def compare_attrs(oa, ob, where, ignore, out):
    ka, kb = set(oa.attrs.keys()), set(ob.attrs.keys())
    for k in sorted(ka | kb):
        rec = {"where": where, "attr": k}
        if k not in ka or k not in kb:
            rec["problem"] = "missing in " + ("reference" if k not in ka else "candidate")
        else:
            if not attrs_equal(oa.attrs[k], ob.attrs[k]):
                rec["problem"] = "value differs"
            elif attr_flavour(oa, k) != attr_flavour(ob, k):
                rec["problem"] = f"type differs {attr_flavour(oa, k)} vs {attr_flavour(ob, k)}"
        if "problem" in rec:
            rec["ignored"] = k in ignore
            out.append(rec)


def ulps(a, b):
    if a.dtype == np.float64:
        ia, ib = a.view(np.int64), b.view(np.int64)
    else:
        ia, ib = a.view(np.int32).astype(np.int64), b.view(np.int32).astype(np.int64)
    ia = np.where(ia < 0, np.iinfo(np.int64).min - ia, ia)
    ib = np.where(ib < 0, np.iinfo(np.int64).min - ib, ib)
    return np.abs(ia - ib)


def compare_dataset(da, db, rows):
    r = {"dtype": [str(da.dtype), str(db.dtype)], "shape": [list(da.shape), list(db.shape)]}
    ok = da.dtype == db.dtype and da.shape == db.shape
    r["layout_equal"] = (da.chunks == db.chunks and da.compression == db.compression
                         and da.compression_opts == db.compression_opts and da.shuffle == db.shuffle)
    if not r["layout_equal"]:
        r["layout"] = [str((da.chunks, da.compression, da.compression_opts, da.shuffle)),
                       str((db.chunks, db.compression, db.compression_opts, db.shuffle))]
    if not ok:
        r["values_equal"] = False
        return r
    n = da.shape[0]
    lo, hi = rows if rows else (0, n)
    hi = min(hi, n)
    nbad = 0
    maxulp = 0
    first = None
    for s in range(lo, hi, STEP):
        e = min(s + STEP, hi)
        a, b = da[s:e], db[s:e]
        ba = a.view(np.uint8).reshape(a.shape[0], -1)
        bb = b.view(np.uint8).reshape(b.shape[0], -1)
        diff = np.any(ba != bb, axis=1)
        k = int(diff.sum())
        if k:
            nbad += k
            if first is None:
                first = s + int(np.flatnonzero(diff)[0])
            if a.dtype.kind == "f":
                both_nan = np.isnan(a[diff]) & np.isnan(b[diff])
                u = ulps(a[diff][~both_nan], b[diff][~both_nan])
                if u.size:
                    maxulp = max(maxulp, int(u.max()))
    r["rows_compared"] = hi - lo
    r["values_equal"] = nbad == 0
    if nbad:
        r["n_mismatch"] = nbad
        r["first_mismatch_row"] = first
        if da.dtype.kind == "f":
            r["max_ulp"] = maxulp
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("candidate")
    ap.add_argument("--rows", default="")
    ap.add_argument("--ignore-attr", nargs="*", default=[])
    ap.add_argument("--only-datasets", default="")
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    rows = tuple(int(x) for x in a.rows.split(":")) if a.rows else None
    only = set(a.only_datasets.split(",")) if a.only_datasets else None
    res = {"reference": a.reference, "candidate": a.candidate, "rows": a.rows or "all",
           "datasets": {}, "attr_problems": []}
    allok = True
    with h5py.File(a.reference, "r") as fa, h5py.File(a.candidate, "r") as fb:
        def walk(g, prefix=""):
            out = []
            for k in g:
                if isinstance(g[k], h5py.Group):
                    out += walk(g[k], prefix + k + "/")
                else:
                    out.append(prefix + k)
            return out
        ka, kb = set(walk(fa)), set(walk(fb))
        if only:
            ka, kb = ka & only, kb & only
        for k in sorted(ka | kb):
            if k not in ka or k not in kb:
                res["datasets"][k] = {"problem": "missing in " + ("reference" if k not in ka else "candidate")}
                allok = False
                continue
            r = compare_dataset(fa[k], fb[k], rows)
            res["datasets"][k] = r
            allok &= r["values_equal"] and r["layout_equal"]
            compare_attrs(fa[k], fb[k], k, set(a.ignore_attr), res["attr_problems"])
        if not only:
            compare_attrs(fa, fb, "/", set(a.ignore_attr), res["attr_problems"])
    allok &= not any(not p["ignored"] for p in res["attr_problems"])
    res["all_equal"] = bool(allok)
    for k, r in res["datasets"].items():
        flag = "OK " if r.get("values_equal") and r.get("layout_equal", False) else "BAD"
        extra = {kk: v for kk, v in r.items() if kk in ("n_mismatch", "max_ulp", "first_mismatch_row", "problem", "layout", "dtype")
                 and not (kk == "dtype" and r.get("values_equal"))}
        print(f"{flag} {k:28s} {extra if extra else ''}")
    for p in res["attr_problems"]:
        print(("IGN " if p["ignored"] else "ATTR") + f" {p['where']} :: {p['attr']} : {p['problem']}")
    print("ALL EQUAL" if allok else "DIFFERENCES FOUND")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
