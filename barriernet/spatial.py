from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Optional, List, Tuple
from collections import deque
from dataclasses import dataclass

HEX_STEPS = [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]

def pick_first_existing(obs: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in obs.columns:
            return c
    return None

def ensure_tls_score(obs: pd.DataFrame) -> str:
    for c in ["TLS_niche_score","TLS_core_score","TLS_score_score","TLS_score"]:
        if c in obs.columns:
            return c
    if all([(f"tg_{x}" in obs.columns) for x in ["Bcell","Plasma","DC"]]):
        obs["TLS_core_score"] = (pd.to_numeric(obs["tg_Bcell"], errors="coerce") +
                                 pd.to_numeric(obs["tg_Plasma"], errors="coerce") +
                                 pd.to_numeric(obs["tg_DC"], errors="coerce"))
        return "TLS_core_score"
    raise KeyError("Cannot find TLS score in obs.")

def build_hex_graph(obs: pd.DataFrame) -> Tuple[list[list[int]], np.ndarray]:
    if ("array_row" not in obs.columns) or ("array_col" not in obs.columns):
        raise KeyError("Need obs['array_row'] and obs['array_col'] for Visium hex-grid adjacency.")
    ar = pd.to_numeric(obs["array_row"], errors="coerce").fillna(-1).astype(int).to_numpy()
    ac = pd.to_numeric(obs["array_col"], errors="coerce").fillna(-1).astype(int).to_numpy()

    keep = np.ones(len(obs), dtype=bool)
    if "in_tissue" in obs.columns:
        keep = (pd.to_numeric(obs["in_tissue"], errors="coerce").fillna(0).astype(int).to_numpy() == 1)

    pos2i = {(int(r), int(c)): i for i,(r,c) in enumerate(zip(ar,ac)) if keep[i] and r>=0 and c>=0}
    nbrs = [[] for _ in range(len(obs))]
    for i,(r,c) in enumerate(zip(ar,ac)):
        if not keep[i] or r < 0 or c < 0:
            continue
        for dr,dc in HEX_STEPS:
            j = pos2i.get((r+dr, c+dc))
            if j is not None:
                nbrs[i].append(j)
    return nbrs, keep

def bfs_dist(nbrs: list[list[int]], keep: np.ndarray, sources_mask: np.ndarray) -> np.ndarray:
    dist = np.full(len(nbrs), np.inf)
    q = deque()
    src = np.where(keep & sources_mask)[0]
    for s in src:
        dist[s] = 0
        q.append(s)
    while q:
        u = q.popleft()
        for v in nbrs[u]:
            if dist[v] > dist[u] + 1:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist

def edge_dist(nbrs: list[list[int]], keep: np.ndarray) -> np.ndarray:
    deg = np.array([len(nbrs[i]) if keep[i] else 0 for i in range(len(nbrs))])
    edge = keep & (deg < 6)
    return bfs_dist(nbrs, keep, edge)

def local_maxima_peaks(xy: np.ndarray, score: np.ndarray, top_q: float=0.95, k: int=9) -> Tuple[np.ndarray, float, int]:
    score = np.asarray(score, dtype=float)
    finite = score[np.isfinite(score)]
    thr = float(np.nanquantile(finite, top_q)) if finite.size else np.nan
    cand = np.isfinite(score) & (score >= thr)

    n = xy.shape[0]
    pk = np.zeros(n, dtype=bool)
    if n == 0 or not cand.any():
        return pk, thr, 0

    # Approx kNN local-max using brute distances (ok for Visium spot counts)
    for i in np.where(cand)[0]:
        d = np.sum((xy - xy[i])**2, axis=1)
        nn = np.argpartition(d, kth=min(k, n-1))[:min(k, n)]
        nn = nn[cand[nn]]
        if nn.size and score[i] >= np.nanmax(score[nn]):
            pk[i] = True
    return pk, thr, int(pk.sum())

@dataclass
class SpatialKeys:
    tls: str
    ecm: Optional[str]
    tam: Optional[str]
    stf3: Optional[str]

def infer_spatial_keys(obs: pd.DataFrame) -> SpatialKeys:
    tls = ensure_tls_score(obs)
    ecm = pick_first_existing(obs, ["state_CAF_ECM","prog_CAF_ECM","score_CAF_ECM"])
    tam = pick_first_existing(obs, ["state_SPP1_TAM","state_TAM_SUP","prog_TAM_SUP","score_TAM_SUP"])
    stf3 = pick_first_existing(obs, ["state_STF3_FAP_INHBA","prog_STF3_FAP_INHBA","prog_CAF_STF3_FAP_INHBA"])
    return SpatialKeys(tls=tls, ecm=ecm, tam=tam, stf3=stf3)

def assign_rings(obs: pd.DataFrame,
                 spatial_xy: np.ndarray,
                 tls_key: str,
                 peak_q: float=0.95,
                 peak_k: int=9,
                 far_min_dist: int=4) -> pd.DataFrame:
    tls = pd.to_numeric(obs[tls_key], errors="coerce").to_numpy(dtype=float)
    tls_pk, tls_thr, tls_n = local_maxima_peaks(spatial_xy, tls, top_q=peak_q, k=peak_k)

    nbrs, keep = build_hex_graph(obs)
    dist_edge = edge_dist(nbrs, keep)
    dist_tls = bfs_dist(nbrs, keep, tls_pk & keep)

    ring = np.full(len(obs), "NA", dtype=object)
    ring[keep] = "Other"
    ring[tls_pk & keep] = "Hotspot"
    ring[(keep) & (dist_tls == 1)] = "Ring1"
    ring[(keep) & (dist_tls == 2)] = "Ring2"

    far = (keep) & (dist_tls >= far_min_dist) & (dist_edge > 1)
    if far.sum() < 30:
        far = (keep) & (dist_tls >= far_min_dist)
    ring[far] = "Far"

    return pd.DataFrame({
        "tls_peak": tls_pk.astype(int),
        "tls_peak_thr": tls_thr,
        "tls_peak_n": tls_n,
        "dist_tls_hops": dist_tls,
        "dist_edge_hops": dist_edge,
        "ring": ring
    }, index=obs.index)

def ring_summary(obs: pd.DataFrame, ring_df: pd.DataFrame, keys: SpatialKeys) -> pd.DataFrame:
    df = pd.concat([obs[[keys.tls]].copy(), ring_df], axis=1)

    def add(key: Optional[str], label: str):
        if key is None:
            df[label] = np.nan
        else:
            df[label] = pd.to_numeric(obs[key], errors="coerce")

    add(keys.ecm, "ECM")
    add(keys.tam, "TAM")
    add(keys.stf3, "STF3")

    ring_order = ["Hotspot","Ring1","Ring2","Far"]
    rows=[]
    for r in ring_order:
        sub = df[df["ring"]==r]
        if sub.shape[0]==0:
            continue
        rows.append({
            "ring": r,
            "n": int(sub.shape[0]),
            "TLS_mean": float(np.nanmean(sub[keys.tls].values)),
            "ECM_mean": float(np.nanmean(sub["ECM"].values)),
            "TAM_mean": float(np.nanmean(sub["TAM"].values)),
            "STF3_mean": float(np.nanmean(sub["STF3"].values)),
        })
    out = pd.DataFrame(rows)

    def get(r, col):
        s = out.loc[out["ring"]==r, col]
        return float(s.iloc[0]) if len(s) else np.nan

    for metric, col in [("ECM","ECM_mean"),("TAM","TAM_mean"),("STF3","STF3_mean")]:
        out[f"{metric}_Ring2_minus_Far"] = get("Ring2", col) - get("Far", col)

    return out
