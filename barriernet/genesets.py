from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List
from pathlib import Path
import yaml

@dataclass
class GeneSet:
    modules: Dict[str, List[str]]
    axes: Dict[str, dict]
    name: str = "geneset"
    version: int = 1

def load_yaml(path: str) -> GeneSet:
    with open(path, "r") as f:
        obj = yaml.safe_load(f)
    modules = {k: [str(x).upper() for x in v] for k, v in obj.get("modules", {}).items()}
    axes = obj.get("axes", {})
    return GeneSet(modules=modules, axes=axes, name=obj.get("name","geneset"), version=int(obj.get("version",1)))

def default_yaml(kind: str) -> str:
    here = Path(__file__).resolve().parent
    if kind.lower() in ("rna","bulk","spatial"):
        return str(here / "data" / "geneset_rna.yaml")
    if kind.lower() in ("ihc","if","protein","proteomics"):
        return str(here / "data" / "geneset_ihc.yaml")
    raise ValueError("kind must be rna|bulk|spatial or ihc|protein")
