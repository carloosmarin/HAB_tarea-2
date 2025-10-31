#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tarea 2: Network Propagation
Autor: Carlos Marín Martínez
Asignatura: Herramientas de Bioinformática Avanzada
Profesor: [nombre del profesor, si quieres incluirlo]

Descripción general
-------------------
Este script implementa un ejemplo de *propagación en redes* utilizando un
algoritmo tipo **GUILD** basado en *Random Walk with Restart (RWR)*.

El objetivo es difundir la influencia de un conjunto de genes semilla
(ENO1, PGK1, HK2) sobre una red de interacción de proteínas y obtener
un ranking de nodos según su relevancia o cercanía funcional a dichas semillas.

El código se inspira en el estilo y estructura de `process_STRING.py`,
siguiendo buenas prácticas de legibilidad, comentarios y trazabilidad.

Entradas
--------
- Archivo de red (`--network`): puede estar en uno de los siguientes formatos:
  1. STRING filtrado (`.tsv`) con columnas `protein1_hugo`, `protein2_hugo`, `combined_score`
  2. Formato “guild” (`.txt`) con tres columnas: nodo1, peso, nodo2
  3. Formato “diamond” (`.txt`) con dos columnas separadas por comas

- Archivo de genes semilla (`--seeds`): lista de genes, uno por línea.

Salidas
--------
- Archivo de resultados (`--out`): tabla TSV con columnas:
    node    score
  donde “score” representa la probabilidad estacionaria del nodo tras la propagación.

Uso
---
python scripts/tu_script.py \
  --algo guild \
  --network data/string_network_filtered_hugo-400.tsv \
  --seeds data/genes_seed.txt \
  --out results/guild_scores.tsv
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm


# ============================================================
# UTILIDADES DE CARGA DE DATOS
# ============================================================

def load_seeds(path: str, graph_nodes: set) -> List[str]:
    """Carga los genes semilla desde un archivo de texto (uno por línea)."""
    seeds_raw = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            gene = line.strip()
            if gene:
                seeds_raw.append(gene)

    seeds = [s for s in seeds_raw if s in graph_nodes]
    missing = [s for s in seeds_raw if s not in graph_nodes]

    if not seeds:
        print("[ADVERTENCIA] Ninguna semilla coincide con nodos de la red.", file=sys.stderr)
    if missing:
        print(f"[INFO] Semillas no encontradas en la red (se ignoran): {', '.join(missing)}", file=sys.stderr)

    return seeds


def _try_load_string_tsv(df: pd.DataFrame) -> Tuple[bool, pd.DataFrame]:
    """Detecta si el DataFrame tiene columnas del formato STRING y las normaliza."""
    cols_lower = {c.lower() for c in df.columns}
    if {"protein1_hugo", "protein2_hugo"}.issubset(cols_lower):
        df.columns = [c.lower() for c in df.columns]
        if "combined_score" not in df.columns:
            df["combined_score"] = 1.0
        df["weight"] = df["combined_score"].astype(float)
        return True, df[["protein1_hugo", "protein2_hugo", "weight"]]
    return False, df


def load_network_auto(path: str) -> nx.Graph:
    """
    Carga automáticamente una red desde distintos formatos (STRING, GUILD o DIAMOnD).
    Devuelve un grafo no dirigido y ponderado.
    """
    p = Path(path)
    G = nx.Graph()

    # --- INTENTO 1: STRING TSV ---
    try:
        df = pd.read_csv(p, sep="\t", dtype=str)
        ok, df2 = _try_load_string_tsv(df)
        if ok:
            for _, r in df2.iterrows():
                u, v, w = str(r["protein1_hugo"]), str(r["protein2_hugo"]), float(r["weight"])
                if u != v:
                    if G.has_edge(u, v):
                        G[u][v]["weight"] = max(G[u][v]["weight"], w)
                    else:
                        G.add_edge(u, v, weight=w)
            return G
    except Exception:
        pass  # si no tiene cabecera, probamos otros formatos

    # --- INTENTO 2: GUILD-like (espacios) ---
    try:
        with open(p, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 3:
                    u, w, v = parts
                    try:
                        w = float(w)
                    except ValueError:
                        w = 1.0
                    if u != v:
                        G.add_edge(u, v, weight=w)
                else:
                    G.clear()
                    break
        if G.number_of_edges() > 0:
            return G
    except Exception:
        pass

    # --- INTENTO 3: DIAMOnD-like (CSV con dos columnas) ---
    try:
        df = pd.read_csv(p, sep=",", header=None, names=["u", "v"], dtype=str)
        for _, r in df.iterrows():
            u, v = str(r["u"]), str(r["v"])
            if u != v:
                G.add_edge(u, v, weight=1.0)
        return G
    except Exception as e:
        raise RuntimeError(f"No se pudo cargar la red desde {path}. Error: {e}")


# ============================================================
# ALGORITMO: Random Walk with Restart (GUILD-like)
# ============================================================

def run_rwr_guild(G: nx.Graph,
                  seeds: List[str],
                  restart: float = 0.5,
                  max_iter: int = 50,
                  tol: float = 1e-6) -> Dict[str, float]:
    """
    Implementación sencilla de Random Walk with Restart (RWR).

    Fórmula:
        p_{t+1} = (1 - α) * W_norm * p_t + α * p0
    donde:
      - α es el parámetro de reinicio ("restart")
      - p0 es el vector inicial (1/n para cada semilla)
      - W_norm es la matriz de transición normalizada por grado

    Devuelve un diccionario {nodo: score}.
    """
    if not seeds:
        raise ValueError("No hay semillas válidas presentes en la red.")

    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    p0 = np.zeros(n)
    for s in seeds:
        p0[idx[s]] = 1.0
    p0 /= p0.sum()

    degrees = np.array([G.degree(n) for n in nodes], dtype=float)
    degrees[degrees == 0.0] = 1.0

    p = p0.copy()
    for _ in tqdm(range(max_iter), desc="Iterando RWR"):
        new_p = np.zeros_like(p)
        for u in nodes:
            pu = p[idx[u]] / degrees[idx[u]]
            if pu == 0.0:
                continue
            for v in G.neighbors(u):
                new_p[idx[v]] += pu
        new_p = (1 - restart) * new_p + restart * p0
        new_p /= new_p.sum()
        if np.linalg.norm(new_p - p, 1) < tol:
            break
        p = new_p

    return {nodes[i]: float(p[i]) for i in range(n)}


# ============================================================
# (OPCIONAL) STUB: DIAMOnD-lite
# ============================================================

def run_diamond_stub(G: nx.Graph, seeds: List[str], k: int = 200) -> pd.DataFrame:
    """Placeholder opcional para DIAMOnD."""
    print("[AVISO] Algoritmo DIAMOnD-lite aún no implementado.", file=sys.stderr)
    return pd.DataFrame(columns=["node", "pvalue", "step"])


# ============================================================
# PROGRAMA PRINCIPAL (CLI)
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Propagación en redes (GUILD-like RWR)")
    parser.add_argument("--algo", choices=["guild", "diamond"], default="guild",
                        help="Algoritmo a ejecutar (por defecto: guild)")
    parser.add_argument("--network", required=True, help="Ruta al archivo de red")
    parser.add_argument("--seeds", default="data/genes_seed.txt", help="Archivo con genes semilla")
    parser.add_argument("--out", required=True, help="Archivo de salida (TSV)")
    parser.add_argument("--restart", type=float, default=0.5, help="Parámetro α del reinicio")
    parser.add_argument("--max-iter", type=int, default=50, help="Número máximo de iteraciones")
    parser.add_argument("--tol", type=float, default=1e-6, help="Tolerancia de convergencia")
    parser.add_argument("--k", type=int, default=200, help="Número de nodos a añadir (DIAMOnD)")
    args = parser.parse_args()

    # --- 1. Cargar red ---
    print(f"[INFO] Cargando red desde: {args.network}")
    G = load_network_auto(args.network)
    print(f"[INFO] Nodos: {G.number_of_nodes()} | Aristas: {G.number_of_edges()}")

    # --- 2. Cargar semillas ---
    seeds = load_seeds(args.seeds, set(G.nodes()))
    print(f"[INFO] Semillas válidas: {len(seeds)} -> {seeds}")

    # --- 3. Ejecutar algoritmo ---
    if args.algo == "guild":
        print("[INFO] Ejecutando propagación tipo GUILD (RWR)...")
        scores = run_rwr_guild(G, seeds, restart=args.restart, max_iter=args.max_iter, tol=args.tol)

        # Guardar resultados
        df = pd.DataFrame(sorted(scores.items(), key=lambda x: x[1], reverse=True),
                          columns=["node", "score"])
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, sep="\t", index=False)
        print(f"[OK] Resultados guardados en {args.out}")
        print(df.head().to_string(index=False))

    elif args.algo == "diamond":
        df = run_diamond_stub(G, seeds, k=args.k)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, sep="\t", index=False)
        print(f"[OK] (Stub) Resultados guardados en {args.out}")

    else:
        raise ValueError("Algoritmo no reconocido.")


if __name__ == "__main__":
    main()
