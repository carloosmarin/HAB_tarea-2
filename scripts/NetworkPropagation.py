#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tarea 2: Network Propagation
Autor: Carlos Marín Martínez

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
import math


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
                  tol: float = 1e-6,
                  use_weights: bool = False) -> Dict[str, float]:
    """
    Random Walk with Restart (RWR).
    Si use_weights=True, el reparto a vecinos se hace proporcional al peso de la arista.
    """
    if not seeds:
        raise ValueError("No hay semillas válidas presentes en la red.")

    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    # p0: prob. uniforme sobre semillas
    p0 = np.zeros(n, dtype=float)
    for s in seeds:
        p0[idx[s]] = 1.0
    p0 /= p0.sum()

    # Denominador por nodo: suma de pesos (o grado si no usamos pesos)
    den = np.zeros(n, dtype=float)
    for u in nodes:
        if use_weights:
            s = 0.0
            for _, data in G[u].items():
                s += float(data.get("weight", 1.0))
            den[idx[u]] = s if s > 0 else 1.0
        else:
            d = G.degree(u)
            den[idx[u]] = float(d) if d > 0 else 1.0

    p = p0.copy()
    for _ in tqdm(range(max_iter), desc="Iterando RWR"):
        new_p = np.zeros_like(p)
        if use_weights:
            # Reparto proporcional a pesos
            for u in nodes:
                pu = p[idx[u]] / den[idx[u]]
                if pu == 0.0:
                    continue
                for v, data in G[u].items():
                    w = float(data.get("weight", 1.0))
                    new_p[idx[v]] += pu * w
        else:
            # Reparto uniforme 1/deg
            for u in nodes:
                pu = p[idx[u]] / den[idx[u]]
                if pu == 0.0:
                    continue
                for v in G.neighbors(u):
                    new_p[idx[v]] += pu

        # Restart
        new_p = (1.0 - restart) * new_p + restart * p0

        s = new_p.sum()
        if s > 0:
            new_p /= s

        if np.linalg.norm(new_p - p, 1) < tol:
            p = new_p
            break
        p = new_p

    return {nodes[i]: float(p[i]) for i in range(n)}

# ============================================================
# (OPCIONAL) STUB: DIAMOnD-lite
# ============================================================

# ---------- Utilidades para DIAMOnD-lite ----------
import math

def _logC(n: int, k: int) -> float:
    """log( combinatoria(n, k) ) usando lgamma para estabilidad numérica."""
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def _hypergeom_sf_nodes(N: int, K: int, n: int, x: int) -> float:
    """
    P(X >= x) para X ~ Hypergeom(N, K, n) con sumatorio estable en log-espacio.
    Interpretación 'sobre nodos':
      - N = población total (número de nodos - 1)
      - K = número de nodos "exitosos" en la población (|S|)
      - n = extracciones (grado del nodo candidato)
      - x = éxitos observados (vecinos del candidato que están en S)
    """
    if x <= 0:
        return 1.0
    xmax = min(n, K)
    if x > xmax:
        return 0.0
    log_den = _logC(N, n)
    s = 0.0
    for i in range(x, xmax + 1):
        term_log = _logC(K, i) + _logC(N - K, n - i) - log_den
        s += math.exp(term_log)
    return min(1.0, max(0.0, s))


# ============================================================
# DIAMOnD-lite: expansión iterativa por significancia hipergeométrica
# ============================================================

def run_diamond_lite(G: nx.Graph, seeds: List[str], k: int = 200) -> pd.DataFrame:
    """
    Implementación ligera de DIAMOnD:
      - S = conjunto actual (inicia en semillas)
      - En cada paso, evalúa para cada candidato v∉S la conectividad k_s(v)
        (nº de vecinos de v que están en S).
      - Calcula p-valor hipergeométrico sobre NODOS:
          N = |V| - 1
          K = |S|
          n = deg(v)
          x = k_s(v)
        p = P[X >= x] con X ~ Hypergeom(N, K, n).
      - Añade el nodo con menor p-valor y repite hasta K pasos.

    Devuelve un DataFrame con columnas: node, pvalue, step
    (solo nodos añadidos; las semillas no aparecen en la tabla).
    """
    if not seeds:
        raise ValueError("No hay semillas válidas presentes en la red.")

    S = set(s for s in seeds if s in G)
    if not S:
        raise ValueError("Ninguna semilla está presente en la red.")

    nodes_all = list(G.nodes())
    Npop = len(nodes_all) - 1  # población para hipergeométrica "sobre nodos"

    added_nodes = []
    current_step = 0

    # Precompute
    neigh = {u: set(G.neighbors(u)) for u in nodes_all}
    deg = {u: len(neigh[u]) for u in nodes_all}

    while current_step < k:
        # Solo consideramos candidatos con al menos 1 vecino en S (rápido)
        candidates = [v for v in nodes_all if v not in S and len(neigh[v] & S) > 0]
        if not candidates:
            break

        best_node, best_p = None, 1.0

        for v in candidates:
            dv = deg[v]
            xs = len(neigh[v] & S)  # éxitos observados
            if dv == 0 or xs <= 0:
                pval = 1.0
            else:
                pval = _hypergeom_sf_nodes(Npop, len(S), dv, xs)
            if pval < best_p:
                best_p, best_node = pval, v

        if best_node is None:
            break

        current_step += 1
        S.add(best_node)
        added_nodes.append((best_node, best_p, current_step))

    df = pd.DataFrame(added_nodes, columns=["node", "pvalue", "step"])
    return df


def run_rwr_nx(G: nx.Graph, seeds: List[str], restart: float = 0.5) -> Dict[str, float]:
    """
    Versión rápida del RWR usando networkx.pagerank con vector de personalización.
    En PageRank: alpha = damping ≈ 1 - restart.
    """
    if not seeds:
        raise ValueError("No hay semillas válidas presentes en la red.")

    damping = 1.0 - restart  # PageRank usa 'alpha' como prob. de continuar
    personalization = {n: 0.0 for n in G}
    w = 1.0 / len(seeds)
    for s in seeds:
        personalization[s] = w

    scores = nx.pagerank(G, alpha=damping,
                         personalization=personalization,
                         weight="weight")
    # Normaliza (por seguridad)
    s = sum(scores.values())
    if s > 0:
        scores = {k: v / s for k, v in scores.items()}
    return scores


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
    parser.add_argument("--use-weights", action="store_true",
                    help="Usar el peso de las aristas (p.ej., combined_score) en la difusión")
    parser.add_argument("--exclude-seeds", action="store_true",
                    help="Excluir semillas del archivo de salida")
    parser.add_argument("--top", type=int, default=0,
                    help="Si >0, limitar la salida a los top-N nodos")
    parser.add_argument("--nx-pagerank", action="store_true",
                    help="Usar networkx.pagerank con personalization (versión rápida de RWR)")


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
        scores = run_rwr_guild(G, seeds, restart=args.restart, max_iter=args.max_iter, tol=args.tol, use_weights=args.use_weights)
        if args.nx_pagerank:
            print("[INFO] Usando versión rápida basada en PageRank de NetworkX.")
            scores = run_rwr_nx(G, seeds, restart=args.restart)
        else:
            scores = run_rwr_guild(G, seeds,
                                   restart=args.restart,
                                   max_iter=args.max_iter,
                                   tol=args.tol,
                                   use_weights=args.use_weights)

        # Guardar resultados
        df = pd.DataFrame(sorted(scores.items(), key=lambda x: x[1], reverse=True),
                  columns=["node", "score"])
        # Añadimos columnas auxiliares
        df["is_seed"] = df["node"].isin(seeds)
        # Si el usuario quiere excluir semillas
        if args.exclude_seeds:
            df = df[~df["is_seed"]]
        # Añadimos ranking
        df["rank"] = np.arange(1, len(df) + 1)
        # Si el usuario pide limitar a top-N
        if args.top and args.top > 0:
            df = df.head(args.top)
        # Cabecera con metadatos para reproducibilidad
        meta = (
            f"# algo=guild use_weights={args.use_weights} restart={args.restart} "
            f"max_iter={args.max_iter} tol={args.tol} exclude_seeds={args.exclude_seeds}\n"
            )
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(meta)
            df.to_csv(f, sep="\t", index=False)
        print(f"[OK] Resultados guardados en: {args.out} (top 5)\n{df.head().to_string(index=False)}")
        

        

    elif args.algo == "diamond":
        print("[INFO] Ejecutando DIAMOnD-lite...")
        df = run_diamond_lite(G, seeds, k=args.k)
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        meta = f"# algo=diamond-lite k={args.k}\n"
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(meta)
            df.to_csv(f, sep="\t", index=False)
        print(f"[OK] Resultados guardados en {args.out} (top 5)\n{df.head().to_string(index=False)}")


    else:
        raise ValueError("Algoritmo no reconocido.")
    
if __name__ == "__main__":
    main()
