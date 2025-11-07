# Tarea 2 — Conclusiones de Network Propagation

Autor: Carlos Marín Martínez
Asignatura: Herramientas y Algoritmos en Bioinformática
Semillas: ENO1, PGK1, HK2
Red usada: data/string_network_filtered_hugo-400.tsv (STRING filtrada, identificadores HUGO)

------------------------------------------------------------
1) Objetivo
------------------------------------------------------------
Aplicar propagación en redes para priorizar genes funcionalmente cercanos a las semillas ENO1, PGK1 y HK2 (glucólisis), utilizando:

- Un enfoque tipo GUILD (Random Walk with Restart, RWR; y su equivalente PageRank con personalización).
- Un enfoque tipo DIAMOnD (Disease Module Detection): expansión iterativa por significancia de conectividad.

------------------------------------------------------------
2) Datos y preprocesado
------------------------------------------------------------
- Red: STRING humana filtrada (~18,889 nodos, 894,129 aristas)
- Pesos: combined_score de STRING (0–1000)
- El script detecta formatos y carga pesos como atributo weight.
- Semillas presentes: 3/3 (ENO1, PGK1, HK2)

------------------------------------------------------------
3) Métodos
------------------------------------------------------------
3.1 GUILD-like (RWR / PageRank personalizado)
Ecuación estacionaria:
p = (1-α)Aᵗp + αp₀

p0: distribución inicial uniforme en semillas.
A: matriz de transición (normalizada por pesos si --use-weights)
Parámetros: restart=0.5, tol=1e-6, max_iter=50.

Variantes:
- Implementación iterativa (RWR)
- Versión rápida con networkx.pagerank (equivalente a RWR con alpha = 1 - restart)

Ejemplos de ejecución:
python scripts/NetworkPropagation.py --algo guild --use-weights --exclude-seeds --top 20 --network data/string_network_filtered_hugo-400.tsv --seeds data/genes_seed.txt --out results/guild_scores_w_top20.tsv
python scripts/NetworkPropagation.py --algo guild --nx-pagerank --use-weights --exclude-seeds --top 20 --network data/string_network_filtered_hugo-400.tsv --seeds data/genes_seed.txt --out results/guild_scores_nx_top20.tsv

3.2 DIAMOnD-lite
Expansión iterativa del conjunto S (semillas + añadidos) con p-valor hipergeométrico:
- En cada paso, para cada candidato v∉S:
  x = nº de vecinos de v en S, n = deg(v), K = |S|, N = |V|-1.
  Se añade el nodo con menor p-valor P(X≥x), X∼Hypergeom(N,K,n).
- Parámetro: k=50 (nº de nodos a añadir)

Comando:
python scripts/NetworkPropagation.py --algo diamond --k 50 --network data/string_network_filtered_hugo-400.tsv --seeds data/genes_seed.txt --out results/diamond_rank_k50.tsv

------------------------------------------------------------
4) Resultados principales
------------------------------------------------------------
GUILD (RWR con pesos) — Top-5 (sin semillas)
node   score     is_seed  rank
GPI    0.002126  False    1
H6PD   0.002004  False    2
GAPDH  0.001928  False    3
PFKM   0.001887  False    4
HK1    0.001881  False    5

GUILD (PageRank personalizado) — Top-5 (sin semillas)
node   score     is_seed  rank
GPI    0.002111  False    1
H6PD   0.001988  False    2
GAPDH  0.001916  False    3
PFKM   0.001876  False    4
HK1    0.001872  False    5

DIAMOnD-lite (k=50) — Top-5 añadidos
node    pvalue          step
ADPGK   1.968134e-08    1
HKDC1   2.552651e-10    2
GALM    1.845415e-12    3
ALDOC   1.292116e-14    4
HK3     4.178568e-16    5

------------------------------------------------------------
5) Comparación entre métodos
------------------------------------------------------------

| Método | Tipo de resultado | Ventajas principales | Coincidencias biológicas |
|---------|------------------|----------------------|---------------------------|
| **RWR (GUILD)** | Probabilidad global continua | Captura relaciones difusas y peso de aristas | GPI, GAPDH, PFKM, HK1 |
| **PageRank (GUILD rápido)** | Igual que RWR | Más eficiente y estable | Mismo top que RWR |
| **DIAMOnD-lite** | Selección iterativa con p-valor | Explicación paso a paso de nodos nuevos | HKDC1, HK3, ALDOC (hexoquinasas y glucólisis) |

Conclusión: RWR y PageRank producen rankings casi idénticos; DIAMOnD complementa con un enfoque discreto e interpretativo.

------------------------------------------------------------
6) Interpretación biológica
------------------------------------------------------------
Ambos enfoques (difusión y expansión) identifican genes relacionados con la glucólisis y el metabolismo energético.

------------------------------------------------------------
7) Sensibilidad de parámetros
------------------------------------------------------------
- Pesos (--use-weights): mejoran resultados biológicos.
- Restart=0.5: equilibrio entre exploración y reinicio.
- PageRank: versión más rápida del RWR.
- DIAMOnD k=50: módulo central robusto.

------------------------------------------------------------
8) Limitaciones
------------------------------------------------------------
- RWR: las semillas dominan el top (por eso se excluyen).
- Dependencia de pesos: ruido puede sesgar.
- DIAMOnD-lite: versión simplificada sin corrección múltiple.
- Identificadores: requiere consistencia HUGO.

------------------------------------------------------------
9) Conclusiones
------------------------------------------------------------
1. Ambos métodos (RWR/GUILD y DIAMOnD) recuperan genes coherentes con glucólisis.
2. Los pesos de STRING mejoran la calidad biológica.
3. PageRank ofrece mismo ranking que RWR con menor coste.
4. DIAMOnD-lite añade nodos relevantes paso a paso.

------------------------------------------------------------
10) Próximos pasos
------------------------------------------------------------
- Análisis GO/KEGG del top.
- Visualización de módulos.
- Corrección FDR para DIAMOnD.
- Opciones CSV/TSV y logging.

------------------------------------------------------------
11) Reproducibilidad
------------------------------------------------------------
python scripts/NetworkPropagation.py --algo guild --use-weights --exclude-seeds --top 20 --network data/string_network_filtered_hugo-400.tsv --seeds data/genes_seed.txt --out results/guild_scores_w_top20.tsv
python scripts/NetworkPropagation.py --algo guild --nx-pagerank --use-weights --exclude-seeds --top 20 --network data/string_network_filtered_hugo-400.tsv --seeds data/genes_seed.txt --out results/guild_scores_nx_top20.tsv
python scripts/NetworkPropagation.py --algo diamond --k 50 --network data/string_network_filtered_hugo-400.tsv --seeds data/genes_seed.txt --out results/diamond_rank_k50.tsv
