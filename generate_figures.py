"""
SheafNMA Figure Generator
Generates 3 publication-quality figures for the manuscript.
Requires: matplotlib, numpy (standard Anaconda/pip install).
Output: figures/ directory with PNG (300dpi) and PDF.

Usage: python generate_figures.py
"""

import sys
import io
import os
import math
import numpy as np

# UTF-8 stdout for Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

# ===== Deterministic PRNG (matching JS xoshiro128**) =====

def splitmix32(seed):
    """Splitmix32 for seed expansion."""
    a = seed & 0xFFFFFFFF
    def next_val():
        nonlocal a
        a = (a + 0x9E3779B9) & 0xFFFFFFFF
        t = a ^ (a >> 16)
        t = (t * 0x21F0AAAD) & 0xFFFFFFFF
        t = t ^ (t >> 15)
        t = (t * 0x735A2D97) & 0xFFFFFFFF
        return (t ^ (t >> 15)) & 0xFFFFFFFF
    return next_val


def xoshiro128ss(seed):
    """xoshiro128** PRNG matching the JS implementation."""
    sm = splitmix32(seed)
    s = [sm(), sm(), sm(), sm()]

    def next_val():
        result = ((s[1] * 5) & 0xFFFFFFFF)
        result = (((result << 7) | (result >> 25)) & 0xFFFFFFFF)
        result = (result * 9) & 0xFFFFFFFF
        t = (s[1] << 9) & 0xFFFFFFFF

        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]

        s[2] ^= t
        s[3] = ((s[3] << 11) | (s[3] >> 21)) & 0xFFFFFFFF

        return result / 4294967296.0

    return next_val


def normal_random(rng):
    """Box-Muller transform."""
    u1 = rng()
    u2 = rng()
    return math.sqrt(-2 * math.log(u1 + 1e-15)) * math.cos(2 * math.pi * u2)


# ===== Simulate data (matching JS exactly) =====

def get_simulated_data():
    """Generate simulated dataset matching JS getSimulatedData() with seed=42."""
    rng = xoshiro128ss(42)
    truth = {'A': 0, 'B': -0.3, 'C': -0.5, 'D': -0.8, 'E': -1.0}

    def make_contrast(study, t1, t2, shift_bias=0):
        true_effect = truth[t1] - truth[t2]
        noise = normal_random(rng) * 0.05
        se = 0.10 + rng() * 0.20
        effect = round((true_effect + noise + shift_bias) * 10000) / 10000
        se = round(se * 10000) / 10000
        return {'study': study, 'treat1': t1, 'treat2': t2, 'effect': effect, 'se': se}

    data = [
        make_contrast("Sim01", "A", "B"),
        make_contrast("Sim02", "A", "C"),
        make_contrast("Sim03", "A", "E"),
        make_contrast("Sim04", "B", "C"),
        make_contrast("Sim05", "B", "D"),
        make_contrast("Sim06", "B", "E"),
        make_contrast("Sim07", "C", "D"),
        make_contrast("Sim08", "C", "E"),
        make_contrast("Sim09", "D", "E"),
        make_contrast("Sim10", "A", "B"),
        make_contrast("Sim11_INC", "A", "D", 0.5),
        make_contrast("Sim12_INC", "A", "D", 0.5),
    ]
    return data


def build_network(contrasts):
    """Build network from contrasts, pool multi-study edges."""
    node_set = set()
    for c in contrasts:
        node_set.add(c['treat1'])
        node_set.add(c['treat2'])
    nodes = sorted(node_set)

    edge_map = {}
    for c in contrasts:
        t1, t2, eff = c['treat1'], c['treat2'], c['effect']
        if t1 > t2:
            t1, t2 = t2, t1
            eff = -eff
        key = f"{t1}-{t2}"
        if key not in edge_map:
            edge_map[key] = {'treat1': t1, 'treat2': t2, 'contrasts': []}
        edge_map[key]['contrasts'].append({'study': c['study'], 'effect': eff, 'se': c['se']})

    edges = []
    for e in edge_map.values():
        cArr = e['contrasts']
        k = len(cArr)
        sum_w = 0
        sum_wd = 0
        studies = []
        for c in cArr:
            w = 1 / (c['se'] ** 2)
            sum_w += w
            sum_wd += w * c['effect']
            studies.append(c['study'])
        d_pooled = sum_wd / sum_w
        se_pooled = 1 / math.sqrt(sum_w)
        edges.append({
            'treat1': e['treat1'],
            'treat2': e['treat2'],
            'effect': d_pooled,
            'se': se_pooled,
            'studies': studies,
            'k': k
        })

    return {'nodes': nodes, 'edges': edges}


def sheaf_analysis(network):
    """Core sheaf analysis matching JS implementation."""
    nodes = network['nodes']
    edges = network['edges']
    n = len(nodes)
    m = len(edges)
    node_idx = {nd: i for i, nd in enumerate(nodes)}

    # Coboundary map F (m x n)
    F = np.zeros((m, n))
    d = np.zeros(m)

    for e_idx, edge in enumerate(edges):
        i = node_idx[edge['treat1']]
        j = node_idx[edge['treat2']]
        w = 1.0 / edge['se']
        F[e_idx, i] = -w
        F[e_idx, j] = w
        d[e_idx] = edge['effect'] * w

    # Sheaf Laplacian
    L = F.T @ F

    # Solve for node estimates (fix node 0 = 0)
    Ftd = F.T @ d
    if n == 1:
        node_estimates = np.array([0.0])
    else:
        L_red = L[1:, 1:]
        b_red = Ftd[1:]
        x_red = np.linalg.solve(L_red, b_red)
        node_estimates = np.zeros(n)
        node_estimates[1:] = x_red

    # Edge residuals and scores
    Fx = F @ node_estimates
    residuals = d - Fx
    edge_scores = residuals ** 2

    max_score = np.max(edge_scores) if np.max(edge_scores) > 1e-15 else 1.0
    normalized_scores = (edge_scores / max_score) * 100

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(L)

    # GII
    total_resid_ss = np.sum(residuals ** 2)
    df_incons = m - (n - 1)
    gii = total_resid_ss / df_incons if df_incons > 0 else 0

    return {
        'nodes': nodes,
        'edges': edges,
        'node_estimates': node_estimates,
        'residuals': residuals,
        'edge_scores': edge_scores,
        'normalized_scores': normalized_scores,
        'laplacian': L,
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'gii': gii,
    }


# ===== Figure Generation =====

def score_to_color(score):
    """Map 0-100 score to green-yellow-red."""
    if score < 50:
        t = score / 50
        r = (34 + t * (234 - 34)) / 255
        g = (197 + t * (179 - 197)) / 255
        b = (94 + t * (8 - 94)) / 255
    else:
        t = (score - 50) / 50
        r = (234 + t * (239 - 234)) / 255
        g = (179 + t * (68 - 179)) / 255
        b = (8 + t * (68 - 8)) / 255
    return (r, g, b)


def generate_figure1(result, network, outdir):
    """Figure 1: Network heatmap with inconsistency coloring."""
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    nodes = network['nodes']
    edges = network['edges']
    n = len(nodes)

    # Circular layout
    positions = {}
    for i, nd in enumerate(nodes):
        angle = (2 * math.pi * i / n) - math.pi / 2
        positions[nd] = (math.cos(angle) * 2.5, math.sin(angle) * 2.5)

    # Draw edges colored by inconsistency score
    for i, edge in enumerate(edges):
        p1 = positions[edge['treat1']]
        p2 = positions[edge['treat2']]
        score = result['normalized_scores'][i]
        color = score_to_color(score)
        precision = 1.0 / edge['se']
        max_prec = max(1.0 / e['se'] for e in edges)
        min_prec = min(1.0 / e['se'] for e in edges)
        prec_range = max_prec - min_prec if max_prec > min_prec else 1
        lw = 2 + ((precision - min_prec) / prec_range) * 6

        ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                color=color, linewidth=lw, solid_capstyle='round', zorder=1)

        # Edge label
        mx = (p1[0] + p2[0]) / 2
        my = (p1[1] + p2[1]) / 2
        # Offset label slightly perpendicular to edge
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        length = math.sqrt(dx**2 + dy**2) + 1e-10
        nx, ny = -dy / length, dx / length
        offset = 0.18
        ax.text(mx + nx * offset, my + ny * offset,
                f"{edge['effect']:.2f}",
                fontsize=7, ha='center', va='center', color='#555555')

    # Draw nodes
    for nd in nodes:
        px, py = positions[nd]
        degree = sum(1 for e in edges if e['treat1'] == nd or e['treat2'] == nd)
        radius = 0.35 + 0.1 * math.sqrt(degree)

        circle = plt.Circle((px, py), radius, facecolor='white',
                            edgecolor='#1e3a5f', linewidth=2.5, zorder=3)
        ax.add_patch(circle)
        ax.text(px, py, nd, fontsize=14, fontweight='bold',
                ha='center', va='center', zorder=4, color='#1a1a2e')

    # Color bar
    cmap = LinearSegmentedColormap.from_list('incons',
        [(34/255, 197/255, 94/255),
         (234/255, 179/255, 8/255),
         (239/255, 68/255, 68/255)])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 100))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, aspect=20, pad=0.05)
    cbar.set_label('Inconsistency Score (%)', fontsize=10)

    ax.set_xlim(-4, 4)
    ax.set_ylim(-4, 4)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Network Inconsistency Heatmap\n(Simulated: 5 treatments, planted A-D inconsistency)',
                 fontsize=12, fontweight='bold', pad=15)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'figure1_network_heatmap.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(outdir, 'figure1_network_heatmap.pdf'), bbox_inches='tight')
    plt.close(fig)
    print("Figure 1: Network heatmap saved.")


def generate_figure2(result, outdir):
    """Figure 2: Eigenvalue scree plot."""
    fig, ax = plt.subplots(1, 1, figsize=(7, 4.5))

    eigenvalues = result['eigenvalues']
    n = len(eigenvalues)

    colors = []
    for i, ev in enumerate(eigenvalues):
        if ev < 1e-8:
            colors.append('#94a3b8')  # zero eigenvalue
        elif i == n - 1:
            colors.append('#d97706')  # largest (highlighted)
        else:
            colors.append('#3b82f6')  # non-zero

    bars = ax.bar(range(1, n + 1), eigenvalues, color=colors,
                  edgecolor=['#b45309' if i == n - 1 else 'none' for i in range(n)],
                  linewidth=[2 if i == n - 1 else 0 for i in range(n)],
                  width=0.6, zorder=3)

    # Value labels on bars
    for i, (bar, ev) in enumerate(zip(bars, eigenvalues)):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f'{ev:.2f}', ha='center', va='bottom', fontsize=9,
                fontfamily='monospace', color='#374151')

    ax.set_xlabel('Eigenvalue Index', fontsize=11)
    ax.set_ylabel('Eigenvalue', fontsize=11)
    ax.set_title('Sheaf Laplacian Eigenvalue Scree Plot\n(Simulated dataset, 5 treatments)',
                 fontsize=12, fontweight='bold')
    ax.set_xticks(range(1, n + 1))
    ax.grid(axis='y', alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#94a3b8', label='Zero (connected component)'),
        mpatches.Patch(facecolor='#3b82f6', label='Non-zero (inconsistency mode)'),
        mpatches.Patch(facecolor='#d97706', edgecolor='#b45309',
                       linewidth=2, label='Dominant mode (highlighted)'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'figure2_scree_plot.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(outdir, 'figure2_scree_plot.pdf'), bbox_inches='tight')
    plt.close(fig)
    print("Figure 2: Eigenvalue scree plot saved.")


def generate_figure3(result, network, outdir):
    """Figure 3: Edge inconsistency scores bar chart."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    edges = network['edges']
    scores = result['normalized_scores']

    # Sort by score descending for clearest visual
    edge_labels = [f"{e['treat1']}-{e['treat2']}" for e in edges]
    sorted_indices = np.argsort(scores)[::-1]
    sorted_labels = [edge_labels[i] for i in sorted_indices]
    sorted_scores = [scores[i] for i in sorted_indices]
    sorted_colors = [score_to_color(s) for s in sorted_scores]

    bars = ax.barh(range(len(sorted_labels)), sorted_scores,
                   color=sorted_colors, edgecolor='#e5e7eb', linewidth=0.5,
                   height=0.6, zorder=3)

    ax.set_yticks(range(len(sorted_labels)))
    ax.set_yticklabels(sorted_labels, fontsize=10, fontfamily='monospace')
    ax.invert_yaxis()

    # Score labels
    for i, (bar, score) in enumerate(zip(bars, sorted_scores)):
        if score > 5:
            ax.text(score - 2, i, f'{score:.1f}%', va='center', ha='right',
                    fontsize=9, fontweight='bold', color='white', fontfamily='monospace')
        else:
            ax.text(score + 1, i, f'{score:.1f}%', va='center', ha='left',
                    fontsize=9, color='#374151', fontfamily='monospace')

    ax.set_xlabel('Normalized Inconsistency Score (%)', fontsize=11)
    ax.set_title('Edge Inconsistency Scores\n(Simulated: A-D = planted inconsistency, all others consistent)',
                 fontsize=12, fontweight='bold')
    ax.set_xlim(0, 110)
    ax.grid(axis='x', alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    # Add ground truth annotation
    ax.axvline(x=50, color='#ef4444', linestyle='--', linewidth=1, alpha=0.5, zorder=2)
    ax.text(52, len(sorted_labels) - 0.5, 'High inconsistency\nthreshold',
            fontsize=8, color='#ef4444', va='bottom')

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, 'figure3_edge_scores.png'), dpi=300, bbox_inches='tight')
    fig.savefig(os.path.join(outdir, 'figure3_edge_scores.pdf'), bbox_inches='tight')
    plt.close(fig)
    print("Figure 3: Edge inconsistency scores saved.")


# ===== Main =====

def main():
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures')
    os.makedirs(outdir, exist_ok=True)

    print("Generating simulated data (seed=42)...")
    data = get_simulated_data()
    network = build_network(data)
    result = sheaf_analysis(network)

    print(f"Network: {len(network['nodes'])} nodes, {len(network['edges'])} edges")
    print(f"GII: {result['gii']:.6f}")
    print(f"Eigenvalues: {[f'{v:.4f}' for v in result['eigenvalues']]}")

    # Print edge scores
    print("\nEdge inconsistency scores:")
    for i, edge in enumerate(network['edges']):
        label = f"{edge['treat1']}-{edge['treat2']}"
        print(f"  {label:6s}: {result['normalized_scores'][i]:6.1f}% (raw={result['edge_scores'][i]:.6f})")

    # Verify A-D is #1
    ad_idx = None
    for i, edge in enumerate(network['edges']):
        if (edge['treat1'] == 'A' and edge['treat2'] == 'D') or \
           (edge['treat1'] == 'D' and edge['treat2'] == 'A'):
            ad_idx = i
            break

    if ad_idx is not None:
        ad_score = result['normalized_scores'][ad_idx]
        print(f"\nA-D edge: normalized score = {ad_score:.1f}% (expected: 100.0%)")
        assert ad_score > 90, f"FAIL: A-D score {ad_score} not highest!"
        print("PASS: A-D correctly identified as most inconsistent edge.")
    else:
        print("WARNING: A-D edge not found!")

    # Generate figures
    print("\nGenerating figures...")
    generate_figure1(result, network, outdir)
    generate_figure2(result, outdir)
    generate_figure3(result, network, outdir)

    print(f"\nAll 3 figures saved to: {outdir}")
    print("Formats: PNG (300 dpi) + PDF")


if __name__ == '__main__':
    main()
