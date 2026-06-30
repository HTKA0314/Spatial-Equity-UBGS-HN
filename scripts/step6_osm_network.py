"""
Step 6: Download OSM walking network for the full study grid + buffer.

Input:
  config.BASE_GRID_GPKG
    Hanoi base grid 250m vector layer

Output:
  config.OSM_NETWORK_GRAPH
    Projected walking network graphml file
  ../outputs/step6_osm_network/network_quality_report.txt
"""

import os
import sys
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import osmnx as ox

import config

warnings.filterwarnings("ignore", category=UserWarning)
sys.stdout.reconfigure(encoding="utf-8")


def _read_base_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(config.BASE_GRID_GPKG):
        raise FileNotFoundError(
            f"BASE_GRID_GPKG not found: {config.BASE_GRID_GPKG}. "
            "Restore or create the base 250 m grid before downloading the OSM network."
        )

    try:
        grid = gpd.read_file(config.BASE_GRID_GPKG, layer="grid_250m")
    except Exception:
        grid = gpd.read_file(config.BASE_GRID_GPKG)

    if grid.empty:
        raise ValueError("Base grid is empty.")

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    if grid.geometry.isna().any():
        grid = grid[grid.geometry.notna()].copy()

    grid = grid[~grid.geometry.is_empty].copy()

    if grid.empty:
        raise ValueError("Base grid has no valid geometries.")

    return grid


def _union_geometry(gdf: gpd.GeoDataFrame):
    try:
        return gdf.geometry.union_all()
    except AttributeError:
        return gdf.geometry.unary_union


def _to_undirected_graph(G):
    """Return an undirected graph with OSMnx-version compatibility."""
    if hasattr(ox, "convert") and hasattr(ox.convert, "to_undirected"):
        return ox.convert.to_undirected(G)
    if hasattr(ox, "utils_graph") and hasattr(ox.utils_graph, "get_undirected"):
        return ox.utils_graph.get_undirected(G)
    return nx.Graph(G)


def _safe_project_graph(G, to_crs):
    """ĐỒNG BỘ 2026: Tránh lỗi loại bỏ hàm ox.project_graph ở cấp rễ"""
    if hasattr(ox, "projection") and hasattr(ox.projection, "project_graph"):
        return ox.projection.project_graph(G, to_crs=to_crs)
    return ox.project_graph(G, to_crs=to_crs)


def _safe_save_graphml(G, filepath):
    """ĐỒNG BỘ 2026: Tránh lỗi loại bỏ hàm ox.save_graphml ở cấp rễ"""
    if hasattr(ox, "io") and hasattr(ox.io, "save_graphml"):
        ox.io.save_graphml(G, filepath=filepath)
    else:
        ox.save_graphml(G, filepath=filepath)


def _safe_graph_from_polygon(polygon, network_type, simplify, retain_all):
    """OSMnx version-safe wrapper for graph_from_polygon."""
    if hasattr(ox, "graph") and hasattr(ox.graph, "graph_from_polygon"):
        func = ox.graph.graph_from_polygon
    else:
        func = ox.graph_from_polygon
    return func(polygon, network_type=network_type, simplify=simplify, retain_all=retain_all)


def _diagnose_graph(G) -> dict:
    n_nodes = len(G.nodes)
    n_edges = len(G.edges)

    if n_nodes == 0:
        raise ValueError("Downloaded OSM graph has zero nodes.")

    edge_node_ratio = n_edges / n_nodes if n_nodes > 0 else float("nan")

    Gu = _to_undirected_graph(G)
    components = list(nx.connected_components(Gu))
    component_sizes = sorted([len(c) for c in components], reverse=True)

    largest_component_nodes = component_sizes[0] if component_sizes else 0
    second_largest_component_nodes = component_sizes[1] if len(component_sizes) > 1 else 0
    largest_component_share = largest_component_nodes / n_nodes if n_nodes > 0 else float("nan")

    lengths = [data.get("length", float("nan")) for _, _, data in G.edges(data=True)]
    lengths = [x for x in lengths if not (x != x)]  # drop NaN

    return {
        "nodes": n_nodes,
        "edges": n_edges,
        "edge_node_ratio": edge_node_ratio,
        "connected_components": len(component_sizes),
        "largest_component_nodes": largest_component_nodes,
        "second_largest_component_nodes": second_largest_component_nodes,
        "largest_component_share": largest_component_share,
        "total_edge_length_m": float(np.sum(lengths)) if lengths else float("nan"),
        "mean_edge_length_m": float(np.mean(lengths)) if lengths else float("nan"),
    }


def _write_quality_report(report: dict, output_path: str):
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("=== STEP 6 OSM WALKING NETWORK QUALITY REPORT ===\n\n")
        for key, value in report.items():
            f.write(f"{key}: {value}\n")


def download_osm_network():
    print("--- STEP 6: DOWNLOAD OSM WALKING NETWORK ---")

    ox.settings.use_cache = True
    ox.settings.timeout = 300

    buffer_m = getattr(config, "NETWORK_BUFFER_M", 2000)
    min_nodes = getattr(config, "MIN_OSM_WALK_NODES", 50000)

    print(f"Loading base grid: {config.BASE_GRID_GPKG}")
    grid = _read_base_grid()

    print(f"Grid cells: {len(grid)}")
    print(f"Grid CRS: {grid.crs}")

    print(f"Creating study-area polygon with {buffer_m} m buffer...")
    area_proj = _union_geometry(grid)

    area_buffered = gpd.GeoSeries([area_proj], crs=config.PROJ_CRS).buffer(buffer_m).iloc[0]
    area_wgs84 = gpd.GeoSeries([area_buffered], crs=config.PROJ_CRS).to_crs("EPSG:4326").iloc[0]

    print("Downloading OSM walking network from polygon...")
    G = _safe_graph_from_polygon(area_wgs84, network_type="walk", simplify=True, retain_all=True)

    raw_report = _diagnose_graph(G)

    print(
        f"Raw graph: "
        f"nodes={raw_report['nodes']:,}, "
        f"edges={raw_report['edges']:,}, "
        f"edge/node ratio={raw_report['edge_node_ratio']:.2f}"
    )

    print(f"Connected components: {raw_report['connected_components']:,}")
    print(f"Largest component nodes: {raw_report['largest_component_nodes']:,}")
    print(f"Largest component share: {raw_report['largest_component_share'] * 100:.2f}%")

    if raw_report["second_largest_component_nodes"] > 0:
        print(f"Second largest component nodes: {raw_report['second_largest_component_nodes']:,}")

    if raw_report["nodes"] < min_nodes:
        print(f"WARNING: Only {raw_report['nodes']:,} nodes. Expected more than {min_nodes:,} for Hanoi.")

    if raw_report["edge_node_ratio"] < 1.5:
        print("WARNING: Low edge/node ratio. Network may be fragmented.")

    if raw_report["largest_component_share"] < 0.90:
        print("WARNING: Largest connected component contains less than 90% of nodes.")

    print(f"Projecting graph to {config.PROJ_CRS}...")
    G_proj = _safe_project_graph(G, to_crs=config.PROJ_CRS)

    projected_report = _diagnose_graph(G_proj)

    report = {
        "input_grid": config.BASE_GRID_GPKG,
        "output_graph": config.OSM_NETWORK_GRAPH,
        "network_type": "walk",
        "buffer_m": buffer_m,
        "graph_crs": str(G_proj.graph.get("crs", "unknown")),
        "download_timestamp": datetime.now().isoformat(timespec="seconds"),
        "osmnx_version": ox.__version__,
        "nodes": projected_report["nodes"],
        "edges": projected_report["edges"],
        "edge_node_ratio": round(projected_report["edge_node_ratio"], 4),
        "total_edge_length_m": round(projected_report["total_edge_length_m"], 2),
        "mean_edge_length_m": round(projected_report["mean_edge_length_m"], 2),
        "connected_components": projected_report["connected_components"],
        "largest_component_nodes": projected_report["largest_component_nodes"],
        "second_largest_component_nodes": projected_report["second_largest_component_nodes"],
        "largest_component_share": round(projected_report["largest_component_share"], 6),
    }

    os.makedirs(os.path.dirname(config.OSM_NETWORK_GRAPH) or ".", exist_ok=True)
    _safe_save_graphml(G_proj, filepath=config.OSM_NETWORK_GRAPH)

    report_dir = os.path.join(config.OUTPUT_DIR, "step6_osm_network")
    report_path = os.path.join(report_dir, "network_quality_report.txt")
    _write_quality_report(report, report_path)

    print("Saved OSM walking network:")
    print(config.OSM_NETWORK_GRAPH)
    print("Saved network quality report:")
    print(report_path)

    return G_proj


if __name__ == "__main__":
    download_osm_network()