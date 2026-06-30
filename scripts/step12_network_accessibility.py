"""
Step 12: Network-based walking accessibility to accessible green-blue spaces.

Input:
  config.GAP_GRID_GPKG
    hanoi_grid_250m_with_gap.gpkg

  config.OSM_NETWORK_GRAPH
    OSM walking network graph from Step 6

  config.ACCESSIBLE_SPACES_GPKG
    accessible_spaces.gpkg from Step 9

Output:
  config.ACCESSIBILITY_GPKG
    hanoi_grid_250m_walking_access.gpkg

Outputs per grid:
  dist_to_accessible_A_m
  dist_to_accessible_AB_m
  within_300m_A / within_500m_A / within_1000m_A
  within_300m_AB / within_500m_AB / within_1000m_AB
  origin_snap_dist_m
"""

import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import osmnx as ox
import config

sys.stdout.reconfigure(encoding="utf-8")


# ---------------------------------------------------------------------
# GRAPH HELPERS
# ---------------------------------------------------------------------

def _get_graph_crs(G):
    crs = G.graph.get("crs", None)
    return crs


def _ensure_graph_projected(G):
    """
    Ensure graph is projected to config.PROJ_CRS.
    Nearest-node queries and distance diagnostics are then in meters.
    """
    graph_crs = _get_graph_crs(G)

    if graph_crs is None:
        print("WARNING: Graph CRS unknown. Assuming it is already in project CRS.")
        G.graph["crs"] = config.PROJ_CRS
        return G

    if str(graph_crs).upper() == str(config.PROJ_CRS).upper():
        return G

    print(f"Projecting graph from {graph_crs} to {config.PROJ_CRS}...")
    return ox.project_graph(G, to_crs=config.PROJ_CRS)


def _to_undirected_walk_graph(G):
    """
    Convert walking graph to undirected graph for accessibility analysis.

    Rationale:
      Step 12 estimates pedestrian access to the nearest green-blue space.
      Unless the study explicitly models one-way pedestrian restrictions,
      undirected walking distance is the more defensible accessibility measure.
    """
    print("Converting walking network to undirected graph for pedestrian accessibility...")

    try:
        # OSMnx >= 2
        Gu = ox.convert.to_undirected(G)
    except Exception:
        try:
            # OSMnx <= 1
            Gu = ox.utils_graph.get_undirected(G)
        except Exception:
            Gu = G.to_undirected()

    Gu.graph["crs"] = G.graph.get("crs", config.PROJ_CRS)

    return Gu


def _nearest_nodes_for_points(G, points_gs: gpd.GeoSeries, return_dist: bool = False):
    """
    Snap points to nearest graph nodes.

    Points are converted to config.PROJ_CRS. Graph is assumed projected.
    """
    points_gs = points_gs.to_crs(config.PROJ_CRS)

    x = points_gs.x.to_numpy()
    y = points_gs.y.to_numpy()

    # OSMnx 2.0+ standard
    nodes = ox.nearest_nodes(G, x, y)
    nodes = np.asarray(nodes)

    if return_dist:
        # Calculate Euclidean distance manually since return_dist is deprecated in OSMnx 2.0+
        node_x = np.array([G.nodes[n]['x'] for n in nodes])
        node_y = np.array([G.nodes[n]['y'] for n in nodes])
        dists = np.sqrt((x - node_x)**2 + (y - node_y)**2)
        return nodes, dists

    return nodes


# ---------------------------------------------------------------------
# READ DATA
# ---------------------------------------------------------------------

def _read_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(config.GAP_GRID_GPKG):
        raise FileNotFoundError(
            f"GAP_GRID_GPKG not found: {config.GAP_GRID_GPKG}. Run Step 10 first."
        )

    try:
        grid = gpd.read_file(config.GAP_GRID_GPKG, layer="grid_gap")
    except Exception:
        grid = gpd.read_file(config.GAP_GRID_GPKG)

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    if "grid_id" not in grid.columns:
        raise ValueError("grid_id missing from grid.")

    if "area_m2" not in grid.columns:
        grid["area_m2"] = grid.geometry.area

    grid["grid_id"] = grid["grid_id"].astype(int)

    return grid


def _read_access_layer(layer_name: str, target_crs) -> gpd.GeoDataFrame:
    if not os.path.exists(config.ACCESSIBLE_SPACES_GPKG):
        raise FileNotFoundError(
            f"ACCESSIBLE_SPACES_GPKG not found: {config.ACCESSIBLE_SPACES_GPKG}. Run Step 9 first."
        )

    gdf = gpd.read_file(config.ACCESSIBLE_SPACES_GPKG, layer=layer_name)

    if gdf.empty:
        return gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs=target_crs)

    if gdf.crs is None:
        gdf = gdf.set_crs(config.PROJ_CRS)
    else:
        gdf = gdf.to_crs(target_crs)

    gdf = gdf[gdf.geometry.notna() & (~gdf.geometry.is_empty)].copy()
    gdf["geometry"] = gdf.geometry.make_valid()
    gdf = gdf[gdf.geometry.notna() & (~gdf.geometry.is_empty)].copy()

    return gdf


# ---------------------------------------------------------------------
# DESTINATION SAMPLING
# ---------------------------------------------------------------------

def _sample_dest_nodes(G, spaces: gpd.GeoDataFrame, label: str) -> list:
    """
    Sample representative point + boundary points per accessible polygon.
    Then snap all sampled points to nearest graph nodes.

    This avoids relying on a single polygon centroid for large parks/lakes.
    """
    if spaces.empty:
        print(f"  No spaces for {label}.")
        return []

    spaces = spaces.to_crs(config.PROJ_CRS)
    
    spaces = spaces.dissolve()
    spaces = gpd.GeoDataFrame(spaces).explode(index_parts=False).reset_index(drop=True)

    pts = []
    skipped = 0

    for geom in spaces.geometry:
        if geom is None or geom.is_empty:
            skipped += 1
            continue

        try:
            # Point inside polygon
            pts.append(geom.representative_point())

            # Boundary samples: sample every 150m (or as configured)
            boundary = geom.boundary
            if boundary is not None and not boundary.is_empty:
                spacing_m = getattr(config, "DESTINATION_SAMPLE_SPACING_M", 150)
                length = boundary.length
                n_samples = max(4, int(np.ceil(length / spacing_m)))

                for i in range(n_samples):
                    pts.append(boundary.interpolate(i / n_samples, normalized=True))

        except Exception:
            skipped += 1

    if skipped:
        print(f"  WARNING: {skipped} space(s) skipped during node sampling for {label}")

    if not pts:
        return []

    pts_gs = gpd.GeoSeries(pts, crs=config.PROJ_CRS)

    print(f"  Querying nearest graph nodes for {len(pts_gs)} sampled points ({label})...")

    nodes, snap_dists = _nearest_nodes_for_points(G, pts_gs, return_dist=True)

    max_dest_snap_m = getattr(config, "DEST_SNAP_MAX_M", 150)
    valid = np.isfinite(snap_dists) & (snap_dists <= max_dest_snap_m)
    
    if valid.sum() == 0:
        print(f"  WARNING: No destination nodes within {max_dest_snap_m} m for {label}.")
        return []

    bad = (~valid).sum()
    if bad > 0:
        print(f"  WARNING: Dropped {bad} sampled destination points > {max_dest_snap_m} m from network.")

    if np.isfinite(snap_dists).any():
        print(
            f"  Destination snap distance {label}: "
            f"mean={np.nanmean(snap_dists):.1f} m, "
            f"p95={np.nanpercentile(snap_dists, 95):.1f} m, "
            f"max={np.nanmax(snap_dists):.1f} m"
        )

    nodes = np.asarray(nodes)[valid]
    return list(set(nodes))


# ---------------------------------------------------------------------
# ACCESSIBILITY CALCULATION
# ---------------------------------------------------------------------

def _dijkstra_distances(G, dest_nodes: list, cutoff_m: int) -> dict:
    if not dest_nodes:
        return {}

    return nx.multi_source_dijkstra_path_length(
        G,
        sources=dest_nodes,
        weight="length",
        cutoff=cutoff_m,
    )


def _add_access_flags(grid: gpd.GeoDataFrame, dist_col: str, suffix: str) -> gpd.GeoDataFrame:
    for threshold in (300, 500, 1000):
        grid[f"within_{threshold}m_{suffix}"] = (grid[dist_col] <= threshold).astype(int)
    return grid


def calculate_walking_access() -> gpd.GeoDataFrame:
    print("--- STEP 12: NETWORK WALKING ACCESSIBILITY ---")

    dijkstra_cutoff = getattr(config, "DIJKSTRA_CUTOFF_M", 1000)
    snap_warning_threshold = getattr(config, "SNAP_WARNING_THRESHOLD_M", 100)

    if not os.path.exists(config.OSM_NETWORK_GRAPH):
        raise FileNotFoundError(
            f"OSM_NETWORK_GRAPH not found: {config.OSM_NETWORK_GRAPH}. Run Step 6 first."
        )

    # ------------------------------------------------------------------
    # 1. Load and prepare graph
    # ------------------------------------------------------------------
    print("Loading OSM walking network...")
    G = ox.load_graphml(config.OSM_NETWORK_GRAPH)
    G = _ensure_graph_projected(G)
    G = _to_undirected_walk_graph(G)

    n, e = len(G.nodes), len(G.edges)
    print(f"Graph nodes: {n:,}, edges: {e:,}")
    print(f"Graph CRS: {G.graph.get('crs', 'unknown')}")

    # ------------------------------------------------------------------
    # 2. Load grid and accessible spaces
    # ------------------------------------------------------------------
    print("Loading grid with gap metrics...")
    grid = _read_grid()
    print(f"Grid cells: {len(grid)}")

    print("Loading accessible-space layers...")
    spaces_A = _read_access_layer("accessible_spaces_A", config.PROJ_CRS)
    spaces_AB = _read_access_layer("accessible_spaces_AB", config.PROJ_CRS)

    print(f"Level A spaces: {len(spaces_A)}")
    print(f"Level A+B spaces: {len(spaces_AB)}")

    if spaces_A.empty:
        raise ValueError("accessible_spaces_A layer is empty. Check Step 9 output.")

    if spaces_AB.empty:
        raise ValueError("accessible_spaces_AB layer is empty. Check Step 9 output.")

    # ------------------------------------------------------------------
    # 3. Origin nodes from grid representative points
    # ------------------------------------------------------------------
    print("Finding nearest graph nodes for grid cells...")

    origin_points = gpd.GeoSeries(
        grid.geometry.representative_point(),
        crs=config.PROJ_CRS,
    )

    origin_nodes, origin_snap_dists = _nearest_nodes_for_points(
        G,
        origin_points,
        return_dist=True,
    )

    grid["origin_node"] = origin_nodes
    grid["origin_snap_dist_m"] = origin_snap_dists

    unique_origin_nodes = len(set(origin_nodes))
    origin_node_ratio = unique_origin_nodes / len(origin_nodes) if len(origin_nodes) > 0 else np.nan

    print(
        f"Origin nodes: {unique_origin_nodes:,} unique / {len(origin_nodes):,} grid cells "
        f"({origin_node_ratio * 100:.1f}% unique)"
    )

    if np.isfinite(origin_node_ratio) and origin_node_ratio < 0.90:
        print(
            "WARNING: <90% unique origin nodes. Multiple grid cells share the same graph node; "
            "network resolution may be coarse relative to the 250m grid."
        )

    if np.isfinite(origin_snap_dists).any():
        print(
            f"Origin snap distance: "
            f"mean={np.nanmean(origin_snap_dists):.1f} m, "
            f"p95={np.nanpercentile(origin_snap_dists, 95):.1f} m, "
            f"max={np.nanmax(origin_snap_dists):.1f} m"
        )

        high_snap = np.sum(origin_snap_dists > snap_warning_threshold)
        if high_snap > 0:
            print(
                f"WARNING: {high_snap} grid cells snapped farther than "
                f"{snap_warning_threshold} m from the walking network."
            )

    # ------------------------------------------------------------------
    # 4. Level A distances
    # ------------------------------------------------------------------
    print("Sampling Level A destination nodes...")
    dest_A = _sample_dest_nodes(G, spaces_A, "A")
    print(f"  Unique destination nodes A: {len(dest_A)}")

    if not dest_A:
        raise ValueError("No Level A destination nodes found.")

    print(f"Running multi-source Dijkstra for Level A, cutoff={dijkstra_cutoff} m...")
    dist_A = _dijkstra_distances(G, dest_A, dijkstra_cutoff)

    grid["dist_to_accessible_A_m"] = [
        dist_A.get(int(node), np.inf) for node in origin_nodes
    ]

    # ------------------------------------------------------------------
    # 5. Level A+B distances
    # ------------------------------------------------------------------
    print("Sampling Level A+B destination nodes...")
    dest_AB = _sample_dest_nodes(G, spaces_AB, "A+B")
    dest_AB = list(set(dest_AB).union(set(dest_A)))
    print(f"  Unique destination nodes A+B after forcing include A: {len(dest_AB)}")

    if not dest_AB:
        raise ValueError("No Level A+B destination nodes found.")

    print(f"Running multi-source Dijkstra for Level A+B, cutoff={dijkstra_cutoff} m...")
    dist_AB = _dijkstra_distances(G, dest_AB, dijkstra_cutoff)

    grid["dist_to_accessible_AB_m"] = [
        dist_AB.get(int(node), np.inf) for node in origin_nodes
    ]

    # ------------------------------------------------------------------
    # ZERO-DISTANCE CORRECTION (Physical Intersection)
    # ------------------------------------------------------------------
    print("Applying Zero-Distance Correction for cells physically overlapping parks...")
    # If a cell physically contains accessible space (calculated in Step 10), distance is 0.
    if "accessible_gb_area_A_final_m2" in grid.columns:
        overlap_A = grid["accessible_gb_area_A_final_m2"] > 0
        grid.loc[overlap_A, "dist_to_accessible_A_m"] = 0.0
        print(f"  Forced 0m distance for {overlap_A.sum()} cells overlapping Level A parks.")

    if "accessible_gb_area_AB_final_m2" in grid.columns:
        overlap_AB = grid["accessible_gb_area_AB_final_m2"] > 0
        grid.loc[overlap_AB, "dist_to_accessible_AB_m"] = 0.0
        # AB inherently includes A, so force AB to 0 if A is 0
        overlap_AB_total = overlap_AB | (grid.get("accessible_gb_area_A_final_m2", pd.Series(False, index=grid.index)) > 0)
        grid.loc[overlap_AB_total, "dist_to_accessible_AB_m"] = 0.0
        print(f"  Forced 0m distance for {overlap_AB_total.sum()} cells overlapping Level A+B parks.")


    # Backward-compatible aliases
    grid["dist_to_accessible_A"] = grid["dist_to_accessible_A_m"]
    grid["dist_to_accessible_AB"] = grid["dist_to_accessible_AB_m"]

    # Access flags
    grid = _add_access_flags(grid, "dist_to_accessible_A_m", "A")
    grid = _add_access_flags(grid, "dist_to_accessible_AB_m", "AB")

    # Backward-compatible conservative Level A flags
    grid["within_300m"] = grid["within_300m_A"]
    grid["within_500m"] = grid["within_500m_A"]
    grid["within_1000m"] = grid["within_1000m_A"]

    # Check monotonicity
    A_fill = grid["dist_to_accessible_A_m"].fillna(np.inf)
    AB_fill = grid["dist_to_accessible_AB_m"].fillna(np.inf)
    violation = (AB_fill > A_fill + 1.0).sum()
    if violation > 0:
        print(f"WARNING: Monotonicity violated in {violation} cells (dist_AB > dist_A).")

    # ------------------------------------------------------------------
    # 6. Validation
    # ------------------------------------------------------------------
    inf_A = int(np.isinf(grid["dist_to_accessible_A_m"]).sum())
    inf_AB = int(np.isinf(grid["dist_to_accessible_AB_m"]).sum())

    reachable_A = grid.loc[
        ~np.isinf(grid["dist_to_accessible_A_m"]),
        "dist_to_accessible_A_m",
    ]

    reachable_AB = grid.loc[
        ~np.isinf(grid["dist_to_accessible_AB_m"]),
        "dist_to_accessible_AB_m",
    ]

    print("\nAccessibility summary:")
    print(f"  Dijkstra cutoff       : {dijkstra_cutoff} m")
    print(f"  Unreachable A         : {inf_A} ({inf_A / len(grid) * 100:.1f}%)")
    print(f"  Unreachable A+B       : {inf_AB} ({inf_AB / len(grid) * 100:.1f}%)")

    if len(reachable_A) > 0:
        print(f"  Mean dist A reachable   : {reachable_A.mean():.0f} m")
        print(f"  Median dist A reachable : {reachable_A.median():.0f} m")

    if len(reachable_AB) > 0:
        print(f"  Mean dist AB reachable  : {reachable_AB.mean():.0f} m")
        print(f"  Median dist AB reachable: {reachable_AB.median():.0f} m")

    for suffix in ("A", "AB"):
        print(f"  Within 300m {suffix}  : {grid[f'within_300m_{suffix}'].mean() * 100:.1f}%")
        print(f"  Within 500m {suffix}  : {grid[f'within_500m_{suffix}'].mean() * 100:.1f}%")
        print(f"  Within 1000m {suffix} : {grid[f'within_1000m_{suffix}'].mean() * 100:.1f}%")

    if grid["within_1000m_A"].mean() < 0.30:
        print(
            "WARNING: <30% grid cells are within 1000m of Level A spaces. "
            "This may reflect conservative access rules or incomplete OSM tagging."
        )

    if grid["within_1000m_A"].mean() > 0.95:
        print(
            "WARNING: >95% grid cells are within 1000m of Level A spaces. "
            "Destination classification may be over-inclusive."
        )

    # Population-weighted validation
    if "population" in grid.columns:
        pop = grid["population"].fillna(0)
        total_pop = pop.sum()
        
        populated = grid[pop > 0]
        if not populated.empty:
            print("\nPopulated-cell accessibility:")
            for suffix in ("A", "AB"):
                for threshold in (300, 500, 1000):
                    col = f"within_{threshold}m_{suffix}"
                    print(f"  Populated cells within {threshold}m {suffix}: {populated[col].mean() * 100:.1f}%")
        
        if total_pop > 0:
            print("\nPopulation-weighted accessibility:")
            for suffix in ("A", "AB"):
                for threshold in (300, 500, 1000):
                    col = f"within_{threshold}m_{suffix}"
                    pct = (grid.loc[grid[col] == 1, "population"].sum() / total_pop) * 100
                    print(f"  Population within {threshold}m {suffix}: {pct:.1f}%")

    # ------------------------------------------------------------------
    # 7. Save
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(config.ACCESSIBILITY_GPKG) or ".", exist_ok=True)

    grid.to_file(
        config.ACCESSIBILITY_GPKG,
        driver="GPKG",
        layer="grid_walking_access",
    )

    print(f"\nSaved -> {config.ACCESSIBILITY_GPKG}")
    print("Next step: run Step 13 spatial regression or Step 14 risk hotspot analysis.")

    return grid


if __name__ == "__main__":
    calculate_walking_access()