"""
Step 9: Extract and classify accessible green-blue spaces from OSM and JRC.

Input:
  config.LANDCOVER_GRID_GPKG
    hanoi_grid_250m_with_landcover.gpkg

  config.WATER_PATH
    JRC permanent water polygon layer

Output:
  config.ACCESSIBLE_SPACES_GPKG
    accessible_spaces.gpkg

Layers written:
  - accessible_spaces_all
  - accessible_spaces_A
  - accessible_spaces_AB
  - jrc_water_reference

Access levels:
  A = clearly public / official public access
  B = semi-public, permissive, or de facto accessible
  C = restricted, private, uncertain, or physical blue reference only

Important fix:
  Overpass is queried using a simple buffered bounding box, not the detailed
  15,077-cell study-area union. Results are then clipped back to the exact
  study area. This avoids the Overpass 16 MB query-length limit.
"""

import os
import sys
import warnings

import geopandas as gpd
import pandas as pd
import osmnx as ox
import config

from shapely.geometry import box

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore", category=UserWarning)


# ---------------------------------------------------------------------
# OSM TAGS
# ---------------------------------------------------------------------

OSM_TAGS = {
    "leisure": [
        "park",
        "garden",
        "recreation_ground",
        "playground",
        "nature_reserve",
        "sports_centre",
        "pitch",
    ],
    "landuse": [
        "grass",
        "recreation_ground",
        "village_green",
        "forest",
        "cemetery",
        "reservoir",
        "basin",
    ],
    "natural": [
        "wood",
        "grassland",
        "scrub",
        "heath",
        "wetland",
        "water",
    ],
    "water": [
        "lake",
        "pond",
        "reservoir",
        "river",
        "basin",
        "canal",
    ],
    "waterway": [
        "river",
        "canal",
        "stream",
        "drain",
        "ditch",
    ],
    "amenity": [
        "grave_yard",
    ],
    "boundary": [
        "protected_area",
    ],
    "tourism": [
        "picnic_site",
    ],
}


PRIVATE_ACCESS_VALUES = {
    "private",
    "no",
    "customers",
    "permit",
    "delivery",
    "destination",
}

PUBLIC_ACCESS_VALUES = {
    "yes",
    "public",
    "designated",
}

PERMISSIVE_ACCESS_VALUES = {
    "permissive",
}


# ---------------------------------------------------------------------
# READ INPUTS
# ---------------------------------------------------------------------

def _read_grid() -> gpd.GeoDataFrame:
    """Read main grid after land-cover merge."""

    if not os.path.exists(config.LANDCOVER_GRID_GPKG):
        raise FileNotFoundError(
            f"LANDCOVER_GRID_GPKG not found: {config.LANDCOVER_GRID_GPKG}. "
            "Run Step 8 first."
        )

    try:
        grid = gpd.read_file(config.LANDCOVER_GRID_GPKG, layer="grid_with_lc")
    except Exception:
        grid = gpd.read_file(config.LANDCOVER_GRID_GPKG)

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    return grid


def _get_study_area_polygon(grid: gpd.GeoDataFrame):
    """Create exact study-area polygon from final grid."""

    try:
        area = grid.geometry.union_all()
    except AttributeError:
        area = grid.geometry.unary_union

    return area


def _get_overpass_query_polygon(grid: gpd.GeoDataFrame):
    """
    Create a simple polygon for Overpass query.

    Do not send the detailed union of 15,586 grid cells to Overpass,
    because it can exceed Overpass query-length limits.

    We query by a buffered bounding box, then clip results back to
    the exact study-area polygon later.
    """

    buffer_m = getattr(config, "OSM_FEATURE_QUERY_BUFFER_M", 500)

    minx, miny, maxx, maxy = grid.total_bounds

    bbox_poly = box(minx, miny, maxx, maxy)

    query_poly_proj = (
        gpd.GeoSeries([bbox_poly], crs=config.PROJ_CRS)
        .buffer(buffer_m)
        .iloc[0]
    )

    query_poly_wgs84 = (
        gpd.GeoSeries([query_poly_proj], crs=config.PROJ_CRS)
        .to_crs("EPSG:4326")
        .iloc[0]
    )

    return query_poly_wgs84


# ---------------------------------------------------------------------
# OSM DOWNLOAD
# ---------------------------------------------------------------------

def _fetch_osm_features(area_wgs84):
    """Fetch OSM features with compatibility across OSMnx versions."""

    print("Downloading OSM green-blue candidate features...")

    if hasattr(ox, "features_from_polygon"):
        osm = ox.features_from_polygon(area_wgs84, tags=OSM_TAGS)
    else:
        osm = ox.geometries_from_polygon(area_wgs84, tags=OSM_TAGS)

    if osm.empty:
        raise ValueError("No OSM features returned for the query area.")

    osm = osm.reset_index()

    print(f"  Raw OSM features returned: {len(osm)}")

    return osm


# ---------------------------------------------------------------------
# CLASSIFICATION HELPERS
# ---------------------------------------------------------------------

def _keep_polygonal_features(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Keep only polygon and multipolygon geometries."""

    gdf = gdf[gdf.geometry.notna() & (~gdf.geometry.is_empty)].copy()

    gdf = gdf[
        gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])
    ].copy()

    if gdf.empty:
        raise ValueError("No polygonal OSM green-blue features after filtering.")

    gdf["geometry"] = gdf.geometry.make_valid()
    gdf = gpd.GeoDataFrame(gdf).explode(index_parts=False).reset_index(drop=True)
    return gdf


def _classify_space_type(row) -> str:
    """Classify candidate into green / blue / mixed / other."""

    natural = str(row.get("natural", "")).lower()
    water = str(row.get("water", "")).lower()
    waterway = str(row.get("waterway", "")).lower()
    landuse = str(row.get("landuse", "")).lower()
    leisure = str(row.get("leisure", "")).lower()
    amenity = str(row.get("amenity", "")).lower()
    boundary = str(row.get("boundary", "")).lower()

    blue_conditions = (
        natural == "water"
        or water in {"lake", "pond", "reservoir", "river", "basin", "canal"}
        or landuse in {"reservoir", "basin"}
        or waterway in {"river", "canal", "stream", "drain", "ditch"}
    )

    green_conditions = (
        leisure in {
            "park",
            "garden",
            "recreation_ground",
            "playground",
            "nature_reserve",
            "sports_centre",
            "pitch",
        }
        or landuse in {
            "grass",
            "recreation_ground",
            "village_green",
            "forest",
            "cemetery",
        }
        or natural in {"wood", "grassland", "scrub", "heath", "wetland"}
        or amenity == "grave_yard"
        or boundary == "protected_area"
    )

    if blue_conditions and green_conditions:
        return "mixed"
    if blue_conditions:
        return "blue"
    if green_conditions:
        return "green"

    return "other"


def _classify_access_level(row) -> str:
    """
    Assign access level:
      A = public
      B = semi-public / de facto accessible
      C = restricted / private / uncertain
    """

    access = str(row.get("access", "")).lower().strip()
    foot = str(row.get("foot", "")).lower().strip()

    leisure = str(row.get("leisure", "")).lower().strip()
    landuse = str(row.get("landuse", "")).lower().strip()
    natural = str(row.get("natural", "")).lower().strip()
    amenity = str(row.get("amenity", "")).lower().strip()
    boundary = str(row.get("boundary", "")).lower().strip()
    tourism = str(row.get("tourism", "")).lower().strip()

    # Explicit restrictions override other rules.
    if access in PRIVATE_ACCESS_VALUES or foot in {"no", "private"}:
        return "C"

    # Explicit public/permissive access.
    if access in PUBLIC_ACCESS_VALUES or foot in PUBLIC_ACCESS_VALUES:
        return "A"

    if access in PERMISSIVE_ACCESS_VALUES or foot in PERMISSIVE_ACCESS_VALUES:
        return "B"

    # Clearly public urban green spaces.
    if leisure in {"park", "playground", "recreation_ground"}:
        return "A"

    if landuse in {"village_green", "recreation_ground"}:
        return "A"

    if tourism == "picnic_site":
        return "A"

    # Gardens can be public or private; without access tag, treat cautiously.
    if leisure == "garden":
        return "B"

    # Cemeteries / graveyards are often physically open but not recreational parks.
    if landuse == "cemetery" or amenity == "grave_yard":
        return "B"

    # Nature reserves / protected areas may be accessible but not always urban public space.
    if leisure == "nature_reserve" or boundary == "protected_area":
        return "B"

    # Sports pitches/centres can be fenced or fee-based.
    if leisure in {"sports_centre", "pitch"}:
        return "B"

    # Natural/grass/wood/wetland areas without access tags: de facto accessible at best.
    if landuse in {"grass", "forest"}:
        return "B"

    if natural in {"wood", "grassland", "scrub", "heath", "wetland"}:
        return "B"

    # Blue polygons alone do not imply public access to water edge.
    if natural == "water" or landuse in {"reservoir", "basin"}:
        return "C"

    # Conservative default.
    return "C"


def _classify_source_tag(row) -> str:
    """Create a compact text field showing the main OSM tag behind the feature."""

    tag_cols = ["leisure", "landuse", "natural", "water", "waterway", "amenity", "boundary", "tourism"]
    parts = []

    for col in tag_cols:
        val = row.get(col, None)
        if pd.notna(val):
            parts.append(f"{col}={val}")

    return ";".join(parts) if parts else "unknown"


# ---------------------------------------------------------------------
# PREPARE OSM / JRC
# ---------------------------------------------------------------------

def _prepare_osm_spaces(osm: gpd.GeoDataFrame, area_proj) -> gpd.GeoDataFrame:
    """Clean, classify, clip, and filter OSM accessible space polygons."""

    print("Preparing OSM accessible-space polygons...")

    if osm.crs is None:
        osm = osm.set_crs("EPSG:4326")

    osm = osm.to_crs(config.PROJ_CRS)

    is_line = osm.geometry.geom_type.isin(["LineString", "MultiLineString"])
    water_cols = [c for c in ["water", "waterway"] if c in osm.columns]
    if water_cols:
        is_water_line = is_line & osm[water_cols].notna().any(axis=1)
        if is_water_line.sum() > 0:
            osm.loc[is_water_line, "geometry"] = osm.loc[is_water_line, "geometry"].buffer(5.0)

    osm = _keep_polygonal_features(osm)

    # Clip to exact study area, not the bbox used for query.
    area_gdf = gpd.GeoDataFrame(geometry=[area_proj], crs=config.PROJ_CRS)

    osm = gpd.clip(osm, area_gdf)

    osm = osm[osm.geometry.notna() & (~osm.geometry.is_empty)].copy()
    osm["geometry"] = osm.geometry.make_valid()

    if osm.empty:
        raise ValueError("No OSM features remain after clipping to study area.")

    osm["area_m2"] = osm.geometry.area

    min_area_m2 = getattr(config, "MIN_ACCESSIBLE_SPACE_AREA_M2", 100)

    osm = osm[osm["area_m2"] >= min_area_m2].copy()

    if osm.empty:
        raise ValueError(
            f"No OSM accessible-space candidates remain after area >= {min_area_m2} m2."
        )

    osm["source"] = "osm"
    osm["space_type"] = osm.apply(_classify_space_type, axis=1)
    osm["access_level"] = osm.apply(_classify_access_level, axis=1)
    osm["source_tag"] = osm.apply(_classify_source_tag, axis=1)

    if "name" not in osm.columns:
        osm["name"] = None

    osm["access_score"] = (
        osm["access_level"]
        .map({"A": 1.0, "B": 0.5, "C": 0.0})
        .fillna(0.0)
    )

    keep_cols = [
        "source",
        "name",
        "space_type",
        "access_level",
        "access_score",
        "source_tag",
        "area_m2",
        "geometry",
    ]

    osm = osm[keep_cols].copy()

    print(f"  OSM spaces after clipping/filtering: {len(osm)}")

    return osm


def _empty_spaces_gdf() -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        columns=[
            "source",
            "name",
            "space_type",
            "access_level",
            "access_score",
            "source_tag",
            "area_m2",
            "geometry",
        ],
        geometry="geometry",
        crs=config.PROJ_CRS,
    )


def _prepare_jrc_reference(area_proj) -> gpd.GeoDataFrame:
    """
    Load JRC permanent water as physical blue reference.

    JRC water is not automatically treated as public access.
    It is assigned access_level = C by default.
    """

    if not os.path.exists(config.WATER_PATH):
        print(f"WARNING: WATER_PATH not found, skipping JRC reference: {config.WATER_PATH}")
        return _empty_spaces_gdf()

    print("Preparing JRC water reference polygons...")

    water = gpd.read_file(config.WATER_PATH)

    if water.crs is None:
        raise ValueError("JRC water layer has no CRS.")

    water = water.to_crs(config.PROJ_CRS)

    water = water[water.geometry.notna() & (~water.geometry.is_empty)].copy()
    water["geometry"] = water.geometry.make_valid()
    water["area_m2"] = water.geometry.area

    min_water_area_m2 = getattr(config, "MIN_WATER_AREA_M2", 1000)

    water = water[water["area_m2"] >= min_water_area_m2].copy()

    if water.empty:
        print("WARNING: No JRC water polygons after filtering.")
        return _empty_spaces_gdf()

    area_gdf = gpd.GeoDataFrame(geometry=[area_proj], crs=config.PROJ_CRS)

    water = gpd.clip(water, area_gdf)

    if water.empty:
        print("WARNING: No JRC water intersects study area.")
        return _empty_spaces_gdf()

    water["area_m2"] = water.geometry.area
    water["source"] = "jrc"
    water["name"] = None
    water["space_type"] = "blue"
    water["access_level"] = "C"
    water["access_score"] = 0.0
    water["source_tag"] = "JRC permanent water reference"

    keep_cols = [
        "source",
        "name",
        "space_type",
        "access_level",
        "access_score",
        "source_tag",
        "area_m2",
        "geometry",
    ]

    water = water[keep_cols].copy()

    print(f"  JRC reference spaces: {len(water)}")

    return water


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def build_accessible_spaces() -> gpd.GeoDataFrame:
    """Main Step 9 function."""

    print("--- STEP 9: ACCESSIBLE GREEN-BLUE SPACES ---")

    ox.settings.use_cache = True
    ox.settings.timeout = 300

    grid = _read_grid()

    print(f"Loaded grid cells: {len(grid)}")
    print(f"Grid CRS: {grid.crs}")

    # Exact study area used for final clipping.
    area_proj = _get_study_area_polygon(grid)

    # Simple query polygon used only for Overpass request.
    query_poly_wgs84 = _get_overpass_query_polygon(grid)

    print("Using simplified buffered bounding box for Overpass query.")
    print(f"OSM query buffer: {getattr(config, 'OSM_FEATURE_QUERY_BUFFER_M', 500)} m")

    osm_raw = _fetch_osm_features(query_poly_wgs84)

    osm_spaces = _prepare_osm_spaces(osm_raw, area_proj)
    jrc_reference = _prepare_jrc_reference(area_proj)

    all_spaces = pd.concat([osm_spaces, jrc_reference], ignore_index=True)

    all_spaces = gpd.GeoDataFrame(
        all_spaces,
        geometry="geometry",
        crs=config.PROJ_CRS,
    )

    all_spaces = all_spaces[
        all_spaces.geometry.notna() & (~all_spaces.geometry.is_empty)
    ].copy()

    all_spaces["geometry"] = all_spaces.geometry.make_valid()
    all_spaces["area_m2"] = all_spaces.geometry.area

    # Stable feature ID
    all_spaces = all_spaces.reset_index(drop=True)
    all_spaces["space_id"] = all_spaces.index.astype(int)

    all_spaces = all_spaces[
        [
            "space_id",
            "source",
            "name",
            "space_type",
            "access_level",
            "access_score",
            "source_tag",
            "area_m2",
            "geometry",
        ]
    ].copy()

    spaces_A = all_spaces[all_spaces["access_level"] == "A"].copy()
    spaces_AB = all_spaces[all_spaces["access_level"].isin(["A", "B"])].copy()

    print("\nAccessible-space summary:")
    print(f"  All candidate spaces : {len(all_spaces)}")
    print(f"  Level A spaces       : {len(spaces_A)}")
    print(f"  Level A+B spaces     : {len(spaces_AB)}")
    print(f"  JRC reference spaces : {(all_spaces['source'] == 'jrc').sum()}")

    print("\nSummary by Source and Access Level:")
    print(all_spaces.groupby(["source", "access_level"]).size().to_string())
    
    print("\nSummary by Source and Space Type:")
    print(all_spaces.groupby(["source", "space_type"]).size().to_string())

    print("\nTop Tags in OSM (by access level):")
    osm_only = all_spaces[all_spaces["source"] == "osm"]
    print(osm_only.groupby(["access_level", "source_tag"]).size().sort_values(ascending=False).head(20).to_string())

    print("\nSpot-check specific spaces:")
    spot_check_names = ["Thống Nhất", "Cầu Giấy", "Hồ Tây", "Yên Sở", "Hoàn Kiếm"]
    for name in spot_check_names:
        matches = all_spaces[all_spaces["name"].str.contains(name, na=False, case=False)]
        if not matches.empty:
            match_str = matches.groupby("access_level").size().to_dict()
            print(f"  {name}: {match_str} (total parts: {len(matches)})")

    print("\nRaw Polygon Area by access level, hectares (May contain overlaps):")
    area_by_access = (
        all_spaces.groupby("access_level")["area_m2"]
        .sum()
        .div(10000)
        .round(2)
    )
    print(area_by_access.to_string())

    print("\nActual Dissolved Area (No double counting), hectares:")
    
    def safe_union(gdf):
        if gdf.empty:
            return None
        try:
            return gdf.geometry.union_all()
        except AttributeError:
            return gdf.geometry.unary_union

    dissolved_A = safe_union(spaces_A)
    dissolved_AB = safe_union(spaces_AB)
    
    if dissolved_A is not None:
        print(f"  Dissolved Level A   : {(dissolved_A.area / 10000):.2f} ha")
    if dissolved_AB is not None:
        print(f"  Dissolved Level A+B : {(dissolved_AB.area / 10000):.2f} ha")

    out_path = config.ACCESSIBLE_SPACES_GPKG

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    print(f"\nSaving accessible spaces to:")
    print(out_path)

    all_spaces.to_file(out_path, driver="GPKG", layer="accessible_spaces_all")
    spaces_A.to_file(out_path, driver="GPKG", layer="accessible_spaces_A")
    spaces_AB.to_file(out_path, driver="GPKG", layer="accessible_spaces_AB")
    jrc_reference.to_file(out_path, driver="GPKG", layer="jrc_water_reference")
    
    if dissolved_A is not None and not dissolved_A.is_empty:
        gpd.GeoDataFrame(geometry=[dissolved_A], crs=config.PROJ_CRS).to_file(
            out_path, driver="GPKG", layer="accessible_spaces_A_dissolved"
        )
    if dissolved_AB is not None and not dissolved_AB.is_empty:
        gpd.GeoDataFrame(geometry=[dissolved_AB], crs=config.PROJ_CRS).to_file(
            out_path, driver="GPKG", layer="accessible_spaces_AB_dissolved"
        )

    print("\nDone.")
    print("Next step: run step10_accessibility_mapping.py")
    print(
        "\nMethod note: accessible_spaces_A and accessible_spaces_AB layers may contain"
        "\noverlapping OSM polygons (e.g. a park containing a pitch sub-polygon)."
        "\nStep 10 (_dissolve_to_single_geometry) unions all polygons before computing"
        "\nintersection area with each grid cell, preventing any double-counting of"
        "\naccessible green-blue area. No pre-dissolve is needed here."
    )

    return all_spaces


if __name__ == "__main__":
    build_accessible_spaces()