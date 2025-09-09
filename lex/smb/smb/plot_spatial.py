#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing routines to create map plots of Surface Mass Balance data.
"""
import os
from typing import OrderedDict
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from scipy.interpolate import griddata
from matplotlib import colors as c
from livvkit import elements as el
import xarray as xr
from pathlib import Path
import smb.preproc as preproc

DESCRIBE_CORE = """
Filled contour of modeled annual surface mass balance of the Greenland ice
sheet, with firn and core field estimates overlaid as filled circles. Data were
compiled from Cogley et al. (2004), Bales et al., 2009, and the PROMICE
database (Machguth et al., 2016), from both ablation and accumulation zones.
"""

DESCRIBE_IB_SPATIAL = """
Spatial plot of annual surface mass balance along IceBridge transects in the
Greenland accumulation zone. SMB is given in kg m^-2 a^-1, and was estimated from
altimetry data.
"""

DESCRIBE_IB_DIFF = """
Spatial plot of differences between modeled and observed annual surface mass
balance along IceBridge altimetry transects in the Greenland ice sheet.
Differences are given in m^-2 a^-1 Modeled SMB values are compared to accumulation
zone data compiled from Lewis et al. (2017). Blue (red) indicates locations at
which the model overestimates (underestimates) surface mass balance when
compared to annual IceBridge estimates.
"""

DESCRIBE_METADATA = """
Map describing the datasets used in model surface mass balance validation, including core/firn
locations, altimetry transects, and drainage basin delineations in the Greenland Ice Sheet.
Drainage basins were denoted by Zwally et al. (2012) are numbered in orange, and their colors
correspond to those used in other histograms and scatterplots that compare modeled and observed
SMB. Altimetry data (white lines) were obtained from the IceBridge project (Lewis et al., 2017)
and SMB values along transects represent the annual average from the temporal domain of each
transect (the oldest resolvable internal reflecting horizon varies for each airborne radar
location, but can reach as far back as the early 1700s). Core/firn estimates are colored based
on their source; dark blue circles are locations compiled by Cogley et al. (2004) and updated
by Bales et al. (2009) PARCA cores, while yellow triangles are estimates supplied by PROMICE
(Machguth et al., 2016). Core data points are sized by the length of their temporal record, with
larger points indicating annual estimates were taken from a greater number of yearly SMB values.
"""

IMG_GROUP = "Spatial"


def mali_to_contour(lon_cell, lat_cell, data_cell):
    """Convert MALI unstructured to gridded data."""
    n_cells = data_cell.squeeze().shape[0]

    # First make a regular grid to interpolate onto
    # Adjust the density of the regular grid (4 * sqrt(nCells))
    numcols = int(n_cells ** 0.5) * 2
    numrows = numcols
    _xc = np.linspace(lon_cell.min(), lon_cell.max(), numcols)
    _yc = np.linspace(lat_cell.min(), lat_cell.max(), numrows)
    x_grid, y_grid = np.meshgrid(_xc, _yc)

    # Interpolate at the points in xi, yi
    z_grid = griddata((lon_cell, lat_cell), data_cell.squeeze(), (x_grid, y_grid))
    return x_grid, y_grid, z_grid


def load_model_data(args, config, regrid=True):
    """Load Model data."""
    nc1 = Dataset(config["climo"].format(m_s=1, m_e=12, clim="ANN"), mode="r")
    nc2 = Dataset(args.latlon, mode="r")
    nc3 = Dataset(args.elevation, mode="r")

    lats_model = nc1.variables[config["latv"]][:]
    lons_model = nc1.variables[config["lonv"]][:]
    if nc2.variables[args.lonv].units == "radians":
        lats_model = np.degrees(lats_model)
        lons_model = np.degrees(lons_model)

    smb_model = nc1.variables[args.smbv][:]

    thk_model = nc1.variables[args.topov][:]
    if args.landfracv in nc1.variables:
        smb_model *= nc1.variables[args.landfracv][:]
    smb_model *= args.smbscale
    nc1.close()
    nc2.close()
    nc3.close()

    mask = thk_model.flatten() < 0.0001
    smb_flat = smb_model.flatten()
    smb_flat[np.where(mask)] = np.nan
    smb_flat[np.where(np.abs(smb_flat) > 1e10)] = np.nan
    smb_flat.shape = smb_model.shape
    msmb = ma.masked_invalid(smb_flat)  # Make sure to convert this to a masked array

    if lons_model.ndim == 1 and lons_model.shape[0] < 1441:
        lons_model, lats_model = np.meshgrid(lons_model, lats_model)
    elif regrid:
        lons_model, lats_model, msmb = mali_to_contour(lons_model, lats_model, msmb)

    return lons_model, lats_model, msmb


def gen_map(lon_0, lat_0):
    """
    Create Basemap on which to plot.

    Parameters
    ----------
    lon_0, lat_0 : float
        Center location of map to be drawn

    Returns
    -------
    pmap : basemap.Basemap
        Lambert conformal conic `Basemap` centered on Greenland
        (future: do a stereographic for Antarctica)

    """
    pmap = Basemap(
        width=2000000,
        height=3000000,
        rsphere=(6378137.00, 6356752.3142),
        resolution="l",
        projection="lcc",
        lat_0=lat_0,
        lon_0=lon_0,
    )
    return pmap


def annotate_map(pmap, axis=None):
    """
    Add common details to basemap.

    Draw grid and coastlines, and fill the continents with colour.

    Parameters
    ----------
    pmap : basemap.Basemap
        Basemap on which to draw annotations

    axis : matplotlib.pyplot.Axis, optional
        Axis where Basemap is drawn, default is None

    """
    pmap.drawcoastlines(ax=axis)
    pmap.drawcountries(ax=axis)
    pmap.fillcontinents(color="gainsboro", lake_color="aqua", zorder=1)
    pmap.drawparallels(
        np.arange(-80.0, 81.0, 5.0),
        labels=[1, 0, 0, 0],
        fontsize=10,
        color="lightgrey",
        zorder=1,
        ax=axis,
    )
    pmap.drawmeridians(
        np.arange(-180.0, 181.0, 10.0),
        labels=[0, 0, 0, 1],
        fontsize=10,
        color="lightgrey",
        zorder=1,
        ax=axis,
    )


def plot_core(args, config):
    """Plot Ice Core data on map with model data."""
    img_list = []
    _, _, _, core_file = preproc.core_file(config)
    smb_avg = pd.read_csv(core_file)

    plt.figure(figsize=(12, 12))

    # lon_0 = lons_model.mean()
    # lat_0 = lats_model.mean()
    lon_0 = -40.591
    lat_0 = 71.308

    lons_model, lats_model, msmb = load_model_data(args, config)
    pmap = gen_map(lon_0, lat_0)
    annotate_map(pmap)
    xi, yi = pmap(lons_model.squeeze(), lats_model.squeeze())
    vabsmax = 2000
    # smb = pmap.scatter(xi, yi, c=msmb.squeeze(), vmin=-1500, vmax=1500,
    #                 cmap='Spectral', zorder=2)
    smb = pmap.pcolormesh(
        xi, yi, msmb.squeeze(), vmin=-vabsmax, vmax=vabsmax, cmap="Spectral", zorder=2
    )
    cbar = pmap.colorbar(smb, location="bottom", pad="5%")
    cbar.set_label("Surface mass balance (kg m$^{-2}$ a$^{-1}$)")

    lat_obs = smb_avg["Y"].values
    lon_obs = smb_avg["X"].values
    xobs, yobs = pmap(lon_obs, lat_obs)
    smbobs = smb_avg.b
    _ = pmap.scatter(
        xobs,
        yobs,
        c=smbobs,
        vmin=-1500,
        vmax=1500,
        marker="o",
        edgecolor="black",
        cmap="Spectral",
        zorder=3,
    )

    plt.tight_layout()
    img_file = os.path.join(args.out, "core_spatial.png")
    plt.savefig(img_file)
    plt.close()

    img_link = os.path.join(
        "imgs", os.path.basename(args.out), os.path.basename(img_file)
    )
    img_elem = el.Image(
        "Modeled SMB with field overlay",
        " ".join(DESCRIBE_CORE.split()),
        img_link,
        height=args.img_height,
        group=IMG_GROUP,
        relative_to="",
    )
    img_list.append(img_elem)

    return img_list


def plot_ib_spatial(args, config):
    """Create map plots of IceBridge transects."""
    img_list = []
    _, _, ib_file = preproc.ib_outfile(config)
    ice_bridge = pd.read_csv(ib_file)
    # lons_model, lats_model, msmb = load_model_data(args, config)
    fig = plt.figure(figsize=(12, 12))
    axis = fig.add_subplot(1, 1, 1)

    # lon_0 = lons_model.mean()
    # lat_0 = lats_model.mean()
    lon_0 = -40.591
    lat_0 = 71.30

    # print(f"LON0: {lon_0} LAT0: {lat_0}")
    pmap = gen_map(lon_0, lat_0)
    annotate_map(pmap, axis)

    lat_obs = ice_bridge["Y"].values
    lon_obs = ice_bridge["X"].values
    xobs, yobs = pmap(lon_obs, lat_obs)
    # xmod, ymod = m(lons_model, lats_model)

    smbobs = ice_bridge["b"].values
    obs = pmap.scatter(
        xobs, yobs, c=smbobs, marker="o", cmap="viridis", zorder=3, lw=0, ax=axis
    )
    # _map = pmap.contour(xmod, ymod, climo_data[0], colors="k")
    # _map = pmap.contour(xmod, ymod, climo_data, color="k")

    cbar = pmap.colorbar(obs, location="bottom", pad="5%", ax=axis)
    cbar.set_label("Surface mass balance (kg m$^-2$ a$^{-1}$)")

    plt.tight_layout()
    img_file = os.path.join(args.out, "IB_spatial.png")
    plt.savefig(img_file)
    plt.close()

    img_link = os.path.join(
        "imgs", os.path.basename(args.out), os.path.basename(img_file)
    )
    img_elem = el.Image(
        "IceBridge transects",
        " ".join(DESCRIBE_IB_SPATIAL.split()),
        img_link,
        height=args.img_height,
        group=IMG_GROUP,
        relative_to="",
    )
    img_list.append(img_elem)

    fig = plt.figure(figsize=(12, 12))
    axis = fig.add_subplot(1, 1, 1)
    annotate_map(pmap, axis)

    lat_obs = ice_bridge["Y"].values
    lon_obs = ice_bridge["X"].values
    xobs, yobs = pmap(lon_obs, lat_obs)
    smbobs = ice_bridge["mod_b"].values - ice_bridge["b"].values
    obs = pmap.scatter(
        xobs,
        yobs,
        c=smbobs,
        marker="o",
        cmap="RdBu",
        zorder=3,
        vmin=-200,
        vmax=200,
        lw=0,
        ax=axis,
    )

    cbar = pmap.colorbar(obs, location="bottom", pad="5%")
    cbar.set_label(
        "Difference in surface mass balance \n (Model - IceBridge) (kg m$^-2$ a$^{-1}$)"
    )

    plt.tight_layout()
    img_file = os.path.join(args.out, "IB_spatial_difference.png")
    plt.savefig(img_file)
    plt.close()

    img_link = os.path.join(
        "imgs", os.path.basename(args.out), os.path.basename(img_file)
    )
    img_elem = el.Image(
        "IceBridge transect differences",
        " ".join(DESCRIBE_IB_DIFF.split()),
        img_link,
        height=args.img_height,
        group=IMG_GROUP,
        relative_to="",
    )
    img_list.append(img_elem)
    return img_list


def plot_metadata(args, config):
    """
    Create plot of basin locations, ice bridge transects, and core locations.
    """
    img_list = []
    _, _, ib_file = preproc.ib_outfile(config)
    _, _, _, core_file = preproc.core_file(config)
    try:
        ice_bridge = pd.read_csv(ib_file)
    except pd.errors.EmptyDataError:
        print(ib_file)
        raise
    smb_avg = pd.read_csv(core_file)
    zwally_data = pd.read_csv(Path(config["preproc_dir"], config["zwally_file"]))

    lons_model, lats_model, smb_model = load_model_data(args, config, regrid=False)
    plt.figure(figsize=(12, 12))

    forcolors = c.ListedColormap(
        [
            "moccasin",
            "steelblue",
            "darkturquoise",
            "green",
            "lightsalmon",
            "mediumpurple",
            "grey",
            "purple",
            "firebrick",
        ]
    )

    # This will be the center of our map
    # lon_0 = lons_model.mean()
    # lat_0 = lats_model.mean()
    lon_0 = -40.591
    lat_0 = 71.308
    pmap = gen_map(lon_0, lat_0)
    annotate_map(pmap)

    # Read in the zwally basins and mask out model cells that are missing a basin designation
    basins = np.floor(zwally_data.zwally_basin.values)
    mask = basins.flatten() < 0.0001
    basins[np.where(mask)] = np.nan
    basins.shape = smb_model.squeeze().shape
    mbasins = ma.masked_invalid(basins)

    # Plotting the basins pseudocolor
    if lons_model.ndim == 1 and lons_model.shape[0] < 1441:
        lons_model, lats_model = np.meshgrid(lons_model, lats_model)
    elif lons_model.shape[0] >= 1441:
        lons_model, lats_model, mbasins = mali_to_contour(
            lons_model, lats_model, mbasins
        )

    xi, yi = pmap(lons_model.squeeze(), lats_model.squeeze())
    _ = pmap.pcolormesh(xi, yi, mbasins, cmap=forcolors, zorder=2)

    basin_labels = {
        "1": (-48, 80.5),
        "2": (-29, 76),
        "3": (-32, 70),
        "4": (-40.5, 66),
        "5": (-47.5, 62),
        "6": (-49, 67.5),
        "7": (-50, 71),
        "8": (-55, 75),
    }
    for lbl, loc in basin_labels.items():
        xi, yi = pmap(loc[0], loc[1])
        plt.text(xi, yi, lbl, fontsize=16, fontweight="bold", color="#FF7900")

    # Adding in firn/core measurement locations. Size by the number of years in the record
    lat_obs = smb_avg["Y"].values
    lon_obs = smb_avg["X"].values
    xobs, yobs = pmap(lon_obs, lat_obs)
    smb_avg.loc[smb_avg["nyears"] > 50] = 50
    smb_avg.loc[smb_avg["nyears"] < 5] = 5
    _ = pmap.scatter(
        xobs, yobs, marker="o", edgecolor="black", s=4 * smb_avg["nyears"], zorder=4
    )

    # Color the ablation zone (PROMICE) measurements yellow
    smb_promice = smb_avg[smb_avg.source == "promice"].copy()
    lat_obs = smb_promice["Y"].values
    lon_obs = smb_promice["X"].values
    xobs, yobs = pmap(lon_obs, lat_obs)
    smb_promice.loc[smb_promice["nyears"] > 50] = 50
    smb_promice.loc[smb_promice["nyears"] < 5] = 5
    _ = pmap.scatter(
        xobs,
        yobs,
        marker="^",
        color="yellow",
        edgecolor="black",
        s=4 * smb_promice["nyears"],
        zorder=4,
    )

    # Add IceBridge transects
    lat_obs = ice_bridge["Y"].values
    lon_obs = ice_bridge["X"].values
    xobs, yobs = pmap(lon_obs, lat_obs)
    _ = pmap.scatter(xobs, yobs, marker="o", lw=0, zorder=3, s=3, color="white")

    plt.tight_layout()
    img_file = os.path.join(args.out, "plot_meta_old.png")
    plt.savefig(img_file)
    plt.close()

    img_link = os.path.join(
        "imgs", os.path.basename(args.out), os.path.basename(img_file)
    )
    img_elem = el.Image(
        "Data locations",
        " ".join(DESCRIBE_METADATA.split()),
        img_link,
        height=args.img_height,
        group=IMG_GROUP,
        relative_to="",
    )
    img_list.append(img_elem)

    return img_list