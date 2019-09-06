#!/usr/bin/env python3

import os
import netCDF4
import matplotlib.pylab as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature


# Location of ECMWF data should be specified in GBCS_NWP_PATH environment variable
nwpbase = os.getenv('GBCS_NWP_PATH', '/badc/ecmwf-era-interim/data')
nwppath = os.path.join(nwpbase, 'gg/as/{Y}/{m}/{d}/ggas{Y}{m}{d}{H}00.nc')


# Increase default font size from 10 to 12
plt.rcParams.update({'font.size': 12})


def cart_fig(figsize=None, projection=ccrs.PlateCarree()):
    fig=plt.figure(figsize=figsize)
    axes_class = GeoAxes, dict(map_projection=projection)
    axes = AxesGrid(fig, 111, (1,1),
                   axes_class=axes_class,
                   axes_pad=0.6,
                   cbar_location='right',
                   cbar_mode='single',
                   cbar_pad=0.2,
                   cbar_size='3%',
                   label_mode='')  # note the empty label_mode
    return fig, axes


def cart_grid(ax):
    gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                    color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator(np.linspace(-180,180,7))
    gl.ylocator = mticker.FixedLocator(np.linspace(-90,90,7))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

# Plot locations of ARC training profiles

dtype = [('date','S10'),
         ('lon','d'),
         ('lat','d'),
         ('t','d'),
         ('p','d'),
         ('u10','d'),
         ('v10','d'),
         ('ice','d'),
         ('tcwv','d')]

dat = np.loadtxt('arc_surface.txt', dtype)
lon = dat['lon']
lon[lon > 180] = lon[lon>180] - 360

fig, axes = cart_fig(figsize=(10,5))
axes[0].set_global()
axes[0].coastlines()
cart_grid(axes[0])
p=axes[0].scatter(dat['lon'], dat['lat'], c=dat['tcwv'],
                  vmin=0, vmax=65,
                  transform=ccrs.PlateCarree())
cb = axes.cbar_axes[0].colorbar(p)
cb.set_label_text('Total Column Water Vapour / kg m$^{-2}$')
#plt.tight_layout()
fig.savefig('ARC_profile_locations.png')
fig.savefig('ARC_profile_locations.pdf')


# Plot ECMWF cloud prior
time=dict(Y='2000',m='12',d='30',H='12')
nwp=netCDF4.Dataset(nwppath.format(**time))
lat = nwp.variables['latitude'][:]
lon = nwp.variables['longitude'][:]
lon[lon >= 180] = lon[lon>=180] - 360
tcc = nwp.variables['TCC'][0,0]

lon = np.roll(lon,256)
tcc = np.roll(tcc,256,1)
pclr = 1 - np.clip(tcc,0.5,0.95)


fig, axes = cart_fig(figsize=(10,5))
axes[0].set_global()
axes[0].coastlines()
cart_grid(axes[0])
axes[0].add_feature(cfeature.LAND)
axes[0].add_feature(cfeature.COASTLINE)
axes[0].set_title('Prior Probability of clear ({Y}-{m}-{d} {H}:00)'.format(**time))
p=axes[0].imshow(pclr, vmin=0, vmax=0.5, transform=ccrs.PlateCarree(),
            extent=(-180,180,-90,90), interpolation='none', origin='upper')
cb = axes.cbar_axes[0].colorbar(p)
cb.set_label_text('Probability')
#plt.tight_layout()
fig.savefig('Bayes_prior_PClear.png')
fig.savefig('Bayes_prior_PClear.pdf')
