#!/usr/bin/env python3
import os
import matplotlib.pylab as plt
import cartopy.crs as ccrs
from matplotlib.gridspec import GridSpec

import netCDF4


gbcsver = 'v2.6.0'
## JASMIN internal path
path1 = '/gws/nopw/j04/cds_c3s_sst/output/{gbcsver}/{lev}/{prod}/{Y}/{m}/{d}'
path2 = '/gws/nopw/j04/cds_c3s_sst/public/data/ICDR_v2/Analysis/L4/v2.0/{Y}/{m}/{d}'

## OpenDAP access
path = 'http://dap.ceda.ac.uk/thredds/dodsC/dap/neodc/esacci/sst/data/CDR_v2/{group}/{lev}/v2.1/{prod}/{Y}/{m}/{d}'

# Increase default font size from 10 to 12
plt.rcParams.update({'font.size': 12})


def fullpath(filename):
    parts = filename.split('-')
    prod = parts[4]
    if 'ATSR' in prod:
        group = 'ATSR'
    elif 'AVHRR' in prod:
        group = 'AVHRR'
    elif 'OSTIA' in prod:
        group = 'Analysis'
        prod = '.'
    subs = {
        'Y': filename[0:4],
        'm': filename[4:6],
        'd': filename[6:8],
        'group': group,
        'lev': parts[2].split('_')[0],
        'prod': prod,
        'gbcsver': gbcsver,
        }
    # path = path2 if prod == 'OSTIA' else path1
    return os.path.join(path.format(**subs), filename)


l2p = '20160107112259-ESACCI-L2P_GHRSST-SSTskin-AVHRR19_G-CDR2.1-v02.0-fv01.0.nc'
l3u = '20160107112259-ESACCI-L3U_GHRSST-SSTskin-AVHRR19_G-CDR2.1-v02.0-fv01.0.nc'
l3c = '20160107120000-ESACCI-L3C_GHRSST-SSTskin-AVHRR19_G-CDR2.1_day-v02.0-fv01.0.nc'
l3n = '20160107120000-ESACCI-L3C_GHRSST-SSTskin-AVHRR19_G-CDR2.1_night-v02.0-fv01.0.nc'
l4 = '20160107120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc'

with netCDF4.Dataset(fullpath(l2p)) as nc:
    st2p = nc.variables['sea_surface_temperature'][0]
    mask = nc.variables['l2p_flags'][0]

with netCDF4.Dataset(fullpath(l3u)) as nc:
    st3u = nc.variables['sea_surface_temperature'][0]

with netCDF4.Dataset(fullpath(l3c)) as nc:
    st3c = nc.variables['sea_surface_temperature'][0]

with netCDF4.Dataset(fullpath(l3n)) as nc:
    st3n = nc.variables['sea_surface_temperature'][0]

with netCDF4.Dataset(fullpath(l4)) as nc:
    st4 = nc.variables['analysed_sst'][0]

lsm = mask & 2

fig = plt.figure(figsize=(10,14))
gs = GridSpec(3,4)
ax0 = plt.subplot(gs[:,0])
ax0.set_title('L2P')
ax0.xaxis.set_ticks([])
ax0.yaxis.set_ticks([])
im0 = ax0.imshow(st2p[2000:6000], vmin=270, vmax=305, interpolation='nearest', aspect='auto')
ax0.contour(lsm[2000:6000], levels=[1], colors='k')

ax1 = plt.subplot(gs[0,1:], projection=ccrs.PlateCarree())
ax1.set_title('L3U')
ax1.coastlines('50m')
im1 = ax1.imshow(st3u, transform=ccrs.PlateCarree(), extent=(-180,180,-90,90), vmin=270, vmax=305)

ax2 = plt.subplot(gs[1,1:], projection=ccrs.PlateCarree())
ax2.set_title('L3C')
ax2.coastlines('50m')
im2 = ax2.imshow(st3c, transform=ccrs.PlateCarree(), extent=(-180,180,-90,90), vmin=270, vmax=305)
im2 = ax2.imshow(st3n, transform=ccrs.PlateCarree(), extent=(-180,180,-90,90), vmin=270, vmax=305)

ax3 = plt.subplot(gs[2,1:], projection=ccrs.PlateCarree())
ax3.set_title('L4')
ax3.coastlines('50m')
im3 = ax3.imshow(st4, transform=ccrs.PlateCarree(), extent=(-180,180,-90,90), vmin=270, vmax=305)
plt.tight_layout()
plt.colorbar(im3, ax=[ax0,ax1,ax2,ax3], orientation='horizontal', label='SST / K', fraction=0.1, pad=0.02, extend='both')
plt.savefig('PUGS_levels.png')
plt.savefig('PUGS_levels.pdf')
