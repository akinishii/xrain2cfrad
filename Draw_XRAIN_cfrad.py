#%%
"""
Draw_XRAIN_cfrad.py ver 1.4 coded by A.NISHII (Nagoya Univ., Japan)
Draw XRAIN PPI data from a cf-radial file converted by Conv_XRAIN2Cfrad.py

USEAGE
1. Setting
*Set parameters in fron around L30-L70

2. Run
python3 Draw_XRAIN_cfrad.py

HISTORY(yyyy/mm/dd)
2024/10/07 ver 1.0 (First created) by A.NISHII
2024/10/12 ver 1.1 Bug fixed & Drawing Doppler width implemented by A.NISHII
2024/11/19 ver 1.2 Bug fixed by A.NISHII
2025/01/10 ver 1.3 Modified order of parameters
2025/06/04 ver 1.4 Modified colormaps

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl
import pyart
import cartopy.crs as ccrs
from sys import argv
from os import makedirs
from os.path import basename
from cftime import num2date
from datetime import timedelta
import warnings
warnings.simplefilter('ignore')

#%%
###Settings###
fname = './testdata/nc/cfrad.TAMURA0000-20220525-1635-EL010000-DEG004.nc'
use_flist = False #True: Read files from flistname (You can process multiple files at once.) False: Read a file set in fname
flistname = 'flist2.txt'

extent=[140.4,37.5,140.6,37.7] ##[min_lon,min_lat,max_lon,max_lat]
lon_ticks = np.arange(140.4,140.61,0.1) #Locations of lon ticks
lat_ticks = np.arange(37.5,37.71,0.1)   #Locations of lat ticks

figdir = './testdata/fig_ppi_ql' #Directory to save figures

#Draw distances from the radar at REF and VEL maps
draw_cirle = True
circle_dis = np.arange(0,81,20)
draw_hair  = True
hairsize = 40 #size of hair (corresponding to plt.scatter's size)

#Set colors for REF
clevsz = np.arange(0,71,5)
cmapz="ChaseSpectral"
ticksz = clevsz[::2] #Locs of ticks in colorbar (can be same as clevs)
##Set colors for VEL
clevsv=np.arange(-20,21,2)
cmapv="balance"
ticksv = clevsv[::2] #Locs of ticks in colorbar (can be same as clevs)
##Set colors for ZDR
clevszdr=np.arange(-2,5.1,0.5)
cmapzdr="HomeyerRainbow"
tickszdr = clevszdr[::2] #Locs of ticks in colorbar (can be same as clevs)
##Set colors for KDP
clevskdp=np.arange(-2,5.1,0.5)
cmapkdp="HomeyerRainbow"
tickskdp = clevskdp[::2] #Locs of ticks in colorbar (can be same as clevs)
##Set colors for RHOV
clevsrhv=np.array([0.7, 0.75, 0.8, 0.85, 0.9,
                   0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0])
cmaprhv="viridis"
ticksrhv = clevsrhv #Locs of ticks in colorbar (can be same as clevs)
##Set colors for W
clevsw=np.arange(0,10.1,1.0)
cmapw="HomeyerRainbow"
ticksw = clevsw #Locs of ticks in colorbar (can be same as clevs)

mask_tuple=None #Conditions for masking (see documents of Pyart for detail)


###End of Settings###

#%%
def parse_flist(flistname):
    with open(flistname,'r') as fp:
        lines = fp.readlines()
        files = []
        for l in lines:
            files.append(l.replace('\n',''))
    return files

def Draw_SupTitle(fig,radar,scnum,utc_to_jst=False,fontsize=16,y=0.98):
    #from cftime import num2date
    #from datetime import timedelta
    sidx = radar.get_start(scnum)
    el = radar.fixed_angle['data'][scnum]
    date_dt = num2date(radar.time['data'][sidx],radar.time['units'])
    if utc_to_jst:
        date_str = (date_dt+timedelta(hours=9)).strftime('%Y-%m-%d %H:%M:%SJST')
    else:
        date_str = date_dt.strftime('%Y-%m-%d %H:%M:%SZ')
    site = radar.metadata['instrument_name']
    title = f'{site} {el:.1f} deg. {date_str}'
    fig.suptitle(title,y=y,fontsize=fontsize)

def Savefig(radar,scnum,figdir,dpi=300,save_jst=False,save_svg=False):
    #from cftime import num2date
    #from datetime import timedelta

    sidx = radar.get_start(scnum)
    el = radar.fixed_angle['data'][scnum]
    date_dt = num2date(radar.time['data'][sidx],radar.time['units'])
    if save_jst:
        date_str = (date_dt+timedelta(hours=9)).strftime('%Y%m%d%H%M%SJST')
    else:
        date_str = date_dt.strftime('%Y%m%d%H%M%SZ')
    site = radar.metadata['instrument_name']

    ext = 'svg' if save_svg else 'jpg'

    figname = f'{figdir}/{site}_{date_str}_{el*10:03.0f}.{ext}'
    plt.savefig(figname,dpi=dpi,bbox_inches='tight')
    print(f'fig saved to {figname}')

def Def_cmap_norm(colormap,clevs,extend='both',under_white=False):
    #Create discrete colormap
    #under_while: not draw (i.e. draw with white below levels[0])
    cmap = plt.get_cmap(colormap).copy()
    if under_white: cmap.set_under('white')
    norm = colors.BoundaryNorm(clevs,cmap.N,extend=extend)

    return cmap,norm

#%%
#Parse flist if use_flist=True
if use_flist:
    fnames = parse_flist(flistname)
else:
    fnames = [fname]

#cmap for REF
cmapz_draw, normz = Def_cmap_norm(cmapz,clevsz,under_white=True)
cmapv_draw, normv = Def_cmap_norm(cmapv,clevsv)
cmapzdr_draw, normzdr = Def_cmap_norm(cmapzdr,clevszdr)
cmapkdp_draw, normkdp = Def_cmap_norm(cmapkdp,clevskdp)
cmaprhv_draw, normrhv = Def_cmap_norm(cmaprhv,clevsrhv)
cmapw_draw, normw = Def_cmap_norm(cmapw,clevsw,under_white=True)

#Slightly extends extent to draw labels at corners
extent[0]-=0.001
extent[1]-=0.001
extent[2]+=0.001
extent[3]+=0.001

makedirs(figdir,exist_ok=True)

#%%
for f in fnames:
    print('Draw file ',f)
    radar = pyart.io.read_cfradial(f)
    display = pyart.graph.RadarMapDisplay(radar)
    #Definition of a plot coodinate.
    rlon = radar.longitude['data'][0]
    rlat = radar.latitude['data'][0]
    projection = ccrs.LambertConformal(central_latitude=rlat, central_longitude=rlon)

    #Define a figure
    fig,axes = plt.subplots(3,2,figsize=(7,8),subplot_kw={'projection':projection})
    #fig,axes = plt.subplots(3,2,figsize=(7,8))

    #REF(DBZ)
    print('DBZ')
    cm = display.plot_ppi_map("DBZ",0,mask_tuple=mask_tuple,cmap=cmapz_draw,norm=normz,resolution='10m',
                            min_lon=extent[0],min_lat=extent[1],max_lon=extent[2],max_lat=extent[3],
                            fig=fig,ax=axes[0][0],raster=True,
                            lon_lines=lon_ticks,lat_lines=lat_ticks,colorbar_label='dBZ',ticks=ticksz)
    axes[0][0].set(xlabel='',ylabel='',title='REF')
    #Draw circles and a hair on REF
    if draw_hair:
        axes[0][0].scatter(rlon,rlat,marker='+',s=hairsize,color='k',lw=1.0,transform=ccrs.PlateCarree())
    if draw_cirle:
        for r in circle_dis:
            if r==0: continue
            display.plot_range_ring(r,ax=axes[0][0],ls='--',col='k',lw=1.0,alpha=0.5)
            
    #VEL
    print('VEL')
    cm = display.plot_ppi_map("VEL",0,mask_tuple=mask_tuple,cmap=cmapv_draw,norm=normv,resolution='10m',
                            min_lon=extent[0],min_lat=extent[1],max_lon=extent[2],max_lat=extent[3],
                            fig=fig,ax=axes[0][1],raster=True,
                            lon_lines=lon_ticks,lat_lines=lat_ticks,colorbar_label='m s$^{-1}$',ticks=ticksv)
    axes[0][1].set(xlabel='',ylabel='',title='VEL')
    if draw_hair:
        axes[0][1].scatter(rlon,rlat,marker='+',s=hairsize,color='k',lw=1.0,transform=ccrs.PlateCarree())
    if draw_cirle:
        for r in circle_dis:
            if r==0: continue
            display.plot_range_ring(r,ax=axes[0][1],ls='--',col='k',lw=1.0,alpha=0.5)
    
    #ZDR
    print('ZDR')
    cm = display.plot_ppi_map("ZDR",0,mask_tuple=mask_tuple,cmap=cmapzdr_draw,norm=normzdr,resolution='10m',
                            min_lon=extent[0],min_lat=extent[1],max_lon=extent[2],max_lat=extent[3],
                            fig=fig,ax=axes[1][0],raster=True,
                            lon_lines=lon_ticks,lat_lines=lat_ticks,colorbar_label='dB',ticks=tickszdr)
    axes[1][0].set(xlabel='',ylabel='',title='ZDR')
    #KDP
    print('KDP')
    cm = display.plot_ppi_map("KDP",0,mask_tuple=mask_tuple,cmap=cmapkdp_draw,norm=normkdp,resolution='10m',
                            min_lon=extent[0],min_lat=extent[1],max_lon=extent[2],max_lat=extent[3],
                            fig=fig,ax=axes[1][1],raster=True,
                            lon_lines=lon_ticks,lat_lines=lat_ticks,colorbar_label='deg. km$^{-1}$',ticks=tickskdp)
    axes[1][1].set(xlabel='',ylabel='',title='KDP')
    #RHOHV
    print('RHOHV')
    cm = display.plot_ppi_map("RHOHV",0,mask_tuple=mask_tuple,cmap=cmaprhv_draw,norm=normrhv,resolution='10m',
                            min_lon=extent[0],min_lat=extent[1],max_lon=extent[2],max_lat=extent[3],
                            fig=fig,ax=axes[2][0],raster=True,
                            lon_lines=lon_ticks,lat_lines=lat_ticks,colorbar_label='dimensionless',ticks=ticksrhv)
    axes[2][0].set(xlabel='',ylabel='',title='RHOHV')
    #Hide the last axis
    print('WIDTH')
    cm = display.plot_ppi_map("WIDTH",0,mask_tuple=mask_tuple,cmap=cmapw_draw,norm=normw,resolution='10m',
                            min_lon=extent[0],min_lat=extent[1],max_lon=extent[2],max_lat=extent[3],
                            fig=fig,ax=axes[2][1],raster=True,
                            lon_lines=lon_ticks,lat_lines=lat_ticks,colorbar_label='m s$^{-1}$',ticks=ticksw)
    axes[2][1].set(xlabel='',ylabel='',title='WIDTH')

    #Add suptitle
    Draw_SupTitle(fig,radar,0,True)
    plt.tight_layout()

    #Save figure
    Savefig(radar,0,figdir,save_jst=True)

    plt.close()

print('Finished')
# %%
