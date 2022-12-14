#%%
"""
check_cfrad_xrain.py ver 0.1 coded by A.NISHII (Nagoya Univ., Japan)
Conv_xrain2cfrad.pyで正常にXRAINデータを変換できたかどうかチェックする
Test code for checking if output from Conv_xrain2cfrad.py is ok.

HISTORY(yyyy/mm/dd)
2022.12.14 First created by A.NISHII

"""


import matplotlib.pyplot as plt
import pyart

# %%
fname = 'cfrad.USHIO00000-20140820-0101-EL010000.nc' #ファイル名を指定(Set Cfrad)

#%%
print('Input filename: '+fname)
radar = pyart.io.read_cfradial(fname)
display = pyart.graph.RadarDisplay(radar)
# %%
#display.scan_type = 'vpt'
fig = plt.figure(figsize=(16,8))
ax_ph = fig.add_subplot(241)

ax_phmti = fig.add_subplot(242)
ax_pv = fig.add_subplot(243)
ax_pvmti = fig.add_subplot(244)
ax_vel = fig.add_subplot(245)
ax_width = fig.add_subplot(246)
ax_pdp = fig.add_subplot(247)
ax_rhv = fig.add_subplot(248)

display.plot_ppi('DBMHC', 0, ax=ax_ph, vmin=-110,vmax=-60, colorbar_label='')
ax_ph.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('PHMTI', 0, ax=ax_phmti, vmin=-110,vmax=-60, colorbar_label='')
ax_phmti.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('DBMVC', 0, ax=ax_pv, vmin=-110,vmax=-60, colorbar_label='')
ax_pv.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('PVMTI', 0, ax=ax_pvmti, vmin=-110,vmax=-60, colorbar_label='')
ax_pvmti.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('VEL', 0, ax=ax_vel, vmin=-20,vmax=20, colorbar_label='')
ax_vel.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('WIDTH', 0, ax=ax_width, vmin=0,vmax=10, colorbar_label='')
ax_width.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('PHIDP', 0, ax=ax_pdp, vmin=0,vmax=120, colorbar_label='')
ax_pdp.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('RHOHV', 0, ax=ax_rhv, vmin=0.8,vmax=1.0, colorbar_label='')
ax_rhv.set(aspect=1,xlabel='',ylabel='')

#plt.show()
plt.savefig('chk_raw.jpg',bbox_inches='tight')

# %%
#display.scan_type = 'vpt'
fig = plt.figure(figsize=(12,8))
ax_zh = fig.add_subplot(231)
ax_zdr = fig.add_subplot(232)
ax_kdp = fig.add_subplot(233)
ax_rain = fig.add_subplot(234)
ax_qf = fig.add_subplot(235)

display.plot_ppi('DBZ', 0, ax=ax_zh, vmin=10,vmax=50, colorbar_label='')
ax_zh.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('ZDR', 0, ax=ax_zdr, vmin=-3,vmax=3, colorbar_label='')
ax_zdr.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('KDP', 0, ax=ax_kdp, vmin=-1,vmax=10.0, colorbar_label='')
ax_kdp.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('RRR', 0, ax=ax_rain, vmin=0,vmax=120, colorbar_label='')
ax_rain.set(aspect=1,xlabel='',ylabel='')
display.plot_ppi('QF', 0, ax=ax_qf, vmin=0,vmax=32, colorbar_label='')
ax_qf.set(aspect=1,xlabel='',ylabel='')

#plt.show()
plt.savefig('chk_intermidiate.jpg',bbox_inches='tight')

print('FINISH')
