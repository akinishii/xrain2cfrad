#%%
"""
Conv_XRAIN2CFrad.py ver 1.2 coded by A.NISHII (Nagoya Univ., Japan)
Convert XRAIN raw and intermediate data to CF-radial ver 1.5

USEAGE
python3 Conv_XRAIN2Cfrad.py path/to/raw(P008)_file
*Archives of raw files (*P008*.tgz) and intermediate files (*R005*.tgz) must be in the same directory
*Converted cfradial file is saved in the current directory.
*Any meta data (e.g., radar coefficient, pulse width) is not output in the current version. 
*You can change the output directory of converted files by changing Converter._outdir

Format of the input file name: cfrad.site_name-yyyymmdd-hhmmss-ELxxxxxx-DEGyyy.nc
                               *xxxxxx: elevation number of input file
                               *yyy   : tenfold of the elevation angle

HISTORY(yyyy/mm/dd)
2022/10/30 ver 0.1 (First created) by A.NISHII
2022/12/14 ver 0.2 Added quality flag in output by A.NISHII
2024/10/07 ver 1.0 Bug fixed by A.NISHII
2024/10/07 ver 1.1 Modified functions for setting output dir
2024/10/14 ver 1.2 Modified the format of instrument namme

MIT License
Copyright (c) 2022 Akira NISHII

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice (including the next paragraph) shall be 
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import netCDF4
import numpy as np
import struct
import datetime
import tarfile
from os.path import basename
from os import makedirs
from sys import argv, exit

#%%
class Converter:
    _FillValueU8  = 255
    _FillValueU16 = 0
    _FillValueF32 = -327.68
    _outdir = '../out_nc' #Save directory for a converted netCDF file.

    def __init__(self,mode,flag_overwrite=False):
        self.mode = mode
        self.flag_ow = False
        if mode == 0:
            #P008 mode
            self.n_var = 8 #Fixed to 8
            #varname_xrain:[varname_cf,standard (long) name,unit]
            self.varinfo = {'PHN0':['DBMHC','log_power_co_polar_h','dbm'],
                            'PHM0':['PHMTI','log_power_co_polar_h_mti','dbm'],
                            'PVN0':['DBMVC','log_power_co_polar_v','dbm'],
                            'PVM0':['PVMTI','log_power_co_polar_v_mti','dbm'],
                            'PPDP':['PHIDP','differential_phase_hv','degrees'],
                            'PRHV':['RHOHV','cross_correlation_ratio_hv','unitless'],
                            'PV00':['VEL',  'radial_velocity_of_scatterers_away_from_instrument','m/s'],
                            'PW00':['WIDTH','doppler_spectrum_width','m/s']}
        else:
            #Z005 mode
            self.n_var = 4 #Fixed to 4 (QF is added later)
                        #varname_xrain:[varname_cf,standard (long) name,unit]
            self.varinfo = {'RZH0':['DBZ','equivalent_reflectivity_factor','DBZ'],
                            'RZDR':['ZDR','log_differential_reflectivity_hv','DB'],
                            'RKDP':['KDP','specific_differential_phase_hv','degrees/km'],
                            'RRR0':['RRR','radar_estimated_rain_rate','mm/hr'],
                            'RQF0':['QF','radar_quality_mask','unitless']}
        

        if flag_overwrite: self.flag_ow = True


    def convert(self,fname,outdir=None):
        self.fname = fname
        #self.ncname = ncname
        if outdir is not None: self._outdir = outdir
        makedirs(self._outdir,exist_ok=True)
        self.unzip_tar()
        self.read_xrain_ppi()
        #Remove -P008 and -R005 from the input file name, and then add the tenfold of the fixed angle to the output netCDF file name.
        self.ncname = self._outdir+'/cfrad.' + basename(fname).split('.')[0].replace('-P008','').replace('-R005','') + f'-DEG{int(self.fixed_angle*10):03d}.nc'
        self.write_cfrad()
    

    def unzip_tar(self):
        with tarfile.open(self.fname,'r:gz') as tf:
            self.fnames = tf.getnames()
            tf.extractall(path='work_conv')


    def read_xrain_ppi(self):
        self.n_read = 0
        #Read header
        self.read_header(self.fnames[0])
        self.vars = np.zeros((self.n_var,self.n_ray,self.n_range),dtype='float32')
        self.pnames = []
        #Read parameters
        for fname in self.fnames:
            pname = basename(fname).split('-')[3]
            if pname == 'RQF0': 
                self.qflg = np.full((self.n_var,self.n_ray,self.n_range),fill_value=255,dtype='uint8')
                self.qflg = self.read_values_uint8(fname)
                continue
            self.pnames.append(pname)
            self.vars[self.n_read] = self.read_values(fname, pname)
            self.n_read += 1


    def read_header(self,fname):
        #Encodes of characters are iso2022_jp.
        enc = 'iso2022_jp'
        with open('./work_conv/' + fname,'rb') as f:
            flag = struct.unpack('>B',f.read(1))[0] #Check whether XRAIN PPI or not
            #print(flag)
            if int(flag) != 253:
                print("ERROR: Input file is not supported in this program!")
                exit(1)
            
            #Read header
            f.seek(0,0)
            header_buf = f.read(512)
            header = struct.unpack('>512B',header_buf)
            self.header = header

            #Read obs date
            f.seek(8,0)
            buf = f.read(16)
            date_str = struct.unpack('>16s',buf)[0].decode(encoding=enc)[0:10]
            self.date_str = date_str

            #Read fixed elevation angle of the PPI scan
            f.seek(48,0)
            buf = f.read(2)
            self.fixed_angle = struct.unpack('>h',buf)[0] / 100.

            #62~173 bytes: Radar informations
            f.seek(62,0)
            rinfo_buf = f.read(112)
            rinfo = struct.unpack('>6H2I22H2B8s8s4L3HL4B',rinfo_buf)
            self.rinfo = rinfo

        #ms4bit: most significant 4 bits, ls4bit: least significant 4 bits
        ms4bit, ls4bit = self.split_byte(header[2])
        self.datacode = ms4bit
        self.paramcode = header[3]

        self.tzone = self.decode_bcd(header[28])
        
        buf = self.decode_bcd(header[40])
        buf2 = self.decode_bcd(header[41])
        self.anntena_rotetion_speed = 100. * buf + buf2 / 10.

        self.location = np.empty(3,dtype='float64')  #List for the anntena location [lon, lat, altitude]
        self.prf = np.empty(2,dtype='float64') #List for pulse repetation frequencies [PRF1, PRF2]

        self.location[1] = rinfo[0] + rinfo[1] / 60. + rinfo[2] / 3600. #Longitude of radar 
        self.location[0] = rinfo[3] + rinfo[4] / 60. + rinfo[5] / 3600. #Latitude of radar 
        self.location[2] = rinfo[6] / 100. #Altitude of radar

        stime_str = rinfo[32].decode(encoding=enc)
        dt_2utc = datetime.timedelta(hours=9)
        self.sdate = datetime.datetime.strptime(date_str+stime_str,"%Y.%m.%d%H.%M.%S") - dt_2utc

        etime_str = rinfo[33].decode(encoding=enc)
        self.edate = datetime.datetime.strptime(date_str+etime_str,"%Y.%m.%d%H.%M.%S") - dt_2utc

        self.rconst_h = (rinfo[12]-32768) / 100. #Radar coefficient (Horizontal wave)
        self.nlev_hs = (rinfo[13]-32768) / 100.  #Noise level of short pulse (Horizontal)
        self.nlev_hl = (rinfo[14]-32768) / 100.  #Noise level of long pulse (Horizontal)

        self.rconst_v = (rinfo[19]-32768) / 100. #Radar coefficient (Vertical wave)
        self.nlev_vs = (rinfo[20]-32768) /  100. #Noise level of short pulse (Horizontal)
        self.nlev_vl = (rinfo[21]-32768) /  100. #Noise level of long pulse (Horizontal)

        self.freq = rinfo[22] #Radar frequency in MHz

        self.pwidth_s = rinfo[23] / 100. #Pulse width (short pulse)
        self.pwidth_l = rinfo[24] / 100. #Pulse width (long pulse)
        
        self.prfmode = rinfo[39] #prfmode (1:Single PRF 2:Dual(staggared) PRF)
        if self.prfmode == 2:
            self.prf[0] = rinfo[26]
            self.prf[1] = rinfo[27]
        else:
            self.prf[0] = rinfo[25]

        self.range_reso = rinfo[36] / 100. #Range resolution in meter
        self.n_range = rinfo[37] #Number of bins for each ray
        self.n_ray = rinfo[38] #Number of rays


    def read_values(self,fname,pname):
        with open('./work_conv/' + fname,'rb') as f:
            buf = f.read(512)
            f.seek(512,0)
            size = (16 + 2*self.n_range) * self.n_ray
            fmt = '>' + '2H2hIi{0}H'.format(self.n_range) * self.n_ray
            buf = f.read(size)
            rdata = np.array(struct.unpack(fmt,buf),dtype='int64').reshape((self.n_ray,self.n_range+6))

        if self.n_read == 0:
            rdata[abs(rdata[:,0]-rdata[:,1])>35000,1] += 36000
            self.az = ((rdata[:,0] + rdata[:,1]) / 200.).astype('float64')
            self.az[self.az>360] -=360
            self.el = ((rdata[:,2] + rdata[:,3]) / 200.).astype('float64')
            self.el = ((rdata[:,2] + rdata[:,3]) / 200.).astype('float64')
            self.nyq = rdata[:,4] * (10.**rdata[:,5])

        if (pname == 'PW00') or (pname == 'RRR0'):
            values = ((rdata[:,6:] - 1) / 100.).astype('float32')
        elif pname == 'PRHV':
            values = ((rdata[:,6:] - 1) / 65533.).astype('float32')
        elif pname == 'PPDP':
            values = (360. * (rdata[:,6:] - 1) / 65534.).astype('float32')
        else:
            values = ((rdata[:,6:] - 32768) / 100.).astype('float32')

        return values

    #Function to decode quality flag (Currently not used)
    def read_values_uint8(self,fname):
        with open('./work_conv/' + fname,'rb') as f:
            buf = f.read(512)

            f.seek(512,0)
            size = (16 + 2*self.n_range) * self.n_ray
            fmt = '>' + '2H2hIi{0}B'.format(self.n_range) * self.n_ray
            buf = f.read(size)
            rdata = np.array(struct.unpack(fmt,buf),dtype='int64').reshape((self.n_ray,self.n_range+6))

        values = rdata[:,6:].astype('uint8')

        return values


    def write_cfrad(self):
        mode = 'w'
        if self.flag_ow : mode = 'a'
        self.cf = netCDF4.Dataset(self.ncname,mode=mode,format='NETCDF4')
        if not self.flag_ow:
            self.write_cf_gl_attr()
            self.write_cf_dims()
            self.write_cf_gl_vars()
            self.write_cf_coord_vars()
            self.write_cf_loc_vars()
            self.write_cf_azelnyq()
            self.write_cf_sweep_vars()
        self.write_cf_data_vars()
        if self.mode == 1:
            self.write_cf_qflg()
        self.cf.close()


    def write_cf_gl_attr(self):
        nc = self.cf

        nc.setncattr("Conventions","Cf/Radial instrument_parameters")
        nc.setncattr("version", "1.5")
        nc.setncattr("title","XRAIN raw and intermediated data")
        nc.setncattr("institution","Ministry of Land, Infrastructure, Transport and Tourism, Japan")
        nc.setncattr("references","Converted by Conv_xrain2cfrad.py coded by A.NISHII (Nagoya Univ.)")
        nc.setncattr("source",self.fname)
        nc.setncattr("history","")
        nc.setncattr("comment","")
        rname = basename(self.fname).split('-')[0].replace('0','')
        nc.setncattr("instrument_name",'XRAIN_'+rname)
        nc.setncattr('site_name',rname)
        nc.setncattr('scan_name','')


    def write_cf_dims(self):
        nc = self.cf

        nc.createDimension('time',self.n_ray)
        nc.createDimension('range',self.n_range)
        nc.createDimension('sweep',1)
        nc.createDimension('string_length',None)


    def write_cf_gl_vars(self):
        nc = self.cf
        volume_number = nc.createVariable('volume_number',np.dtype('int32').char)
        volume_number[:] = 0

        time_start = nc.createVariable('time_converge_start','S1',('string_length'))
        time_start[:] = np.array(self.sdate.strftime("%Y-%m-%dT%H:%M:%SZ"),dtype="S20")
        time_start.long_name = "data_volume_start_time_utc"
        time_start.units = "unitless"

        time_end = nc.createVariable('time_converge_end','S1',('string_length'))
        time_end[:] = np.array(self.edate.strftime("%Y-%m-%dT%H:%M:%SZ"),dtype="S20")
        time_end.long_name = "data_volume_end_time_utc"
        time_start.units = "unitless"

        time_ref = nc.createVariable('time_reference','S1',('string_length'))
        time_ref[:] = np.array(self.sdate.strftime("%Y-%m-%dT%H:%M:%SZ"),dtype="S20")
        time_ref.long_name = "time_reference"
        time_ref.units = "unitless"


    def write_cf_coord_vars(self):
        nc = self.cf

        time = nc.createVariable('time',np.dtype('double').char,('time'))
        obsdur = (self.edate - self.sdate).total_seconds()
        seconds = np.linspace(0,obsdur,self.n_ray,dtype='double')
        time[:] = seconds
        time.standard_name = 'time'
        time.long_name = 'time_in_seconds_since_volume_start'
        time.units = "seconds since {0}Z".format(self.sdate.strftime("%Y-%m-%dT%H:%M:%SZ"))
        time.calendar = "gregorian"
        
        radar_range = nc.createVariable('range',np.dtype('float32').char,('range'))
        range_array = (np.arange(0,self.n_range)+0.5) * self.range_reso
        radar_range[:] = range_array
        radar_range.standard_name = 'projection_range_coordinate'
        radar_range.long_name = 'range_to_measurement_volume'
        radar_range.units ='meters'
        radar_range.spacing_is_constant = 'true'
        radar_range.meters_to_center_of_first_gate = range_array[0]
        radar_range.axis = 'radial_range_coordinate'
    

    def write_cf_loc_vars(self):
        nc = self.cf

        lat = nc.createVariable('latitude',np.dtype('double').char)
        lat[:] = self.location[1]
        lat.units = 'degrees_north'

        lon = nc.createVariable('longitude',np.dtype('double').char)
        lon[:] = self.location[0]
        lon.units = 'degrees_east'

        alt = nc.createVariable('altitude',np.dtype('double').char)
        alt[:] = self.location[2]
        alt.units = 'meters'


    def write_cf_sweep_vars(self):
        nc = self.cf

        sweep_n = nc.createVariable('sweep_number',np.dtype('int32').char,('sweep'))
        sweep_n[:] = 0

        sweep_mode = nc.createVariable('sweep_mode','S1',('sweep','string_length'))
        sweep_mode[:] = 'sector'

        fx_angle = nc.createVariable('fixed_angle',np.dtype('float32').char,('sweep'))
        fx_angle[:] = self.fixed_angle
        fx_angle.units = 'degrees'

        sweep_start_ray = nc.createVariable('sweep_start_ray_index',np.dtype('int32').char,('sweep'))
        sweep_start_ray[:] = 0

        sweep_end_ray = nc.createVariable('sweep_end_ray_index',np.dtype('int32').char,('sweep'))
        sweep_end_ray[:] = len(self.el) - 1


    def write_cf_azelnyq(self):
        nc = self.cf

        az = nc.createVariable('azimuth',np.dtype('float32').char,('time'))
        az[:] = np.round(self.az,1)
        az.standard_name = 'ray_azimuth_angle'
        az.long_name = 'azimuth_angle_from_true_north'
        az.units = 'degrees'
        az.axis = 'radial_azimuth_coordinate'
        
        el = nc.createVariable('elevation',np.dtype('float32').char,('time'))
        el[:] = np.round(self.el,1)
        el.standard_name = 'ray_elevation_angle'
        el.long_name = 'elavation_angle_from_horizontal_plane'
        el.units = 'degrees'
        el.axis = 'radial_elevation_coordinate'

        nyq = nc.createVariable('nyquist_velocity',np.dtype('float32').char,('time'))
        #Set the extended nyquist velocity for dual-PRF mode because velocity unfolding is already applied.
        if self.prfmode == 2:
            self.nyq_unfolded = self.calc_nyq_unfolded()
        else:
            self.nyq_unfolded = self.nyq[0]
        nyq[:] = round(self.nyq_unfolded,2)
        nyq.standard_name = 'nyquist_velocity'
        nyq.long_name = 'unambiguous_doppler_velocity'
        nyq.units = 'm/s'
        nyq.axis = 'radial_azimuth_coordinate'
    

    def write_cf_data_vars(self):
        nc = self.cf
        
        for n in range(self.n_var):
            varinfo = self.varinfo[self.pnames[n]]
            self.Define_variable(nc,self.vars[n],varinfo[0],np.dtype('float32').char, varinfo[1],varinfo[2])
    

    def Define_variable(self,nc,var,varname,dtype_char,standard_name,units):
        ncvar = nc.createVariable(varname,dtype_char,('time','range'),fill_value=self._FillValueF32)
        ncvar[:] = var
        ncvar.long_name = standard_name
        ncvar.standard_name = standard_name
        ncvar.units = units
        ncvar.coodinates = 'evevation azimuth range'

        return ncvar


    def write_cf_qflg(self):
        nc = self.cf
        
        varinfo = self.varinfo['RQF0']
        #print('Save QF')
        ncvar = nc.createVariable(varinfo[0],np.dtype('uint8').char,('time','range'),fill_value=self._FillValueU8)
        ncvar[:] = self.qflg
        ncvar.long_name = varinfo[0]
        ncvar.standard_name = varinfo[1]
        ncvar.units = varinfo[2]
        ncvar.coodinates = 'evevation azimuth range'


    def calc_nyq_unfolded(self):
        prf1 = self.prf[0]
        prf2 = self.prf[1]
        prf_ref = max(prf1,prf2)
        vel1 = self.nyq[0]
        vel2 = self.nyq[1]
        vel_ref = min(vel1,vel2)

        gcd = self.calc_gcd(int(prf1),int(prf2))

        return vel_ref * (prf_ref/gcd)


    def split_byte(self,inbyte):
        ms4bit = bin(inbyte>>4)
        ls4bit = bin(inbyte  & 0b00001111)

        return ms4bit, ls4bit


    def decode_bcd(self,bcd_byte):
        ord_10, ord_1 = self.split_byte(bcd_byte)
        
        dec = 10 * int(ord_10,2) + int(ord_1,2)
        return dec


    def calc_gcd(self,a,b):
        while b != 0:
            rem = a%b
            a = b
            b = rem
        return a


#%%
def main():
    fname_P008_tar = argv[1]
    fname_R005_tar = fname_P008_tar.replace('P008','R005')
    print('Input file name(Raw, P008): '+fname_P008_tar)
    print('Input file name(Intermediated, R005): '+fname_R005_tar)
    conv_P = Converter(mode=0,flag_overwrite=False)
    conv_P.convert(fname_P008_tar)
    conv_Z = Converter(mode=1,flag_overwrite=True)
    conv_Z.convert(fname_R005_tar)
    print('Convert success: Saved to '+conv_Z.ncname)

#%%
if __name__ == '__main__':
    if len(argv) != 2:
        print(f'ERROR: Invalid number of arguments: N_of_args = {len(argv)}')
        exit(1)
    main()
