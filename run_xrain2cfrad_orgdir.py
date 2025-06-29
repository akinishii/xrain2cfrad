#%%
"""
run_xrain2cfrad_orgdir.py ver 1.0 coded by A.NISHII
This script converts XRAIN raw (P008) and intermediate (R005) files to CF-Radial format.
Original data is loaded from the untered XRAIN directory downloaded from DIAS.

Dependencies:
    Code:
    - Conv_xrain2cfrad.py: A module to convert XRAIN files must be in the same directory.

    Libraries:
    - pandas: For date range generation.
    - numpy: For numerical operations.
    - netcdf4: For writing CF-Radial format files.

Useage:
    python3 run_xrain2cfrad_orgdir.py [options] SITE [SITE ...
    e.g.) python3 run_xrain2cfrad_orgdir.py FUNABASHI SHINYOKO -d 202406010000 202406010300 -i /data/original -o ./output/nc

    positional arguments:
    SITE                  Name of Sites (e.g., KITAHIRO, FUNABASHI, KUSENBU, etc.) to be processed. Multiple sites can be specified by space-separated values.

    options:
    -h, --help            show this help message and exit
    -d DATERANGE DATERANGE, --daterange DATERANGE DATERANGE
                            REQUIRED: Date range in JST (first: start_date, second: end_date) in the format YYYYmmddHHMM. If not specified, processing all dates in XRAIN directory.
    -i INPUTPATH, --inputpath INPUTPATH
                            Path of untared XRAIN directory (untared from files.tar). If not specified, the program will search in the current directory.
    -o OUTDIR, --outdir OUTDIR
                            Output directory of the CF-Radial files. If not specified, the program will save converted files in ./nc

NOTES:
    *Both P008 (kyoku1) and R005 (kyoku2) files must be in the XRAIN directory.
    *Some meta data (e.g., radar coefficient, pulse width) is not output in the current version.
    *Format of the output file name: cfrad.site_name-yyyymmdd-hhmmss-ELxxxxxx-DEGyyy.nc
                                     xxxxxx: elevation number of the oroginal file.
                                     yyy   : tenfold of the elevation angle e.g. 1.5 deg -> DEG015
                                     Datetime in filename is in JST (Japan Standard Time, UTC+9), but time recorded in the file is in UTC.

HISTORY(yyyy/mm/dd):
    2025/06/29 ver 1.0 (First created) by A.NISHII

LICENSE: MIT License
"""

import argparse
import datetime
import glob
import os

from pandas import date_range

from Conv_xrain2cfrad import Converter

#%%
# Function to analyze the directory name of the location and region
def display_arguments(args):
    """
    Displays the arguments parsed by argparse in a formatted way.
    """
    print("--- Input Arguments ---")
    
    # Convert the argparse.Namespace object to a dictionary
    args_dict = vars(args)

    # Calculate the maximum key length for alignment
    max_key_length = 0
    if args_dict: # Only proceed if the dictionary is not empty
        max_key_length = max(len(key) for key in args_dict.keys())
    
    # Calculate the maximum key length for alignment

    for key, value in args_dict.items():
        # Format and print each key-value pair
        print(f"  {key.ljust(max_key_length)} : {value}")
    
    print("------------------------------------")


def get_dirname_of_location_region(site: str) -> tuple:
    """Analyze the directory name of the location and region.
       The directory names are 10-lengths strings whose spaces are filled with 0 (location) or 0*1 (region).
       e.g.) DATE -> DATE000000, SAPPORO -> SAPPORO001
       The function returns a tuple of (location, region).

    Args:
        location (str): Location name.
    
    Returns:
        tuple: directory names of (location, region).
    """
    # A dictionary to find regions corresponding to locations
    sites_areas = {
        'KITAHIRO':'SAPPORO', 'ISHIKARI':'SAPPORO', 
        'IWANUMA':'KURIKOMA', 'WAKUYA':'KURIKOMA', 'MORIOKA':'KURIKOMA', 'TAKANOSU':'KURIKOMA', 'ICHIHASAMA':'KURIKOMA', 'ICHINOSEKI':'KURIKOMA', 
        'DATE':'FUKUSHIMA', 'TAMURA':'FUKUSHIMA', 
        'FUNABASHI':'TOUKYOU', 'KANTOU':'TOUKYOU', 'SHINYOKO':'TOUKYOU', 'UJIIE':'TOUKYOU', 'YATTAJIMA':'TOUKYOU', 
        'KYOUGASE':'NIIGATA', 'NAKANOKUTI':'NIIGATA', 
        'MIZUHASHI':'HOKURIKU', 'NOUMI':'HOKURIKU', 
        'FUJINOMIYA':'SHIZUOKA', 'KANUKI':'SHIZUOKA', 'SHIZUKITA':'SHIZUOKA', 'HAMAMATSU':'SHIZUOKA',
        'ANJOU':'NAGOYA', 'BISAI':'NAGOYA', 'SUZUKA':'NAGOYA', 
        'JUUBUSAN':'KINKI', 'KATSURAGI':'KINKI', 'ROKKO':'KINKI', 'TANOKUCHI':'KINKI',
        'KUMAYAMA':'OKAYAMA', 'TSUNEYAMA':'OKAYAMA', 
        'NOGAIBARA':'HIROSHIMA', 'USHIO':'HIROSHIMA', 
        'FURUTSUKI':'FUKUOKA', 'KAZASI':'FUKUOKA', 'KUSENBU':'FUKUOKA', 'SUGADAKE':'FUKUOKA', 
        'UKI':'KUMAMOTO', 'YAMAGA':'KUMAMOTO', 
        'SAKURAJIMA':'OOSUMI'
    }

    site_list = list(sites_areas.keys())
    if site in site_list:
        area = sites_areas[site]
        # Create directory names
        area = area.ljust(10, '0')  # Fill with '0' to make it 10 characters
        area = area[:-1] + '1'
        site = site.ljust(10, '0')  # Fill with '0' to make it 10 characters
    else:
        raise ValueError()
    
    return (site, area)


def get_date_range(start_date: str, end_date: str) -> list:
    """Generate a list of dates from start_date to end_date in the format YYYYmmddHHMM.
    
    Args:
        start_date (str): Start date in the format YYYYmmddHHMM.
        end_date (str): End date in the format YYYYmmddHHMM.
    
    Returns:
        list: List of dates in the format YYYYmmddHHMM.
    """
    start = datetime.datetime.strptime(start_date, '%Y%m%d%H%M')
    end = datetime.datetime.strptime(end_date, '%Y%m%d%H%M')
    return date_range(start, end, freq='1min')


def skip_message(code: int, input_file: str):
    """Show massages to skip the convertsion because of the issues in the input file.
    """
    if code == 1:
        print(f'ERROR: Input file name must be P008, not R005.\nInput file name: {input_file}')
    elif code == 2:
        print(f'ERROR: {input_file} Not Found!')
    elif code == 3:
        (f'ERROR: {input_file} Not Found! P008 & R005 files must in the same directory')

    print('Skip this file.\n')


def main():
    """Main function.
    """
    parser = argparse.ArgumentParser(description='Convert XRAIN raw (P008) and intermediate (R005) files to CF-Radial format. The original data is loaded from the untered XRAIN directory.\n'\
                                     'Both P008 (kyoku1) and R005 (kyoku2) files must be in the XRAIN directory. \n')
    parser.add_argument('SITE', nargs='+', type=str, help='Name of Sites (e.g., KITAHIRO, FUNABASHI, KUSENBU, etc.) to be processed. Multiple sites can be specified by space-separated values.')
    parser.add_argument('-d', '--daterange', nargs = 2, type=str, required=True, help='REQUIRED: Date range in JST (first: start_date, second: end_date) in the format YYYYmmddHHMM. If not specified, processing all dates in XRAIN directory.')
    parser.add_argument('-i', '--inputpath', type=str, help='Path of untared XRAIN directory (untared from files.tar). If not specified, the program will search in the current directory.')
    parser.add_argument('-o', '--outdir', type=str, help='Output directory of the CF-Radial files. If not specified, the program will save converted files in ./nc')
    args = parser.parse_args()

    display_arguments(args)

    err_flag = False

    # Get directory names of locations and regions
    site_dirs = []
    for site in args.SITE:
        site = site.upper()
        try:
            site_dir, area_dir = get_dirname_of_location_region(site)
            print(f'Site: {site}, Site(formatted): {site_dir} Area(formatted): {area_dir}')
        except ValueError as e:
            print(f'Warning: Site {site} not found in the dictionary. Skipping this site.')
            continue
        site_dirs.append((site_dir, area_dir))

    # Parse and get date_range
    start_date, end_date = args.daterange
    date_time_list = get_date_range(start_date, end_date)

    # Parse the XRAIN path if specified
    if args.inputpath:
        xrain_path = args.inputpath + '/XRAIN'
        if not os.path.exists(xrain_path):
            raise FileNotFoundError(f'Error: Input XRAIN path {xrain_path} does not exist. Please check the path.')
    else:
        xrain_path = './XRAIN'
    
    # Parse output directory if specified
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = None

    # Convert each site
    print('Start converting XRAIN data to CF-Radial format.\n')
    for site_dir, area_dir in site_dirs:
        for date_time in date_time_list:
            date_str = date_time.strftime('%Y%m%d')
            time_str = date_time.strftime('%H%M')
            print(f'Processing site: {site_dir}, Area: {area_dir}, Date: {date_str}, Time: {time_str}JST')
            
            # Construct the input file name
            base_dir_P008 = f'{xrain_path}/kyoku1/{area_dir}/{site_dir}/{date_str[0:4]}/{date_str[4:6]}/{date_str[6:8]}/{time_str[0:2]}'
            search_pattern = os.path.join(base_dir_P008, f'{site_dir}-{date_str}-{time_str}-P008-EL*0000.tgz')
            files_P008 = glob.glob(search_pattern)
            if not files_P008:
                print(f'Warning: No P008 file found for {site_dir} on {date_str} at {time_str} JST. Skipping this datetime.\n')
                err_flag = True
                continue

            for fname_P008 in files_P008:
                print(f'Input file name (Raw, P008): {fname_P008}')
                if '-R005-' in fname_P008:
                    skip_message(1, fname_P008)
                    err_flag = True
                    continue
                
                # Construct the corresponding R005 file name
                fname_R005 = fname_P008.replace('P008', 'R005')
                fname_R005 = fname_R005.replace('kyoku1', 'kyoku2')
                if not os.path.exists(fname_R005):
                    skip_message(3, fname_R005)
                    err_flag = True
                    continue
            
                # Convert files
                conv_P = Converter(mode=0, flag_overwrite=False)
                conv_P.convert(fname_P008, outdir)
                conv_Z = Converter(mode=1, flag_overwrite=True)
                conv_Z.convert(fname_R005, outdir)
                print(f'Convert success. Saved to: {conv_Z.ncname}\n')

                del conv_P, conv_Z # Clear memory to avoid memory issues

        print(f'Finished processing site: {site_dir}\n')
    
    print('All sites processed.')
    if err_flag:
        print('Some files and datetimes skipped. Please check the messages above.')

if __name__=='__main__':
    main()



