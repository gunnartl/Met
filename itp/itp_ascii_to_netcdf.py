import pandas as pd 
import glob
import xarray as xr
import time
import datetime as dt
import sys
import getpass

def itp_ascii_to_netcdf(in_path=None, out_file,existing_netcdf=None,min_length=4):
    """
    Function that reads in all profiles from a WHOI-ITP from the directory in_path and dumps the data in a single
    NETcdf-file: out_file.
    For new data from a ITP that is still active, existing_netcdf can be specified, and the new data will be added to this file.

    Renames the varaibles, and metadata to follow CF-1.8 standard
    and adds discovermetadata to follow the ACDD-1.3 standard

    min_length is set to 4 to skip very short profiles. can be changed 
    """
    if in_path == None:
        sys.exit("In-path to ipt-.dat files must be set")
    if out_file == None:
        sys.exit("Out-file name must be set")

    files = sorted(glob.glob(in_path + "/*.dat"))[:60] #just 40 files to be able to test-run on crappy laptop
    if existing_netcdf == None:
        first = True
        
    else:
        first = False
        buoy  = xr.open_dataset(existing_netcdf,engine="netcdf4")
        buoy.close() #lukker bare NETcdf-fila så man kan skrive til den etterpå
        changes = False

    start = time.time()
    for i in files:
        if (existing_netcdf!= None) and (int(i[-8:-4]) in buoy.profile.values): # checks that the given profile has not already been rea
            continue

        meta = pd.read_table(i,skiprows=None,sep="\s+",nrows=1,engine="python")

        if(meta.values[0,4]<min_length): # hopper over de korteste profilene
            if i == files[-1]:
                if "buoy" not in locals(): #sier ifra hvis netcdfen blir tom
                    sys.exit("No profiles of desired lenght in target directory")
                if changes == False: #sier ifra hvis det ikke blir noen endringer i eksisterende fil
                    sys.exit("No new profiles, or no new profiles of desired lenght in target directory. No changes made to {file}".format(file=existing_netcdf), )
            continue


        df = pd.read_table(i,skiprows=2, delim_whitespace=True,skipfooter=1,engine="python")

        measurement_time = pd.to_datetime(float(meta.values[0,1]),origin=str(int(meta.values[0,0])),unit="D").timestamp()
        measurement_lat  = float(meta.values[0,2])
        measurement_lon  = float(meta.values[0,3])

        #removing useless columns nobs and nacm, and combines year and day to "times"
        if "%year" in df.columns:
            df["times"] = 0.0     #makes a new column to keep trak of individual measuremnt times, if included
            for i in range(len(df["%year"])):
                df.times.values[i] = pd.to_datetime(float(df.day[i]),origin=str(int(df["%year"][i])),unit="D").timestamp()
            df = df.drop(["%year","day"],axis=1)
        if "nobs" in df.columns:
            df = df.drop("nobs",axis=1)
        if "nacm" in df.columns:
            df = df.drop("nacm",axis=1)

        #standard names: 
        df.rename(columns={"%pressure(dbar)":"sea_water_pressure",
                           "pressure(dbar)":"sea_water_pressure",
                           "temperature(C)":"sea_water_temperature",
                           "salinity":"sea_water_salinity",
                           "dissolved_oxygen":"moles_of_oxygen_per_unit_mass_in_sea_water",
                           "oxygen(umol/kg)":"moles_of_oxygen_per_unit_mass_in_sea_water",
                           "CDOM(ppb)":"concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent_mass_fraction_of_quinine_sulfate_dihydrate",
                           "turbidity(/m/sr)x10^4":"sea_water_turbidity",
                           "chlorophyll-a(ug/l)":"mass_concentration_of_chlorophyll_a_in_sea_water",
                           "PAR(uE/m^2/s)":"downwelling_photosynthetic_radiative_flux_in_sea_water",
                           "east(cm/s)":"eastward_sea_water_velocity",
                           "north(cm/s)":"northward_sea_water_velocity",
                           "vert(cm/s)":"upward_sea_water_velocity"
                           }, inplace=True)
        
        #setter trykket som koordinat
        df = df.set_index("sea_water_pressure")

        ds = xr.Dataset.from_dataframe(df)

        ds["time"] = measurement_time
        ds["latitude"]  = measurement_lat
        ds["longitude"]  = measurement_lon

        profile_nr = int(str(meta.head().columns[3])[:-1])

        #setter profil som koordinat
        ds = ds.assign_coords(profile=profile_nr) 
        ds = ds.expand_dims("profile")


        # joining files
        if first==True:
            buoy= ds
            first=False
        else:
            buoy=xr.concat([buoy,ds],dim = "profile")
            changes = True

    #lager metadata:
    units = {"time":"Seconds since 1970-01-01 00:00:00+0",
             "latitude":"degree_north",
             "longitude":"degree_east",
             "sea_water_pressure":"dBar",
             "sea_water_salinity":"1e-3",
             "sea_water_temperature":"celsius",
             "sea_water_turbidity":"(m-1 sr-1) x 10e-4",
             "moles_of_oxygen_per_unit_mass_in_sea_water":"umol/kg",
             "mass_concentration_of_chlorophyll_a_in_sea_water":"ug/l",
             "downwelling_photosynthetic_radiative_flux_in_sea_water":"uE/m^2/s",
             "concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent_mass_fraction_of_quinine_sulfate_dihydrate":"ppb",
             "eastward_sea_water_velocity":"cm/s",
             "northward_sea_water_velocity":"cm/s",
             "upward_sea_water_velocity":"cm/s"
            }
    
    #søksmetadata
    for i in buoy:
        if i == "times": #behandler times for seg selv da dette er tidspunkt for individuelle målinger 
            buoy[i].attrs["long_name"] = "individual time for each measurement in a profile"
            buoy[i].attrs["unit"] = units["time"]
            continue
        buoy[i].attrs["standard_name"] = i
        buoy[i].attrs["unit"] = units[i]
    buoy["time"].attrs["long_name"] = "starting time for each profile"

    #global attributes
    #list of affiliated projects to itp numbers
    project_names= {"1" : "Beaufort Gyre Observing System (BGOS)",
                    "2" : "Beaufort Gyre Freshwater Experiment (BGFE)",
                    "3" : "Beaufort Gyre Observing System (BGOS)",
                    "4" : "Beaufort Gyre Observing System (BGOS)",
                    "5" : "Beaufort Gyre Observing System (BGOS)",
                    "6" : "Beaufort Gyre Observing System (BGOS)",
                    "7" : "North Pole Environmental Observatory (NPEO)",
                    "8" : "Beaufort Gyre Observing System (BGOS)",
                    "9" : "Damocles",
                    "10" : "Damocles",
                    "11" : "Damocles",
                    "12" : "Damocles",
                    "13" : "Beaufort Gyre Observing System (BGOS)",
                    "14" : "Damocles",
                    "15" : "Damocles",
                    "16" : "Damocles",
                    "17" : "Damocles",
                    "18" : "Beaufort Gyre Observing System (BGOS)",
                    "19" : "North Pole Environmental Observatory (NPEO)",
                    "20" : "Beaufort Gyre Observing System (BGOS)",
                    "21" : "Beaufort Gyre Observing System (BGOS)",
                    "22" : "Beaufort Gyre Observing System (BGOS)",
                    "23" : "Beaufort Gyre Observing System (BGOS)",
                    "24" : "Damocles",
                    "25" : "Damocles",
                    "26" : "Damocles",
                    "27" : "Damocles",
                    "28" : "Damocles",
                    "29" : "Damocles",
                    "30" : "Beaufort Gyre Observing System (BGOS)",
                    "31" : "Not available",
                    "32" : "Beaufort Gyre Observing System (BGOS)",
                    "33" : "Beaufort Gyre Observing System (BGOS)",
                    "34" : "Beaufort Gyre Observing System (BGOS)",
                    "36" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "37" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "38" : "North Pole Environmental Observatory (NPEO)",
                    "40" : "National Institute of Water and Atmospheric Research (NIWA)",
                    "41" : "Beaufort Gyre Observing System (BGOS)",
                    "42" : "Beaufort Gyre Observing System (BGOS)",
                    "43" : "Beaufort Gyre Observing System (BGOS)",
                    "44" : "Beaufort Gyre Observing System (BGOS)",
                    "47" : "North Pole Environmental Observatory (NPEO)",
                    "48" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "49" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "50" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "51" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "52" : "Beaufort Gyre Observing System (BGOS)",
                    "53" : "Beaufort Gyre Observing System (BGOS)",
                    "54" : "Beaufort Gyre Observing System (BGOS)",
                    "55" : "Beaufort Gyre Observing System (BGOS)",
                    "56" : "North Pole Environmental Observatory (NPEO)",
                    "57" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "58" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "60" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "62" : "Beaufort Gyre Observing System (BGOS)",
                    "63" : "NP-39 drifting ice station (no project specified)",
                    "64" : "Beaufort Gyre Observing System (BGOS)",
                    "65" : "Beaufort Gyre Observing System (BGOS)",
                    "66" : "Beaufort Gyre Observing System (BGOS)",
                    "68" : "Beaufort Gyre Observing System (BGOS)",
                    "69" : "Beaufort Gyre Observing System (BGOS)",
                    "70" : "Beaufort Gyre Observing System (BGOS)",
                    "71" : "Beaufort Gyre Observing System (BGOS)",
                    "72" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "73" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "74" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "75" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "76" : "North Pole Environmental Observatory (NPEO)",
                    "77" : "ONR Marginal Ice Zone (MIZ)",
                    "78" : "ONR Marginal Ice Zone (MIZ)",
                    "79" : "ONR Marginal Ice Zone (MIZ)",
                    "80" : "ONR Marginal Ice Zone (MIZ)",
                    "81" : "CHINARE 2014 Expedition",
                    "82" : "CHINARE 2014 Expedition",
                    "83" : "North Pole Environmental Observatory (NPEO)",
                    "84" : "Beaufort Gyre Observing System (BGOS)",
                    "85" : "Beaufort Gyre Observing System (BGOS)",
                    "86" : "", # står bare: from the Korean Research vessel Araon .
                    "87" : "CHINARE 2014 Expedition",
                    "88" : "Beaufort Gyre Observing System (BGOS)",
                    "89" : "Beaufort Gyre Observing System (BGOS)",
                    "90" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "91" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "92" : "Nansen and Amundsen Basins Observational System (NABOS)",
                    "93" : "Frontiers in Arctic Marine Monitoring (FRAM)",
                    "94" : "Hybrid Arctic/Antarctic Float Observation System (HAFOS)",
                    "95" : "",#from the Russian ice camp Barneo
                    "96" : "", #fins ikke? 
                    "97" : "Beaufort Gyre Observing System (BGOS)",
                    "98" : "Beaufort Gyre Observing System (BGOS)",
                    "99" : "Beaufort Gyre Observing System (BGOS)",
                    "100" : "Beaufort Gyre Observing System (BGOS)",
                    "101" : "Beaufort Gyre Observing System (BGOS)",
                    "102" : "Multidisciplinary drifting Observatory for the Study of Arctic Climate (MOSAiC)",
                    "103" : "Stratified Ocean Dynamics of the Arctic (SODA)",
                    "104" : "Stratified Ocean Dynamics of the Arctic (SODA)",
                    "105" : "Stratified Ocean Dynamics of the Arctic (SODA)",
                    "106" : "",#fins ikke?
                    "107" : "Beaufort Gyre Observing System (BGOS)",
                    "108" : "Beaufort Gyre Observing System (BGOS)",
                    "109" : "Beaufort Gyre Observing System (BGOS)",
                    "110" : "Beaufort Gyre Observing System (BGOS)",
                    "111" : "Multidisciplinary drifting Observatory for the Study of Arctic Climate (MOSAiC)",
                    "112" : "Beaufort Gyre Observing System (BGOS)",
                    "113" : "Stratified Ocean Dynamics of the Arctic (SODA)",
                    "114" : "Stratified Ocean Dynamics of the Arctic (SODA)",
                    "115" : "",#fins ikke?
                    "116" : "Coordinated Arctic Acoustic Thermometry Experiment (CAATEX)",
                    "117" : "Beaufort Gyre Observing System (BGOS)",
                    "118" : "Beaufort Gyre Observing System (BGOS)",
                    "119" : "Beaufort Gyre Observing System (BGOS)",
                    "120" : "Beaufort Gyre Observing System (BGOS)",
                    "121" : "Beaufort Gyre Observing System (BGOS)" #siste per 22/6-2021
                  }
    itp_nr = str(meta.head().columns[1][:-1])

    buoy.attrs["title"] = ("Trajectory of profiles from WHOI-ITP " + itp_nr) #change a0 to the meta-indexing
    #summary for normal grd-files, Level 2


    if "times" in df.columns:
        buoy.attrs["summary"] = "Trajectory of ITP (Ice-Tethered Profiler) profiles, that use pressure in dbar as vertical coordinate "+ \
                                "All profiles contain measurement times, temperature and salinity, and may include dissolved oxygen, "+ \
                                "chromophoric dissolved organic matter (CDOM), turbidity, mass concentration of chlorophyll, " +\
                                "photosynthetically active radiation (PAR) and velocities. Metadata include time of initialization, "+ \
                                "coordinates and profile data points (ndepths)."
    else: #summary for final files, averaged
        buoy.attrs["summary"] = "Trajectory of ITP (Ice-Tethered Profiler) profiles, that use pressure in dbar as vertical coordinate. "+\
                                "All profiles contain averaged measurements of temperature and salinity, and may include dissolved oxygen ,"+\
                                "chromophoric dissolved organic matter (CDOM), turbidity, mass concentration of chlorophyll, "+\
                                "photosynthetically active radiation (PAR) and velocities. Metadata include time of initialization, "+\
                                "coordinates and profile data points (ndepths)."
    buoy.attrs["keywords"] = "EARTH SCIENCE > OCEANS > SALINITY/DENSITY > DENSITY,\n"+\
                             "EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > WATER TEMPERATURE,\n"+\
                             "EARTH SCIENCE > OCEANS > SALINITY/DENSITY > SALINITY,"

    if "moles_of_oxygen_per_unit_mass_in_sea_water" in df.columns:
        buoy.attrs["keywords"] += "\nEARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > OXYGEN,"
    if "mass_concentration_of_chlorophyll_a_in_sea_water" in df.columns:
        buoy.attrs["keywords"] += "\nEARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > CHLOROPHYLL,"
    if "concentration_of_colored_dissolved_organic_matter_in_sea_water_expressed_as_equivalent_mass_fraction_of_quinine_sulfate_dihydrate" in df.columns:
        buoy.attrs["keywords"] += "\nEARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > CHLOROPHYLL,"
    if "sea_water_turbidity" in df.columns:
        buoy.attrs["keywords"] += "\nEARTH SCIENCE > OCEANS > OCEAN OPTICS > TURBIDITY,"
    if "eastward_sea_water_velocity" in df.columns:
        buoy.attrs["keywords"] += "\nEARTH SCIENCE > OCEANS > OCEAN CIRCULATION > ADVECTION,"

    buoy.attrs["keywords"] = buoy.attrs["keywords"][:-1] #bare fjerner siste komma
    buoy.attrs["keywords_vocabulary"] = "GCMD"
    buoy.attrs["featureType"] = "trajectoryProfile"

    buoy.attrs["geospatial_lat_min"] = min(buoy.latitude.values)
    buoy.attrs["geospatial_lat_max"] = max(buoy.latitude.values)
    buoy.attrs["geospatial_lon_min"] = min(buoy.longitude.values)
    buoy.attrs["geospatial_lon_max"] = max(buoy.longitude.values)

    buoy.attrs["time_coverage_start"] = str(min(buoy.time.values))
    buoy.attrs["time_coverage_end"] = str(max(buoy.time.values))

    buoy.attrs["Conventions"] = "ACDD-1.3, CF-1.8"
    if existing_netcdf==None:
        buoy.attrs["history"] = "{time}: user: {user}, program:{program}".format(time=dt.datetime.now(), user=getpass.getuser(), program=sys.argv)#str([str(dt.datetime.now()),getpass.getuser(), "program name:",sys.argv])
    else:
        buoy.attrs["history"] += "\n{time}: user: {user}, program:{program}".format(time=dt.datetime.now(), user=getpass.getuser(), program=sys.argv)
    buoy.attrs["date_created"] = str(dt.date.today())
    buoy.attrs["creator_type"] = "Institution"
    buoy.attrs["creator_institution"] = "Woods Hole Oceanographic Institute (WHOI)"
    buoy.attrs["creator_name"] = "Woods Hole Oceanographic Institution"
    buoy.attrs["creator_email"] = "information@whoi.edu" #?
    buoy.attrs["creator_url"] = "https://www2.whoi.edu/site/itp/"
    buoy.attrs["project"] = project_names[itp_nr]
    buoy.attrs["license"] = "Free"
    buoy.attrs["metadata_author"] = "Magnus Dyrmose Ryseth and Gunnar Thorsen Liahjell for MET Norway"

    buoy = buoy.sortby("profile")
    buoy.to_netcdf(out_file)
    print("Det tok", time.time()-start)
    return(0)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inPath", help="Path to directory containing ITP .dat files")
    parser.add_argument("-o","--outFile", help="Filename fot the outputfile, with path")
    parser.add_argument("-e","--existingFile",help="Filename of existing NETcdf,to be updated")
    args = parser.parse_args()
    
    
    #path = "data/114"
    itp_ascii_to_netcdf(args.inPath,args.outFile,args.existingFile)



    """
    eksamle run for a new itp: 
    gunnartl@Abel:~/Documents/met/itp (master)$ python3 itp_ascii_to_netcdf.py -i data/114 -o itp144.nc

    eksample run for an already existing itp file that is to be updated with new data 
    (data is here stored to the same file as the original, this is optional): 
    gunnartl@Abel:~/Documents/met/itp (master)$ python3 itp_ascii_to_netcdf.py -i newdata/114 -o itp144.nc -e itp144.nc
    """