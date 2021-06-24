import pandas as pd 
import glob
import xarray as xr
import time
import datetime as dt
import sys
import getpass

def itp_ascii_to_netcdf(in_path, out_file,min_length=4):
    """
    Function that reads in all profiles from a WHOI-ITP from the directory in_path and dumps the data in a single
    NETcdf-file out_file.
    Renames the varaibles, and metadata to follow CF-1.8 standard
    and adds discovermetadata to follow the ACDD-1.3 standard
    min_length is set to 4 to not use very short profiles. 
    """
    
    files = sorted(glob.glob(in_path + "/*.dat"))[:40] #dont increase it will be craaaazy
    first = True
    start = time.time()
    for i in files:


        meta = pd.read_table(i,skiprows=None,sep="\s+",nrows=1,engine="python")

        if(meta.values[0,4]<min_length): # hopper over de korteste profilene
            if i == files[-1] and "buoy" not in globals(): #sier ifra hvis netcdfen blir tom
                sys.exit("No profiles of desired lenght in directory")
            continue


        #print(df.values.shape[0], int(str(meta.head().columns[3])[:-1]))

        df = pd.read_table(i,skiprows=2, delim_whitespace=True,skipfooter=1,engine="python")

        measurement_time = pd.to_datetime(float(meta.values[0,1]),origin=str(int(meta.values[0,0])),unit="D").timestamp()
        measurement_lat  = float(meta.values[0,2])
        measurement_lon  = float(meta.values[0,3])


        if "%year" in df.columns:
            df["%year"] = df["%year"].astype(int)
            df["times"] = pd.to_datetime(df["day"], unit = 'D', 
                                         origin = str(df["%year"][0]))
            df = df.drop(["%year","day"],axis=1)
        if "nobs" in df.columns:
            df = df.drop("nobs",axis=1)
        if "nacm" in df.columns:
            df = df.drop("nacm",axis=1)


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
        df = df.set_index("sea_water_pressure")

        ds = xr.Dataset.from_dataframe(df)

        ds["time"] = measurement_time
        ds["latitude"]  = measurement_lat
        ds["longitude"]  = measurement_lon

        #profile_nr = int(list(pd.read_table(i,sep="[:, ]",nrows=0,engine="python"))[4])
        profile_nr = int(str(meta.head().columns[3])[:-1])


        ds = ds.assign_coords(profile=profile_nr)
        ds = ds.expand_dims("profile")


        # joining files
        if first==True:
            buoy= ds
            first=False
        else:
            buoy=xr.concat([buoy,ds],dim = "profile")
            #concat gjør at det blir mye nan i temp og salinity, men det følger cf. hør med Steingod

    


    #La oss fikse litt metadata da
    #legg til alle de her greiene: https://adc.met.no/node/4
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
            continue
        buoy[i].attrs["standard_name"] = i
        buoy[i].attrs["units"] = units[i]

    #global attributes
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
    buoy.attrs["summary"] = ("""Trajectory of ITP (Ice-Tethered Profiler) profiles, that use pressure in dbar as vertical coordinate
                             All profiles contain measurement times, temperature and salinity, and may include dissolved oxygen,
                             chromophoric dissolved organic matter (CDOM), turbidity, mass concentration of chlorophyll,
                             photosynthetically active radiation (PAR) and velocities. Metadata include time of initialization,
                             coordinates and profile data points (ndepths).""")
    #summary for final files, averaged
    buoy.attrs["summary"] = ("""Trajectory of ITP (Ice-Tethered Profiler) profiles, that use pressure in dbar as vertical coordinate.
                             All profiles contain averaged measurements of temperature and salinity, and may include dissolved oxygen,
                             chromophoric dissolved organic matter (CDOM), turbidity, mass concentration of chlorophyll,
                             photosynthetically active radiation (PAR) and velocities. Metadata include time of initialization,
                             coordinates and profile data points (ndepths).""")
    buoy.attrs["keywords"] = ["EARTH SCIENCE > OCEANS > SALINITY/DENSITY > DENSITY",
                          "EARTH SCIENCE > OCEANS > OCEAN TEMPERATURE > WATER TEMPERATURE",
                          "EARTH SCIENCE > OCEANS > SALINITY/DENSITY > SALINITY",
                          "EARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > OXYGEN",
                          "EARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > ORGANIC MATTER",
                          "EARTH SCIENCE > OCEANS > OCEAN OPTICS > TURBIDITY",
                          "EARTH SCIENCE > OCEANS > OCEAN CHEMISTRY > CHLOROPHYLL",
                          "EARTH SCIENCE > OCEANS > OCEAN CIRCULATION > ADVECTION"]
    buoy.attrs["keywords_vocabulary"] = "GCMD"
    buoy.attrs["featureType"] = "trajectoryProfile"

    buoy.attrs["geospatial_lat_min"] = min(buoy.latitude.values)
    buoy.attrs["geospatial_lat_max"] = max(buoy.latitude.values)
    buoy.attrs["geospatial_lon_min"] = min(buoy.longitude.values)
    buoy.attrs["geospatial_lon_max"] = max(buoy.longitude.values)

    buoy.attrs["time_coverage_start"] = min(buoy.time.values)
    buoy.attrs["time_coverage_end"] = max(buoy.time.values)

    buoy.attrs["Conventions"] = "ACDD-1.3, CF-1.8"
    buoy.attrs["history"] = str([str(dt.datetime.now()),getpass.getuser(), "program name:",sys.argv])
    buoy.attrs["date_created"] = str(dt.date.today())
    buoy.attrs["creator_type"] = "Institution"
    buoy.attrs["creator_institution"] = "Woods Hole Oceanographic Institute (WHOI)"
    buoy.attrs["creator_name"] = "Woods Hole Oceanographic Institution"
    buoy.attrs["creator_email"] = "information@whoi.edu" #?
    buoy.attrs["creator_url"] = "https://www2.whoi.edu/site/itp/"
    buoy.attrs["project"] = project_names[itp_nr]
    buoy.attrs["license"] = "Free"
    buoy.attrs["metadata_author"] = "Magnus Dyrmose Ryseth and Gunnar Thorsen Liahjell for MET Norway"


    
    buoy.to_netcdf(out_file)
    print("Det tok", time.time()-start)
    return(0)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inPath", help="Path to directory containing ITP .dat files")
    parser.add_argument("-o","--outFile", help="Filename fot the outputfile, with path")
    args = parser.parse_args()
    
    in_path  = args.inPath
    out_file = args.outFile

    
    #path = "data/114"
    itp_ascii_to_netcdf(in_path, out_file)

    
    