import cdsapi

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],             
    "year": [
        "1990", "1991", "1992",
        "1993", "1994", "1995",
        "1996", "1997", "1998",
        "1999", "2000", "2001",
        "2002", "2003", "2004",
        "2005", "2006", "2007",
        "2008", "2009", "2010",
        "2011", "2012", "2013",
        "2014", "2015", "2016",
        "2017", "2018", "2019",
        "2020"
    ],
    "month": ["02"],
    "time": ["00:00"],
    "data_format": "netcdf"
}


req1 = {**request, "variable": ["sea_surface_temperature", 
                 "mean_sea_level_pressure",
                 "boundary_layer_height", 
                 "10m_u_component_of_wind",
                 "10m_v_component_of_wind",
                 "2m_dewpoint_temperature",
                 "2m_temperature",]
                 }
req1["download_format"] = "unarchived"

req2 = {**request, 
        "variable": [
        "surface_latent_heat_flux",
        "surface_sensible_heat_flux"
        ]}
req2["download_format"] = "zip"

client = cdsapi.Client()
client.retrieve(dataset, req1).download("/work/mh1498/m301248/TCO_data/U_era5_climatology/E5_1MSfc_02.nc")

client.retrieve(dataset, req2).download("/work/mh1498/m301248/TCO_data/U_era5_climatology/E5_1MFluxes_02.zip")
import zipfile
with zipfile.ZipFile("/work/mh1498/m301248/TCO_data/U_era5_climatology/E5_1MFluxes_02.zip") as z:
    z.extractall("/work/mh1498/m301248/TCO_data/U_era5_climatology/")