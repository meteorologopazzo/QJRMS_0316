import cdsapi

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],             
    "year": ["2020"],
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


client = cdsapi.Client()
client.retrieve(dataset, req1).download("/work/mh1498/m301248/TCO_data/E5_1M_sst_202002.nc")