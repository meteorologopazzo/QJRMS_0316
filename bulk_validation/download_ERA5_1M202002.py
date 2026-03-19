import cdsapi

dataset = "reanalysis-era5-pressure-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "variable": [
        "divergence",
        "specific_humidity",
        "temperature",
        "u_component_of_wind",
        "v_component_of_wind",
        "vertical_velocity"
    ],
    "pressure_level": [
        "200", "225", "250",
        "300", "350", "400",
        "450", "500", "550",
        "600", "650", "700",
        "750", "775", "800",
        "825", "850", "875",
        "900", "925", "950",
        "1000"
    ],
    "year": ["2020"],
    "month": ["02"],
    "time": ["00:00"],
    "data_format": "netcdf",
    "download_format": "unarchived"
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()