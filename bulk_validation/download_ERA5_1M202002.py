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
    "data_format": "netcdf",
    "download_format": "unarchived"
}

client = cdsapi.Client()
client.retrieve(dataset, request).download("/work/mh1498/m301248/TCO_data/U_era5_climatology/E5_1MAir_02.nc")