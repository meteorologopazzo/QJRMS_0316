!#/bin/bash

# /work/bk1099/data/E5/sf/an/1M/034/E5sf00_1M_2020_034.grb   # sst
# /work/bk1099/data/E5/sf/an/1M/134/E5sf00_1M_2020_134.grb   # sfc pressure

# /work/bk1099/data/E5/pl/an/1M/130/E5pl00_1M_2020_130.grb   # temperature
# /work/bk1099/data/E5/pl/an/1M/133/E5pl00_1M_2020_133.grb   # specific humidity
# /work/bk1099/data/E5/pl/an/1M/155/E5pl00_1M_2020_155.grb   # divergence
# /work/bk1099/data/E5/pl/an/1M/130/E5pl00_1M_2020_131.grb   # U
# /work/bk1099/data/E5/pl/an/1M/130/E5pl00_1M_2020_132.grb   # V

module load cdo

paths=(
  "/work/bk1099/data/E5/sf/an/1M/034/E5sf00_1M_2020_034.grb"
  "/work/bk1099/data/E5/sf/an/1M/134/E5sf00_1M_2020_134.grb"
  "/work/bk1099/data/E5/pl/an/1M/130/E5pl00_1M_2020_130.grb"
  "/work/bk1099/data/E5/pl/an/1M/133/E5pl00_1M_2020_133.grb"
  "/work/bk1099/data/E5/pl/an/1M/155/E5pl00_1M_2020_155.grb"
  "/work/bk1099/data/E5/pl/an/1M/131/E5pl00_1M_2020_131.grb"
  "/work/bk1099/data/E5/pl/an/1M/132/E5pl00_1M_2020_132.grb"
)


names_in=(
  "E5sf00_1M_2020_034.grb"
  "E5sf00_1M_2020_134.grb"
  "E5pl00_1M_2020_130.grb"
  "E5pl00_1M_2020_133.grb"
  "E5pl00_1M_2020_155.grb"
  "E5pl00_1M_2020_131.grb"
  "E5pl00_1M_2020_132.grb"
)

names_out=(
  "E5sf00_1M_2020_sst.nc"
  "E5sf00_1M_2020_sp.nc"
  "E5pl00_1M_2020_ta.nc"
  "E5pl00_1M_2020_q.nc"
  "E5pl00_1M_2020_D.nc"
  "E5pl00_1M_2020_u.nc"
  "E5pl00_1M_2020_v.nc"
)

for i in "${!paths[@]}"; do
    # Only copy if the source file isn't already in the working directory
    if [[ ! -f "${names_in[$i]}" ]]; then
        cp "${paths[$i]}" .
    fi

    # Only convert if the output .nc file doesn't already exist
    if [[ ! -f "${names_out[$i]}" ]]; then
        cdo -f nc copy "${names_in[$i]}" "${names_out[$i]}"
    fi
done

if [[ ! -f "E5_2020_02_1M_selVars_merged.nc" ]]; then
    cdo merge "${names_out[@]}" E5_2020_02_1M_selVars_merged.nc
fi