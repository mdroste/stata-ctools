clear
quietly set obs 200000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen double v1 = (_n+1)/7
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen double v2 = (_n+2)/7
quietly replace v2 = . if mod(_n,101)==0
quietly replace v2 = .z if mod(_n,103)==0
quietly gen double v3 = (_n+3)/7
quietly replace v3 = . if mod(_n,101)==0
quietly replace v3 = .z if mod(_n,103)==0
quietly gen double v4 = (_n+4)/7
quietly replace v4 = . if mod(_n,101)==0
quietly replace v4 = .z if mod(_n,103)==0
quietly gen double v5 = (_n+5)/7
quietly replace v5 = . if mod(_n,101)==0
quietly replace v5 = .z if mod(_n,103)==0
quietly gen double v6 = (_n+6)/7
quietly replace v6 = . if mod(_n,101)==0
quietly replace v6 = .z if mod(_n,103)==0
quietly gen double v7 = (_n+7)/7
quietly replace v7 = . if mod(_n,101)==0
quietly replace v7 = .z if mod(_n,103)==0
quietly gen double v8 = (_n+8)/7
quietly replace v8 = . if mod(_n,101)==0
quietly replace v8 = .z if mod(_n,103)==0
quietly gen double v9 = (_n+9)/7
quietly replace v9 = . if mod(_n,101)==0
quietly replace v9 = .z if mod(_n,103)==0
quietly gen double v10 = (_n+10)/7
quietly replace v10 = . if mod(_n,101)==0
quietly replace v10 = .z if mod(_n,103)==0
quietly gen double v11 = (_n+11)/7
quietly replace v11 = . if mod(_n,101)==0
quietly replace v11 = .z if mod(_n,103)==0
quietly gen double v12 = (_n+12)/7
quietly replace v12 = . if mod(_n,101)==0
quietly replace v12 = .z if mod(_n,103)==0
quietly gen double v13 = (_n+13)/7
quietly replace v13 = . if mod(_n,101)==0
quietly replace v13 = .z if mod(_n,103)==0
quietly gen double v14 = (_n+14)/7
quietly replace v14 = . if mod(_n,101)==0
quietly replace v14 = .z if mod(_n,103)==0
quietly gen double v15 = (_n+15)/7
quietly replace v15 = . if mod(_n,101)==0
quietly replace v15 = .z if mod(_n,103)==0
quietly gen double v16 = (_n+16)/7
quietly replace v16 = . if mod(_n,101)==0
quietly replace v16 = .z if mod(_n,103)==0
quietly gen double v17 = (_n+17)/7
quietly replace v17 = . if mod(_n,101)==0
quietly replace v17 = .z if mod(_n,103)==0
quietly gen double v18 = (_n+18)/7
quietly replace v18 = . if mod(_n,101)==0
quietly replace v18 = .z if mod(_n,103)==0
quietly gen double v19 = (_n+19)/7
quietly replace v19 = . if mod(_n,101)==0
quietly replace v19 = .z if mod(_n,103)==0
quietly gen double v20 = (_n+20)/7
quietly replace v20 = . if mod(_n,101)==0
quietly replace v20 = .z if mod(_n,103)==0
quietly gen double v21 = (_n+21)/7
quietly replace v21 = . if mod(_n,101)==0
quietly replace v21 = .z if mod(_n,103)==0
quietly gen double v22 = (_n+22)/7
quietly replace v22 = . if mod(_n,101)==0
quietly replace v22 = .z if mod(_n,103)==0
quietly gen double v23 = (_n+23)/7
quietly replace v23 = . if mod(_n,101)==0
quietly replace v23 = .z if mod(_n,103)==0
quietly gen double v24 = (_n+24)/7
quietly replace v24 = . if mod(_n,101)==0
quietly replace v24 = .z if mod(_n,103)==0
quietly gen double v25 = (_n+25)/7
quietly replace v25 = . if mod(_n,101)==0
quietly replace v25 = .z if mod(_n,103)==0
quietly gen double v26 = (_n+26)/7
quietly replace v26 = . if mod(_n,101)==0
quietly replace v26 = .z if mod(_n,103)==0
quietly gen double v27 = (_n+27)/7
quietly replace v27 = . if mod(_n,101)==0
quietly replace v27 = .z if mod(_n,103)==0
quietly gen double v28 = (_n+28)/7
quietly replace v28 = . if mod(_n,101)==0
quietly replace v28 = .z if mod(_n,103)==0
quietly gen double v29 = (_n+29)/7
quietly replace v29 = . if mod(_n,101)==0
quietly replace v29 = .z if mod(_n,103)==0
quietly gen double v30 = (_n+30)/7
quietly replace v30 = . if mod(_n,101)==0
quietly replace v30 = .z if mod(_n,103)==0
quietly gen double v31 = (_n+31)/7
quietly replace v31 = . if mod(_n,101)==0
quietly replace v31 = .z if mod(_n,103)==0
quietly gen double v32 = (_n+32)/7
quietly replace v32 = . if mod(_n,101)==0
quietly replace v32 = .z if mod(_n,103)==0
quietly gen double v33 = (_n+33)/7
quietly replace v33 = . if mod(_n,101)==0
quietly replace v33 = .z if mod(_n,103)==0
quietly gen double v34 = (_n+34)/7
quietly replace v34 = . if mod(_n,101)==0
quietly replace v34 = .z if mod(_n,103)==0
quietly gen double v35 = (_n+35)/7
quietly replace v35 = . if mod(_n,101)==0
quietly replace v35 = .z if mod(_n,103)==0
quietly gen double v36 = (_n+36)/7
quietly replace v36 = . if mod(_n,101)==0
quietly replace v36 = .z if mod(_n,103)==0
quietly gen double v37 = (_n+37)/7
quietly replace v37 = . if mod(_n,101)==0
quietly replace v37 = .z if mod(_n,103)==0
quietly gen double v38 = (_n+38)/7
quietly replace v38 = . if mod(_n,101)==0
quietly replace v38 = .z if mod(_n,103)==0
quietly gen double v39 = (_n+39)/7
quietly replace v39 = . if mod(_n,101)==0
quietly replace v39 = .z if mod(_n,103)==0
quietly gen double v40 = (_n+40)/7
quietly replace v40 = . if mod(_n,101)==0
quietly replace v40 = .z if mod(_n,103)==0
quietly gen double v41 = (_n+41)/7
quietly replace v41 = . if mod(_n,101)==0
quietly replace v41 = .z if mod(_n,103)==0
quietly gen double v42 = (_n+42)/7
quietly replace v42 = . if mod(_n,101)==0
quietly replace v42 = .z if mod(_n,103)==0
quietly gen double v43 = (_n+43)/7
quietly replace v43 = . if mod(_n,101)==0
quietly replace v43 = .z if mod(_n,103)==0
quietly gen double v44 = (_n+44)/7
quietly replace v44 = . if mod(_n,101)==0
quietly replace v44 = .z if mod(_n,103)==0
quietly gen double v45 = (_n+45)/7
quietly replace v45 = . if mod(_n,101)==0
quietly replace v45 = .z if mod(_n,103)==0
quietly gen double v46 = (_n+46)/7
quietly replace v46 = . if mod(_n,101)==0
quietly replace v46 = .z if mod(_n,103)==0
quietly gen double v47 = (_n+47)/7
quietly replace v47 = . if mod(_n,101)==0
quietly replace v47 = .z if mod(_n,103)==0
quietly gen double v48 = (_n+48)/7
quietly replace v48 = . if mod(_n,101)==0
quietly replace v48 = .z if mod(_n,103)==0
quietly gen double v49 = (_n+49)/7
quietly replace v49 = . if mod(_n,101)==0
quietly replace v49 = .z if mod(_n,103)==0
quietly gen double v50 = (_n+50)/7
quietly replace v50 = . if mod(_n,101)==0
quietly replace v50 = .z if mod(_n,103)==0
quietly gen double v51 = (_n+51)/7
quietly replace v51 = . if mod(_n,101)==0
quietly replace v51 = .z if mod(_n,103)==0
quietly gen double v52 = (_n+52)/7
quietly replace v52 = . if mod(_n,101)==0
quietly replace v52 = .z if mod(_n,103)==0
quietly gen double v53 = (_n+53)/7
quietly replace v53 = . if mod(_n,101)==0
quietly replace v53 = .z if mod(_n,103)==0
quietly gen double v54 = (_n+54)/7
quietly replace v54 = . if mod(_n,101)==0
quietly replace v54 = .z if mod(_n,103)==0
quietly gen double v55 = (_n+55)/7
quietly replace v55 = . if mod(_n,101)==0
quietly replace v55 = .z if mod(_n,103)==0
quietly gen double v56 = (_n+56)/7
quietly replace v56 = . if mod(_n,101)==0
quietly replace v56 = .z if mod(_n,103)==0
quietly gen double v57 = (_n+57)/7
quietly replace v57 = . if mod(_n,101)==0
quietly replace v57 = .z if mod(_n,103)==0
quietly gen double v58 = (_n+58)/7
quietly replace v58 = . if mod(_n,101)==0
quietly replace v58 = .z if mod(_n,103)==0
quietly gen double v59 = (_n+59)/7
quietly replace v59 = . if mod(_n,101)==0
quietly replace v59 = .z if mod(_n,103)==0
quietly gen double v60 = (_n+60)/7
quietly replace v60 = . if mod(_n,101)==0
quietly replace v60 = .z if mod(_n,103)==0
quietly gen double v61 = (_n+61)/7
quietly replace v61 = . if mod(_n,101)==0
quietly replace v61 = .z if mod(_n,103)==0
quietly gen double v62 = (_n+62)/7
quietly replace v62 = . if mod(_n,101)==0
quietly replace v62 = .z if mod(_n,103)==0
quietly gen double v63 = (_n+63)/7
quietly replace v63 = . if mod(_n,101)==0
quietly replace v63 = .z if mod(_n,103)==0
quietly gen double v64 = (_n+64)/7
quietly replace v64 = . if mod(_n,101)==0
quietly replace v64 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "baseline_n200k_k64_double_r6" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
