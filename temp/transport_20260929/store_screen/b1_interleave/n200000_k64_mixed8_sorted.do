clear
quietly set obs 200000
local transport_pad "xxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen str8 v2 = cond(mod(_n,17)==0,"",substr(string(mod(_n+2,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v3 = _n+3
quietly replace v3 = . if mod(_n,101)==0
quietly replace v3 = .z if mod(_n,103)==0
quietly gen str8 v4 = cond(mod(_n,17)==0,"",substr(string(mod(_n+4,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen double v5 = (_n+5)/7
quietly replace v5 = . if mod(_n,101)==0
quietly replace v5 = .z if mod(_n,103)==0
quietly gen str8 v6 = cond(mod(_n,17)==0,"",substr(string(mod(_n+6,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen int v7 = mod(_n+7,101)-50
quietly replace v7 = . if mod(_n,101)==0
quietly replace v7 = .z if mod(_n,103)==0
quietly gen str8 v8 = cond(mod(_n,17)==0,"",substr(string(mod(_n+8,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen float v9 = (_n+9)/7
quietly replace v9 = . if mod(_n,101)==0
quietly replace v9 = .z if mod(_n,103)==0
quietly gen str8 v10 = cond(mod(_n,17)==0,"",substr(string(mod(_n+10,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen byte v11 = mod(_n+11,101)-50
quietly replace v11 = . if mod(_n,101)==0
quietly replace v11 = .z if mod(_n,103)==0
quietly gen str8 v12 = cond(mod(_n,17)==0,"",substr(string(mod(_n+12,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v13 = _n+13
quietly replace v13 = . if mod(_n,101)==0
quietly replace v13 = .z if mod(_n,103)==0
quietly gen str8 v14 = cond(mod(_n,17)==0,"",substr(string(mod(_n+14,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen double v15 = (_n+15)/7
quietly replace v15 = . if mod(_n,101)==0
quietly replace v15 = .z if mod(_n,103)==0
quietly gen str8 v16 = cond(mod(_n,17)==0,"",substr(string(mod(_n+16,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen int v17 = mod(_n+17,101)-50
quietly replace v17 = . if mod(_n,101)==0
quietly replace v17 = .z if mod(_n,103)==0
quietly gen str8 v18 = cond(mod(_n,17)==0,"",substr(string(mod(_n+18,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen float v19 = (_n+19)/7
quietly replace v19 = . if mod(_n,101)==0
quietly replace v19 = .z if mod(_n,103)==0
quietly gen str8 v20 = cond(mod(_n,17)==0,"",substr(string(mod(_n+20,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen byte v21 = mod(_n+21,101)-50
quietly replace v21 = . if mod(_n,101)==0
quietly replace v21 = .z if mod(_n,103)==0
quietly gen str8 v22 = cond(mod(_n,17)==0,"",substr(string(mod(_n+22,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v23 = _n+23
quietly replace v23 = . if mod(_n,101)==0
quietly replace v23 = .z if mod(_n,103)==0
quietly gen str8 v24 = cond(mod(_n,17)==0,"",substr(string(mod(_n+24,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen double v25 = (_n+25)/7
quietly replace v25 = . if mod(_n,101)==0
quietly replace v25 = .z if mod(_n,103)==0
quietly gen str8 v26 = cond(mod(_n,17)==0,"",substr(string(mod(_n+26,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen int v27 = mod(_n+27,101)-50
quietly replace v27 = . if mod(_n,101)==0
quietly replace v27 = .z if mod(_n,103)==0
quietly gen str8 v28 = cond(mod(_n,17)==0,"",substr(string(mod(_n+28,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen float v29 = (_n+29)/7
quietly replace v29 = . if mod(_n,101)==0
quietly replace v29 = .z if mod(_n,103)==0
quietly gen str8 v30 = cond(mod(_n,17)==0,"",substr(string(mod(_n+30,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen byte v31 = mod(_n+31,101)-50
quietly replace v31 = . if mod(_n,101)==0
quietly replace v31 = .z if mod(_n,103)==0
quietly gen str8 v32 = cond(mod(_n,17)==0,"",substr(string(mod(_n+32,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v33 = _n+33
quietly replace v33 = . if mod(_n,101)==0
quietly replace v33 = .z if mod(_n,103)==0
quietly gen str8 v34 = cond(mod(_n,17)==0,"",substr(string(mod(_n+34,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen double v35 = (_n+35)/7
quietly replace v35 = . if mod(_n,101)==0
quietly replace v35 = .z if mod(_n,103)==0
quietly gen str8 v36 = cond(mod(_n,17)==0,"",substr(string(mod(_n+36,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen int v37 = mod(_n+37,101)-50
quietly replace v37 = . if mod(_n,101)==0
quietly replace v37 = .z if mod(_n,103)==0
quietly gen str8 v38 = cond(mod(_n,17)==0,"",substr(string(mod(_n+38,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen float v39 = (_n+39)/7
quietly replace v39 = . if mod(_n,101)==0
quietly replace v39 = .z if mod(_n,103)==0
quietly gen str8 v40 = cond(mod(_n,17)==0,"",substr(string(mod(_n+40,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen byte v41 = mod(_n+41,101)-50
quietly replace v41 = . if mod(_n,101)==0
quietly replace v41 = .z if mod(_n,103)==0
quietly gen str8 v42 = cond(mod(_n,17)==0,"",substr(string(mod(_n+42,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v43 = _n+43
quietly replace v43 = . if mod(_n,101)==0
quietly replace v43 = .z if mod(_n,103)==0
quietly gen str8 v44 = cond(mod(_n,17)==0,"",substr(string(mod(_n+44,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen double v45 = (_n+45)/7
quietly replace v45 = . if mod(_n,101)==0
quietly replace v45 = .z if mod(_n,103)==0
quietly gen str8 v46 = cond(mod(_n,17)==0,"",substr(string(mod(_n+46,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen int v47 = mod(_n+47,101)-50
quietly replace v47 = . if mod(_n,101)==0
quietly replace v47 = .z if mod(_n,103)==0
quietly gen str8 v48 = cond(mod(_n,17)==0,"",substr(string(mod(_n+48,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen float v49 = (_n+49)/7
quietly replace v49 = . if mod(_n,101)==0
quietly replace v49 = .z if mod(_n,103)==0
quietly gen str8 v50 = cond(mod(_n,17)==0,"",substr(string(mod(_n+50,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen byte v51 = mod(_n+51,101)-50
quietly replace v51 = . if mod(_n,101)==0
quietly replace v51 = .z if mod(_n,103)==0
quietly gen str8 v52 = cond(mod(_n,17)==0,"",substr(string(mod(_n+52,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v53 = _n+53
quietly replace v53 = . if mod(_n,101)==0
quietly replace v53 = .z if mod(_n,103)==0
quietly gen str8 v54 = cond(mod(_n,17)==0,"",substr(string(mod(_n+54,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen double v55 = (_n+55)/7
quietly replace v55 = . if mod(_n,101)==0
quietly replace v55 = .z if mod(_n,103)==0
quietly gen str8 v56 = cond(mod(_n,17)==0,"",substr(string(mod(_n+56,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen int v57 = mod(_n+57,101)-50
quietly replace v57 = . if mod(_n,101)==0
quietly replace v57 = .z if mod(_n,103)==0
quietly gen str8 v58 = cond(mod(_n,17)==0,"",substr(string(mod(_n+58,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen float v59 = (_n+59)/7
quietly replace v59 = . if mod(_n,101)==0
quietly replace v59 = .z if mod(_n,103)==0
quietly gen str8 v60 = cond(mod(_n,17)==0,"",substr(string(mod(_n+60,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen byte v61 = mod(_n+61,101)-50
quietly replace v61 = . if mod(_n,101)==0
quietly replace v61 = .z if mod(_n,103)==0
quietly gen str8 v62 = cond(mod(_n,17)==0,"",substr(string(mod(_n+62,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly gen long v63 = _n+63
quietly replace v63 = . if mod(_n,101)==0
quietly replace v63 = .z if mod(_n,103)==0
quietly gen str8 v64 = cond(mod(_n,17)==0,"",substr(string(mod(_n+64,1000000),"%06.0f")+"`transport_pad'",1,8))
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r0" "8" "sorted" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r1" "8" "sorted" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r2" "8" "sorted" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r3" "8" "sorted" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r4" "8" "sorted" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r5" "8" "sorted" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64, "interleave_n200000_k64_mixed8_sorted_r6" "8" "sorted" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
