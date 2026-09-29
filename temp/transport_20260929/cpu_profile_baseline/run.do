clear all
set more off
set linesize 255
log using "/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/cpu_profile_baseline/transport.log", text replace
adopath ++ "/Users/Mike/Documents/GitHub/stata-ctools/build"
program io0, plugin using("/Users/Mike/Documents/GitHub/stata-ctools/temp/transport_20260929/cpu_profile_baseline/baseline.plugin")
program run_benchmark
version 16
clear
quietly set obs 200000
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen int v2 = mod(_n+2,101)-50
quietly replace v2 = . if mod(_n,101)==0
quietly replace v2 = .z if mod(_n,103)==0
quietly gen long v3 = _n+3
quietly replace v3 = . if mod(_n,101)==0
quietly replace v3 = .z if mod(_n,103)==0
quietly gen float v4 = (_n+4)/7
quietly replace v4 = . if mod(_n,101)==0
quietly replace v4 = .z if mod(_n,103)==0
quietly gen double v5 = (_n+5)/7
quietly replace v5 = . if mod(_n,101)==0
quietly replace v5 = .z if mod(_n,103)==0
quietly gen byte v6 = mod(_n+6,101)-50
quietly replace v6 = . if mod(_n,101)==0
quietly replace v6 = .z if mod(_n,103)==0
quietly gen int v7 = mod(_n+7,101)-50
quietly replace v7 = . if mod(_n,101)==0
quietly replace v7 = .z if mod(_n,103)==0
quietly gen long v8 = _n+8
quietly replace v8 = . if mod(_n,101)==0
quietly replace v8 = .z if mod(_n,103)==0
quietly gen float v9 = (_n+9)/7
quietly replace v9 = . if mod(_n,101)==0
quietly replace v9 = .z if mod(_n,103)==0
quietly gen double v10 = (_n+10)/7
quietly replace v10 = . if mod(_n,101)==0
quietly replace v10 = .z if mod(_n,103)==0
quietly gen byte v11 = mod(_n+11,101)-50
quietly replace v11 = . if mod(_n,101)==0
quietly replace v11 = .z if mod(_n,103)==0
quietly gen int v12 = mod(_n+12,101)-50
quietly replace v12 = . if mod(_n,101)==0
quietly replace v12 = .z if mod(_n,103)==0
quietly gen long v13 = _n+13
quietly replace v13 = . if mod(_n,101)==0
quietly replace v13 = .z if mod(_n,103)==0
quietly gen float v14 = (_n+14)/7
quietly replace v14 = . if mod(_n,101)==0
quietly replace v14 = .z if mod(_n,103)==0
quietly gen double v15 = (_n+15)/7
quietly replace v15 = . if mod(_n,101)==0
quietly replace v15 = .z if mod(_n,103)==0
quietly gen byte v16 = mod(_n+16,101)-50
quietly replace v16 = . if mod(_n,101)==0
quietly replace v16 = .z if mod(_n,103)==0
quietly gen int v17 = mod(_n+17,101)-50
quietly replace v17 = . if mod(_n,101)==0
quietly replace v17 = .z if mod(_n,103)==0
quietly gen long v18 = _n+18
quietly replace v18 = . if mod(_n,101)==0
quietly replace v18 = .z if mod(_n,103)==0
quietly gen float v19 = (_n+19)/7
quietly replace v19 = . if mod(_n,101)==0
quietly replace v19 = .z if mod(_n,103)==0
quietly gen double v20 = (_n+20)/7
quietly replace v20 = . if mod(_n,101)==0
quietly replace v20 = .z if mod(_n,103)==0
quietly gen byte v21 = mod(_n+21,101)-50
quietly replace v21 = . if mod(_n,101)==0
quietly replace v21 = .z if mod(_n,103)==0
quietly gen int v22 = mod(_n+22,101)-50
quietly replace v22 = . if mod(_n,101)==0
quietly replace v22 = .z if mod(_n,103)==0
quietly gen long v23 = _n+23
quietly replace v23 = . if mod(_n,101)==0
quietly replace v23 = .z if mod(_n,103)==0
quietly gen float v24 = (_n+24)/7
quietly replace v24 = . if mod(_n,101)==0
quietly replace v24 = .z if mod(_n,103)==0
quietly gen double v25 = (_n+25)/7
quietly replace v25 = . if mod(_n,101)==0
quietly replace v25 = .z if mod(_n,103)==0
quietly gen byte v26 = mod(_n+26,101)-50
quietly replace v26 = . if mod(_n,101)==0
quietly replace v26 = .z if mod(_n,103)==0
quietly gen int v27 = mod(_n+27,101)-50
quietly replace v27 = . if mod(_n,101)==0
quietly replace v27 = .z if mod(_n,103)==0
quietly gen long v28 = _n+28
quietly replace v28 = . if mod(_n,101)==0
quietly replace v28 = .z if mod(_n,103)==0
quietly gen float v29 = (_n+29)/7
quietly replace v29 = . if mod(_n,101)==0
quietly replace v29 = .z if mod(_n,103)==0
quietly gen double v30 = (_n+30)/7
quietly replace v30 = . if mod(_n,101)==0
quietly replace v30 = .z if mod(_n,103)==0
quietly gen byte v31 = mod(_n+31,101)-50
quietly replace v31 = . if mod(_n,101)==0
quietly replace v31 = .z if mod(_n,103)==0
quietly gen int v32 = mod(_n+32,101)-50
quietly replace v32 = . if mod(_n,101)==0
quietly replace v32 = .z if mod(_n,103)==0
quietly gen long v33 = _n+33
quietly replace v33 = . if mod(_n,101)==0
quietly replace v33 = .z if mod(_n,103)==0
quietly gen float v34 = (_n+34)/7
quietly replace v34 = . if mod(_n,101)==0
quietly replace v34 = .z if mod(_n,103)==0
quietly gen double v35 = (_n+35)/7
quietly replace v35 = . if mod(_n,101)==0
quietly replace v35 = .z if mod(_n,103)==0
quietly gen byte v36 = mod(_n+36,101)-50
quietly replace v36 = . if mod(_n,101)==0
quietly replace v36 = .z if mod(_n,103)==0
quietly gen int v37 = mod(_n+37,101)-50
quietly replace v37 = . if mod(_n,101)==0
quietly replace v37 = .z if mod(_n,103)==0
quietly gen long v38 = _n+38
quietly replace v38 = . if mod(_n,101)==0
quietly replace v38 = .z if mod(_n,103)==0
quietly gen float v39 = (_n+39)/7
quietly replace v39 = . if mod(_n,101)==0
quietly replace v39 = .z if mod(_n,103)==0
quietly gen double v40 = (_n+40)/7
quietly replace v40 = . if mod(_n,101)==0
quietly replace v40 = .z if mod(_n,103)==0
quietly gen byte v41 = mod(_n+41,101)-50
quietly replace v41 = . if mod(_n,101)==0
quietly replace v41 = .z if mod(_n,103)==0
quietly gen int v42 = mod(_n+42,101)-50
quietly replace v42 = . if mod(_n,101)==0
quietly replace v42 = .z if mod(_n,103)==0
quietly gen long v43 = _n+43
quietly replace v43 = . if mod(_n,101)==0
quietly replace v43 = .z if mod(_n,103)==0
quietly gen float v44 = (_n+44)/7
quietly replace v44 = . if mod(_n,101)==0
quietly replace v44 = .z if mod(_n,103)==0
quietly gen double v45 = (_n+45)/7
quietly replace v45 = . if mod(_n,101)==0
quietly replace v45 = .z if mod(_n,103)==0
quietly gen byte v46 = mod(_n+46,101)-50
quietly replace v46 = . if mod(_n,101)==0
quietly replace v46 = .z if mod(_n,103)==0
quietly gen int v47 = mod(_n+47,101)-50
quietly replace v47 = . if mod(_n,101)==0
quietly replace v47 = .z if mod(_n,103)==0
quietly gen long v48 = _n+48
quietly replace v48 = . if mod(_n,101)==0
quietly replace v48 = .z if mod(_n,103)==0
quietly gen float v49 = (_n+49)/7
quietly replace v49 = . if mod(_n,101)==0
quietly replace v49 = .z if mod(_n,103)==0
quietly gen double v50 = (_n+50)/7
quietly replace v50 = . if mod(_n,101)==0
quietly replace v50 = .z if mod(_n,103)==0
quietly gen byte v51 = mod(_n+51,101)-50
quietly replace v51 = . if mod(_n,101)==0
quietly replace v51 = .z if mod(_n,103)==0
quietly gen int v52 = mod(_n+52,101)-50
quietly replace v52 = . if mod(_n,101)==0
quietly replace v52 = .z if mod(_n,103)==0
quietly gen long v53 = _n+53
quietly replace v53 = . if mod(_n,101)==0
quietly replace v53 = .z if mod(_n,103)==0
quietly gen float v54 = (_n+54)/7
quietly replace v54 = . if mod(_n,101)==0
quietly replace v54 = .z if mod(_n,103)==0
quietly gen double v55 = (_n+55)/7
quietly replace v55 = . if mod(_n,101)==0
quietly replace v55 = .z if mod(_n,103)==0
quietly gen byte v56 = mod(_n+56,101)-50
quietly replace v56 = . if mod(_n,101)==0
quietly replace v56 = .z if mod(_n,103)==0
quietly gen int v57 = mod(_n+57,101)-50
quietly replace v57 = . if mod(_n,101)==0
quietly replace v57 = .z if mod(_n,103)==0
quietly gen long v58 = _n+58
quietly replace v58 = . if mod(_n,101)==0
quietly replace v58 = .z if mod(_n,103)==0
quietly gen float v59 = (_n+59)/7
quietly replace v59 = . if mod(_n,101)==0
quietly replace v59 = .z if mod(_n,103)==0
quietly gen double v60 = (_n+60)/7
quietly replace v60 = . if mod(_n,101)==0
quietly replace v60 = .z if mod(_n,103)==0
quietly gen byte v61 = mod(_n+61,101)-50
quietly replace v61 = . if mod(_n,101)==0
quietly replace v61 = .z if mod(_n,103)==0
quietly gen int v62 = mod(_n+62,101)-50
quietly replace v62 = . if mod(_n,101)==0
quietly replace v62 = .z if mod(_n,103)==0
quietly gen long v63 = _n+63
quietly replace v63 = . if mod(_n,101)==0
quietly replace v63 = .z if mod(_n,103)==0
quietly gen float v64 = (_n+64)/7
quietly replace v64 = . if mod(_n,101)==0
quietly replace v64 = .z if mod(_n,103)==0
quietly gen double v65 = (_n+65)/7
quietly replace v65 = . if mod(_n,101)==0
quietly replace v65 = .z if mod(_n,103)==0
quietly gen byte v66 = mod(_n+66,101)-50
quietly replace v66 = . if mod(_n,101)==0
quietly replace v66 = .z if mod(_n,103)==0
quietly gen int v67 = mod(_n+67,101)-50
quietly replace v67 = . if mod(_n,101)==0
quietly replace v67 = .z if mod(_n,103)==0
quietly gen long v68 = _n+68
quietly replace v68 = . if mod(_n,101)==0
quietly replace v68 = .z if mod(_n,103)==0
quietly gen float v69 = (_n+69)/7
quietly replace v69 = . if mod(_n,101)==0
quietly replace v69 = .z if mod(_n,103)==0
quietly gen double v70 = (_n+70)/7
quietly replace v70 = . if mod(_n,101)==0
quietly replace v70 = .z if mod(_n,103)==0
quietly gen byte v71 = mod(_n+71,101)-50
quietly replace v71 = . if mod(_n,101)==0
quietly replace v71 = .z if mod(_n,103)==0
quietly gen int v72 = mod(_n+72,101)-50
quietly replace v72 = . if mod(_n,101)==0
quietly replace v72 = .z if mod(_n,103)==0
quietly gen long v73 = _n+73
quietly replace v73 = . if mod(_n,101)==0
quietly replace v73 = .z if mod(_n,103)==0
quietly gen float v74 = (_n+74)/7
quietly replace v74 = . if mod(_n,101)==0
quietly replace v74 = .z if mod(_n,103)==0
quietly gen double v75 = (_n+75)/7
quietly replace v75 = . if mod(_n,101)==0
quietly replace v75 = .z if mod(_n,103)==0
quietly gen byte v76 = mod(_n+76,101)-50
quietly replace v76 = . if mod(_n,101)==0
quietly replace v76 = .z if mod(_n,103)==0
quietly gen int v77 = mod(_n+77,101)-50
quietly replace v77 = . if mod(_n,101)==0
quietly replace v77 = .z if mod(_n,103)==0
quietly gen long v78 = _n+78
quietly replace v78 = . if mod(_n,101)==0
quietly replace v78 = .z if mod(_n,103)==0
quietly gen float v79 = (_n+79)/7
quietly replace v79 = . if mod(_n,101)==0
quietly replace v79 = .z if mod(_n,103)==0
quietly gen double v80 = (_n+80)/7
quietly replace v80 = . if mod(_n,101)==0
quietly replace v80 = .z if mod(_n,103)==0
quietly gen byte v81 = mod(_n+81,101)-50
quietly replace v81 = . if mod(_n,101)==0
quietly replace v81 = .z if mod(_n,103)==0
quietly gen int v82 = mod(_n+82,101)-50
quietly replace v82 = . if mod(_n,101)==0
quietly replace v82 = .z if mod(_n,103)==0
quietly gen long v83 = _n+83
quietly replace v83 = . if mod(_n,101)==0
quietly replace v83 = .z if mod(_n,103)==0
quietly gen float v84 = (_n+84)/7
quietly replace v84 = . if mod(_n,101)==0
quietly replace v84 = .z if mod(_n,103)==0
quietly gen double v85 = (_n+85)/7
quietly replace v85 = . if mod(_n,101)==0
quietly replace v85 = .z if mod(_n,103)==0
quietly gen byte v86 = mod(_n+86,101)-50
quietly replace v86 = . if mod(_n,101)==0
quietly replace v86 = .z if mod(_n,103)==0
quietly gen int v87 = mod(_n+87,101)-50
quietly replace v87 = . if mod(_n,101)==0
quietly replace v87 = .z if mod(_n,103)==0
quietly gen long v88 = _n+88
quietly replace v88 = . if mod(_n,101)==0
quietly replace v88 = .z if mod(_n,103)==0
quietly gen float v89 = (_n+89)/7
quietly replace v89 = . if mod(_n,101)==0
quietly replace v89 = .z if mod(_n,103)==0
quietly gen double v90 = (_n+90)/7
quietly replace v90 = . if mod(_n,101)==0
quietly replace v90 = .z if mod(_n,103)==0
quietly gen byte v91 = mod(_n+91,101)-50
quietly replace v91 = . if mod(_n,101)==0
quietly replace v91 = .z if mod(_n,103)==0
quietly gen int v92 = mod(_n+92,101)-50
quietly replace v92 = . if mod(_n,101)==0
quietly replace v92 = .z if mod(_n,103)==0
quietly gen long v93 = _n+93
quietly replace v93 = . if mod(_n,101)==0
quietly replace v93 = .z if mod(_n,103)==0
quietly gen float v94 = (_n+94)/7
quietly replace v94 = . if mod(_n,101)==0
quietly replace v94 = .z if mod(_n,103)==0
quietly gen double v95 = (_n+95)/7
quietly replace v95 = . if mod(_n,101)==0
quietly replace v95 = .z if mod(_n,103)==0
quietly gen byte v96 = mod(_n+96,101)-50
quietly replace v96 = . if mod(_n,101)==0
quietly replace v96 = .z if mod(_n,103)==0
quietly gen int v97 = mod(_n+97,101)-50
quietly replace v97 = . if mod(_n,101)==0
quietly replace v97 = .z if mod(_n,103)==0
quietly gen long v98 = _n+98
quietly replace v98 = . if mod(_n,101)==0
quietly replace v98 = .z if mod(_n,103)==0
quietly gen float v99 = (_n+99)/7
quietly replace v99 = . if mod(_n,101)==0
quietly replace v99 = .z if mod(_n,103)==0
quietly gen double v100 = (_n+100)/7
quietly replace v100 = . if mod(_n,101)==0
quietly replace v100 = .z if mod(_n,103)==0
quietly gen byte v101 = mod(_n+101,101)-50
quietly replace v101 = . if mod(_n,101)==0
quietly replace v101 = .z if mod(_n,103)==0
quietly gen int v102 = mod(_n+102,101)-50
quietly replace v102 = . if mod(_n,101)==0
quietly replace v102 = .z if mod(_n,103)==0
quietly gen long v103 = _n+103
quietly replace v103 = . if mod(_n,101)==0
quietly replace v103 = .z if mod(_n,103)==0
quietly gen float v104 = (_n+104)/7
quietly replace v104 = . if mod(_n,101)==0
quietly replace v104 = .z if mod(_n,103)==0
quietly gen double v105 = (_n+105)/7
quietly replace v105 = . if mod(_n,101)==0
quietly replace v105 = .z if mod(_n,103)==0
quietly gen byte v106 = mod(_n+106,101)-50
quietly replace v106 = . if mod(_n,101)==0
quietly replace v106 = .z if mod(_n,103)==0
quietly gen int v107 = mod(_n+107,101)-50
quietly replace v107 = . if mod(_n,101)==0
quietly replace v107 = .z if mod(_n,103)==0
quietly gen long v108 = _n+108
quietly replace v108 = . if mod(_n,101)==0
quietly replace v108 = .z if mod(_n,103)==0
quietly gen float v109 = (_n+109)/7
quietly replace v109 = . if mod(_n,101)==0
quietly replace v109 = .z if mod(_n,103)==0
quietly gen double v110 = (_n+110)/7
quietly replace v110 = . if mod(_n,101)==0
quietly replace v110 = .z if mod(_n,103)==0
quietly gen byte v111 = mod(_n+111,101)-50
quietly replace v111 = . if mod(_n,101)==0
quietly replace v111 = .z if mod(_n,103)==0
quietly gen int v112 = mod(_n+112,101)-50
quietly replace v112 = . if mod(_n,101)==0
quietly replace v112 = .z if mod(_n,103)==0
quietly gen long v113 = _n+113
quietly replace v113 = . if mod(_n,101)==0
quietly replace v113 = .z if mod(_n,103)==0
quietly gen float v114 = (_n+114)/7
quietly replace v114 = . if mod(_n,101)==0
quietly replace v114 = .z if mod(_n,103)==0
quietly gen double v115 = (_n+115)/7
quietly replace v115 = . if mod(_n,101)==0
quietly replace v115 = .z if mod(_n,103)==0
quietly gen byte v116 = mod(_n+116,101)-50
quietly replace v116 = . if mod(_n,101)==0
quietly replace v116 = .z if mod(_n,103)==0
quietly gen int v117 = mod(_n+117,101)-50
quietly replace v117 = . if mod(_n,101)==0
quietly replace v117 = .z if mod(_n,103)==0
quietly gen long v118 = _n+118
quietly replace v118 = . if mod(_n,101)==0
quietly replace v118 = .z if mod(_n,103)==0
quietly gen float v119 = (_n+119)/7
quietly replace v119 = . if mod(_n,101)==0
quietly replace v119 = .z if mod(_n,103)==0
quietly gen double v120 = (_n+120)/7
quietly replace v120 = . if mod(_n,101)==0
quietly replace v120 = .z if mod(_n,103)==0
quietly gen byte v121 = mod(_n+121,101)-50
quietly replace v121 = . if mod(_n,101)==0
quietly replace v121 = .z if mod(_n,103)==0
quietly gen int v122 = mod(_n+122,101)-50
quietly replace v122 = . if mod(_n,101)==0
quietly replace v122 = .z if mod(_n,103)==0
quietly gen long v123 = _n+123
quietly replace v123 = . if mod(_n,101)==0
quietly replace v123 = .z if mod(_n,103)==0
quietly gen float v124 = (_n+124)/7
quietly replace v124 = . if mod(_n,101)==0
quietly replace v124 = .z if mod(_n,103)==0
quietly gen double v125 = (_n+125)/7
quietly replace v125 = . if mod(_n,101)==0
quietly replace v125 = .z if mod(_n,103)==0
quietly gen byte v126 = mod(_n+126,101)-50
quietly replace v126 = . if mod(_n,101)==0
quietly replace v126 = .z if mod(_n,103)==0
quietly gen int v127 = mod(_n+127,101)-50
quietly replace v127 = . if mod(_n,101)==0
quietly replace v127 = .z if mod(_n,103)==0
quietly gen long v128 = _n+128
quietly replace v128 = . if mod(_n,101)==0
quietly replace v128 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
local __ctools_strw ""
plugin call io0 v*, "baseline_numeric_w8_h0_r0" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r1" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r2" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r3" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r4" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r5" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r6" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r7" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r8" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r9" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r10" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r11" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r12" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r13" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r14" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r15" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r16" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r17" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r18" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r19" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r20" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r21" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r22" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r23" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r24" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r25" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r26" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r27" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r28" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r29" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r30" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r31" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r32" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r33" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r34" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r35" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r36" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r37" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r38" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r39" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r40" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r41" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r42" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r43" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r44" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r45" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r46" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r47" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r48" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r49" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r50" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r51" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r52" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r53" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r54" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r55" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r56" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r57" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r58" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r59" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r60" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r61" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r62" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r63" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r64" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r65" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r66" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r67" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r68" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r69" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r70" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r71" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r72" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r73" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r74" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r75" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r76" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r77" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r78" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r79" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r80" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r81" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r82" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r83" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r84" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r85" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r86" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r87" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r88" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r89" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r90" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r91" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r92" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r93" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r94" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r95" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r96" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r97" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r98" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r99" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r100" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r101" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r102" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r103" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r104" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r105" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r106" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r107" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r108" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r109" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r110" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r111" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r112" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r113" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r114" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r115" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r116" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r117" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r118" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r119" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r120" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r121" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r122" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r123" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r124" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r125" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r126" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r127" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r128" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r129" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r130" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r131" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r132" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r133" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r134" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r135" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r136" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r137" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r138" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r139" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r140" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r141" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r142" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r143" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r144" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r145" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r146" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r147" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r148" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r149" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r150" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r151" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r152" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r153" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r154" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r155" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r156" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r157" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r158" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r159" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r160" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r161" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r162" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r163" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r164" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r165" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r166" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r167" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r168" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r169" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r170" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r171" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r172" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r173" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r174" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r175" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r176" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r177" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r178" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r179" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r180" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r181" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r182" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r183" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r184" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r185" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r186" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r187" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r188" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r189" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r190" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r191" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r192" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r193" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r194" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r195" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r196" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r197" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r198" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r199" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r200" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r201" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r202" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r203" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r204" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r205" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r206" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r207" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r208" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r209" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r210" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r211" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r212" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r213" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r214" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r215" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r216" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r217" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r218" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r219" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r220" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r221" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r222" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r223" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r224" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r225" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r226" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r227" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r228" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r229" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r230" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r231" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r232" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r233" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r234" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r235" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r236" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r237" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r238" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r239" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r240" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r241" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r242" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r243" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r244" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r245" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r246" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r247" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r248" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r249" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r250" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r251" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r252" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r253" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r254" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r255" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r256" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r257" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r258" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r259" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r260" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r261" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r262" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r263" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r264" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r265" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r266" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r267" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r268" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r269" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r270" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r271" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r272" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r273" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r274" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r275" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r276" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r277" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r278" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r279" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r280" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r281" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r282" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r283" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r284" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r285" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r286" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r287" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r288" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r289" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r290" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r291" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r292" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r293" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r294" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r295" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r296" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r297" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r298" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r299" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r300" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r301" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r302" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r303" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r304" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r305" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r306" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r307" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r308" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r309" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r310" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r311" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r312" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r313" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r314" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r315" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r316" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r317" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r318" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r319" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r320" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r321" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r322" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r323" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r324" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r325" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r326" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r327" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r328" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r329" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r330" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r331" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r332" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r333" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r334" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r335" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r336" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r337" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r338" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r339" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r340" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r341" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r342" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r343" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r344" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r345" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r346" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r347" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r348" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r349" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r350" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r351" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r352" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r353" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r354" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r355" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r356" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r357" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r358" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r359" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r360" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r361" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r362" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r363" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r364" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r365" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r366" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r367" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r368" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r369" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r370" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r371" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r372" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r373" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r374" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r375" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r376" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r377" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r378" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r379" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r380" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r381" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r382" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r383" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r384" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r385" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r386" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r387" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r388" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r389" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r390" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r391" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r392" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r393" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r394" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r395" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r396" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r397" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r398" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r399" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r400" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r401" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r402" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r403" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r404" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r405" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r406" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r407" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r408" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r409" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r410" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r411" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r412" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r413" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r414" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r415" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r416" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r417" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r418" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r419" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r420" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r421" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r422" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r423" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r424" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r425" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r426" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r427" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r428" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r429" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r430" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r431" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r432" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r433" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r434" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r435" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r436" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r437" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r438" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r439" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r440" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r441" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r442" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r443" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r444" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r445" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r446" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r447" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r448" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r449" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r450" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r451" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r452" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r453" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r454" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r455" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r456" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r457" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r458" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r459" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r460" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r461" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r462" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r463" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r464" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r465" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r466" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r467" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r468" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r469" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r470" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r471" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r472" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r473" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r474" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r475" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r476" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r477" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r478" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r479" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r480" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r481" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r482" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r483" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r484" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r485" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r486" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r487" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r488" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r489" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r490" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r491" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r492" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r493" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r494" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r495" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r496" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r497" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r498" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r499" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r500" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r501" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r502" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r503" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r504" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r505" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r506" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r507" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r508" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r509" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r510" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r511" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r512" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r513" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r514" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r515" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r516" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r517" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r518" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r519" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r520" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r521" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r522" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r523" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r524" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r525" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r526" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r527" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r528" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r529" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r530" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r531" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r532" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r533" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r534" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r535" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r536" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r537" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r538" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r539" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r540" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r541" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r542" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r543" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r544" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r545" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r546" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r547" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r548" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r549" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r550" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r551" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r552" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r553" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r554" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r555" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r556" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r557" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r558" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r559" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r560" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r561" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r562" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r563" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r564" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r565" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r566" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r567" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r568" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r569" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r570" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r571" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r572" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r573" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r574" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r575" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r576" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r577" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r578" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r579" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r580" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r581" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r582" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r583" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r584" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r585" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r586" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r587" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r588" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r589" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r590" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r591" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r592" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r593" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r594" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r595" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r596" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r597" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r598" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r599" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r600" "8" "identity" "0"
plugin call io0 v*, "baseline_numeric_w8_h0_r601" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
end
capture noisily run_benchmark
local rc = _rc
di "TRANSPORT_COMPLETE RC=`rc'"
log close
exit, clear
