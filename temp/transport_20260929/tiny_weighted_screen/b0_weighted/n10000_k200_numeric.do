clear
quietly set obs 10000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
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
quietly gen float v129 = (_n+129)/7
quietly replace v129 = . if mod(_n,101)==0
quietly replace v129 = .z if mod(_n,103)==0
quietly gen double v130 = (_n+130)/7
quietly replace v130 = . if mod(_n,101)==0
quietly replace v130 = .z if mod(_n,103)==0
quietly gen byte v131 = mod(_n+131,101)-50
quietly replace v131 = . if mod(_n,101)==0
quietly replace v131 = .z if mod(_n,103)==0
quietly gen int v132 = mod(_n+132,101)-50
quietly replace v132 = . if mod(_n,101)==0
quietly replace v132 = .z if mod(_n,103)==0
quietly gen long v133 = _n+133
quietly replace v133 = . if mod(_n,101)==0
quietly replace v133 = .z if mod(_n,103)==0
quietly gen float v134 = (_n+134)/7
quietly replace v134 = . if mod(_n,101)==0
quietly replace v134 = .z if mod(_n,103)==0
quietly gen double v135 = (_n+135)/7
quietly replace v135 = . if mod(_n,101)==0
quietly replace v135 = .z if mod(_n,103)==0
quietly gen byte v136 = mod(_n+136,101)-50
quietly replace v136 = . if mod(_n,101)==0
quietly replace v136 = .z if mod(_n,103)==0
quietly gen int v137 = mod(_n+137,101)-50
quietly replace v137 = . if mod(_n,101)==0
quietly replace v137 = .z if mod(_n,103)==0
quietly gen long v138 = _n+138
quietly replace v138 = . if mod(_n,101)==0
quietly replace v138 = .z if mod(_n,103)==0
quietly gen float v139 = (_n+139)/7
quietly replace v139 = . if mod(_n,101)==0
quietly replace v139 = .z if mod(_n,103)==0
quietly gen double v140 = (_n+140)/7
quietly replace v140 = . if mod(_n,101)==0
quietly replace v140 = .z if mod(_n,103)==0
quietly gen byte v141 = mod(_n+141,101)-50
quietly replace v141 = . if mod(_n,101)==0
quietly replace v141 = .z if mod(_n,103)==0
quietly gen int v142 = mod(_n+142,101)-50
quietly replace v142 = . if mod(_n,101)==0
quietly replace v142 = .z if mod(_n,103)==0
quietly gen long v143 = _n+143
quietly replace v143 = . if mod(_n,101)==0
quietly replace v143 = .z if mod(_n,103)==0
quietly gen float v144 = (_n+144)/7
quietly replace v144 = . if mod(_n,101)==0
quietly replace v144 = .z if mod(_n,103)==0
quietly gen double v145 = (_n+145)/7
quietly replace v145 = . if mod(_n,101)==0
quietly replace v145 = .z if mod(_n,103)==0
quietly gen byte v146 = mod(_n+146,101)-50
quietly replace v146 = . if mod(_n,101)==0
quietly replace v146 = .z if mod(_n,103)==0
quietly gen int v147 = mod(_n+147,101)-50
quietly replace v147 = . if mod(_n,101)==0
quietly replace v147 = .z if mod(_n,103)==0
quietly gen long v148 = _n+148
quietly replace v148 = . if mod(_n,101)==0
quietly replace v148 = .z if mod(_n,103)==0
quietly gen float v149 = (_n+149)/7
quietly replace v149 = . if mod(_n,101)==0
quietly replace v149 = .z if mod(_n,103)==0
quietly gen double v150 = (_n+150)/7
quietly replace v150 = . if mod(_n,101)==0
quietly replace v150 = .z if mod(_n,103)==0
quietly gen byte v151 = mod(_n+151,101)-50
quietly replace v151 = . if mod(_n,101)==0
quietly replace v151 = .z if mod(_n,103)==0
quietly gen int v152 = mod(_n+152,101)-50
quietly replace v152 = . if mod(_n,101)==0
quietly replace v152 = .z if mod(_n,103)==0
quietly gen long v153 = _n+153
quietly replace v153 = . if mod(_n,101)==0
quietly replace v153 = .z if mod(_n,103)==0
quietly gen float v154 = (_n+154)/7
quietly replace v154 = . if mod(_n,101)==0
quietly replace v154 = .z if mod(_n,103)==0
quietly gen double v155 = (_n+155)/7
quietly replace v155 = . if mod(_n,101)==0
quietly replace v155 = .z if mod(_n,103)==0
quietly gen byte v156 = mod(_n+156,101)-50
quietly replace v156 = . if mod(_n,101)==0
quietly replace v156 = .z if mod(_n,103)==0
quietly gen int v157 = mod(_n+157,101)-50
quietly replace v157 = . if mod(_n,101)==0
quietly replace v157 = .z if mod(_n,103)==0
quietly gen long v158 = _n+158
quietly replace v158 = . if mod(_n,101)==0
quietly replace v158 = .z if mod(_n,103)==0
quietly gen float v159 = (_n+159)/7
quietly replace v159 = . if mod(_n,101)==0
quietly replace v159 = .z if mod(_n,103)==0
quietly gen double v160 = (_n+160)/7
quietly replace v160 = . if mod(_n,101)==0
quietly replace v160 = .z if mod(_n,103)==0
quietly gen byte v161 = mod(_n+161,101)-50
quietly replace v161 = . if mod(_n,101)==0
quietly replace v161 = .z if mod(_n,103)==0
quietly gen int v162 = mod(_n+162,101)-50
quietly replace v162 = . if mod(_n,101)==0
quietly replace v162 = .z if mod(_n,103)==0
quietly gen long v163 = _n+163
quietly replace v163 = . if mod(_n,101)==0
quietly replace v163 = .z if mod(_n,103)==0
quietly gen float v164 = (_n+164)/7
quietly replace v164 = . if mod(_n,101)==0
quietly replace v164 = .z if mod(_n,103)==0
quietly gen double v165 = (_n+165)/7
quietly replace v165 = . if mod(_n,101)==0
quietly replace v165 = .z if mod(_n,103)==0
quietly gen byte v166 = mod(_n+166,101)-50
quietly replace v166 = . if mod(_n,101)==0
quietly replace v166 = .z if mod(_n,103)==0
quietly gen int v167 = mod(_n+167,101)-50
quietly replace v167 = . if mod(_n,101)==0
quietly replace v167 = .z if mod(_n,103)==0
quietly gen long v168 = _n+168
quietly replace v168 = . if mod(_n,101)==0
quietly replace v168 = .z if mod(_n,103)==0
quietly gen float v169 = (_n+169)/7
quietly replace v169 = . if mod(_n,101)==0
quietly replace v169 = .z if mod(_n,103)==0
quietly gen double v170 = (_n+170)/7
quietly replace v170 = . if mod(_n,101)==0
quietly replace v170 = .z if mod(_n,103)==0
quietly gen byte v171 = mod(_n+171,101)-50
quietly replace v171 = . if mod(_n,101)==0
quietly replace v171 = .z if mod(_n,103)==0
quietly gen int v172 = mod(_n+172,101)-50
quietly replace v172 = . if mod(_n,101)==0
quietly replace v172 = .z if mod(_n,103)==0
quietly gen long v173 = _n+173
quietly replace v173 = . if mod(_n,101)==0
quietly replace v173 = .z if mod(_n,103)==0
quietly gen float v174 = (_n+174)/7
quietly replace v174 = . if mod(_n,101)==0
quietly replace v174 = .z if mod(_n,103)==0
quietly gen double v175 = (_n+175)/7
quietly replace v175 = . if mod(_n,101)==0
quietly replace v175 = .z if mod(_n,103)==0
quietly gen byte v176 = mod(_n+176,101)-50
quietly replace v176 = . if mod(_n,101)==0
quietly replace v176 = .z if mod(_n,103)==0
quietly gen int v177 = mod(_n+177,101)-50
quietly replace v177 = . if mod(_n,101)==0
quietly replace v177 = .z if mod(_n,103)==0
quietly gen long v178 = _n+178
quietly replace v178 = . if mod(_n,101)==0
quietly replace v178 = .z if mod(_n,103)==0
quietly gen float v179 = (_n+179)/7
quietly replace v179 = . if mod(_n,101)==0
quietly replace v179 = .z if mod(_n,103)==0
quietly gen double v180 = (_n+180)/7
quietly replace v180 = . if mod(_n,101)==0
quietly replace v180 = .z if mod(_n,103)==0
quietly gen byte v181 = mod(_n+181,101)-50
quietly replace v181 = . if mod(_n,101)==0
quietly replace v181 = .z if mod(_n,103)==0
quietly gen int v182 = mod(_n+182,101)-50
quietly replace v182 = . if mod(_n,101)==0
quietly replace v182 = .z if mod(_n,103)==0
quietly gen long v183 = _n+183
quietly replace v183 = . if mod(_n,101)==0
quietly replace v183 = .z if mod(_n,103)==0
quietly gen float v184 = (_n+184)/7
quietly replace v184 = . if mod(_n,101)==0
quietly replace v184 = .z if mod(_n,103)==0
quietly gen double v185 = (_n+185)/7
quietly replace v185 = . if mod(_n,101)==0
quietly replace v185 = .z if mod(_n,103)==0
quietly gen byte v186 = mod(_n+186,101)-50
quietly replace v186 = . if mod(_n,101)==0
quietly replace v186 = .z if mod(_n,103)==0
quietly gen int v187 = mod(_n+187,101)-50
quietly replace v187 = . if mod(_n,101)==0
quietly replace v187 = .z if mod(_n,103)==0
quietly gen long v188 = _n+188
quietly replace v188 = . if mod(_n,101)==0
quietly replace v188 = .z if mod(_n,103)==0
quietly gen float v189 = (_n+189)/7
quietly replace v189 = . if mod(_n,101)==0
quietly replace v189 = .z if mod(_n,103)==0
quietly gen double v190 = (_n+190)/7
quietly replace v190 = . if mod(_n,101)==0
quietly replace v190 = .z if mod(_n,103)==0
quietly gen byte v191 = mod(_n+191,101)-50
quietly replace v191 = . if mod(_n,101)==0
quietly replace v191 = .z if mod(_n,103)==0
quietly gen int v192 = mod(_n+192,101)-50
quietly replace v192 = . if mod(_n,101)==0
quietly replace v192 = .z if mod(_n,103)==0
quietly gen long v193 = _n+193
quietly replace v193 = . if mod(_n,101)==0
quietly replace v193 = .z if mod(_n,103)==0
quietly gen float v194 = (_n+194)/7
quietly replace v194 = . if mod(_n,101)==0
quietly replace v194 = .z if mod(_n,103)==0
quietly gen double v195 = (_n+195)/7
quietly replace v195 = . if mod(_n,101)==0
quietly replace v195 = .z if mod(_n,103)==0
quietly gen byte v196 = mod(_n+196,101)-50
quietly replace v196 = . if mod(_n,101)==0
quietly replace v196 = .z if mod(_n,103)==0
quietly gen int v197 = mod(_n+197,101)-50
quietly replace v197 = . if mod(_n,101)==0
quietly replace v197 = .z if mod(_n,103)==0
quietly gen long v198 = _n+198
quietly replace v198 = . if mod(_n,101)==0
quietly replace v198 = .z if mod(_n,103)==0
quietly gen float v199 = (_n+199)/7
quietly replace v199 = . if mod(_n,101)==0
quietly replace v199 = .z if mod(_n,103)==0
quietly gen double v200 = (_n+200)/7
quietly replace v200 = . if mod(_n,101)==0
quietly replace v200 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r0" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r1" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r2" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r3" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r4" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r5" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r6" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r7" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r8" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r9" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r10" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r11" "8" "identity" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200, "weighted_n10000_k200_numeric_r12" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
