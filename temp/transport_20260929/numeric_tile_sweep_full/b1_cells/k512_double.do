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
quietly gen double v65 = (_n+65)/7
quietly replace v65 = . if mod(_n,101)==0
quietly replace v65 = .z if mod(_n,103)==0
quietly gen double v66 = (_n+66)/7
quietly replace v66 = . if mod(_n,101)==0
quietly replace v66 = .z if mod(_n,103)==0
quietly gen double v67 = (_n+67)/7
quietly replace v67 = . if mod(_n,101)==0
quietly replace v67 = .z if mod(_n,103)==0
quietly gen double v68 = (_n+68)/7
quietly replace v68 = . if mod(_n,101)==0
quietly replace v68 = .z if mod(_n,103)==0
quietly gen double v69 = (_n+69)/7
quietly replace v69 = . if mod(_n,101)==0
quietly replace v69 = .z if mod(_n,103)==0
quietly gen double v70 = (_n+70)/7
quietly replace v70 = . if mod(_n,101)==0
quietly replace v70 = .z if mod(_n,103)==0
quietly gen double v71 = (_n+71)/7
quietly replace v71 = . if mod(_n,101)==0
quietly replace v71 = .z if mod(_n,103)==0
quietly gen double v72 = (_n+72)/7
quietly replace v72 = . if mod(_n,101)==0
quietly replace v72 = .z if mod(_n,103)==0
quietly gen double v73 = (_n+73)/7
quietly replace v73 = . if mod(_n,101)==0
quietly replace v73 = .z if mod(_n,103)==0
quietly gen double v74 = (_n+74)/7
quietly replace v74 = . if mod(_n,101)==0
quietly replace v74 = .z if mod(_n,103)==0
quietly gen double v75 = (_n+75)/7
quietly replace v75 = . if mod(_n,101)==0
quietly replace v75 = .z if mod(_n,103)==0
quietly gen double v76 = (_n+76)/7
quietly replace v76 = . if mod(_n,101)==0
quietly replace v76 = .z if mod(_n,103)==0
quietly gen double v77 = (_n+77)/7
quietly replace v77 = . if mod(_n,101)==0
quietly replace v77 = .z if mod(_n,103)==0
quietly gen double v78 = (_n+78)/7
quietly replace v78 = . if mod(_n,101)==0
quietly replace v78 = .z if mod(_n,103)==0
quietly gen double v79 = (_n+79)/7
quietly replace v79 = . if mod(_n,101)==0
quietly replace v79 = .z if mod(_n,103)==0
quietly gen double v80 = (_n+80)/7
quietly replace v80 = . if mod(_n,101)==0
quietly replace v80 = .z if mod(_n,103)==0
quietly gen double v81 = (_n+81)/7
quietly replace v81 = . if mod(_n,101)==0
quietly replace v81 = .z if mod(_n,103)==0
quietly gen double v82 = (_n+82)/7
quietly replace v82 = . if mod(_n,101)==0
quietly replace v82 = .z if mod(_n,103)==0
quietly gen double v83 = (_n+83)/7
quietly replace v83 = . if mod(_n,101)==0
quietly replace v83 = .z if mod(_n,103)==0
quietly gen double v84 = (_n+84)/7
quietly replace v84 = . if mod(_n,101)==0
quietly replace v84 = .z if mod(_n,103)==0
quietly gen double v85 = (_n+85)/7
quietly replace v85 = . if mod(_n,101)==0
quietly replace v85 = .z if mod(_n,103)==0
quietly gen double v86 = (_n+86)/7
quietly replace v86 = . if mod(_n,101)==0
quietly replace v86 = .z if mod(_n,103)==0
quietly gen double v87 = (_n+87)/7
quietly replace v87 = . if mod(_n,101)==0
quietly replace v87 = .z if mod(_n,103)==0
quietly gen double v88 = (_n+88)/7
quietly replace v88 = . if mod(_n,101)==0
quietly replace v88 = .z if mod(_n,103)==0
quietly gen double v89 = (_n+89)/7
quietly replace v89 = . if mod(_n,101)==0
quietly replace v89 = .z if mod(_n,103)==0
quietly gen double v90 = (_n+90)/7
quietly replace v90 = . if mod(_n,101)==0
quietly replace v90 = .z if mod(_n,103)==0
quietly gen double v91 = (_n+91)/7
quietly replace v91 = . if mod(_n,101)==0
quietly replace v91 = .z if mod(_n,103)==0
quietly gen double v92 = (_n+92)/7
quietly replace v92 = . if mod(_n,101)==0
quietly replace v92 = .z if mod(_n,103)==0
quietly gen double v93 = (_n+93)/7
quietly replace v93 = . if mod(_n,101)==0
quietly replace v93 = .z if mod(_n,103)==0
quietly gen double v94 = (_n+94)/7
quietly replace v94 = . if mod(_n,101)==0
quietly replace v94 = .z if mod(_n,103)==0
quietly gen double v95 = (_n+95)/7
quietly replace v95 = . if mod(_n,101)==0
quietly replace v95 = .z if mod(_n,103)==0
quietly gen double v96 = (_n+96)/7
quietly replace v96 = . if mod(_n,101)==0
quietly replace v96 = .z if mod(_n,103)==0
quietly gen double v97 = (_n+97)/7
quietly replace v97 = . if mod(_n,101)==0
quietly replace v97 = .z if mod(_n,103)==0
quietly gen double v98 = (_n+98)/7
quietly replace v98 = . if mod(_n,101)==0
quietly replace v98 = .z if mod(_n,103)==0
quietly gen double v99 = (_n+99)/7
quietly replace v99 = . if mod(_n,101)==0
quietly replace v99 = .z if mod(_n,103)==0
quietly gen double v100 = (_n+100)/7
quietly replace v100 = . if mod(_n,101)==0
quietly replace v100 = .z if mod(_n,103)==0
quietly gen double v101 = (_n+101)/7
quietly replace v101 = . if mod(_n,101)==0
quietly replace v101 = .z if mod(_n,103)==0
quietly gen double v102 = (_n+102)/7
quietly replace v102 = . if mod(_n,101)==0
quietly replace v102 = .z if mod(_n,103)==0
quietly gen double v103 = (_n+103)/7
quietly replace v103 = . if mod(_n,101)==0
quietly replace v103 = .z if mod(_n,103)==0
quietly gen double v104 = (_n+104)/7
quietly replace v104 = . if mod(_n,101)==0
quietly replace v104 = .z if mod(_n,103)==0
quietly gen double v105 = (_n+105)/7
quietly replace v105 = . if mod(_n,101)==0
quietly replace v105 = .z if mod(_n,103)==0
quietly gen double v106 = (_n+106)/7
quietly replace v106 = . if mod(_n,101)==0
quietly replace v106 = .z if mod(_n,103)==0
quietly gen double v107 = (_n+107)/7
quietly replace v107 = . if mod(_n,101)==0
quietly replace v107 = .z if mod(_n,103)==0
quietly gen double v108 = (_n+108)/7
quietly replace v108 = . if mod(_n,101)==0
quietly replace v108 = .z if mod(_n,103)==0
quietly gen double v109 = (_n+109)/7
quietly replace v109 = . if mod(_n,101)==0
quietly replace v109 = .z if mod(_n,103)==0
quietly gen double v110 = (_n+110)/7
quietly replace v110 = . if mod(_n,101)==0
quietly replace v110 = .z if mod(_n,103)==0
quietly gen double v111 = (_n+111)/7
quietly replace v111 = . if mod(_n,101)==0
quietly replace v111 = .z if mod(_n,103)==0
quietly gen double v112 = (_n+112)/7
quietly replace v112 = . if mod(_n,101)==0
quietly replace v112 = .z if mod(_n,103)==0
quietly gen double v113 = (_n+113)/7
quietly replace v113 = . if mod(_n,101)==0
quietly replace v113 = .z if mod(_n,103)==0
quietly gen double v114 = (_n+114)/7
quietly replace v114 = . if mod(_n,101)==0
quietly replace v114 = .z if mod(_n,103)==0
quietly gen double v115 = (_n+115)/7
quietly replace v115 = . if mod(_n,101)==0
quietly replace v115 = .z if mod(_n,103)==0
quietly gen double v116 = (_n+116)/7
quietly replace v116 = . if mod(_n,101)==0
quietly replace v116 = .z if mod(_n,103)==0
quietly gen double v117 = (_n+117)/7
quietly replace v117 = . if mod(_n,101)==0
quietly replace v117 = .z if mod(_n,103)==0
quietly gen double v118 = (_n+118)/7
quietly replace v118 = . if mod(_n,101)==0
quietly replace v118 = .z if mod(_n,103)==0
quietly gen double v119 = (_n+119)/7
quietly replace v119 = . if mod(_n,101)==0
quietly replace v119 = .z if mod(_n,103)==0
quietly gen double v120 = (_n+120)/7
quietly replace v120 = . if mod(_n,101)==0
quietly replace v120 = .z if mod(_n,103)==0
quietly gen double v121 = (_n+121)/7
quietly replace v121 = . if mod(_n,101)==0
quietly replace v121 = .z if mod(_n,103)==0
quietly gen double v122 = (_n+122)/7
quietly replace v122 = . if mod(_n,101)==0
quietly replace v122 = .z if mod(_n,103)==0
quietly gen double v123 = (_n+123)/7
quietly replace v123 = . if mod(_n,101)==0
quietly replace v123 = .z if mod(_n,103)==0
quietly gen double v124 = (_n+124)/7
quietly replace v124 = . if mod(_n,101)==0
quietly replace v124 = .z if mod(_n,103)==0
quietly gen double v125 = (_n+125)/7
quietly replace v125 = . if mod(_n,101)==0
quietly replace v125 = .z if mod(_n,103)==0
quietly gen double v126 = (_n+126)/7
quietly replace v126 = . if mod(_n,101)==0
quietly replace v126 = .z if mod(_n,103)==0
quietly gen double v127 = (_n+127)/7
quietly replace v127 = . if mod(_n,101)==0
quietly replace v127 = .z if mod(_n,103)==0
quietly gen double v128 = (_n+128)/7
quietly replace v128 = . if mod(_n,101)==0
quietly replace v128 = .z if mod(_n,103)==0
quietly gen double v129 = (_n+129)/7
quietly replace v129 = . if mod(_n,101)==0
quietly replace v129 = .z if mod(_n,103)==0
quietly gen double v130 = (_n+130)/7
quietly replace v130 = . if mod(_n,101)==0
quietly replace v130 = .z if mod(_n,103)==0
quietly gen double v131 = (_n+131)/7
quietly replace v131 = . if mod(_n,101)==0
quietly replace v131 = .z if mod(_n,103)==0
quietly gen double v132 = (_n+132)/7
quietly replace v132 = . if mod(_n,101)==0
quietly replace v132 = .z if mod(_n,103)==0
quietly gen double v133 = (_n+133)/7
quietly replace v133 = . if mod(_n,101)==0
quietly replace v133 = .z if mod(_n,103)==0
quietly gen double v134 = (_n+134)/7
quietly replace v134 = . if mod(_n,101)==0
quietly replace v134 = .z if mod(_n,103)==0
quietly gen double v135 = (_n+135)/7
quietly replace v135 = . if mod(_n,101)==0
quietly replace v135 = .z if mod(_n,103)==0
quietly gen double v136 = (_n+136)/7
quietly replace v136 = . if mod(_n,101)==0
quietly replace v136 = .z if mod(_n,103)==0
quietly gen double v137 = (_n+137)/7
quietly replace v137 = . if mod(_n,101)==0
quietly replace v137 = .z if mod(_n,103)==0
quietly gen double v138 = (_n+138)/7
quietly replace v138 = . if mod(_n,101)==0
quietly replace v138 = .z if mod(_n,103)==0
quietly gen double v139 = (_n+139)/7
quietly replace v139 = . if mod(_n,101)==0
quietly replace v139 = .z if mod(_n,103)==0
quietly gen double v140 = (_n+140)/7
quietly replace v140 = . if mod(_n,101)==0
quietly replace v140 = .z if mod(_n,103)==0
quietly gen double v141 = (_n+141)/7
quietly replace v141 = . if mod(_n,101)==0
quietly replace v141 = .z if mod(_n,103)==0
quietly gen double v142 = (_n+142)/7
quietly replace v142 = . if mod(_n,101)==0
quietly replace v142 = .z if mod(_n,103)==0
quietly gen double v143 = (_n+143)/7
quietly replace v143 = . if mod(_n,101)==0
quietly replace v143 = .z if mod(_n,103)==0
quietly gen double v144 = (_n+144)/7
quietly replace v144 = . if mod(_n,101)==0
quietly replace v144 = .z if mod(_n,103)==0
quietly gen double v145 = (_n+145)/7
quietly replace v145 = . if mod(_n,101)==0
quietly replace v145 = .z if mod(_n,103)==0
quietly gen double v146 = (_n+146)/7
quietly replace v146 = . if mod(_n,101)==0
quietly replace v146 = .z if mod(_n,103)==0
quietly gen double v147 = (_n+147)/7
quietly replace v147 = . if mod(_n,101)==0
quietly replace v147 = .z if mod(_n,103)==0
quietly gen double v148 = (_n+148)/7
quietly replace v148 = . if mod(_n,101)==0
quietly replace v148 = .z if mod(_n,103)==0
quietly gen double v149 = (_n+149)/7
quietly replace v149 = . if mod(_n,101)==0
quietly replace v149 = .z if mod(_n,103)==0
quietly gen double v150 = (_n+150)/7
quietly replace v150 = . if mod(_n,101)==0
quietly replace v150 = .z if mod(_n,103)==0
quietly gen double v151 = (_n+151)/7
quietly replace v151 = . if mod(_n,101)==0
quietly replace v151 = .z if mod(_n,103)==0
quietly gen double v152 = (_n+152)/7
quietly replace v152 = . if mod(_n,101)==0
quietly replace v152 = .z if mod(_n,103)==0
quietly gen double v153 = (_n+153)/7
quietly replace v153 = . if mod(_n,101)==0
quietly replace v153 = .z if mod(_n,103)==0
quietly gen double v154 = (_n+154)/7
quietly replace v154 = . if mod(_n,101)==0
quietly replace v154 = .z if mod(_n,103)==0
quietly gen double v155 = (_n+155)/7
quietly replace v155 = . if mod(_n,101)==0
quietly replace v155 = .z if mod(_n,103)==0
quietly gen double v156 = (_n+156)/7
quietly replace v156 = . if mod(_n,101)==0
quietly replace v156 = .z if mod(_n,103)==0
quietly gen double v157 = (_n+157)/7
quietly replace v157 = . if mod(_n,101)==0
quietly replace v157 = .z if mod(_n,103)==0
quietly gen double v158 = (_n+158)/7
quietly replace v158 = . if mod(_n,101)==0
quietly replace v158 = .z if mod(_n,103)==0
quietly gen double v159 = (_n+159)/7
quietly replace v159 = . if mod(_n,101)==0
quietly replace v159 = .z if mod(_n,103)==0
quietly gen double v160 = (_n+160)/7
quietly replace v160 = . if mod(_n,101)==0
quietly replace v160 = .z if mod(_n,103)==0
quietly gen double v161 = (_n+161)/7
quietly replace v161 = . if mod(_n,101)==0
quietly replace v161 = .z if mod(_n,103)==0
quietly gen double v162 = (_n+162)/7
quietly replace v162 = . if mod(_n,101)==0
quietly replace v162 = .z if mod(_n,103)==0
quietly gen double v163 = (_n+163)/7
quietly replace v163 = . if mod(_n,101)==0
quietly replace v163 = .z if mod(_n,103)==0
quietly gen double v164 = (_n+164)/7
quietly replace v164 = . if mod(_n,101)==0
quietly replace v164 = .z if mod(_n,103)==0
quietly gen double v165 = (_n+165)/7
quietly replace v165 = . if mod(_n,101)==0
quietly replace v165 = .z if mod(_n,103)==0
quietly gen double v166 = (_n+166)/7
quietly replace v166 = . if mod(_n,101)==0
quietly replace v166 = .z if mod(_n,103)==0
quietly gen double v167 = (_n+167)/7
quietly replace v167 = . if mod(_n,101)==0
quietly replace v167 = .z if mod(_n,103)==0
quietly gen double v168 = (_n+168)/7
quietly replace v168 = . if mod(_n,101)==0
quietly replace v168 = .z if mod(_n,103)==0
quietly gen double v169 = (_n+169)/7
quietly replace v169 = . if mod(_n,101)==0
quietly replace v169 = .z if mod(_n,103)==0
quietly gen double v170 = (_n+170)/7
quietly replace v170 = . if mod(_n,101)==0
quietly replace v170 = .z if mod(_n,103)==0
quietly gen double v171 = (_n+171)/7
quietly replace v171 = . if mod(_n,101)==0
quietly replace v171 = .z if mod(_n,103)==0
quietly gen double v172 = (_n+172)/7
quietly replace v172 = . if mod(_n,101)==0
quietly replace v172 = .z if mod(_n,103)==0
quietly gen double v173 = (_n+173)/7
quietly replace v173 = . if mod(_n,101)==0
quietly replace v173 = .z if mod(_n,103)==0
quietly gen double v174 = (_n+174)/7
quietly replace v174 = . if mod(_n,101)==0
quietly replace v174 = .z if mod(_n,103)==0
quietly gen double v175 = (_n+175)/7
quietly replace v175 = . if mod(_n,101)==0
quietly replace v175 = .z if mod(_n,103)==0
quietly gen double v176 = (_n+176)/7
quietly replace v176 = . if mod(_n,101)==0
quietly replace v176 = .z if mod(_n,103)==0
quietly gen double v177 = (_n+177)/7
quietly replace v177 = . if mod(_n,101)==0
quietly replace v177 = .z if mod(_n,103)==0
quietly gen double v178 = (_n+178)/7
quietly replace v178 = . if mod(_n,101)==0
quietly replace v178 = .z if mod(_n,103)==0
quietly gen double v179 = (_n+179)/7
quietly replace v179 = . if mod(_n,101)==0
quietly replace v179 = .z if mod(_n,103)==0
quietly gen double v180 = (_n+180)/7
quietly replace v180 = . if mod(_n,101)==0
quietly replace v180 = .z if mod(_n,103)==0
quietly gen double v181 = (_n+181)/7
quietly replace v181 = . if mod(_n,101)==0
quietly replace v181 = .z if mod(_n,103)==0
quietly gen double v182 = (_n+182)/7
quietly replace v182 = . if mod(_n,101)==0
quietly replace v182 = .z if mod(_n,103)==0
quietly gen double v183 = (_n+183)/7
quietly replace v183 = . if mod(_n,101)==0
quietly replace v183 = .z if mod(_n,103)==0
quietly gen double v184 = (_n+184)/7
quietly replace v184 = . if mod(_n,101)==0
quietly replace v184 = .z if mod(_n,103)==0
quietly gen double v185 = (_n+185)/7
quietly replace v185 = . if mod(_n,101)==0
quietly replace v185 = .z if mod(_n,103)==0
quietly gen double v186 = (_n+186)/7
quietly replace v186 = . if mod(_n,101)==0
quietly replace v186 = .z if mod(_n,103)==0
quietly gen double v187 = (_n+187)/7
quietly replace v187 = . if mod(_n,101)==0
quietly replace v187 = .z if mod(_n,103)==0
quietly gen double v188 = (_n+188)/7
quietly replace v188 = . if mod(_n,101)==0
quietly replace v188 = .z if mod(_n,103)==0
quietly gen double v189 = (_n+189)/7
quietly replace v189 = . if mod(_n,101)==0
quietly replace v189 = .z if mod(_n,103)==0
quietly gen double v190 = (_n+190)/7
quietly replace v190 = . if mod(_n,101)==0
quietly replace v190 = .z if mod(_n,103)==0
quietly gen double v191 = (_n+191)/7
quietly replace v191 = . if mod(_n,101)==0
quietly replace v191 = .z if mod(_n,103)==0
quietly gen double v192 = (_n+192)/7
quietly replace v192 = . if mod(_n,101)==0
quietly replace v192 = .z if mod(_n,103)==0
quietly gen double v193 = (_n+193)/7
quietly replace v193 = . if mod(_n,101)==0
quietly replace v193 = .z if mod(_n,103)==0
quietly gen double v194 = (_n+194)/7
quietly replace v194 = . if mod(_n,101)==0
quietly replace v194 = .z if mod(_n,103)==0
quietly gen double v195 = (_n+195)/7
quietly replace v195 = . if mod(_n,101)==0
quietly replace v195 = .z if mod(_n,103)==0
quietly gen double v196 = (_n+196)/7
quietly replace v196 = . if mod(_n,101)==0
quietly replace v196 = .z if mod(_n,103)==0
quietly gen double v197 = (_n+197)/7
quietly replace v197 = . if mod(_n,101)==0
quietly replace v197 = .z if mod(_n,103)==0
quietly gen double v198 = (_n+198)/7
quietly replace v198 = . if mod(_n,101)==0
quietly replace v198 = .z if mod(_n,103)==0
quietly gen double v199 = (_n+199)/7
quietly replace v199 = . if mod(_n,101)==0
quietly replace v199 = .z if mod(_n,103)==0
quietly gen double v200 = (_n+200)/7
quietly replace v200 = . if mod(_n,101)==0
quietly replace v200 = .z if mod(_n,103)==0
quietly gen double v201 = (_n+201)/7
quietly replace v201 = . if mod(_n,101)==0
quietly replace v201 = .z if mod(_n,103)==0
quietly gen double v202 = (_n+202)/7
quietly replace v202 = . if mod(_n,101)==0
quietly replace v202 = .z if mod(_n,103)==0
quietly gen double v203 = (_n+203)/7
quietly replace v203 = . if mod(_n,101)==0
quietly replace v203 = .z if mod(_n,103)==0
quietly gen double v204 = (_n+204)/7
quietly replace v204 = . if mod(_n,101)==0
quietly replace v204 = .z if mod(_n,103)==0
quietly gen double v205 = (_n+205)/7
quietly replace v205 = . if mod(_n,101)==0
quietly replace v205 = .z if mod(_n,103)==0
quietly gen double v206 = (_n+206)/7
quietly replace v206 = . if mod(_n,101)==0
quietly replace v206 = .z if mod(_n,103)==0
quietly gen double v207 = (_n+207)/7
quietly replace v207 = . if mod(_n,101)==0
quietly replace v207 = .z if mod(_n,103)==0
quietly gen double v208 = (_n+208)/7
quietly replace v208 = . if mod(_n,101)==0
quietly replace v208 = .z if mod(_n,103)==0
quietly gen double v209 = (_n+209)/7
quietly replace v209 = . if mod(_n,101)==0
quietly replace v209 = .z if mod(_n,103)==0
quietly gen double v210 = (_n+210)/7
quietly replace v210 = . if mod(_n,101)==0
quietly replace v210 = .z if mod(_n,103)==0
quietly gen double v211 = (_n+211)/7
quietly replace v211 = . if mod(_n,101)==0
quietly replace v211 = .z if mod(_n,103)==0
quietly gen double v212 = (_n+212)/7
quietly replace v212 = . if mod(_n,101)==0
quietly replace v212 = .z if mod(_n,103)==0
quietly gen double v213 = (_n+213)/7
quietly replace v213 = . if mod(_n,101)==0
quietly replace v213 = .z if mod(_n,103)==0
quietly gen double v214 = (_n+214)/7
quietly replace v214 = . if mod(_n,101)==0
quietly replace v214 = .z if mod(_n,103)==0
quietly gen double v215 = (_n+215)/7
quietly replace v215 = . if mod(_n,101)==0
quietly replace v215 = .z if mod(_n,103)==0
quietly gen double v216 = (_n+216)/7
quietly replace v216 = . if mod(_n,101)==0
quietly replace v216 = .z if mod(_n,103)==0
quietly gen double v217 = (_n+217)/7
quietly replace v217 = . if mod(_n,101)==0
quietly replace v217 = .z if mod(_n,103)==0
quietly gen double v218 = (_n+218)/7
quietly replace v218 = . if mod(_n,101)==0
quietly replace v218 = .z if mod(_n,103)==0
quietly gen double v219 = (_n+219)/7
quietly replace v219 = . if mod(_n,101)==0
quietly replace v219 = .z if mod(_n,103)==0
quietly gen double v220 = (_n+220)/7
quietly replace v220 = . if mod(_n,101)==0
quietly replace v220 = .z if mod(_n,103)==0
quietly gen double v221 = (_n+221)/7
quietly replace v221 = . if mod(_n,101)==0
quietly replace v221 = .z if mod(_n,103)==0
quietly gen double v222 = (_n+222)/7
quietly replace v222 = . if mod(_n,101)==0
quietly replace v222 = .z if mod(_n,103)==0
quietly gen double v223 = (_n+223)/7
quietly replace v223 = . if mod(_n,101)==0
quietly replace v223 = .z if mod(_n,103)==0
quietly gen double v224 = (_n+224)/7
quietly replace v224 = . if mod(_n,101)==0
quietly replace v224 = .z if mod(_n,103)==0
quietly gen double v225 = (_n+225)/7
quietly replace v225 = . if mod(_n,101)==0
quietly replace v225 = .z if mod(_n,103)==0
quietly gen double v226 = (_n+226)/7
quietly replace v226 = . if mod(_n,101)==0
quietly replace v226 = .z if mod(_n,103)==0
quietly gen double v227 = (_n+227)/7
quietly replace v227 = . if mod(_n,101)==0
quietly replace v227 = .z if mod(_n,103)==0
quietly gen double v228 = (_n+228)/7
quietly replace v228 = . if mod(_n,101)==0
quietly replace v228 = .z if mod(_n,103)==0
quietly gen double v229 = (_n+229)/7
quietly replace v229 = . if mod(_n,101)==0
quietly replace v229 = .z if mod(_n,103)==0
quietly gen double v230 = (_n+230)/7
quietly replace v230 = . if mod(_n,101)==0
quietly replace v230 = .z if mod(_n,103)==0
quietly gen double v231 = (_n+231)/7
quietly replace v231 = . if mod(_n,101)==0
quietly replace v231 = .z if mod(_n,103)==0
quietly gen double v232 = (_n+232)/7
quietly replace v232 = . if mod(_n,101)==0
quietly replace v232 = .z if mod(_n,103)==0
quietly gen double v233 = (_n+233)/7
quietly replace v233 = . if mod(_n,101)==0
quietly replace v233 = .z if mod(_n,103)==0
quietly gen double v234 = (_n+234)/7
quietly replace v234 = . if mod(_n,101)==0
quietly replace v234 = .z if mod(_n,103)==0
quietly gen double v235 = (_n+235)/7
quietly replace v235 = . if mod(_n,101)==0
quietly replace v235 = .z if mod(_n,103)==0
quietly gen double v236 = (_n+236)/7
quietly replace v236 = . if mod(_n,101)==0
quietly replace v236 = .z if mod(_n,103)==0
quietly gen double v237 = (_n+237)/7
quietly replace v237 = . if mod(_n,101)==0
quietly replace v237 = .z if mod(_n,103)==0
quietly gen double v238 = (_n+238)/7
quietly replace v238 = . if mod(_n,101)==0
quietly replace v238 = .z if mod(_n,103)==0
quietly gen double v239 = (_n+239)/7
quietly replace v239 = . if mod(_n,101)==0
quietly replace v239 = .z if mod(_n,103)==0
quietly gen double v240 = (_n+240)/7
quietly replace v240 = . if mod(_n,101)==0
quietly replace v240 = .z if mod(_n,103)==0
quietly gen double v241 = (_n+241)/7
quietly replace v241 = . if mod(_n,101)==0
quietly replace v241 = .z if mod(_n,103)==0
quietly gen double v242 = (_n+242)/7
quietly replace v242 = . if mod(_n,101)==0
quietly replace v242 = .z if mod(_n,103)==0
quietly gen double v243 = (_n+243)/7
quietly replace v243 = . if mod(_n,101)==0
quietly replace v243 = .z if mod(_n,103)==0
quietly gen double v244 = (_n+244)/7
quietly replace v244 = . if mod(_n,101)==0
quietly replace v244 = .z if mod(_n,103)==0
quietly gen double v245 = (_n+245)/7
quietly replace v245 = . if mod(_n,101)==0
quietly replace v245 = .z if mod(_n,103)==0
quietly gen double v246 = (_n+246)/7
quietly replace v246 = . if mod(_n,101)==0
quietly replace v246 = .z if mod(_n,103)==0
quietly gen double v247 = (_n+247)/7
quietly replace v247 = . if mod(_n,101)==0
quietly replace v247 = .z if mod(_n,103)==0
quietly gen double v248 = (_n+248)/7
quietly replace v248 = . if mod(_n,101)==0
quietly replace v248 = .z if mod(_n,103)==0
quietly gen double v249 = (_n+249)/7
quietly replace v249 = . if mod(_n,101)==0
quietly replace v249 = .z if mod(_n,103)==0
quietly gen double v250 = (_n+250)/7
quietly replace v250 = . if mod(_n,101)==0
quietly replace v250 = .z if mod(_n,103)==0
quietly gen double v251 = (_n+251)/7
quietly replace v251 = . if mod(_n,101)==0
quietly replace v251 = .z if mod(_n,103)==0
quietly gen double v252 = (_n+252)/7
quietly replace v252 = . if mod(_n,101)==0
quietly replace v252 = .z if mod(_n,103)==0
quietly gen double v253 = (_n+253)/7
quietly replace v253 = . if mod(_n,101)==0
quietly replace v253 = .z if mod(_n,103)==0
quietly gen double v254 = (_n+254)/7
quietly replace v254 = . if mod(_n,101)==0
quietly replace v254 = .z if mod(_n,103)==0
quietly gen double v255 = (_n+255)/7
quietly replace v255 = . if mod(_n,101)==0
quietly replace v255 = .z if mod(_n,103)==0
quietly gen double v256 = (_n+256)/7
quietly replace v256 = . if mod(_n,101)==0
quietly replace v256 = .z if mod(_n,103)==0
quietly gen double v257 = (_n+257)/7
quietly replace v257 = . if mod(_n,101)==0
quietly replace v257 = .z if mod(_n,103)==0
quietly gen double v258 = (_n+258)/7
quietly replace v258 = . if mod(_n,101)==0
quietly replace v258 = .z if mod(_n,103)==0
quietly gen double v259 = (_n+259)/7
quietly replace v259 = . if mod(_n,101)==0
quietly replace v259 = .z if mod(_n,103)==0
quietly gen double v260 = (_n+260)/7
quietly replace v260 = . if mod(_n,101)==0
quietly replace v260 = .z if mod(_n,103)==0
quietly gen double v261 = (_n+261)/7
quietly replace v261 = . if mod(_n,101)==0
quietly replace v261 = .z if mod(_n,103)==0
quietly gen double v262 = (_n+262)/7
quietly replace v262 = . if mod(_n,101)==0
quietly replace v262 = .z if mod(_n,103)==0
quietly gen double v263 = (_n+263)/7
quietly replace v263 = . if mod(_n,101)==0
quietly replace v263 = .z if mod(_n,103)==0
quietly gen double v264 = (_n+264)/7
quietly replace v264 = . if mod(_n,101)==0
quietly replace v264 = .z if mod(_n,103)==0
quietly gen double v265 = (_n+265)/7
quietly replace v265 = . if mod(_n,101)==0
quietly replace v265 = .z if mod(_n,103)==0
quietly gen double v266 = (_n+266)/7
quietly replace v266 = . if mod(_n,101)==0
quietly replace v266 = .z if mod(_n,103)==0
quietly gen double v267 = (_n+267)/7
quietly replace v267 = . if mod(_n,101)==0
quietly replace v267 = .z if mod(_n,103)==0
quietly gen double v268 = (_n+268)/7
quietly replace v268 = . if mod(_n,101)==0
quietly replace v268 = .z if mod(_n,103)==0
quietly gen double v269 = (_n+269)/7
quietly replace v269 = . if mod(_n,101)==0
quietly replace v269 = .z if mod(_n,103)==0
quietly gen double v270 = (_n+270)/7
quietly replace v270 = . if mod(_n,101)==0
quietly replace v270 = .z if mod(_n,103)==0
quietly gen double v271 = (_n+271)/7
quietly replace v271 = . if mod(_n,101)==0
quietly replace v271 = .z if mod(_n,103)==0
quietly gen double v272 = (_n+272)/7
quietly replace v272 = . if mod(_n,101)==0
quietly replace v272 = .z if mod(_n,103)==0
quietly gen double v273 = (_n+273)/7
quietly replace v273 = . if mod(_n,101)==0
quietly replace v273 = .z if mod(_n,103)==0
quietly gen double v274 = (_n+274)/7
quietly replace v274 = . if mod(_n,101)==0
quietly replace v274 = .z if mod(_n,103)==0
quietly gen double v275 = (_n+275)/7
quietly replace v275 = . if mod(_n,101)==0
quietly replace v275 = .z if mod(_n,103)==0
quietly gen double v276 = (_n+276)/7
quietly replace v276 = . if mod(_n,101)==0
quietly replace v276 = .z if mod(_n,103)==0
quietly gen double v277 = (_n+277)/7
quietly replace v277 = . if mod(_n,101)==0
quietly replace v277 = .z if mod(_n,103)==0
quietly gen double v278 = (_n+278)/7
quietly replace v278 = . if mod(_n,101)==0
quietly replace v278 = .z if mod(_n,103)==0
quietly gen double v279 = (_n+279)/7
quietly replace v279 = . if mod(_n,101)==0
quietly replace v279 = .z if mod(_n,103)==0
quietly gen double v280 = (_n+280)/7
quietly replace v280 = . if mod(_n,101)==0
quietly replace v280 = .z if mod(_n,103)==0
quietly gen double v281 = (_n+281)/7
quietly replace v281 = . if mod(_n,101)==0
quietly replace v281 = .z if mod(_n,103)==0
quietly gen double v282 = (_n+282)/7
quietly replace v282 = . if mod(_n,101)==0
quietly replace v282 = .z if mod(_n,103)==0
quietly gen double v283 = (_n+283)/7
quietly replace v283 = . if mod(_n,101)==0
quietly replace v283 = .z if mod(_n,103)==0
quietly gen double v284 = (_n+284)/7
quietly replace v284 = . if mod(_n,101)==0
quietly replace v284 = .z if mod(_n,103)==0
quietly gen double v285 = (_n+285)/7
quietly replace v285 = . if mod(_n,101)==0
quietly replace v285 = .z if mod(_n,103)==0
quietly gen double v286 = (_n+286)/7
quietly replace v286 = . if mod(_n,101)==0
quietly replace v286 = .z if mod(_n,103)==0
quietly gen double v287 = (_n+287)/7
quietly replace v287 = . if mod(_n,101)==0
quietly replace v287 = .z if mod(_n,103)==0
quietly gen double v288 = (_n+288)/7
quietly replace v288 = . if mod(_n,101)==0
quietly replace v288 = .z if mod(_n,103)==0
quietly gen double v289 = (_n+289)/7
quietly replace v289 = . if mod(_n,101)==0
quietly replace v289 = .z if mod(_n,103)==0
quietly gen double v290 = (_n+290)/7
quietly replace v290 = . if mod(_n,101)==0
quietly replace v290 = .z if mod(_n,103)==0
quietly gen double v291 = (_n+291)/7
quietly replace v291 = . if mod(_n,101)==0
quietly replace v291 = .z if mod(_n,103)==0
quietly gen double v292 = (_n+292)/7
quietly replace v292 = . if mod(_n,101)==0
quietly replace v292 = .z if mod(_n,103)==0
quietly gen double v293 = (_n+293)/7
quietly replace v293 = . if mod(_n,101)==0
quietly replace v293 = .z if mod(_n,103)==0
quietly gen double v294 = (_n+294)/7
quietly replace v294 = . if mod(_n,101)==0
quietly replace v294 = .z if mod(_n,103)==0
quietly gen double v295 = (_n+295)/7
quietly replace v295 = . if mod(_n,101)==0
quietly replace v295 = .z if mod(_n,103)==0
quietly gen double v296 = (_n+296)/7
quietly replace v296 = . if mod(_n,101)==0
quietly replace v296 = .z if mod(_n,103)==0
quietly gen double v297 = (_n+297)/7
quietly replace v297 = . if mod(_n,101)==0
quietly replace v297 = .z if mod(_n,103)==0
quietly gen double v298 = (_n+298)/7
quietly replace v298 = . if mod(_n,101)==0
quietly replace v298 = .z if mod(_n,103)==0
quietly gen double v299 = (_n+299)/7
quietly replace v299 = . if mod(_n,101)==0
quietly replace v299 = .z if mod(_n,103)==0
quietly gen double v300 = (_n+300)/7
quietly replace v300 = . if mod(_n,101)==0
quietly replace v300 = .z if mod(_n,103)==0
quietly gen double v301 = (_n+301)/7
quietly replace v301 = . if mod(_n,101)==0
quietly replace v301 = .z if mod(_n,103)==0
quietly gen double v302 = (_n+302)/7
quietly replace v302 = . if mod(_n,101)==0
quietly replace v302 = .z if mod(_n,103)==0
quietly gen double v303 = (_n+303)/7
quietly replace v303 = . if mod(_n,101)==0
quietly replace v303 = .z if mod(_n,103)==0
quietly gen double v304 = (_n+304)/7
quietly replace v304 = . if mod(_n,101)==0
quietly replace v304 = .z if mod(_n,103)==0
quietly gen double v305 = (_n+305)/7
quietly replace v305 = . if mod(_n,101)==0
quietly replace v305 = .z if mod(_n,103)==0
quietly gen double v306 = (_n+306)/7
quietly replace v306 = . if mod(_n,101)==0
quietly replace v306 = .z if mod(_n,103)==0
quietly gen double v307 = (_n+307)/7
quietly replace v307 = . if mod(_n,101)==0
quietly replace v307 = .z if mod(_n,103)==0
quietly gen double v308 = (_n+308)/7
quietly replace v308 = . if mod(_n,101)==0
quietly replace v308 = .z if mod(_n,103)==0
quietly gen double v309 = (_n+309)/7
quietly replace v309 = . if mod(_n,101)==0
quietly replace v309 = .z if mod(_n,103)==0
quietly gen double v310 = (_n+310)/7
quietly replace v310 = . if mod(_n,101)==0
quietly replace v310 = .z if mod(_n,103)==0
quietly gen double v311 = (_n+311)/7
quietly replace v311 = . if mod(_n,101)==0
quietly replace v311 = .z if mod(_n,103)==0
quietly gen double v312 = (_n+312)/7
quietly replace v312 = . if mod(_n,101)==0
quietly replace v312 = .z if mod(_n,103)==0
quietly gen double v313 = (_n+313)/7
quietly replace v313 = . if mod(_n,101)==0
quietly replace v313 = .z if mod(_n,103)==0
quietly gen double v314 = (_n+314)/7
quietly replace v314 = . if mod(_n,101)==0
quietly replace v314 = .z if mod(_n,103)==0
quietly gen double v315 = (_n+315)/7
quietly replace v315 = . if mod(_n,101)==0
quietly replace v315 = .z if mod(_n,103)==0
quietly gen double v316 = (_n+316)/7
quietly replace v316 = . if mod(_n,101)==0
quietly replace v316 = .z if mod(_n,103)==0
quietly gen double v317 = (_n+317)/7
quietly replace v317 = . if mod(_n,101)==0
quietly replace v317 = .z if mod(_n,103)==0
quietly gen double v318 = (_n+318)/7
quietly replace v318 = . if mod(_n,101)==0
quietly replace v318 = .z if mod(_n,103)==0
quietly gen double v319 = (_n+319)/7
quietly replace v319 = . if mod(_n,101)==0
quietly replace v319 = .z if mod(_n,103)==0
quietly gen double v320 = (_n+320)/7
quietly replace v320 = . if mod(_n,101)==0
quietly replace v320 = .z if mod(_n,103)==0
quietly gen double v321 = (_n+321)/7
quietly replace v321 = . if mod(_n,101)==0
quietly replace v321 = .z if mod(_n,103)==0
quietly gen double v322 = (_n+322)/7
quietly replace v322 = . if mod(_n,101)==0
quietly replace v322 = .z if mod(_n,103)==0
quietly gen double v323 = (_n+323)/7
quietly replace v323 = . if mod(_n,101)==0
quietly replace v323 = .z if mod(_n,103)==0
quietly gen double v324 = (_n+324)/7
quietly replace v324 = . if mod(_n,101)==0
quietly replace v324 = .z if mod(_n,103)==0
quietly gen double v325 = (_n+325)/7
quietly replace v325 = . if mod(_n,101)==0
quietly replace v325 = .z if mod(_n,103)==0
quietly gen double v326 = (_n+326)/7
quietly replace v326 = . if mod(_n,101)==0
quietly replace v326 = .z if mod(_n,103)==0
quietly gen double v327 = (_n+327)/7
quietly replace v327 = . if mod(_n,101)==0
quietly replace v327 = .z if mod(_n,103)==0
quietly gen double v328 = (_n+328)/7
quietly replace v328 = . if mod(_n,101)==0
quietly replace v328 = .z if mod(_n,103)==0
quietly gen double v329 = (_n+329)/7
quietly replace v329 = . if mod(_n,101)==0
quietly replace v329 = .z if mod(_n,103)==0
quietly gen double v330 = (_n+330)/7
quietly replace v330 = . if mod(_n,101)==0
quietly replace v330 = .z if mod(_n,103)==0
quietly gen double v331 = (_n+331)/7
quietly replace v331 = . if mod(_n,101)==0
quietly replace v331 = .z if mod(_n,103)==0
quietly gen double v332 = (_n+332)/7
quietly replace v332 = . if mod(_n,101)==0
quietly replace v332 = .z if mod(_n,103)==0
quietly gen double v333 = (_n+333)/7
quietly replace v333 = . if mod(_n,101)==0
quietly replace v333 = .z if mod(_n,103)==0
quietly gen double v334 = (_n+334)/7
quietly replace v334 = . if mod(_n,101)==0
quietly replace v334 = .z if mod(_n,103)==0
quietly gen double v335 = (_n+335)/7
quietly replace v335 = . if mod(_n,101)==0
quietly replace v335 = .z if mod(_n,103)==0
quietly gen double v336 = (_n+336)/7
quietly replace v336 = . if mod(_n,101)==0
quietly replace v336 = .z if mod(_n,103)==0
quietly gen double v337 = (_n+337)/7
quietly replace v337 = . if mod(_n,101)==0
quietly replace v337 = .z if mod(_n,103)==0
quietly gen double v338 = (_n+338)/7
quietly replace v338 = . if mod(_n,101)==0
quietly replace v338 = .z if mod(_n,103)==0
quietly gen double v339 = (_n+339)/7
quietly replace v339 = . if mod(_n,101)==0
quietly replace v339 = .z if mod(_n,103)==0
quietly gen double v340 = (_n+340)/7
quietly replace v340 = . if mod(_n,101)==0
quietly replace v340 = .z if mod(_n,103)==0
quietly gen double v341 = (_n+341)/7
quietly replace v341 = . if mod(_n,101)==0
quietly replace v341 = .z if mod(_n,103)==0
quietly gen double v342 = (_n+342)/7
quietly replace v342 = . if mod(_n,101)==0
quietly replace v342 = .z if mod(_n,103)==0
quietly gen double v343 = (_n+343)/7
quietly replace v343 = . if mod(_n,101)==0
quietly replace v343 = .z if mod(_n,103)==0
quietly gen double v344 = (_n+344)/7
quietly replace v344 = . if mod(_n,101)==0
quietly replace v344 = .z if mod(_n,103)==0
quietly gen double v345 = (_n+345)/7
quietly replace v345 = . if mod(_n,101)==0
quietly replace v345 = .z if mod(_n,103)==0
quietly gen double v346 = (_n+346)/7
quietly replace v346 = . if mod(_n,101)==0
quietly replace v346 = .z if mod(_n,103)==0
quietly gen double v347 = (_n+347)/7
quietly replace v347 = . if mod(_n,101)==0
quietly replace v347 = .z if mod(_n,103)==0
quietly gen double v348 = (_n+348)/7
quietly replace v348 = . if mod(_n,101)==0
quietly replace v348 = .z if mod(_n,103)==0
quietly gen double v349 = (_n+349)/7
quietly replace v349 = . if mod(_n,101)==0
quietly replace v349 = .z if mod(_n,103)==0
quietly gen double v350 = (_n+350)/7
quietly replace v350 = . if mod(_n,101)==0
quietly replace v350 = .z if mod(_n,103)==0
quietly gen double v351 = (_n+351)/7
quietly replace v351 = . if mod(_n,101)==0
quietly replace v351 = .z if mod(_n,103)==0
quietly gen double v352 = (_n+352)/7
quietly replace v352 = . if mod(_n,101)==0
quietly replace v352 = .z if mod(_n,103)==0
quietly gen double v353 = (_n+353)/7
quietly replace v353 = . if mod(_n,101)==0
quietly replace v353 = .z if mod(_n,103)==0
quietly gen double v354 = (_n+354)/7
quietly replace v354 = . if mod(_n,101)==0
quietly replace v354 = .z if mod(_n,103)==0
quietly gen double v355 = (_n+355)/7
quietly replace v355 = . if mod(_n,101)==0
quietly replace v355 = .z if mod(_n,103)==0
quietly gen double v356 = (_n+356)/7
quietly replace v356 = . if mod(_n,101)==0
quietly replace v356 = .z if mod(_n,103)==0
quietly gen double v357 = (_n+357)/7
quietly replace v357 = . if mod(_n,101)==0
quietly replace v357 = .z if mod(_n,103)==0
quietly gen double v358 = (_n+358)/7
quietly replace v358 = . if mod(_n,101)==0
quietly replace v358 = .z if mod(_n,103)==0
quietly gen double v359 = (_n+359)/7
quietly replace v359 = . if mod(_n,101)==0
quietly replace v359 = .z if mod(_n,103)==0
quietly gen double v360 = (_n+360)/7
quietly replace v360 = . if mod(_n,101)==0
quietly replace v360 = .z if mod(_n,103)==0
quietly gen double v361 = (_n+361)/7
quietly replace v361 = . if mod(_n,101)==0
quietly replace v361 = .z if mod(_n,103)==0
quietly gen double v362 = (_n+362)/7
quietly replace v362 = . if mod(_n,101)==0
quietly replace v362 = .z if mod(_n,103)==0
quietly gen double v363 = (_n+363)/7
quietly replace v363 = . if mod(_n,101)==0
quietly replace v363 = .z if mod(_n,103)==0
quietly gen double v364 = (_n+364)/7
quietly replace v364 = . if mod(_n,101)==0
quietly replace v364 = .z if mod(_n,103)==0
quietly gen double v365 = (_n+365)/7
quietly replace v365 = . if mod(_n,101)==0
quietly replace v365 = .z if mod(_n,103)==0
quietly gen double v366 = (_n+366)/7
quietly replace v366 = . if mod(_n,101)==0
quietly replace v366 = .z if mod(_n,103)==0
quietly gen double v367 = (_n+367)/7
quietly replace v367 = . if mod(_n,101)==0
quietly replace v367 = .z if mod(_n,103)==0
quietly gen double v368 = (_n+368)/7
quietly replace v368 = . if mod(_n,101)==0
quietly replace v368 = .z if mod(_n,103)==0
quietly gen double v369 = (_n+369)/7
quietly replace v369 = . if mod(_n,101)==0
quietly replace v369 = .z if mod(_n,103)==0
quietly gen double v370 = (_n+370)/7
quietly replace v370 = . if mod(_n,101)==0
quietly replace v370 = .z if mod(_n,103)==0
quietly gen double v371 = (_n+371)/7
quietly replace v371 = . if mod(_n,101)==0
quietly replace v371 = .z if mod(_n,103)==0
quietly gen double v372 = (_n+372)/7
quietly replace v372 = . if mod(_n,101)==0
quietly replace v372 = .z if mod(_n,103)==0
quietly gen double v373 = (_n+373)/7
quietly replace v373 = . if mod(_n,101)==0
quietly replace v373 = .z if mod(_n,103)==0
quietly gen double v374 = (_n+374)/7
quietly replace v374 = . if mod(_n,101)==0
quietly replace v374 = .z if mod(_n,103)==0
quietly gen double v375 = (_n+375)/7
quietly replace v375 = . if mod(_n,101)==0
quietly replace v375 = .z if mod(_n,103)==0
quietly gen double v376 = (_n+376)/7
quietly replace v376 = . if mod(_n,101)==0
quietly replace v376 = .z if mod(_n,103)==0
quietly gen double v377 = (_n+377)/7
quietly replace v377 = . if mod(_n,101)==0
quietly replace v377 = .z if mod(_n,103)==0
quietly gen double v378 = (_n+378)/7
quietly replace v378 = . if mod(_n,101)==0
quietly replace v378 = .z if mod(_n,103)==0
quietly gen double v379 = (_n+379)/7
quietly replace v379 = . if mod(_n,101)==0
quietly replace v379 = .z if mod(_n,103)==0
quietly gen double v380 = (_n+380)/7
quietly replace v380 = . if mod(_n,101)==0
quietly replace v380 = .z if mod(_n,103)==0
quietly gen double v381 = (_n+381)/7
quietly replace v381 = . if mod(_n,101)==0
quietly replace v381 = .z if mod(_n,103)==0
quietly gen double v382 = (_n+382)/7
quietly replace v382 = . if mod(_n,101)==0
quietly replace v382 = .z if mod(_n,103)==0
quietly gen double v383 = (_n+383)/7
quietly replace v383 = . if mod(_n,101)==0
quietly replace v383 = .z if mod(_n,103)==0
quietly gen double v384 = (_n+384)/7
quietly replace v384 = . if mod(_n,101)==0
quietly replace v384 = .z if mod(_n,103)==0
quietly gen double v385 = (_n+385)/7
quietly replace v385 = . if mod(_n,101)==0
quietly replace v385 = .z if mod(_n,103)==0
quietly gen double v386 = (_n+386)/7
quietly replace v386 = . if mod(_n,101)==0
quietly replace v386 = .z if mod(_n,103)==0
quietly gen double v387 = (_n+387)/7
quietly replace v387 = . if mod(_n,101)==0
quietly replace v387 = .z if mod(_n,103)==0
quietly gen double v388 = (_n+388)/7
quietly replace v388 = . if mod(_n,101)==0
quietly replace v388 = .z if mod(_n,103)==0
quietly gen double v389 = (_n+389)/7
quietly replace v389 = . if mod(_n,101)==0
quietly replace v389 = .z if mod(_n,103)==0
quietly gen double v390 = (_n+390)/7
quietly replace v390 = . if mod(_n,101)==0
quietly replace v390 = .z if mod(_n,103)==0
quietly gen double v391 = (_n+391)/7
quietly replace v391 = . if mod(_n,101)==0
quietly replace v391 = .z if mod(_n,103)==0
quietly gen double v392 = (_n+392)/7
quietly replace v392 = . if mod(_n,101)==0
quietly replace v392 = .z if mod(_n,103)==0
quietly gen double v393 = (_n+393)/7
quietly replace v393 = . if mod(_n,101)==0
quietly replace v393 = .z if mod(_n,103)==0
quietly gen double v394 = (_n+394)/7
quietly replace v394 = . if mod(_n,101)==0
quietly replace v394 = .z if mod(_n,103)==0
quietly gen double v395 = (_n+395)/7
quietly replace v395 = . if mod(_n,101)==0
quietly replace v395 = .z if mod(_n,103)==0
quietly gen double v396 = (_n+396)/7
quietly replace v396 = . if mod(_n,101)==0
quietly replace v396 = .z if mod(_n,103)==0
quietly gen double v397 = (_n+397)/7
quietly replace v397 = . if mod(_n,101)==0
quietly replace v397 = .z if mod(_n,103)==0
quietly gen double v398 = (_n+398)/7
quietly replace v398 = . if mod(_n,101)==0
quietly replace v398 = .z if mod(_n,103)==0
quietly gen double v399 = (_n+399)/7
quietly replace v399 = . if mod(_n,101)==0
quietly replace v399 = .z if mod(_n,103)==0
quietly gen double v400 = (_n+400)/7
quietly replace v400 = . if mod(_n,101)==0
quietly replace v400 = .z if mod(_n,103)==0
quietly gen double v401 = (_n+401)/7
quietly replace v401 = . if mod(_n,101)==0
quietly replace v401 = .z if mod(_n,103)==0
quietly gen double v402 = (_n+402)/7
quietly replace v402 = . if mod(_n,101)==0
quietly replace v402 = .z if mod(_n,103)==0
quietly gen double v403 = (_n+403)/7
quietly replace v403 = . if mod(_n,101)==0
quietly replace v403 = .z if mod(_n,103)==0
quietly gen double v404 = (_n+404)/7
quietly replace v404 = . if mod(_n,101)==0
quietly replace v404 = .z if mod(_n,103)==0
quietly gen double v405 = (_n+405)/7
quietly replace v405 = . if mod(_n,101)==0
quietly replace v405 = .z if mod(_n,103)==0
quietly gen double v406 = (_n+406)/7
quietly replace v406 = . if mod(_n,101)==0
quietly replace v406 = .z if mod(_n,103)==0
quietly gen double v407 = (_n+407)/7
quietly replace v407 = . if mod(_n,101)==0
quietly replace v407 = .z if mod(_n,103)==0
quietly gen double v408 = (_n+408)/7
quietly replace v408 = . if mod(_n,101)==0
quietly replace v408 = .z if mod(_n,103)==0
quietly gen double v409 = (_n+409)/7
quietly replace v409 = . if mod(_n,101)==0
quietly replace v409 = .z if mod(_n,103)==0
quietly gen double v410 = (_n+410)/7
quietly replace v410 = . if mod(_n,101)==0
quietly replace v410 = .z if mod(_n,103)==0
quietly gen double v411 = (_n+411)/7
quietly replace v411 = . if mod(_n,101)==0
quietly replace v411 = .z if mod(_n,103)==0
quietly gen double v412 = (_n+412)/7
quietly replace v412 = . if mod(_n,101)==0
quietly replace v412 = .z if mod(_n,103)==0
quietly gen double v413 = (_n+413)/7
quietly replace v413 = . if mod(_n,101)==0
quietly replace v413 = .z if mod(_n,103)==0
quietly gen double v414 = (_n+414)/7
quietly replace v414 = . if mod(_n,101)==0
quietly replace v414 = .z if mod(_n,103)==0
quietly gen double v415 = (_n+415)/7
quietly replace v415 = . if mod(_n,101)==0
quietly replace v415 = .z if mod(_n,103)==0
quietly gen double v416 = (_n+416)/7
quietly replace v416 = . if mod(_n,101)==0
quietly replace v416 = .z if mod(_n,103)==0
quietly gen double v417 = (_n+417)/7
quietly replace v417 = . if mod(_n,101)==0
quietly replace v417 = .z if mod(_n,103)==0
quietly gen double v418 = (_n+418)/7
quietly replace v418 = . if mod(_n,101)==0
quietly replace v418 = .z if mod(_n,103)==0
quietly gen double v419 = (_n+419)/7
quietly replace v419 = . if mod(_n,101)==0
quietly replace v419 = .z if mod(_n,103)==0
quietly gen double v420 = (_n+420)/7
quietly replace v420 = . if mod(_n,101)==0
quietly replace v420 = .z if mod(_n,103)==0
quietly gen double v421 = (_n+421)/7
quietly replace v421 = . if mod(_n,101)==0
quietly replace v421 = .z if mod(_n,103)==0
quietly gen double v422 = (_n+422)/7
quietly replace v422 = . if mod(_n,101)==0
quietly replace v422 = .z if mod(_n,103)==0
quietly gen double v423 = (_n+423)/7
quietly replace v423 = . if mod(_n,101)==0
quietly replace v423 = .z if mod(_n,103)==0
quietly gen double v424 = (_n+424)/7
quietly replace v424 = . if mod(_n,101)==0
quietly replace v424 = .z if mod(_n,103)==0
quietly gen double v425 = (_n+425)/7
quietly replace v425 = . if mod(_n,101)==0
quietly replace v425 = .z if mod(_n,103)==0
quietly gen double v426 = (_n+426)/7
quietly replace v426 = . if mod(_n,101)==0
quietly replace v426 = .z if mod(_n,103)==0
quietly gen double v427 = (_n+427)/7
quietly replace v427 = . if mod(_n,101)==0
quietly replace v427 = .z if mod(_n,103)==0
quietly gen double v428 = (_n+428)/7
quietly replace v428 = . if mod(_n,101)==0
quietly replace v428 = .z if mod(_n,103)==0
quietly gen double v429 = (_n+429)/7
quietly replace v429 = . if mod(_n,101)==0
quietly replace v429 = .z if mod(_n,103)==0
quietly gen double v430 = (_n+430)/7
quietly replace v430 = . if mod(_n,101)==0
quietly replace v430 = .z if mod(_n,103)==0
quietly gen double v431 = (_n+431)/7
quietly replace v431 = . if mod(_n,101)==0
quietly replace v431 = .z if mod(_n,103)==0
quietly gen double v432 = (_n+432)/7
quietly replace v432 = . if mod(_n,101)==0
quietly replace v432 = .z if mod(_n,103)==0
quietly gen double v433 = (_n+433)/7
quietly replace v433 = . if mod(_n,101)==0
quietly replace v433 = .z if mod(_n,103)==0
quietly gen double v434 = (_n+434)/7
quietly replace v434 = . if mod(_n,101)==0
quietly replace v434 = .z if mod(_n,103)==0
quietly gen double v435 = (_n+435)/7
quietly replace v435 = . if mod(_n,101)==0
quietly replace v435 = .z if mod(_n,103)==0
quietly gen double v436 = (_n+436)/7
quietly replace v436 = . if mod(_n,101)==0
quietly replace v436 = .z if mod(_n,103)==0
quietly gen double v437 = (_n+437)/7
quietly replace v437 = . if mod(_n,101)==0
quietly replace v437 = .z if mod(_n,103)==0
quietly gen double v438 = (_n+438)/7
quietly replace v438 = . if mod(_n,101)==0
quietly replace v438 = .z if mod(_n,103)==0
quietly gen double v439 = (_n+439)/7
quietly replace v439 = . if mod(_n,101)==0
quietly replace v439 = .z if mod(_n,103)==0
quietly gen double v440 = (_n+440)/7
quietly replace v440 = . if mod(_n,101)==0
quietly replace v440 = .z if mod(_n,103)==0
quietly gen double v441 = (_n+441)/7
quietly replace v441 = . if mod(_n,101)==0
quietly replace v441 = .z if mod(_n,103)==0
quietly gen double v442 = (_n+442)/7
quietly replace v442 = . if mod(_n,101)==0
quietly replace v442 = .z if mod(_n,103)==0
quietly gen double v443 = (_n+443)/7
quietly replace v443 = . if mod(_n,101)==0
quietly replace v443 = .z if mod(_n,103)==0
quietly gen double v444 = (_n+444)/7
quietly replace v444 = . if mod(_n,101)==0
quietly replace v444 = .z if mod(_n,103)==0
quietly gen double v445 = (_n+445)/7
quietly replace v445 = . if mod(_n,101)==0
quietly replace v445 = .z if mod(_n,103)==0
quietly gen double v446 = (_n+446)/7
quietly replace v446 = . if mod(_n,101)==0
quietly replace v446 = .z if mod(_n,103)==0
quietly gen double v447 = (_n+447)/7
quietly replace v447 = . if mod(_n,101)==0
quietly replace v447 = .z if mod(_n,103)==0
quietly gen double v448 = (_n+448)/7
quietly replace v448 = . if mod(_n,101)==0
quietly replace v448 = .z if mod(_n,103)==0
quietly gen double v449 = (_n+449)/7
quietly replace v449 = . if mod(_n,101)==0
quietly replace v449 = .z if mod(_n,103)==0
quietly gen double v450 = (_n+450)/7
quietly replace v450 = . if mod(_n,101)==0
quietly replace v450 = .z if mod(_n,103)==0
quietly gen double v451 = (_n+451)/7
quietly replace v451 = . if mod(_n,101)==0
quietly replace v451 = .z if mod(_n,103)==0
quietly gen double v452 = (_n+452)/7
quietly replace v452 = . if mod(_n,101)==0
quietly replace v452 = .z if mod(_n,103)==0
quietly gen double v453 = (_n+453)/7
quietly replace v453 = . if mod(_n,101)==0
quietly replace v453 = .z if mod(_n,103)==0
quietly gen double v454 = (_n+454)/7
quietly replace v454 = . if mod(_n,101)==0
quietly replace v454 = .z if mod(_n,103)==0
quietly gen double v455 = (_n+455)/7
quietly replace v455 = . if mod(_n,101)==0
quietly replace v455 = .z if mod(_n,103)==0
quietly gen double v456 = (_n+456)/7
quietly replace v456 = . if mod(_n,101)==0
quietly replace v456 = .z if mod(_n,103)==0
quietly gen double v457 = (_n+457)/7
quietly replace v457 = . if mod(_n,101)==0
quietly replace v457 = .z if mod(_n,103)==0
quietly gen double v458 = (_n+458)/7
quietly replace v458 = . if mod(_n,101)==0
quietly replace v458 = .z if mod(_n,103)==0
quietly gen double v459 = (_n+459)/7
quietly replace v459 = . if mod(_n,101)==0
quietly replace v459 = .z if mod(_n,103)==0
quietly gen double v460 = (_n+460)/7
quietly replace v460 = . if mod(_n,101)==0
quietly replace v460 = .z if mod(_n,103)==0
quietly gen double v461 = (_n+461)/7
quietly replace v461 = . if mod(_n,101)==0
quietly replace v461 = .z if mod(_n,103)==0
quietly gen double v462 = (_n+462)/7
quietly replace v462 = . if mod(_n,101)==0
quietly replace v462 = .z if mod(_n,103)==0
quietly gen double v463 = (_n+463)/7
quietly replace v463 = . if mod(_n,101)==0
quietly replace v463 = .z if mod(_n,103)==0
quietly gen double v464 = (_n+464)/7
quietly replace v464 = . if mod(_n,101)==0
quietly replace v464 = .z if mod(_n,103)==0
quietly gen double v465 = (_n+465)/7
quietly replace v465 = . if mod(_n,101)==0
quietly replace v465 = .z if mod(_n,103)==0
quietly gen double v466 = (_n+466)/7
quietly replace v466 = . if mod(_n,101)==0
quietly replace v466 = .z if mod(_n,103)==0
quietly gen double v467 = (_n+467)/7
quietly replace v467 = . if mod(_n,101)==0
quietly replace v467 = .z if mod(_n,103)==0
quietly gen double v468 = (_n+468)/7
quietly replace v468 = . if mod(_n,101)==0
quietly replace v468 = .z if mod(_n,103)==0
quietly gen double v469 = (_n+469)/7
quietly replace v469 = . if mod(_n,101)==0
quietly replace v469 = .z if mod(_n,103)==0
quietly gen double v470 = (_n+470)/7
quietly replace v470 = . if mod(_n,101)==0
quietly replace v470 = .z if mod(_n,103)==0
quietly gen double v471 = (_n+471)/7
quietly replace v471 = . if mod(_n,101)==0
quietly replace v471 = .z if mod(_n,103)==0
quietly gen double v472 = (_n+472)/7
quietly replace v472 = . if mod(_n,101)==0
quietly replace v472 = .z if mod(_n,103)==0
quietly gen double v473 = (_n+473)/7
quietly replace v473 = . if mod(_n,101)==0
quietly replace v473 = .z if mod(_n,103)==0
quietly gen double v474 = (_n+474)/7
quietly replace v474 = . if mod(_n,101)==0
quietly replace v474 = .z if mod(_n,103)==0
quietly gen double v475 = (_n+475)/7
quietly replace v475 = . if mod(_n,101)==0
quietly replace v475 = .z if mod(_n,103)==0
quietly gen double v476 = (_n+476)/7
quietly replace v476 = . if mod(_n,101)==0
quietly replace v476 = .z if mod(_n,103)==0
quietly gen double v477 = (_n+477)/7
quietly replace v477 = . if mod(_n,101)==0
quietly replace v477 = .z if mod(_n,103)==0
quietly gen double v478 = (_n+478)/7
quietly replace v478 = . if mod(_n,101)==0
quietly replace v478 = .z if mod(_n,103)==0
quietly gen double v479 = (_n+479)/7
quietly replace v479 = . if mod(_n,101)==0
quietly replace v479 = .z if mod(_n,103)==0
quietly gen double v480 = (_n+480)/7
quietly replace v480 = . if mod(_n,101)==0
quietly replace v480 = .z if mod(_n,103)==0
quietly gen double v481 = (_n+481)/7
quietly replace v481 = . if mod(_n,101)==0
quietly replace v481 = .z if mod(_n,103)==0
quietly gen double v482 = (_n+482)/7
quietly replace v482 = . if mod(_n,101)==0
quietly replace v482 = .z if mod(_n,103)==0
quietly gen double v483 = (_n+483)/7
quietly replace v483 = . if mod(_n,101)==0
quietly replace v483 = .z if mod(_n,103)==0
quietly gen double v484 = (_n+484)/7
quietly replace v484 = . if mod(_n,101)==0
quietly replace v484 = .z if mod(_n,103)==0
quietly gen double v485 = (_n+485)/7
quietly replace v485 = . if mod(_n,101)==0
quietly replace v485 = .z if mod(_n,103)==0
quietly gen double v486 = (_n+486)/7
quietly replace v486 = . if mod(_n,101)==0
quietly replace v486 = .z if mod(_n,103)==0
quietly gen double v487 = (_n+487)/7
quietly replace v487 = . if mod(_n,101)==0
quietly replace v487 = .z if mod(_n,103)==0
quietly gen double v488 = (_n+488)/7
quietly replace v488 = . if mod(_n,101)==0
quietly replace v488 = .z if mod(_n,103)==0
quietly gen double v489 = (_n+489)/7
quietly replace v489 = . if mod(_n,101)==0
quietly replace v489 = .z if mod(_n,103)==0
quietly gen double v490 = (_n+490)/7
quietly replace v490 = . if mod(_n,101)==0
quietly replace v490 = .z if mod(_n,103)==0
quietly gen double v491 = (_n+491)/7
quietly replace v491 = . if mod(_n,101)==0
quietly replace v491 = .z if mod(_n,103)==0
quietly gen double v492 = (_n+492)/7
quietly replace v492 = . if mod(_n,101)==0
quietly replace v492 = .z if mod(_n,103)==0
quietly gen double v493 = (_n+493)/7
quietly replace v493 = . if mod(_n,101)==0
quietly replace v493 = .z if mod(_n,103)==0
quietly gen double v494 = (_n+494)/7
quietly replace v494 = . if mod(_n,101)==0
quietly replace v494 = .z if mod(_n,103)==0
quietly gen double v495 = (_n+495)/7
quietly replace v495 = . if mod(_n,101)==0
quietly replace v495 = .z if mod(_n,103)==0
quietly gen double v496 = (_n+496)/7
quietly replace v496 = . if mod(_n,101)==0
quietly replace v496 = .z if mod(_n,103)==0
quietly gen double v497 = (_n+497)/7
quietly replace v497 = . if mod(_n,101)==0
quietly replace v497 = .z if mod(_n,103)==0
quietly gen double v498 = (_n+498)/7
quietly replace v498 = . if mod(_n,101)==0
quietly replace v498 = .z if mod(_n,103)==0
quietly gen double v499 = (_n+499)/7
quietly replace v499 = . if mod(_n,101)==0
quietly replace v499 = .z if mod(_n,103)==0
quietly gen double v500 = (_n+500)/7
quietly replace v500 = . if mod(_n,101)==0
quietly replace v500 = .z if mod(_n,103)==0
quietly gen double v501 = (_n+501)/7
quietly replace v501 = . if mod(_n,101)==0
quietly replace v501 = .z if mod(_n,103)==0
quietly gen double v502 = (_n+502)/7
quietly replace v502 = . if mod(_n,101)==0
quietly replace v502 = .z if mod(_n,103)==0
quietly gen double v503 = (_n+503)/7
quietly replace v503 = . if mod(_n,101)==0
quietly replace v503 = .z if mod(_n,103)==0
quietly gen double v504 = (_n+504)/7
quietly replace v504 = . if mod(_n,101)==0
quietly replace v504 = .z if mod(_n,103)==0
quietly gen double v505 = (_n+505)/7
quietly replace v505 = . if mod(_n,101)==0
quietly replace v505 = .z if mod(_n,103)==0
quietly gen double v506 = (_n+506)/7
quietly replace v506 = . if mod(_n,101)==0
quietly replace v506 = .z if mod(_n,103)==0
quietly gen double v507 = (_n+507)/7
quietly replace v507 = . if mod(_n,101)==0
quietly replace v507 = .z if mod(_n,103)==0
quietly gen double v508 = (_n+508)/7
quietly replace v508 = . if mod(_n,101)==0
quietly replace v508 = .z if mod(_n,103)==0
quietly gen double v509 = (_n+509)/7
quietly replace v509 = . if mod(_n,101)==0
quietly replace v509 = .z if mod(_n,103)==0
quietly gen double v510 = (_n+510)/7
quietly replace v510 = . if mod(_n,101)==0
quietly replace v510 = .z if mod(_n,103)==0
quietly gen double v511 = (_n+511)/7
quietly replace v511 = . if mod(_n,101)==0
quietly replace v511 = .z if mod(_n,103)==0
quietly gen double v512 = (_n+512)/7
quietly replace v512 = . if mod(_n,101)==0
quietly replace v512 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r0" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r1" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r2" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r3" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r4" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r5" "8" "read" "0"
plugin call transport_io v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 v15 v16 v17 v18 v19 v20 v21 v22 v23 v24 v25 v26 v27 v28 v29 v30 v31 v32 v33 v34 v35 v36 v37 v38 v39 v40 v41 v42 v43 v44 v45 v46 v47 v48 v49 v50 v51 v52 v53 v54 v55 v56 v57 v58 v59 v60 v61 v62 v63 v64 v65 v66 v67 v68 v69 v70 v71 v72 v73 v74 v75 v76 v77 v78 v79 v80 v81 v82 v83 v84 v85 v86 v87 v88 v89 v90 v91 v92 v93 v94 v95 v96 v97 v98 v99 v100 v101 v102 v103 v104 v105 v106 v107 v108 v109 v110 v111 v112 v113 v114 v115 v116 v117 v118 v119 v120 v121 v122 v123 v124 v125 v126 v127 v128 v129 v130 v131 v132 v133 v134 v135 v136 v137 v138 v139 v140 v141 v142 v143 v144 v145 v146 v147 v148 v149 v150 v151 v152 v153 v154 v155 v156 v157 v158 v159 v160 v161 v162 v163 v164 v165 v166 v167 v168 v169 v170 v171 v172 v173 v174 v175 v176 v177 v178 v179 v180 v181 v182 v183 v184 v185 v186 v187 v188 v189 v190 v191 v192 v193 v194 v195 v196 v197 v198 v199 v200 v201 v202 v203 v204 v205 v206 v207 v208 v209 v210 v211 v212 v213 v214 v215 v216 v217 v218 v219 v220 v221 v222 v223 v224 v225 v226 v227 v228 v229 v230 v231 v232 v233 v234 v235 v236 v237 v238 v239 v240 v241 v242 v243 v244 v245 v246 v247 v248 v249 v250 v251 v252 v253 v254 v255 v256 v257 v258 v259 v260 v261 v262 v263 v264 v265 v266 v267 v268 v269 v270 v271 v272 v273 v274 v275 v276 v277 v278 v279 v280 v281 v282 v283 v284 v285 v286 v287 v288 v289 v290 v291 v292 v293 v294 v295 v296 v297 v298 v299 v300 v301 v302 v303 v304 v305 v306 v307 v308 v309 v310 v311 v312 v313 v314 v315 v316 v317 v318 v319 v320 v321 v322 v323 v324 v325 v326 v327 v328 v329 v330 v331 v332 v333 v334 v335 v336 v337 v338 v339 v340 v341 v342 v343 v344 v345 v346 v347 v348 v349 v350 v351 v352 v353 v354 v355 v356 v357 v358 v359 v360 v361 v362 v363 v364 v365 v366 v367 v368 v369 v370 v371 v372 v373 v374 v375 v376 v377 v378 v379 v380 v381 v382 v383 v384 v385 v386 v387 v388 v389 v390 v391 v392 v393 v394 v395 v396 v397 v398 v399 v400 v401 v402 v403 v404 v405 v406 v407 v408 v409 v410 v411 v412 v413 v414 v415 v416 v417 v418 v419 v420 v421 v422 v423 v424 v425 v426 v427 v428 v429 v430 v431 v432 v433 v434 v435 v436 v437 v438 v439 v440 v441 v442 v443 v444 v445 v446 v447 v448 v449 v450 v451 v452 v453 v454 v455 v456 v457 v458 v459 v460 v461 v462 v463 v464 v465 v466 v467 v468 v469 v470 v471 v472 v473 v474 v475 v476 v477 v478 v479 v480 v481 v482 v483 v484 v485 v486 v487 v488 v489 v490 v491 v492 v493 v494 v495 v496 v497 v498 v499 v500 v501 v502 v503 v504 v505 v506 v507 v508 v509 v510 v511 v512, "cells_k512_double_r6" "8" "read" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
