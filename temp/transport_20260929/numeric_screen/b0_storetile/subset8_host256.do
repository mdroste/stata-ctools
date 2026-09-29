clear
quietly set obs 200000
local transport_pad "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
quietly gen byte v1 = mod(_n+1,101)-50
quietly replace v1 = . if mod(_n,101)==0
quietly replace v1 = .z if mod(_n,103)==0
quietly gen double v2 = (_n+2)/7
quietly gen double v3 = (_n+3)/7
quietly gen double v4 = (_n+4)/7
quietly gen double v5 = (_n+5)/7
quietly gen double v6 = (_n+6)/7
quietly gen double v7 = (_n+7)/7
quietly gen double v8 = (_n+8)/7
quietly gen double v9 = (_n+9)/7
quietly gen double v10 = (_n+10)/7
quietly gen double v11 = (_n+11)/7
quietly gen double v12 = (_n+12)/7
quietly gen double v13 = (_n+13)/7
quietly gen double v14 = (_n+14)/7
quietly gen double v15 = (_n+15)/7
quietly gen double v16 = (_n+16)/7
quietly gen double v17 = (_n+17)/7
quietly gen double v18 = (_n+18)/7
quietly gen double v19 = (_n+19)/7
quietly gen double v20 = (_n+20)/7
quietly gen double v21 = (_n+21)/7
quietly gen double v22 = (_n+22)/7
quietly gen double v23 = (_n+23)/7
quietly gen double v24 = (_n+24)/7
quietly gen double v25 = (_n+25)/7
quietly gen double v26 = (_n+26)/7
quietly gen double v27 = (_n+27)/7
quietly gen double v28 = (_n+28)/7
quietly gen double v29 = (_n+29)/7
quietly gen double v30 = (_n+30)/7
quietly gen double v31 = (_n+31)/7
quietly gen double v32 = (_n+32)/7
quietly gen double v33 = (_n+33)/7
quietly gen double v34 = (_n+34)/7
quietly gen double v35 = (_n+35)/7
quietly gen double v36 = (_n+36)/7
quietly gen int v37 = mod(_n+37,101)-50
quietly replace v37 = . if mod(_n,101)==0
quietly replace v37 = .z if mod(_n,103)==0
quietly gen double v38 = (_n+38)/7
quietly gen double v39 = (_n+39)/7
quietly gen double v40 = (_n+40)/7
quietly gen double v41 = (_n+41)/7
quietly gen double v42 = (_n+42)/7
quietly gen double v43 = (_n+43)/7
quietly gen double v44 = (_n+44)/7
quietly gen double v45 = (_n+45)/7
quietly gen double v46 = (_n+46)/7
quietly gen double v47 = (_n+47)/7
quietly gen double v48 = (_n+48)/7
quietly gen double v49 = (_n+49)/7
quietly gen double v50 = (_n+50)/7
quietly gen double v51 = (_n+51)/7
quietly gen double v52 = (_n+52)/7
quietly gen double v53 = (_n+53)/7
quietly gen double v54 = (_n+54)/7
quietly gen double v55 = (_n+55)/7
quietly gen double v56 = (_n+56)/7
quietly gen double v57 = (_n+57)/7
quietly gen double v58 = (_n+58)/7
quietly gen double v59 = (_n+59)/7
quietly gen double v60 = (_n+60)/7
quietly gen double v61 = (_n+61)/7
quietly gen double v62 = (_n+62)/7
quietly gen double v63 = (_n+63)/7
quietly gen double v64 = (_n+64)/7
quietly gen double v65 = (_n+65)/7
quietly gen double v66 = (_n+66)/7
quietly gen double v67 = (_n+67)/7
quietly gen double v68 = (_n+68)/7
quietly gen double v69 = (_n+69)/7
quietly gen double v70 = (_n+70)/7
quietly gen double v71 = (_n+71)/7
quietly gen double v72 = (_n+72)/7
quietly gen long v73 = _n+73
quietly replace v73 = . if mod(_n,101)==0
quietly replace v73 = .z if mod(_n,103)==0
quietly gen double v74 = (_n+74)/7
quietly gen double v75 = (_n+75)/7
quietly gen double v76 = (_n+76)/7
quietly gen double v77 = (_n+77)/7
quietly gen double v78 = (_n+78)/7
quietly gen double v79 = (_n+79)/7
quietly gen double v80 = (_n+80)/7
quietly gen double v81 = (_n+81)/7
quietly gen double v82 = (_n+82)/7
quietly gen double v83 = (_n+83)/7
quietly gen double v84 = (_n+84)/7
quietly gen double v85 = (_n+85)/7
quietly gen double v86 = (_n+86)/7
quietly gen double v87 = (_n+87)/7
quietly gen double v88 = (_n+88)/7
quietly gen double v89 = (_n+89)/7
quietly gen double v90 = (_n+90)/7
quietly gen double v91 = (_n+91)/7
quietly gen double v92 = (_n+92)/7
quietly gen double v93 = (_n+93)/7
quietly gen double v94 = (_n+94)/7
quietly gen double v95 = (_n+95)/7
quietly gen double v96 = (_n+96)/7
quietly gen double v97 = (_n+97)/7
quietly gen double v98 = (_n+98)/7
quietly gen double v99 = (_n+99)/7
quietly gen double v100 = (_n+100)/7
quietly gen double v101 = (_n+101)/7
quietly gen double v102 = (_n+102)/7
quietly gen double v103 = (_n+103)/7
quietly gen double v104 = (_n+104)/7
quietly gen double v105 = (_n+105)/7
quietly gen double v106 = (_n+106)/7
quietly gen double v107 = (_n+107)/7
quietly gen double v108 = (_n+108)/7
quietly gen double v109 = (_n+109)/7
quietly gen float v110 = (_n+110)/7
quietly replace v110 = . if mod(_n,101)==0
quietly replace v110 = .z if mod(_n,103)==0
quietly gen double v111 = (_n+111)/7
quietly gen double v112 = (_n+112)/7
quietly gen double v113 = (_n+113)/7
quietly gen double v114 = (_n+114)/7
quietly gen double v115 = (_n+115)/7
quietly gen double v116 = (_n+116)/7
quietly gen double v117 = (_n+117)/7
quietly gen double v118 = (_n+118)/7
quietly gen double v119 = (_n+119)/7
quietly gen double v120 = (_n+120)/7
quietly gen double v121 = (_n+121)/7
quietly gen double v122 = (_n+122)/7
quietly gen double v123 = (_n+123)/7
quietly gen double v124 = (_n+124)/7
quietly gen double v125 = (_n+125)/7
quietly gen double v126 = (_n+126)/7
quietly gen double v127 = (_n+127)/7
quietly gen double v128 = (_n+128)/7
quietly gen double v129 = (_n+129)/7
quietly gen double v130 = (_n+130)/7
quietly gen double v131 = (_n+131)/7
quietly gen double v132 = (_n+132)/7
quietly gen double v133 = (_n+133)/7
quietly gen double v134 = (_n+134)/7
quietly gen double v135 = (_n+135)/7
quietly gen double v136 = (_n+136)/7
quietly gen double v137 = (_n+137)/7
quietly gen double v138 = (_n+138)/7
quietly gen double v139 = (_n+139)/7
quietly gen double v140 = (_n+140)/7
quietly gen double v141 = (_n+141)/7
quietly gen double v142 = (_n+142)/7
quietly gen double v143 = (_n+143)/7
quietly gen double v144 = (_n+144)/7
quietly gen double v145 = (_n+145)/7
quietly gen double v146 = (_n+146)/7
quietly replace v146 = . if mod(_n,101)==0
quietly replace v146 = .z if mod(_n,103)==0
quietly gen double v147 = (_n+147)/7
quietly gen double v148 = (_n+148)/7
quietly gen double v149 = (_n+149)/7
quietly gen double v150 = (_n+150)/7
quietly gen double v151 = (_n+151)/7
quietly gen double v152 = (_n+152)/7
quietly gen double v153 = (_n+153)/7
quietly gen double v154 = (_n+154)/7
quietly gen double v155 = (_n+155)/7
quietly gen double v156 = (_n+156)/7
quietly gen double v157 = (_n+157)/7
quietly gen double v158 = (_n+158)/7
quietly gen double v159 = (_n+159)/7
quietly gen double v160 = (_n+160)/7
quietly gen double v161 = (_n+161)/7
quietly gen double v162 = (_n+162)/7
quietly gen double v163 = (_n+163)/7
quietly gen double v164 = (_n+164)/7
quietly gen double v165 = (_n+165)/7
quietly gen double v166 = (_n+166)/7
quietly gen double v167 = (_n+167)/7
quietly gen double v168 = (_n+168)/7
quietly gen double v169 = (_n+169)/7
quietly gen double v170 = (_n+170)/7
quietly gen double v171 = (_n+171)/7
quietly gen double v172 = (_n+172)/7
quietly gen double v173 = (_n+173)/7
quietly gen double v174 = (_n+174)/7
quietly gen double v175 = (_n+175)/7
quietly gen double v176 = (_n+176)/7
quietly gen double v177 = (_n+177)/7
quietly gen double v178 = (_n+178)/7
quietly gen double v179 = (_n+179)/7
quietly gen double v180 = (_n+180)/7
quietly gen double v181 = (_n+181)/7
quietly gen double v182 = (_n+182)/7
quietly gen byte v183 = mod(_n+183,101)-50
quietly replace v183 = . if mod(_n,101)==0
quietly replace v183 = .z if mod(_n,103)==0
quietly gen double v184 = (_n+184)/7
quietly gen double v185 = (_n+185)/7
quietly gen double v186 = (_n+186)/7
quietly gen double v187 = (_n+187)/7
quietly gen double v188 = (_n+188)/7
quietly gen double v189 = (_n+189)/7
quietly gen double v190 = (_n+190)/7
quietly gen double v191 = (_n+191)/7
quietly gen double v192 = (_n+192)/7
quietly gen double v193 = (_n+193)/7
quietly gen double v194 = (_n+194)/7
quietly gen double v195 = (_n+195)/7
quietly gen double v196 = (_n+196)/7
quietly gen double v197 = (_n+197)/7
quietly gen double v198 = (_n+198)/7
quietly gen double v199 = (_n+199)/7
quietly gen double v200 = (_n+200)/7
quietly gen double v201 = (_n+201)/7
quietly gen double v202 = (_n+202)/7
quietly gen double v203 = (_n+203)/7
quietly gen double v204 = (_n+204)/7
quietly gen double v205 = (_n+205)/7
quietly gen double v206 = (_n+206)/7
quietly gen double v207 = (_n+207)/7
quietly gen double v208 = (_n+208)/7
quietly gen double v209 = (_n+209)/7
quietly gen double v210 = (_n+210)/7
quietly gen double v211 = (_n+211)/7
quietly gen double v212 = (_n+212)/7
quietly gen double v213 = (_n+213)/7
quietly gen double v214 = (_n+214)/7
quietly gen double v215 = (_n+215)/7
quietly gen double v216 = (_n+216)/7
quietly gen double v217 = (_n+217)/7
quietly gen double v218 = (_n+218)/7
quietly gen int v219 = mod(_n+219,101)-50
quietly replace v219 = . if mod(_n,101)==0
quietly replace v219 = .z if mod(_n,103)==0
quietly gen double v220 = (_n+220)/7
quietly gen double v221 = (_n+221)/7
quietly gen double v222 = (_n+222)/7
quietly gen double v223 = (_n+223)/7
quietly gen double v224 = (_n+224)/7
quietly gen double v225 = (_n+225)/7
quietly gen double v226 = (_n+226)/7
quietly gen double v227 = (_n+227)/7
quietly gen double v228 = (_n+228)/7
quietly gen double v229 = (_n+229)/7
quietly gen double v230 = (_n+230)/7
quietly gen double v231 = (_n+231)/7
quietly gen double v232 = (_n+232)/7
quietly gen double v233 = (_n+233)/7
quietly gen double v234 = (_n+234)/7
quietly gen double v235 = (_n+235)/7
quietly gen double v236 = (_n+236)/7
quietly gen double v237 = (_n+237)/7
quietly gen double v238 = (_n+238)/7
quietly gen double v239 = (_n+239)/7
quietly gen double v240 = (_n+240)/7
quietly gen double v241 = (_n+241)/7
quietly gen double v242 = (_n+242)/7
quietly gen double v243 = (_n+243)/7
quietly gen double v244 = (_n+244)/7
quietly gen double v245 = (_n+245)/7
quietly gen double v246 = (_n+246)/7
quietly gen double v247 = (_n+247)/7
quietly gen double v248 = (_n+248)/7
quietly gen double v249 = (_n+249)/7
quietly gen double v250 = (_n+250)/7
quietly gen double v251 = (_n+251)/7
quietly gen double v252 = (_n+252)/7
quietly gen double v253 = (_n+253)/7
quietly gen double v254 = (_n+254)/7
quietly gen double v255 = (_n+255)/7
quietly gen long v256 = _n+256
quietly replace v256 = . if mod(_n,101)==0
quietly replace v256 = .z if mod(_n,103)==0
quietly datasignature
local signature "`r(datasignature)'"
_ctools_strw v1 v37 v73 v110 v146 v183 v219 v256
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r0" "8" "identity" "0"
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r1" "8" "identity" "0"
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r2" "8" "identity" "0"
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r3" "8" "identity" "0"
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r4" "8" "identity" "0"
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r5" "8" "identity" "0"
plugin call transport_io v1 v37 v73 v110 v146 v183 v219 v256, "storetile_subset8_host256_r6" "8" "identity" "1"
quietly datasignature
assert "`r(datasignature)'" == "`signature'"
