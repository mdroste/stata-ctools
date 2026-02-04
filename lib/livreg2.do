*! livreg2 1.2.01  28july2015
*! authors cfb & mes
*! compiled in Stata 11.2
* Mata library for ivreg2 and ranktest.
* Introduced with ivreg2 v 3.1.01 and ranktest v 1.3.01.

* Version comments:
* 1.1.01     First version of library.
*            Contains struct ms_vcvorthog, m_omega, m_calckw, s_vkernel.
*            Compiled in Stata 9.2 for compatibility with ranktest 1.3.01 (a 9.2 program).
* 1.1.02     Add routine cdsy. Standardized spelling/caps/etc. of QS as "Quadratic Spectral"
* 1.1.03     Corrected spelling of "Danielle" kernel in m_omega()
* 1.1.04     Fixed weighting bugs in robust and cluster code of m_omega where K>1
* 1.1.05     Added whichlivreg2(.) to aid in version control
* 1.1.06     Fixed remaining weighting bug (see 1.1.04) in 2-way clustering when interection
*            of clustering levels is groups
* 1.1.07     Fixed HAC bug that crashed m_omega(.) when there were no obs for a particular lag
* 1.2.01     Promotion to compilation in Stata 11.2.  Added support for -center- option.

* Locals used in whichlivreg2; set when compiled
local stata_version `c(stata_version)'
local born_date `c(born_date)'

version 11.2
mata:
mata clear

void whichlivreg2()
{
""
"livreg2 01.2.01  28july2015"
"compiled under Stata " + "`stata_version'" + " born " + "`born_date'"
"Mata library for ivreg2 and related programs"
"authors CFB/MS"
st_sclear()
st_global("s(stata_born_date)","`born_date'")
st_global("s(stata_version)","`stata_version'")
st_global("s(compiled_date)","28july2015")
st_global("s(ver)", "01.2.01")
}


// ********* struct ms_vcvorthog - shared by ivreg2 and ranktest ******************* //
struct ms_vcvorthog {
	string scalar	ename, Znames, touse, weight, wvarname
	string scalar	robust, clustvarname, clustvarname2, clustvarname3, kernel
	string scalar	sw, psd, ivarname, tvarname, tindexname
	real scalar		wf, N, bw, tdelta, dofminus
	real scalar		center
	real matrix		ZZ
	pointer matrix	e
	pointer matrix	Z
	pointer matrix	wvar
}

// ********* s_vkernel - shared by ivreg2 and ranktest ******************* //
// Program checks whether kernel and bw choices are valid.
// s_vkernel is called from Stata.
// Arguments are the kernel name (req), bandwidth (req) and ivar name (opt).
// All 3 are strings.
// Returns results in r() macros.
// r(kernel) - name of kernel (string)
// r(bw) - bandwidth (scalar)

void s_vkernel(	string scalar kernel,
			string scalar bwstring,
			string scalar ivar
				)
{

// Check bandwidth
	if (bwstring=="auto") {
		bw=-1
	}
	else {
		bw=strtoreal(bwstring)
		if (bw==.) {
			printf("{err}bandwidth option bw() required for HAC-robust estimation\n")
			exit(102)
		}
		if (bw<=0) {
			printf("{err}invalid bandwidth in option bw() - must be real > 0\n")
			exit(198)
		}
	}
	
// Check ivar
	if (bwstring=="auto" & ivar~="") {
			printf("{err}Automatic bandwidth selection not available for panel data\n")
			exit(198)
	}

// Check kernel
// Valid kernel list is abbrev, full name, whether special case if bw=1
// First in list is default kernel = Barlett
	vklist = 	(	("", "bartlett", "0")
				\	("bar", "bartlett", "0")
				\	("bartlett", "bartlett", "0")
				\	("par", "parzen", "0")
				\	("parzen", "parzen", "0")
				\	("tru", "truncated", "1")
				\	("truncated", "truncated", "1")
				\	("thann", "tukey-hanning", "0")
				\	("tukey-hanning", "tukey-hanning", "0")
				\	("thamm", "tukey-hamming", "0")
				\	("tukey-hamming", "tukey-hamming", "0")
				\	("qua", "quadratic spectral", "1")
				\	("qs", "quadratic spectral", "1")
				\	("quadratic-spectral", "quadratic spectral", "1")
				\	("quadratic spectral", "quadratic spectral", "1")
				\	("dan", "danielle", "1")
				\	("danielle", "danielle", "1")
				\	("ten", "tent", "1")
				\	("tent", "tent", "1")
			)
	kname=strltrim(strlower(kernel))
	pos = (vklist[.,1] :== kname)

// Exit with error if not in list
	if (sum(pos)==0) {
		printf("{err}invalid kernel\n")
		exit(198)
		}

	vkname=strproper(select(vklist[.,2],pos))
	st_global("r(kernel)", vkname)
	st_numscalar("r(bw)",bw)

// Warn if kernel is type where bw=1 means no lags are used
	if (bw==1 & select(vklist[.,3],pos)=="0") {
		printf("{result}Note: kernel=%s", vkname)
		printf("{result} and bw=1 implies zero lags used.  Standard errors and\n")
		printf("{result}      test statistics are not autocorrelation-consistent.\n")
	}
}  // end of program s_vkernel

// ********* m_omega - shared by ivreg2 and ranktest ********************* //

// NB: ivreg2 always calls m_omega with e as column vector, i.e., K=1      //
//     ranktest can call m_omega with e as matrix, i.e., K>=1              //

real matrix m_omega(struct ms_vcvorthog scalar vcvo) 
{
	if (vcvo.clustvarname~="") {
		st_view(clustvar, ., vcvo.clustvarname, vcvo.touse)
		info = panelsetup(clustvar, 1)
		N_clust=rows(info)
		if (vcvo.clustvarname2~="") {
			st_view(clustvar2, ., vcvo.clustvarname2, vcvo.touse)
			if (vcvo.kernel=="") {
				st_view(clustvar3, ., vcvo.clustvarname3, vcvo.touse) // needed only if not panel tsset
			}
		}
	}

	if (vcvo.kernel~="") {
		st_view(t,    ., st_tsrevar(vcvo.tvarname),  vcvo.touse)
		T=max(t)-min(t)+1
	}

	if ((vcvo.kernel=="Bartlett") | (vcvo.kernel=="Parzen") | (vcvo.kernel=="Truncated") ///
		 | (vcvo.kernel=="Tukey-Hanning")| (vcvo.kernel=="Tukey-Hamming")) {
		window="lag"
	}
	else if ((vcvo.kernel=="Quadratic Spectral") | (vcvo.kernel=="Danielle") | (vcvo.kernel=="Tent")) {
		window="spectral"
	}
	else if (vcvo.kernel~="") {
// Should never reach this point
printf("\n{error:Error: invalid kernel}\n")
		exit(error(3351))
	}

	L=cols(*vcvo.Z)
	K=cols(*vcvo.e)		// ivreg2 always calls with K=1; ranktest may call with K>=1.

// If centered moments specified, need rowvector of mean organized as (e1'Z e2'Z ...).
// Needed only if non-homoskedastic
	if ((vcvo.center==1) & ((vcvo.robust~="") | (vcvo.clustvarname==""))) {
		if ((vcvo.weight=="fweight") | (vcvo.weight=="iweight")) {
			wv = *vcvo.wvar
		}
		else {
			wv = (*vcvo.wvar * vcvo.wf):^2		// wf needed for aweights and pweights
		}
		eZmean=quadcross(*vcvo.e, wv, *vcvo.Z) / vcvo.N
		eZmean=vec(eZmean')'
	}

// Covariance matrices
// shat * 1/N is same as estimated S matrix of orthog conditions

// Block for homoskedastic and AC.  dof correction if any incorporated into sigma estimates.
	if ((vcvo.robust=="") & (vcvo.clustvarname=="")) {
// ZZ is already calculated as an external
		ee = quadcross(*vcvo.e, vcvo.wf*(*vcvo.wvar), *vcvo.e)
		sigma2=ee/(vcvo.N-vcvo.dofminus)
		shat=sigma2#vcvo.ZZ
		if (vcvo.kernel~="") {
			if (window=="spectral") {
				TAU=T/vcvo.tdelta-1
			}
			else {
				TAU=vcvo.bw
			}
			tnow=st_data(., vcvo.tindexname)
			for (tau=1; tau<=TAU; tau++) {
				kw = m_calckw(tau, vcvo.bw, vcvo.kernel)
				if (kw~=0) {						// zero weight possible with some kernels
													// save an unnecessary loop if kw=0
													// remember, kw<0 possible with some kernels!
					lstau = "L"+strofreal(tau)
					tlag=st_data(., lstau+"."+vcvo.tindexname)
					tmatrix = tnow, tlag
					svar=(tnow:<.):*(tlag:<.)		// multiply column vectors of 1s and 0s
					tmatrix=select(tmatrix,svar)	// to get intersection, and replace tmatrix
// if no lags exist, tmatrix has zero rows.
					if (rows(tmatrix)>0) {
// col 1 of tmatrix has row numbers of all rows of data with this time period that have a corresponding lag
// col 2 of tmatrix has row numbers of all rows of data with lag tau that have a corresponding ob this time period.
// Should never happen that fweights or iweights make it here,
// but if they did the next line would be sqrt(wvari)*sqrt(wvari1) [with no wf since not needed for fw or iw]
						wv = (*vcvo.wvar)[tmatrix[.,1]]		///
									:* (*vcvo.wvar)[tmatrix[.,2]]*(vcvo.wf^2)	// inner weighting matrix for quadcross
						sigmahat = quadcross((*vcvo.e)[tmatrix[.,1],.],   wv ,(*vcvo.e)[tmatrix[.,2],.])	///
									/ (vcvo.N-vcvo.dofminus)					// large dof correction
						ZZhat    = quadcross((*vcvo.Z)[tmatrix[.,1],.], wv, (*vcvo.Z)[tmatrix[.,2],.])
						ghat = sigmahat#ZZhat
						shat=shat+kw*(ghat+ghat')
					}
				}	// end non-zero kernel weight block
			}	// end tau loop
		}  // end kernel code
// Note large dof correction (if there is one) has already been incorporated
	shat=shat/vcvo.N
	}  // end homoskedastic, AC code

// Block for robust HC and HAC but not Stock-Watson and single clustering.
// Need to enter for double-clustering if one cluster is time.
	if ( (vcvo.robust~="") & (vcvo.sw=="") & ((vcvo.clustvarname=="")		///
			| ((vcvo.clustvarname2~="") & (vcvo.kernel~="")))  ) {
		if ((K==1) & (vcvo.center~=1)) {				// simple/fast where e is a column vector
														// and no centering
			if ((vcvo.weight=="fweight") | (vcvo.weight=="iweight")) {
				wv = (*vcvo.e:^2) :* *vcvo.wvar
			}
			else {
				wv = (*vcvo.e :* *vcvo.wvar * vcvo.wf):^2		// wf needed for aweights and pweights
			}
			shat=quadcross(*vcvo.Z, wv, *vcvo.Z)		// basic Eicker-Huber-White-sandwich-robust vce
		}
		else {											// e is a matrix so must loop
			shat=J(L*K,L*K,0)
			for (i=1; i<=rows(*vcvo.e); i++) {
				eZi=((*vcvo.e)[i,.])#((*vcvo.Z)[i,.])
				if (vcvo.center==1) {
					eZi=eZi-eZmean
				}
				if ((vcvo.weight=="fweight") | (vcvo.weight=="iweight")) {
// wvar is a column vector. wf not needed for fw and iw (=1 by dfn so redundant).
					shat=shat+quadcross(eZi,eZi)*((*vcvo.wvar)[i])
				}
				else {
					shat=shat+quadcross(eZi,eZi)*((*vcvo.wvar)[i] * vcvo.wf)^2
				}
			}
		}
		if (vcvo.kernel~="") {
// Spectral windows require looping through all T-1 autocovariances
			if (window=="spectral") {
				TAU=T/vcvo.tdelta-1
			}
			else {
				TAU=vcvo.bw
			}
			tnow=st_data(., vcvo.tindexname)
			for (tau=1; tau<=TAU; tau++) {
				kw = m_calckw(tau, vcvo.bw, vcvo.kernel)
				if (kw~=0) {						// zero weight possible with some kernels
													// save an unnecessary loop if kw=0
													// remember, kw<0 possible with some kernels!
					lstau = "L"+strofreal(tau)
					tlag=st_data(., lstau+"."+vcvo.tindexname)
					tmatrix = tnow, tlag
					svar=(tnow:<.):*(tlag:<.)		// multiply column vectors of 1s and 0s
					tmatrix=select(tmatrix,svar)	// to get intersection, and replace tmatrix

// col 1 of tmatrix has row numbers of all rows of data with this time period that have a corresponding lag
// col 2 of tmatrix has row numbers of all rows of data with lag tau that have a corresponding ob this time period.
// if no lags exist, tmatrix has zero rows
					if (rows(tmatrix)>0) {
						if ((K==1) & (vcvo.center~=1)) {			// simple/fast where e is a column vector
																	// and no centering
// wv is inner weighting matrix for quadcross
							wv   = (*vcvo.e)[tmatrix[.,1]] :* (*vcvo.e)[tmatrix[.,2]]		///
								:* (*vcvo.wvar)[tmatrix[.,1]] :* (*vcvo.wvar)[tmatrix[.,2]] * (vcvo.wf^2)
							ghat = quadcross((*vcvo.Z)[tmatrix[.,1],.], wv, (*vcvo.Z)[tmatrix[.,2],.])
						}
						else {										// e is a matrix so must loop
							ghat=J(L*K,L*K,0)
							for (i=1; i<=rows(tmatrix); i++) {
								wvari =(*vcvo.wvar)[tmatrix[i,1]]
								wvari1=(*vcvo.wvar)[tmatrix[i,2]]
								ei    =(*vcvo.e)[tmatrix[i,1],.]
								ei1   =(*vcvo.e)[tmatrix[i,2],.]
								Zi    =(*vcvo.Z)[tmatrix[i,1],.]
								Zi1   =(*vcvo.Z)[tmatrix[i,2],.]
								eZi =ei#Zi
								eZi1=ei1#Zi1
								if (vcvo.center==1) {
									eZi=eZi-eZmean
									eZi1=eZi1-eZmean
								}
// Should never happen that fweights or iweights make it here, but if they did
// the next line would be ghat=ghat+eZi'*eZi1*sqrt(wvari)*sqrt(wvari1)
// [without *vcvo.wf since wf=1 for fw and iw]
								ghat=ghat+quadcross(eZi,eZi1)*wvari*wvari1 * (vcvo.wf^2)
							}
						}
						shat=shat+kw*(ghat+ghat')
					}	// end non-zero-obs accumulation block
				}	// end non-zero kernel weight block
			}	// end tau loop
		}  // end kernel code
// Incorporate large dof correction if there is one
	shat=shat/(vcvo.N-vcvo.dofminus)
	}  // end HC/HAC code

	if (vcvo.clustvarname~="") {
// Block for cluster-robust
// 2-level clustering: S = S(level 1) + S(level 2) - S(level 3 = intersection of levels 1 & 2)
// Prepare shat3 if 2-level clustering
		if (vcvo.clustvarname2~="") {
			if (vcvo.kernel~="") {	// second cluster variable is time
									// shat3 was already calculated above as shat
				shat3=shat*(vcvo.N-vcvo.dofminus)
			}
			else {					// calculate shat3
									// data were sorted on clustvar3-clustvar1 so
									// clustvar3 is nested in clustvar1 and Mata panel functions
									// work for both.
				info3 = panelsetup(clustvar3, 1)
				if (rows(info3)==rows(*vcvo.e)) {	// intersection of levels 1 & 2 are all single obs
													// so no need to loop through row by row
				if ((K==1) & (vcvo.center~=1)) {			// simple/fast where e is a column vector
															// and no centering
						wv = (*vcvo.e :* *vcvo.wvar * vcvo.wf):^2
						shat3=quadcross(*vcvo.Z, wv, *vcvo.Z)		// basic Eicker-Huber-White-sandwich-robust vce
					}
					else {											// e is a matrix so must loop
						shat3=J(L*K,L*K,0)
						for (i=1; i<=rows(*vcvo.e); i++) {
							eZi=((*vcvo.e)[i,.])#((*vcvo.Z)[i,.])
							if (vcvo.center==1) {
								eZi=eZi-eZmean
							}
							shat3=shat3+quadcross(eZi,eZi)*((*vcvo.wvar)[i] * vcvo.wf)^2
							}
						}
				}
				else {								// intersection of levels 1 & 2 includes some groups of obs
					N_clust3=rows(info3)
					shat3=J(L*K,L*K,0)
					for (i=1; i<=N_clust3; i++) {
						esub=panelsubmatrix(*vcvo.e,i,info3)
						Zsub=panelsubmatrix(*vcvo.Z,i,info3)
						wsub=panelsubmatrix(*vcvo.wvar,i,info3)
						wv = esub :* wsub * vcvo.wf
						if ((K==1) & (vcvo.center~=1)) {	// simple/fast where e is a column vector
															// and no centering
							eZ = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
						}
						else {
							eZ = J(1,L*K,0)
							for (j=1; j<=rows(esub); j++) {
								eZ=eZ+(esub[j,.]#Zsub[j,.])*wsub[j,.] * vcvo.wf
							}
							if (vcvo.center==1) {
								eZ=eZ-eZmean
							}
						}
						shat3=shat3+quadcross(eZ,eZ)
					}
				}
			}
		}

// 1st level of clustering, no kernel-robust
// Entered unless 1-level clustering and kernel-robust
		if (!((vcvo.kernel~="") & (vcvo.clustvarname2==""))) {
			shat=J(L*K,L*K,0)
			for (i=1; i<=N_clust; i++) {		// loop through clusters, adding Z'ee'Z
												// for indiv cluster in each loop
				esub=panelsubmatrix(*vcvo.e,i,info)
				Zsub=panelsubmatrix(*vcvo.Z,i,info)
				wsub=panelsubmatrix(*vcvo.wvar,i,info)
				if ((K==1) & (vcvo.center~=1)) {	// simple/fast where e is a column vector
													// and no centering
					wv = esub :* wsub * vcvo.wf
					eZ = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
				}
				else {
					eZ=J(1,L*K,0)
					for (j=1; j<=rows(esub); j++) {
						eZj=(esub[j,.]#Zsub[j,.])*wsub[j,.]*vcvo.wf
							if (vcvo.center==1) {
								eZj=eZj-eZmean
							}
						eZ=eZ+eZj
					}
				}
				shat=shat+quadcross(eZ,eZ)
			}	// end loop through clusters
		}

// 2-level clustering, no kernel-robust
		if ((vcvo.clustvarname2~="") & (vcvo.kernel=="")) {
			imax=max(clustvar2)					// clustvar2 is numbered 1..N_clust2
			shat2=J(L*K,L*K,0)
			for (i=1; i<=imax; i++) {			// loop through clusters, adding Z'ee'Z
												// for indiv cluster in each loop
				svar=(clustvar2:==i)			// mimics panelsubmatrix but doesn't require sorted data
				esub=select(*vcvo.e,svar)		// it is, however, noticably slower.
				Zsub=select(*vcvo.Z,svar)
				wsub=select(*vcvo.wvar,svar)
				if ((K==1) & (vcvo.center~=1)) {	// simple/fast where e is a column vector
													// and no centering
					wv = esub :* wsub * vcvo.wf
					eZ = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
				}
				else {
					eZ=J(1,L*K,0)
					for (j=1; j<=rows(esub); j++) {
						eZ=eZ+(esub[j,.]#Zsub[j,.])*wsub[j,.]*vcvo.wf
					}
				}
				if (vcvo.center==1) {
					eZ=eZ-eZmean
				}
				shat2=shat2+quadcross(eZ,eZ)
			}
		}

// 1st level of cluster, kernel-robust OR
// 2-level clustering, kernel-robust and time is 2nd cluster variable
		if (vcvo.kernel~="") {
			shat2=J(L*K,L*K,0)
// First, standard cluster-robust, i.e., no lags.
			i=min(t)
			while (i<=max(t)) {  				// loop through all T clusters, adding Z'ee'Z
												// for indiv cluster in each loop
				eZ=J(1,L*K,0)
				svar=(t:==i)					// select obs with t=i
				if (colsum(svar)>0) {			// there are obs with t=i
					esub=select(*vcvo.e,svar)
					Zsub=select(*vcvo.Z,svar)
					wsub=select(*vcvo.wvar,svar)
					if ((K==1) & (vcvo.center~=1)) {	// simple/fast where e is a column vector
														// and no centering
						wv = esub :* wsub * vcvo.wf
						eZ = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
					}
					else {
						eZ=J(1,L*K,0)
						for (j=1; j<=rows(esub); j++) {
							eZ=eZ+(esub[j,.]#Zsub[j,.])*wsub[j,.]*vcvo.wf
						}
					}
					if (vcvo.center==1) {
						eZ=eZ-eZmean
					}
					shat2=shat2+quadcross(eZ,eZ)
				}
				i=i+vcvo.tdelta
			} // end i loop through all T clusters

// Spectral windows require looping through all T-1 autocovariances
			if (window=="spectral") {
				TAU=T/vcvo.tdelta-1
			}
			else {
				TAU=vcvo.bw
			}

			for (tau=1; tau<=TAU; tau++) {
				kw = m_calckw(tau, vcvo.bw, vcvo.kernel)	// zero weight possible with some kernels
															// save an unnecessary loop if kw=0
															// remember, kw<0 possible with some kernels!
				if (kw~=0) {
					i=min(t)+tau*vcvo.tdelta				// Loop through all possible ts (time clusters)
					while (i<=max(t)) {						// Start at earliest possible t
						svar=t[.,]:==i						// svar is current, svar1 is tau-th lag
						svar1=t[.,]:==(i-tau*vcvo.tdelta)	// tau*vcvo.tdelta is usually just tau
						if ((colsum(svar)>0)				// there are current & lagged obs
								& (colsum(svar1)>0))	 {
							wv  = select((*vcvo.e),svar)  :* select((*vcvo.wvar),svar)  * vcvo.wf
							wv1 = select((*vcvo.e),svar1) :* select((*vcvo.wvar),svar1) * vcvo.wf
							Zsub =select((*vcvo.Z),svar)
							Zsub1=select((*vcvo.Z),svar1)
						if ((K==1) & (vcvo.center~=1)) {			// simple/fast where e is a column vector
																	// and no centering
								eZsub = quadcross(1, wv, Zsub)		// equivalent to colsum(wv :* Zsub)
								eZsub1= quadcross(1, wv1, Zsub1)	// equivalent to colsum(wv :* Zsub)
							}
							else {
								eZsub=J(1,L*K,0)
								for (j=1; j<=rows(Zsub); j++) {
									wvj =wv[j,.]
									Zj  =Zsub[j,.]
									eZsub=eZsub+(wvj#Zj)
								}
								eZsub1=J(1,L*K,0)
								for (j=1; j<=rows(Zsub1); j++) {
									wv1j =wv1[j,.]
									Z1j  =Zsub1[j,.]
									eZsub1=eZsub1+(wv1j#Z1j)
								}
							}
							if (vcvo.center==1) {
								eZsub=eZsub-eZmean
								eZsub1=eZsub1-eZmean
							}
							ghat=quadcross(eZsub,eZsub1)
							shat2=shat2+kw*(ghat+ghat')
						}
						i=i+vcvo.tdelta
					}
				}	// end non-zero kernel weight block
			}	// end tau loop

// If 1-level clustering, shat2 just calculated above is actually the desired shat
			if (vcvo.clustvarname2=="") {
				shat=shat2
			}
		}

// 2-level clustering, completion
// Cameron-Gelbach-Miller/Thompson method:
// Add 2 cluster variance matrices and subtract 3rd
		if (vcvo.clustvarname2~="") {
			shat = shat+shat2-shat3
		}		

// Note no dof correction required for cluster-robust
	shat=shat/vcvo.N
	} // end cluster-robust code

	if (vcvo.sw~="") {
// Stock-Watson adjustment.  Calculate Bhat in their equation (6).  Also need T=panel length.
// They define for balanced panels.  Since T is not constant for unbalanced panels, need
// to incorporate panel-varying 1/T, 1/(T-1) and 1/(T-2) as weights in summation.

		st_view(ivar, ., st_tsrevar(vcvo.ivarname), vcvo.touse)
		info_ivar = panelsetup(ivar, 1)

		shat=J(L*K,L*K,0)
		bhat=J(L*K,L*K,0)
		N_panels=0
		for (i=1; i<=rows(info_ivar); i++) {
			esub=panelsubmatrix(*vcvo.e,i,info_ivar)
			Zsub=panelsubmatrix(*vcvo.Z,i,info_ivar)
			wsub=panelsubmatrix(*vcvo.wvar,i,info_ivar)
			Tsub=rows(esub)
			if (Tsub>2) {			// SW cov estimator defined only for T>2
				N_panels=N_panels+1
				sigmahatsub=J(K,K,0)
				ZZsub=J(L*K,L*K,0)
				shatsub=J(L*K,L*K,0)
				for (j=1; j<=rows(esub); j++) {
					eZi=esub[j,1]#Zsub[j,.]
					if (vcvo.center==1) {
						eZi=eZi-eZmean
					}
					if ((vcvo.weight=="fweight") | (vcvo.weight=="iweight")) {
						shatsub=shatsub+quadcross(eZi,eZi)*wsub[j]*vcvo.wf
						sigmahatsub=sigmahatsub + quadcross(esub[j,1],esub[j,1])*wsub[j]*vcvo.wf
						ZZsub=ZZsub+quadcross(Zsub[j,.],Zsub[j,.])*wsub[j]*vcvo.wf
					}
					else {
						shatsub=shatsub+quadcross(eZi,eZi)*((wsub[j]*vcvo.wf)^2)
						sigmahatsub=sigmahatsub + quadcross(esub[j,1],esub[j,1])*((wsub[j]*vcvo.wf)^2)
						ZZsub=ZZsub+quadcross(Zsub[j,.],Zsub[j,.])*((wsub[j]*vcvo.wf)^2)
					}
				} // end loop through j obs of panel i
				shat=shat + shatsub*(Tsub-1)/(Tsub-2)
				bhat=bhat + ZZsub/Tsub#sigmahatsub/(Tsub-1)/(Tsub-2)
			}
		} // end loop through i panels

// Note that Stock-Watson incorporate an N-n-k degrees of freedom correction in their eqn 4
// for what we call shat.  We use only an N-n degrees of freedom correction, i.e., we ignore
// the k regressors.  This is because this is an estimate of S, the VCV of orthogonality conditions,
// independently of its use to obtain an estimate of the variance of beta.  Makes no diff aysmptotically.
// Ignore dofminus correction since this is explicitly handled here.
// Use number of valid panels in denominator (SW cov estimator defined only for panels with T>2).
		shat=shat/(vcvo.N-N_panels)
		bhat=bhat/N_panels
		shat=shat-bhat
	} // end Stock-Watson block

	_makesymmetric(shat)

// shat may not be positive-definite.  Use spectral decomposition to obtain an invertable version.
// Extract Eigenvector and Eigenvalues, replace EVs, and reassemble shat.
// psda option: Stock-Watson 2008 Econometrica, Remark 8, say replace neg EVs with abs(EVs).
// psd0 option: Politis (2007) says replace neg EVs with zeros.
	if (vcvo.psd~="") {
		symeigensystem(shat,Evec,Eval)
		if (vcvo.psd=="psda") {
			Eval = abs(Eval)
		}
		else {
			Eval = Eval + (abs(Eval) - Eval)/2
		}
		shat = Evec*diag(Eval)*Evec'
		_makesymmetric(shat)
	}

	return(shat)

} // end of program m_omega

// *********************************************************************** //
// ********* m_calckw - shared by ivreg2 and ranktest ********************* //
// *********************************************************************** //

real scalar m_calckw(	real scalar tau,
						real scalar bw,
						string scalar kernel) 
	{
				karg = tau / bw
				if (kernel=="Truncated") {
					kw=1
				}
				if (kernel=="Bartlett") {
					kw=(1-karg)
				}
				if (kernel=="Parzen") {
					if (karg <= 0.5) {
						kw = 1-6*karg^2+6*karg^3
					}
					else {
						kw = 2*(1-karg)^3
					}
				}
				if (kernel=="Tukey-Hanning") {
					kw=0.5+0.5*cos(pi()*karg)
				}
				if (kernel=="Tukey-Hamming") {
					kw=0.54+0.46*cos(pi()*karg)
				}
				if (kernel=="Tent") {
					kw=2*(1-cos(tau*karg)) / (karg^2)
				}
				if (kernel=="Danielle") {
					kw=sin(pi()*karg) / (pi()*karg)
				}
				if (kernel=="Quadratic Spectral") {
					kw=25/(12*pi()^2*karg^2) /*
						*/ * ( sin(6*pi()*karg/5)/(6*pi()*karg/5) /*
						*/     - cos(6*pi()*karg/5) )
				}
				return(kw)
	}  // end kw

// *********************************************************************** //
// ********* END CODE SHARED BY ivreg2 AND ranktest ******************** //
// *********************************************************************** //

// cdsy: used by ivreg2

void s_cdsy( string scalar temp, scalar choice)
{
string scalar s_ivbias5, s_ivbias10, s_ivbias20, s_ivbias30
string scalar s_ivsize10, s_ivsize15, s_ivsize20, s_ivsize25
string scalar s_fullrel5, s_fullrel10, s_fullrel20, s_fullrel30
string scalar s_fullmax5, s_fullmax10, s_fullmax20, s_fullmax30
string scalar s_limlsize10, s_limlsize15, s_limlsize20, s_limlsize25

s_ivbias5  =  
". , . , . \  . , . , . \  13.91 , . , . \  16.85 , 11.04 , . \  18.37 , 13.97 , 9.53  \ 19.28 , 15.72 , 12.20 \  19.86 , 16.88 , 13.95 \  20.25 , 17.70 , 15.18 \  20.53 , 18.30 , 16.10 \  20.74 , 18.76 , 16.80 \  20.90 , 19.12 , 17.35 \  21.01 , 19.40 , 17.80 \  21.10 , 19.64 , 18.17 \  21.18 , 19.83 , 18.47 \  21.23 , 19.98 , 18.73 \  21.28 , 20.12 , 18.94 \  21.31 , 20.23 , 19.13 \  21.34 , 20.33 , 19.29 \  21.36 , 20.41 , 19.44 \  21.38 , 20.48 , 19.56 \  21.39 , 20.54 , 19.67 \  21.40 , 20.60 , 19.77 \  21.41 , 20.65 , 19.86 \  21.41 , 20.69 , 19.94 \  21.42 , 20.73 , 20.01 \  21.42 , 20.76 , 20.07 \  21.42 , 20.79 , 20.13 \  21.42 , 20.82 , 20.18 \  21.42 , 20.84 , 20.23 \  21.42 , 20.86 , 20.27 \  21.41 , 20.88 , 20.31 \  21.41 , 20.90 , 20.35 \  21.41 , 20.91 , 20.38 \  21.40 , 20.93 , 20.41 \  21.40 , 20.94 , 20.44 \  21.39 , 20.95 , 20.47 \  21.39 , 20.96 , 20.49 \  21.38 , 20.97 , 20.51 \  21.38 , 20.98 , 20.54 \  21.37 , 20.99 , 20.56 \  21.37 , 20.99 , 20.57 \  21.36 , 21.00 , 20.59 \  21.35 , 21.00 , 20.61 \  21.35 , 21.01 , 20.62 \  21.34 , 21.01 , 20.64 \  21.34 , 21.02 , 20.65 \  21.33 , 21.02 , 20.66 \  21.32 , 21.02 , 20.67 \  21.32 , 21.03 , 20.68 \  21.31 , 21.03 , 20.69 \  21.31 , 21.03 , 20.70 \  21.30 , 21.03 , 20.71 \  21.30 , 21.03 , 20.72 \  21.29 , 21.03 , 20.73 \  21.28 , 21.03 , 20.73 \  21.28 , 21.04 , 20.74 \  21.27 , 21.04 , 20.75 \  21.27 , 21.04 , 20.75 \  21.26 , 21.04 , 20.76 \  21.26 , 21.04 , 20.76 \  21.25 , 21.04 , 20.77 \  21.24 , 21.04 , 20.77 \  21.24 , 21.04 , 20.78 \  21.23 , 21.04 , 20.78 \  21.23 , 21.03 , 20.79 \  21.22 , 21.03 , 20.79 \  21.22 , 21.03 , 20.79 \  21.21 , 21.03 , 20.80 \  21.21 , 21.03 , 20.80 \  21.20 , 21.03 , 20.80 \  21.20 , 21.03 , 20.80 \  21.19 , 21.03 , 20.81 \  21.19 , 21.03 , 20.81 \  21.18 , 21.03 , 20.81 \  21.18 , 21.02 , 20.81 \  21.17 , 21.02 , 20.82 \  21.17 , 21.02 , 20.82 \  21.16 , 21.02 , 20.82 \  21.16 , 21.02 , 20.82 \  21.15 , 21.02 , 20.82 \  21.15 , 21.02 , 20.82 \  21.15 , 21.02 , 20.83 \  21.14 , 21.01 , 20.83 \  21.14 , 21.01 , 20.83 \  21.13 , 21.01 , 20.83 \  21.13 , 21.01 , 20.83 \  21.12 , 21.01 , 20.84 \  21.12 , 21.01 , 20.84 \  21.11 , 21.01 , 20.84 \  21.11 , 21.01 , 20.84 \  21.10 , 21.00 , 20.84 \  21.10 , 21.00 , 20.84 \  21.09 , 21.00 , 20.85 \  21.09 , 21.00 , 20.85 \  21.08 , 21.00 , 20.85 \  21.08 , 21.00 , 20.85 \  21.07 , 21.00 , 20.85 \  21.07 , 20.99 , 20.86 \  21.06 , 20.99 , 20.86 \  21.06 , 20.99 , 20.86 \" 
ivbias5 = strtoreal(colshape(colshape(tokens(s_ivbias5), 2)[.,1], 3))

s_ivbias10 = 
". , . , .			\	 	 	. , . , .			\	 	 	9.08 , . , .		\	 	 	10.27 , 7.56 , .		\	 	 	10.83 , 8.78 , 6.61		\	 	 	11.12 , 9.48 , 7.77		\	 	 	11.29 , 9.92 , 8.5		\	 	 	11.39 , 10.22 , 9.01	\	 	 	11.46 , 10.43 , 9.37	\	 	 	11.49 , 10.58 , 9.64	\	 	 	11.51 , 10.69 , 9.85	\	 	 	11.52 , 10.78 , 10.01	\	 	 	11.52 , 10.84 , 10.14	\	 	 	11.52 , 10.89 , 10.25	\	 	 	11.51 , 10.93 , 10.33	\	 	 	11.5 , 10.96 , 10.41	\	 	 	11.49 , 10.99 , 10.47	\	 	 	11.48 , 11 , 10.52		\	 	 	11.46 , 11.02 , 10.56	\	 	 	11.45 , 11.03 , 10.6	\	 	 	11.44 , 11.04 , 10.63	\	 	 	11.42 , 11.05 , 10.65	\	 	 	11.41 , 11.05 , 10.68	\	 	 	11.4 , 11.05 , 10.7		\	 	 	11.38 , 11.06 , 10.71	\	 	 	11.37 , 11.06 , 10.73	\	 	 	11.36 , 11.06 , 10.74	\	 	 	11.34 , 11.05 , 10.75	\	 	 	11.33 , 11.05 , 10.76	\	 	 	11.32 , 11.05 , 10.77	\	 	 	11.3 , 11.05 , 10.78	\	 	 	11.29 , 11.05 , 10.79	\	 	 	11.28 , 11.04 , 10.79	\	 	 	11.27 , 11.04 , 10.8	\	 	 	11.26 , 11.04 , 10.8	\	 	 	11.25 , 11.03 , 10.8	\	 	 	11.24 , 11.03 , 10.81	\	 	 	11.23 , 11.02 , 10.81	\	 	 	11.22 , 11.02 , 10.81	\	 	 	11.21 , 11.02 , 10.81	\	 	 	11.2 , 11.01 , 10.81	\	 	 	11.19 , 11.01 , 10.81	\	 	 	11.18 , 11 , 10.81		\	 	 	11.17 , 11 , 10.81		\	 	 	11.16 , 10.99 , 10.81	\	 	 	11.15 , 10.99 , 10.81	\	 	 	11.14 , 10.98 , 10.81	\	 	 	11.13 , 10.98 , 10.81	\	 	 	11.13 , 10.98 , 10.81	\	 	 	11.12 , 10.97 , 10.81	\	 	 	11.11 , 10.97 , 10.81	\	 	 	11.1 , 10.96 , 10.81	\	 	 	11.1 , 10.96 , 10.81	\	 	 	11.09 , 10.95 , 10.81	\	 	 	11.08 , 10.95 , 10.81	\	 	 	11.07 , 10.94 , 10.8	\	 	 	11.07 , 10.94 , 10.8	\	 	 	11.06 , 10.94 , 10.8	\	 	 	11.05 , 10.93 , 10.8	\	 	 	11.05 , 10.93 , 10.8	\	 	 	11.04 , 10.92 , 10.8	\	 	 	11.03 , 10.92 , 10.79	\	 	 	11.03 , 10.92 , 10.79	\	 	 	11.02 , 10.91 , 10.79	\	 	 	11.02 , 10.91 , 10.79	\	 	 	11.01 , 10.9 , 10.79	\	 	 	11 , 10.9 , 10.79		\	 	 	11 , 10.9 , 10.78		\	 	 	10.99 , 10.89 , 10.78	\	 	 	10.99 , 10.89 , 10.78	\	 	 	10.98 , 10.89 , 10.78	\	 	 	10.98 , 10.88 , 10.78	\	 	 	10.97 , 10.88 , 10.77	\	 	 	10.97 , 10.88 , 10.77	\	 	 	10.96 , 10.87 , 10.77	\	 	 	10.96 , 10.87 , 10.77	\	 	 	10.95 , 10.86 , 10.77	\	 	 	10.95 , 10.86 , 10.76	\	 	 	10.94 , 10.86 , 10.76	\	 	 	10.94 , 10.85 , 10.76	\	 	 	10.93 , 10.85 , 10.76	\	 	 	10.93 , 10.85 , 10.76	\	 	 	10.92 , 10.84 , 10.75	\	 	 	10.92 , 10.84 , 10.75	\	 	 	10.91 , 10.84 , 10.75	\	 	 	10.91 , 10.84 , 10.75	\	 	 	10.91 , 10.83 , 10.75	\	 	 	10.9 , 10.83 , 10.74	\	 	 	10.9 , 10.83 , 10.74	\	 	 	10.89 , 10.82 , 10.74	\	 	 	10.89 , 10.82 , 10.74	\	 	 	10.89 , 10.82 , 10.74	\	 	 	10.88 , 10.81 , 10.74	\	 	 	10.88 , 10.81 , 10.73	\	 	 	10.87 , 10.81 , 10.73	\	 	 	10.87 , 10.81 , 10.73	\	 	 	10.87 , 10.8 , 10.73	\	 	 	10.86 , 10.8 , 10.73	\	 	 	10.86 , 10.8 , 10.72	\	 	 	10.86 , 10.8 , 10.72 \"
ivbias10 = strtoreal(colshape(colshape(tokens(s_ivbias10), 2)[.,1], 3))

s_ivbias20 = 
"	.	,	.	,	.	\  	.	,	.	,	.	\  	6.46	,	.	,	.	\  	6.71	,	5.57	,	.	\  	6.77	,	5.91	,	4.99	\  	6.76	,	6.08	,	5.35	\  	6.73	,	6.16	,	5.56	\  	6.69	,	6.20	,	5.69	\  	6.65	,	6.22	,	5.78	\  	6.61	,	6.23	,	5.83	\  	6.56	,	6.23	,	5.87	\  	6.53	,	6.22	,	5.90	\  	6.49	,	6.21	,	5.92	\  	6.45	,	6.20	,	5.93	\  	6.42	,	6.19	,	5.94	\  	6.39	,	6.17	,	5.94	\  	6.36	,	6.16	,	5.94	\  	6.33	,	6.14	,	5.94	\  	6.31	,	6.13	,	5.94	\  	6.28	,	6.11	,	5.93	\  	6.26	,	6.10	,	5.93	\  	6.24	,	6.08	,	5.92	\  	6.22	,	6.07	,	5.92	\  	6.20	,	6.06	,	5.91	\  	6.18	,	6.05	,	5.90	\  	6.16	,	6.03	,	5.90	\  	6.14	,	6.02	,	5.89	\  	6.13	,	6.01	,	5.88	\  	6.11	,	6.00	,	5.88	\  	6.09	,	5.99	,	5.87	\  	6.08	,	5.98	,	5.87	\  	6.07	,	5.97	,	5.86	\  	6.05	,	5.96	,	5.85	\  	6.04	,	5.95	,	5.85	\  	6.03	,	5.94	,	5.84	\  	6.01	,	5.93	,	5.83	\  	6.00	,	5.92	,	5.83	\  	5.99	,	5.91	,	5.82	\  	5.98	,	5.90	,	5.82	\  	5.97	,	5.89	,	5.81	\  	5.96	,	5.89	,	5.80	\  	5.95	,	5.88	,	5.80	\  	5.94	,	5.87	,	5.79	\  	5.93	,	5.86	,	5.79	\  	5.92	,	5.86	,	5.78	\  	5.91	,	5.85	,	5.78	\  	5.91	,	5.84	,	5.77	\  	5.90	,	5.83	,	5.77	\  	5.89	,	5.83	,	5.76	\  	5.88	,	5.82	,	5.76	\  	5.87	,	5.82	,	5.75	\  	5.87	,	5.81	,	5.75	\  	5.86	,	5.80	,	5.74	\  	5.85	,	5.80	,	5.74	\  	5.85	,	5.79	,	5.73	\  	5.84	,	5.79	,	5.73	\  	5.83	,	5.78	,	5.72	\  	5.83	,	5.78	,	5.72	\  	5.82	,	5.77	,	5.72	\  	5.81	,	5.77	,	5.71	\  	5.81	,	5.76	,	5.71	\  	5.80	,	5.76	,	5.70	\  	5.80	,	5.75	,	5.70	\  	5.79	,	5.75	,	5.70	\  	5.78	,	5.74	,	5.69	\  	5.78	,	5.74	,	5.69	\  	5.77	,	5.73	,	5.68	\  	5.77	,	5.73	,	5.68	\  	5.76	,	5.72	,	5.68	\  	5.76	,	5.72	,	5.67	\  	5.75	,	5.72	,	5.67	\  	5.75	,	5.71	,	5.67	\  	5.75	,	5.71	,	5.66	\  	5.74	,	5.70	,	5.66	\  	5.74	,	5.70	,	5.66	\  	5.73	,	5.70	,	5.65	\  	5.73	,	5.69	,	5.65	\  	5.72	,	5.69	,	5.65	\  	5.72	,	5.68	,	5.65	\  	5.71	,	5.68	,	5.64	\  	5.71	,	5.68	,	5.64	\  	5.71	,	5.67	,	5.64	\  	5.70	,	5.67	,	5.63	\  	5.70	,	5.67	,	5.63	\  	5.70	,	5.66	,	5.63	\  	5.69	,	5.66	,	5.62	\  	5.69	,	5.66	,	5.62	\  	5.68	,	5.65	,	5.62	\  	5.68	,	5.65	,	5.62	\  	5.68	,	5.65	,	5.61	\  	5.67	,	5.65	,	5.61	\  	5.67	,	5.64	,	5.61	\  	5.67	,	5.64	,	5.61	\  	5.66	,	5.64	,	5.60	\  	5.66	,	5.63	,	5.60	\  	5.66	,	5.63	,	5.60	\  	5.65	,	5.63	,	5.60	\  	5.65	,	5.63	,	5.59	\  	5.65	,	5.62	,	5.59	\  	5.65	,	5.62	,	5.59	\"
ivbias20 = strtoreal(colshape(colshape(tokens(s_ivbias20), 2)[.,1], 3))

s_ivbias30 = 
"	.	,	.	,	.	\   	.	,	.	,	.	\   	5.39	,	.	,	.	\   	5.34	,	4.73	,	.	\   	5.25	,	4.79	,	4.30	\   	5.15	,	4.78	,	4.40	\   	5.07	,	4.76	,	4.44	\   	4.99	,	4.73	,	4.46	\   	4.92	,	4.69	,	4.46	\   	4.86	,	4.66	,	4.45	\   	4.80	,	4.62	,	4.44	\   	4.75	,	4.59	,	4.42	\   	4.71	,	4.56	,	4.41	\   	4.67	,	4.53	,	4.39	\   	4.63	,	4.50	,	4.37	\   	4.59	,	4.48	,	4.36	\   	4.56	,	4.45	,	4.34	\   	4.53	,	4.43	,	4.32	\   	4.51	,	4.41	,	4.31	\   	4.48	,	4.39	,	4.29	\   	4.46	,	4.37	,	4.28	\   	4.43	,	4.35	,	4.27	\   	4.41	,	4.33	,	4.25	\   	4.39	,	4.32	,	4.24	\   	4.37	,	4.30	,	4.23	\   	4.35	,	4.29	,	4.21	\   	4.34	,	4.27	,	4.20	\   	4.32	,	4.26	,	4.19	\   	4.31	,	4.24	,	4.18	\   	4.29	,	4.23	,	4.17	\   	4.28	,	4.22	,	4.16	\   	4.26	,	4.21	,	4.15	\   	4.25	,	4.20	,	4.14	\   	4.24	,	4.19	,	4.13	\   	4.23	,	4.18	,	4.13	\   	4.22	,	4.17	,	4.12	\   	4.20	,	4.16	,	4.11	\   	4.19	,	4.15	,	4.10	\   	4.18	,	4.14	,	4.09	\   	4.17	,	4.13	,	4.09	\   	4.16	,	4.12	,	4.08	\   	4.15	,	4.11	,	4.07	\   	4.15	,	4.11	,	4.07	\   	4.14	,	4.10	,	4.06	\   	4.13	,	4.09	,	4.05	\   	4.12	,	4.08	,	4.05	\   	4.11	,	4.08	,	4.04	\   	4.11	,	4.07	,	4.03	\   	4.10	,	4.06	,	4.03	\   	4.09	,	4.06	,	4.02	\   	4.08	,	4.05	,	4.02	\   	4.08	,	4.05	,	4.01	\   	4.07	,	4.04	,	4.01	\   	4.06	,	4.03	,	4.00	\   	4.06	,	4.03	,	4.00	\   	4.05	,	4.02	,	3.99	\   	4.05	,	4.02	,	3.99	\   	4.04	,	4.01	,	3.98	\   	4.04	,	4.01	,	3.98	\   	4.03	,	4.00	,	3.97	\   	4.02	,	4.00	,	3.97	\   	4.02	,	3.99	,	3.96	\   	4.01	,	3.99	,	3.96	\   	4.01	,	3.98	,	3.96	\   	4.00	,	3.98	,	3.95	\   	4.00	,	3.97	,	3.95	\   	3.99	,	3.97	,	3.94	\   	3.99	,	3.97	,	3.94	\   	3.99	,	3.96	,	3.94	\   	3.98	,	3.96	,	3.93	\   	3.98	,	3.95	,	3.93	\   	3.97	,	3.95	,	3.93	\   	3.97	,	3.95	,	3.92	\   	3.96	,	3.94	,	3.92	\   	3.96	,	3.94	,	3.92	\   	3.96	,	3.93	,	3.91	\   	3.95	,	3.93	,	3.91	\   	3.95	,	3.93	,	3.91	\   	3.95	,	3.92	,	3.90	\   	3.94	,	3.92	,	3.90	\   	3.94	,	3.92	,	3.90	\   	3.93	,	3.91	,	3.89	\   	3.93	,	3.91	,	3.89	\   	3.93	,	3.91	,	3.89	\   	3.92	,	3.91	,	3.89	\   	3.92	,	3.90	,	3.88	\   	3.92	,	3.90	,	3.88	\   	3.91	,	3.90	,	3.88	\   	3.91	,	3.89	,	3.87	\   	3.91	,	3.89	,	3.87	\   	3.91	,	3.89	,	3.87	\   	3.90	,	3.89	,	3.87	\   	3.90	,	3.88	,	3.86	\   	3.90	,	3.88	,	3.86	\   	3.89	,	3.88	,	3.86	\   	3.89	,	3.87	,	3.86	\   	3.89	,	3.87	,	3.85	\   	3.89	,	3.87	,	3.85	\   	3.88	,	3.87	,	3.85	\   	3.88	,	3.86	,	3.85	\"
ivbias30 = strtoreal(colshape(colshape(tokens(s_ivbias30), 2)[.,1], 3))


s_ivsize10 = 
"16.38 , .	\	  	19.93 , 7.03	\	  	22.3 , 13.43	\	  	24.58 , 16.87	\	  	26.87 , 19.45	\	  	29.18 , 21.68	\	  	31.5 , 23.72	\	  	33.84 , 25.64	\	  	36.19 , 27.51	\	  	38.54 , 29.32	\	  	40.9 , 31.11	\	  	43.27 , 32.88	\	  	45.64 , 34.62	\	  	48.01 , 36.36	\	  	50.39 , 38.08	\	  	52.77 , 39.8	\	  	55.15 , 41.51	\	  	57.53 , 43.22	\	  	59.92 , 44.92	\	  	62.3 , 46.62	\	  	64.69 , 48.31	\	  	67.07 , 50.01	\	  	69.46 , 51.7	\	  	71.85 , 53.39	\	  	74.24 , 55.07	\	  	76.62 , 56.76	\	  	79.01 , 58.45	\	  	81.4 , 60.13	\	  	83.79 , 61.82	\	  	86.17 , 63.51	\	  	88.56 , 65.19	\	  	90.95 , 66.88	\	  	93.33 , 68.56	\	  	95.72 , 70.25	\	  	98.11 , 71.94	\	  	100.5 , 73.62	\	  	102.88 , 75.31	\	  	105.27 , 76.99	\	  	107.66 , 78.68	\	  	110.04 , 80.37	\	  	112.43 , 82.05	\	  	114.82 , 83.74	\	  	117.21 , 85.42	\	  	119.59 , 87.11	\	  	121.98 , 88.8	\	  	124.37 , 90.48	\	  	126.75 , 92.17	\	  	129.14 , 93.85	\	  	131.53 , 95.54	\	  	133.92 , 97.23	\	  	136.3 , 98.91	\	  	138.69 , 100.6	\	  	141.08 , 102.29	\	  	143.47 , 103.97	\	  	145.85 , 105.66	\	  	148.24 , 107.34	\	  	150.63 , 109.03	\	  	153.01 , 110.72	\	  	155.4 , 112.4	\	  	157.79 , 114.09	\	  	160.18 , 115.77	\	  	162.56 , 117.46	\	  	164.95 , 119.15	\	  	167.34 , 120.83	\	  	169.72 , 122.52	\	  	172.11 , 124.2	\	  	174.5 , 125.89	\	  	176.89 , 127.58	\	  	179.27 , 129.26	\	  	181.66 , 130.95	\	  	184.05 , 132.63	\	  	186.44 , 134.32	\	  	188.82 , 136.01	\	  	191.21 , 137.69	\	  	193.6 , 139.38	\	  	195.98 , 141.07	\	  	198.37 , 142.75	\	  	200.76 , 144.44	\	  	203.15 , 146.12	\	  	205.53 , 147.81	\	  	207.92 , 149.5	\	  	210.31 , 151.18	\	  	212.69 , 152.87	\	  	215.08 , 154.55	\	  	217.47 , 156.24	\	  	219.86 , 157.93	\	  	222.24 , 159.61	\	  	224.63 , 161.3	\	  	227.02 , 162.98	\	  	229.41 , 164.67	\	  	231.79 , 166.36	\	  	234.18 , 168.04	\	  	236.57 , 169.73	\	  	238.95 , 171.41	\	  	241.34 , 173.1	\	  	243.73 , 174.79	\	  	246.12 , 176.47	\	  	248.5 , 178.16	\	  	250.89 , 179.84	\	  	253.28 , 181.53 \"
ivsize10 = strtoreal(colshape(colshape(tokens(s_ivsize10), 2)[.,1], 2))

s_ivsize15 = 
 	"8.96	,	.	\   	11.59	,	4.58	\   	12.83	,	8.18	\   	13.96	,	9.93	\   	15.09	,	11.22	\   	16.23	,	12.33	\   	17.38	,	13.34	\   	18.54	,	14.31	\   	19.71	,	15.24	\   	20.88	,	16.16	\   	22.06	,	17.06	\   	23.24	,	17.95	\   	24.42	,	18.84	\   	25.61	,	19.72	\   	26.80	,	20.60	\   	27.99	,	21.48	\   	29.19	,	22.35	\   	30.38	,	23.22	\   	31.58	,	24.09	\   	32.77	,	24.96	\   	33.97	,	25.82	\   	35.17	,	26.69	\   	36.37	,	27.56	\   	37.57	,	28.42	\   	38.77	,	29.29	\   	39.97	,	30.15	\   	41.17	,	31.02	\   	42.37	,	31.88	\   	43.57	,	32.74	\   	44.78	,	33.61	\   	45.98	,	34.47	\   	47.18	,	35.33	\   	48.38	,	36.19	\   	49.59	,	37.06	\   	50.79	,	37.92	\   	51.99	,	38.78	\   	53.19	,	39.64	\   	54.40	,	40.50	\   	55.60	,	41.37	\   	56.80	,	42.23	\   	58.01	,	43.09	\   	59.21	,	43.95	\   	60.41	,	44.81	\   	61.61	,	45.68	\   	62.82	,	46.54	\   	64.02	,	47.40	\   	65.22	,	48.26	\   	66.42	,	49.12	\   	67.63	,	49.99	\   	68.83	,	50.85	\   	70.03	,	51.71	\   	71.24	,	52.57	\   	72.44	,	53.43	\   	73.64	,	54.30	\   	74.84	,	55.16	\   	76.05	,	56.02	\   	77.25	,	56.88	\   	78.45	,	57.74	\   	79.66	,	58.61	\   	80.86	,	59.47	\   	82.06	,	60.33	\   	83.26	,	61.19	\   	84.47	,	62.05	\   	85.67	,	62.92	\   	86.87	,	63.78	\   	88.07	,	64.64	\   	89.28	,	65.50	\   	90.48	,	66.36	\   	91.68	,	67.22	\   	92.89	,	68.09	\   	94.09	,	68.95	\   	95.29	,	69.81	\   	96.49	,	70.67	\   	97.70	,	71.53	\   	98.90	,	72.40	\   	100.10	,	73.26	\   	101.30	,	74.12	\   	102.51	,	74.98	\   	103.71	,	75.84	\   	104.91	,	76.71	\   	106.12	,	77.57	\   	107.32	,	78.43	\   	108.52	,	79.29	\   	109.72	,	80.15	\   	110.93	,	81.02	\   	112.13	,	81.88	\   	113.33	,	82.74	\   	114.53	,	83.60	\   	115.74	,	84.46	\   	116.94	,	85.33	\   	118.14	,	86.19	\   	119.35	,	87.05	\   	120.55	,	87.91	\   	121.75	,	88.77	\   	122.95	,	89.64	\   	124.16	,	90.50	\   	125.36	,	91.36	\   	126.56	,	92.22	\   	127.76	,	93.08	\   	128.97	,	93.95	\"
ivsize15 = strtoreal(colshape(colshape(tokens(s_ivsize15), 2)[.,1], 2))

s_ivsize20 = 
 "	6.66	,	.	\   	8.75	,	3.95	\   	9.54	,	6.40	\   	10.26	,	7.54	\   	10.98	,	8.38	\   	11.72	,	9.10	\   	12.48	,	9.77	\   	13.24	,	10.41	\   	14.01	,	11.03	\   	14.78	,	11.65	\   	15.56	,	12.25	\   	16.35	,	12.86	\   	17.14	,	13.45	\   	17.93	,	14.05	\   	18.72	,	14.65	\   	19.51	,	15.24	\   	20.31	,	15.83	\   	21.10	,	16.42	\   	21.90	,	17.02	\   	22.70	,	17.61	\   	23.50	,	18.20	\   	24.30	,	18.79	\   	25.10	,	19.38	\   	25.90	,	19.97	\   	26.71	,	20.56	\   	27.51	,	21.15	\   	28.31	,	21.74	\   	29.12	,	22.33	\   	29.92	,	22.92	\   	30.72	,	23.51	\   	31.53	,	24.10	\   	32.33	,	24.69	\   	33.14	,	25.28	\   	33.94	,	25.87	\   	34.75	,	26.46	\   	35.55	,	27.05	\   	36.36	,	27.64	\   	37.17	,	28.23	\   	37.97	,	28.82	\   	38.78	,	29.41	\   	39.58	,	30.00	\   	40.39	,	30.59	\   	41.20	,	31.18	\   	42.00	,	31.77	\   	42.81	,	32.36	\   	43.62	,	32.95	\   	44.42	,	33.54	\   	45.23	,	34.13	\   	46.03	,	34.72	\   	46.84	,	35.31	\   	47.65	,	35.90	\   	48.45	,	36.49	\   	49.26	,	37.08	\   	50.06	,	37.67	\   	50.87	,	38.26	\   	51.68	,	38.85	\   	52.48	,	39.44	\   	53.29	,	40.02	\   	54.09	,	40.61	\   	54.90	,	41.20	\   	55.71	,	41.79	\   	56.51	,	42.38	\   	57.32	,	42.97	\   	58.13	,	43.56	\   	58.93	,	44.15	\   	59.74	,	44.74	\   	60.54	,	45.33	\   	61.35	,	45.92	\   	62.16	,	46.51	\   	62.96	,	47.10	\   	63.77	,	47.69	\   	64.57	,	48.28	\   	65.38	,	48.87	\   	66.19	,	49.46	\   	66.99	,	50.05	\   	67.80	,	50.64	\   	68.60	,	51.23	\   	69.41	,	51.82	\   	70.22	,	52.41	\   	71.02	,	53.00	\   	71.83	,	53.59	\   	72.64	,	54.18	\   	73.44	,	54.77	\   	74.25	,	55.36	\   	75.05	,	55.95	\   	75.86	,	56.54	\   	76.67	,	57.13	\   	77.47	,	57.72	\   	78.28	,	58.31	\   	79.08	,	58.90	\   	79.89	,	59.49	\   	80.70	,	60.08	\   	81.50	,	60.67	\   	82.31	,	61.26	\   	83.12	,	61.85	\   	83.92	,	62.44	\   	84.73	,	63.03	\   	85.53	,	63.62	\   	86.34	,	64.21	\   	87.15	,	64.80	\"
ivsize20 = strtoreal(colshape(colshape(tokens(s_ivsize20), 2)[.,1], 2))

s_ivsize25 = 
 "	5.53	,	.	\   	7.25	,	3.63	\   	7.80	,	5.45	\   	8.31	,	6.28	\   	8.84	,	6.89	\   	9.38	,	7.42	\   	9.93	,	7.91	\   	10.50	,	8.39	\   	11.07	,	8.85	\   	11.65	,	9.31	\   	12.23	,	9.77	\   	12.82	,	10.22	\   	13.41	,	10.68	\   	14.00	,	11.13	\   	14.60	,	11.58	\   	15.19	,	12.03	\   	15.79	,	12.49	\   	16.39	,	12.94	\   	16.99	,	13.39	\   	17.60	,	13.84	\   	18.20	,	14.29	\   	18.80	,	14.74	\   	19.41	,	15.19	\   	20.01	,	15.64	\   	20.61	,	16.10	\   	21.22	,	16.55	\   	21.83	,	17.00	\   	22.43	,	17.45	\   	23.04	,	17.90	\   	23.65	,	18.35	\   	24.25	,	18.81	\   	24.86	,	19.26	\   	25.47	,	19.71	\   	26.08	,	20.16	\   	26.68	,	20.61	\   	27.29	,	21.06	\   	27.90	,	21.52	\   	28.51	,	21.97	\   	29.12	,	22.42	\   	29.73	,	22.87	\   	30.33	,	23.32	\   	30.94	,	23.78	\   	31.55	,	24.23	\   	32.16	,	24.68	\   	32.77	,	25.13	\   	33.38	,	25.58	\   	33.99	,	26.04	\   	34.60	,	26.49	\   	35.21	,	26.94	\   	35.82	,	27.39	\   	36.43	,	27.85	\   	37.04	,	28.30	\   	37.65	,	28.75	\   	38.25	,	29.20	\   	38.86	,	29.66	\   	39.47	,	30.11	\   	40.08	,	30.56	\   	40.69	,	31.01	\   	41.30	,	31.47	\   	41.91	,	31.92	\   	42.52	,	32.37	\   	43.13	,	32.82	\   	43.74	,	33.27	\   	44.35	,	33.73	\   	44.96	,	34.18	\   	45.57	,	34.63	\   	46.18	,	35.08	\   	46.78	,	35.54	\   	47.39	,	35.99	\   	48.00	,	36.44	\   	48.61	,	36.89	\   	49.22	,	37.35	\   	49.83	,	37.80	\   	50.44	,	38.25	\   	51.05	,	38.70	\   	51.66	,	39.16	\   	52.27	,	39.61	\   	52.88	,	40.06	\   	53.49	,	40.51	\   	54.10	,	40.96	\   	54.71	,	41.42	\   	55.32	,	41.87	\   	55.92	,	42.32	\   	56.53	,	42.77	\   	57.14	,	43.23	\   	57.75	,	43.68	\   	58.36	,	44.13	\   	58.97	,	44.58	\   	59.58	,	45.04	\   	60.19	,	45.49	\   	60.80	,	45.94	\   	61.41	,	46.39	\   	62.02	,	46.85	\   	62.63	,	47.30	\   	63.24	,	47.75	\   	63.85	,	48.20	\   	64.45	,	48.65	\   	65.06	,	49.11	\   	65.67	,	49.56	\   	66.28	,	50.01	\"
ivsize25 = strtoreal(colshape(colshape(tokens(s_ivsize25), 2)[.,1], 2))


s_fullrel5 = 
" 	24.09	,	.	\   	13.46	,	15.50	\   	9.61	,	10.83	\   	7.63	,	8.53	\   	6.42	,	7.16	\   	5.61	,	6.24	\   	5.02	,	5.59	\   	4.58	,	5.10	\   	4.23	,	4.71	\   	3.96	,	4.41	\   	3.73	,	4.15	\   	3.54	,	3.94	\   	3.38	,	3.76	\   	3.24	,	3.60	\   	3.12	,	3.47	\   	3.01	,	3.35	\   	2.92	,	3.24	\   	2.84	,	3.15	\   	2.76	,	3.06	\   	2.69	,	2.98	\   	2.63	,	2.91	\   	2.58	,	2.85	\   	2.52	,	2.79	\   	2.48	,	2.73	\   	2.43	,	2.68	\   	2.39	,	2.63	\   	2.36	,	2.59	\   	2.32	,	2.55	\   	2.29	,	2.51	\   	2.26	,	2.47	\   	2.23	,	2.44	\   	2.20	,	2.41	\   	2.18	,	2.37	\   	2.16	,	2.35	\   	2.13	,	2.32	\   	2.11	,	2.29	\   	2.09	,	2.27	\   	2.07	,	2.24	\   	2.05	,	2.22	\   	2.04	,	2.20	\   	2.02	,	2.18	\   	2.00	,	2.16	\   	1.99	,	2.14	\   	1.97	,	2.12	\   	1.96	,	2.10	\   	1.94	,	2.09	\   	1.93	,	2.07	\   	1.92	,	2.05	\   	1.91	,	2.04	\   	1.89	,	2.02	\   	1.88	,	2.01	\   	1.87	,	2.00	\   	1.86	,	1.98	\   	1.85	,	1.97	\   	1.84	,	1.96	\   	1.83	,	1.95	\   	1.82	,	1.94	\   	1.81	,	1.92	\   	1.80	,	1.91	\   	1.79	,	1.90	\   	1.79	,	1.89	\   	1.78	,	1.88	\   	1.77	,	1.87	\   	1.76	,	1.87	\   	1.75	,	1.86	\   	1.75	,	1.85	\   	1.74	,	1.84	\   	1.73	,	1.83	\   	1.72	,	1.83	\   	1.72	,	1.82	\   	1.71	,	1.81	\   	1.70	,	1.80	\   	1.70	,	1.80	\   	1.69	,	1.79	\   	1.68	,	1.79	\   	1.68	,	1.78	\   	1.67	,	1.77	\   	1.67	,	1.77	\   	1.66	,	1.76	\   	1.65	,	1.76	\   	1.65	,	1.75	\   	1.64	,	1.75	\   	1.64	,	1.74	\   	1.63	,	1.74	\   	1.63	,	1.73	\   	1.62	,	1.73	\   	1.61	,	1.73	\   	1.61	,	1.72	\   	1.60	,	1.72	\   	1.60	,	1.71	\   	1.59	,	1.71	\   	1.59	,	1.71	\   	1.58	,	1.71	\   	1.58	,	1.70	\   	1.57	,	1.70	\   	1.57	,	1.70	\   	1.56	,	1.69	\   	1.56	,	1.69	\   	1.55	,	1.69	\   	1.55	,	1.69	)"
fullrel5 = strtoreal(colshape(colshape(tokens(s_fullrel5), 2)[.,1], 2))

s_fullrel10 = 
 "	19.36	,	.	\   	10.89	,	12.55	\   	7.90	,	8.96	\   	6.37	,	7.15	\   	5.44	,	6.07	\   	4.81	,	5.34	\   	4.35	,	4.82	\   	4.01	,	4.43	\   	3.74	,	4.12	\   	3.52	,	3.87	\   	3.34	,	3.67	\   	3.19	,	3.49	\   	3.06	,	3.35	\   	2.95	,	3.22	\   	2.85	,	3.11	\   	2.76	,	3.01	\   	2.69	,	2.92	\   	2.62	,	2.84	\   	2.56	,	2.77	\   	2.50	,	2.71	\   	2.45	,	2.65	\   	2.40	,	2.60	\   	2.36	,	2.55	\   	2.32	,	2.50	\   	2.28	,	2.46	\   	2.24	,	2.42	\   	2.21	,	2.38	\   	2.18	,	2.35	\   	2.15	,	2.31	\   	2.12	,	2.28	\   	2.10	,	2.25	\   	2.07	,	2.23	\   	2.05	,	2.20	\   	2.03	,	2.17	\   	2.01	,	2.15	\   	1.99	,	2.13	\   	1.97	,	2.11	\   	1.95	,	2.09	\   	1.93	,	2.07	\   	1.92	,	2.05	\   	1.90	,	2.03	\   	1.88	,	2.01	\   	1.87	,	2.00	\   	1.86	,	1.98	\   	1.84	,	1.96	\   	1.83	,	1.95	\   	1.82	,	1.93	\   	1.81	,	1.92	\   	1.79	,	1.91	\   	1.78	,	1.89	\   	1.77	,	1.88	\   	1.76	,	1.87	\   	1.75	,	1.86	\   	1.74	,	1.85	\   	1.73	,	1.84	\   	1.72	,	1.83	\   	1.71	,	1.82	\   	1.70	,	1.81	\   	1.70	,	1.80	\   	1.69	,	1.79	\   	1.68	,	1.78	\   	1.67	,	1.77	\   	1.67	,	1.76	\   	1.66	,	1.75	\   	1.65	,	1.75	\   	1.64	,	1.74	\   	1.64	,	1.73	\   	1.63	,	1.72	\   	1.63	,	1.72	\   	1.62	,	1.71	\   	1.61	,	1.70	\   	1.61	,	1.70	\   	1.60	,	1.69	\   	1.60	,	1.68	\   	1.59	,	1.68	\   	1.59	,	1.67	\   	1.58	,	1.67	\   	1.58	,	1.66	\   	1.57	,	1.66	\   	1.57	,	1.65	\   	1.56	,	1.65	\   	1.56	,	1.64	\   	1.56	,	1.64	\   	1.55	,	1.63	\   	1.55	,	1.63	\   	1.54	,	1.62	\   	1.54	,	1.62	\   	1.54	,	1.62	\   	1.53	,	1.61	\   	1.53	,	1.61	\   	1.53	,	1.61	\   	1.52	,	1.60	\   	1.52	,	1.60	\   	1.52	,	1.60	\   	1.52	,	1.59	\   	1.51	,	1.59	\   	1.51	,	1.59	\   	1.51	,	1.59	\   	1.51	,	1.58	\   	1.50	,	1.58	)"
fullrel10 = strtoreal(colshape(colshape(tokens(s_fullrel10), 2)[.,1], 2))

s_fullrel20 = 
" 	15.64	,	.	\   	9.00	,	9.72	\   	6.61	,	7.18	\   	5.38	,	5.85	\   	4.62	,	5.04	\   	4.11	,	4.48	\   	3.75	,	4.08	\   	3.47	,	3.77	\   	3.25	,	3.53	\   	3.07	,	3.33	\   	2.92	,	3.17	\   	2.80	,	3.04	\   	2.70	,	2.92	\   	2.61	,	2.82	\   	2.53	,	2.73	\   	2.46	,	2.65	\   	2.39	,	2.58	\   	2.34	,	2.52	\   	2.29	,	2.46	\   	2.24	,	2.41	\   	2.20	,	2.36	\   	2.16	,	2.32	\   	2.13	,	2.28	\   	2.10	,	2.24	\   	2.06	,	2.21	\   	2.04	,	2.18	\   	2.01	,	2.15	\   	1.99	,	2.12	\   	1.96	,	2.09	\   	1.94	,	2.07	\   	1.92	,	2.04	\   	1.90	,	2.02	\   	1.88	,	2.00	\   	1.87	,	1.98	\   	1.85	,	1.96	\   	1.83	,	1.94	\   	1.82	,	1.93	\   	1.80	,	1.91	\   	1.79	,	1.89	\   	1.78	,	1.88	\   	1.76	,	1.86	\   	1.75	,	1.85	\   	1.74	,	1.84	\   	1.73	,	1.82	\   	1.72	,	1.81	\   	1.71	,	1.80	\   	1.70	,	1.79	\   	1.69	,	1.78	\   	1.68	,	1.77	\   	1.67	,	1.76	\   	1.66	,	1.75	\   	1.65	,	1.74	\   	1.65	,	1.73	\   	1.64	,	1.72	\   	1.63	,	1.71	\   	1.62	,	1.70	\   	1.62	,	1.69	\   	1.61	,	1.68	\   	1.60	,	1.68	\   	1.60	,	1.67	\   	1.59	,	1.66	\   	1.58	,	1.65	\   	1.58	,	1.65	\   	1.57	,	1.64	\   	1.57	,	1.63	\   	1.56	,	1.63	\   	1.56	,	1.62	\   	1.55	,	1.62	\   	1.55	,	1.61	\   	1.54	,	1.60	\   	1.54	,	1.60	\   	1.53	,	1.59	\   	1.53	,	1.59	\   	1.52	,	1.58	\   	1.52	,	1.58	\   	1.51	,	1.57	\   	1.51	,	1.57	\   	1.51	,	1.56	\   	1.50	,	1.56	\   	1.50	,	1.56	\   	1.49	,	1.55	\   	1.49	,	1.55	\   	1.49	,	1.54	\   	1.48	,	1.54	\   	1.48	,	1.54	\   	1.48	,	1.53	\   	1.47	,	1.53	\   	1.47	,	1.53	\   	1.47	,	1.52	\   	1.46	,	1.52	\   	1.46	,	1.52	\   	1.46	,	1.51	\   	1.46	,	1.51	\   	1.45	,	1.51	\   	1.45	,	1.50	\   	1.45	,	1.50	\   	1.45	,	1.50	\   	1.44	,	1.50	\   	1.44	,	1.49	\   	1.44	,	1.49	)"
fullrel20 = strtoreal(colshape(colshape(tokens(s_fullrel20), 2)[.,1], 2))

s_fullrel30 = 
 "	12.71	,	.	\   	7.49	,	8.03	\   	5.60	,	6.15	\   	4.63	,	5.10	\   	4.03	,	4.44	\   	3.63	,	3.98	\   	3.33	,	3.65	\   	3.11	,	3.39	\   	2.93	,	3.19	\   	2.79	,	3.02	\   	2.67	,	2.88	\   	2.57	,	2.77	\   	2.48	,	2.67	\   	2.41	,	2.58	\   	2.34	,	2.51	\   	2.28	,	2.44	\   	2.23	,	2.38	\   	2.18	,	2.33	\   	2.14	,	2.28	\   	2.10	,	2.23	\   	2.07	,	2.19	\   	2.04	,	2.16	\   	2.01	,	2.12	\   	1.98	,	2.09	\   	1.95	,	2.06	\   	1.93	,	2.03	\   	1.90	,	2.01	\   	1.88	,	1.98	\   	1.86	,	1.96	\   	1.84	,	1.94	\   	1.83	,	1.92	\   	1.81	,	1.90	\   	1.79	,	1.88	\   	1.78	,	1.87	\   	1.76	,	1.85	\   	1.75	,	1.83	\   	1.74	,	1.82	\   	1.72	,	1.80	\   	1.71	,	1.79	\   	1.70	,	1.78	\   	1.69	,	1.77	\   	1.68	,	1.75	\   	1.67	,	1.74	\   	1.66	,	1.73	\   	1.65	,	1.72	\   	1.64	,	1.71	\   	1.63	,	1.70	\   	1.62	,	1.69	\   	1.61	,	1.68	\   	1.60	,	1.67	\   	1.60	,	1.66	\   	1.59	,	1.66	\   	1.58	,	1.65	\   	1.57	,	1.64	\   	1.57	,	1.63	\   	1.56	,	1.63	\   	1.55	,	1.62	\   	1.55	,	1.61	\   	1.54	,	1.61	\   	1.54	,	1.60	\   	1.53	,	1.59	\   	1.53	,	1.59	\   	1.52	,	1.58	\   	1.51	,	1.57	\   	1.51	,	1.57	\   	1.50	,	1.56	\   	1.50	,	1.56	\   	1.50	,	1.55	\   	1.49	,	1.55	\   	1.49	,	1.54	\   	1.48	,	1.54	\   	1.48	,	1.53	\   	1.47	,	1.53	\   	1.47	,	1.52	\   	1.47	,	1.52	\   	1.46	,	1.52	\   	1.46	,	1.51	\   	1.46	,	1.51	\   	1.45	,	1.50	\   	1.45	,	1.50	\   	1.45	,	1.50	\   	1.44	,	1.49	\   	1.44	,	1.49	\   	1.44	,	1.48	\   	1.43	,	1.48	\   	1.43	,	1.48	\   	1.43	,	1.47	\   	1.43	,	1.47	\   	1.42	,	1.47	\   	1.42	,	1.47	\   	1.42	,	1.46	\   	1.42	,	1.46	\   	1.41	,	1.46	\   	1.41	,	1.45	\   	1.41	,	1.45	\   	1.41	,	1.45	\   	1.41	,	1.45	\   	1.40	,	1.44	\   	1.40	,	1.44	\   	1.40	,	1.44	\"
fullrel30 = strtoreal(colshape(colshape(tokens(s_fullrel30), 2)[.,1], 2))


s_fullmax5 = 
 "	23.81	,	.	\   	12.38	,	14.19	\   	8.66	,	10.00	\   	6.81	,	7.88	\   	5.71	,	6.60	\   	4.98	,	5.74	\   	4.45	,	5.13	\   	4.06	,	4.66	\   	3.76	,	4.30	\   	3.51	,	4.01	\   	3.31	,	3.77	\   	3.15	,	3.57	\   	3.00	,	3.41	\   	2.88	,	3.26	\   	2.78	,	3.13	\   	2.69	,	3.02	\   	2.61	,	2.92	\   	2.53	,	2.84	\   	2.47	,	2.76	\   	2.41	,	2.69	\   	2.36	,	2.62	\   	2.31	,	2.56	\   	2.27	,	2.51	\   	2.23	,	2.46	\   	2.19	,	2.42	\   	2.15	,	2.37	\   	2.12	,	2.33	\   	2.09	,	2.30	\   	2.07	,	2.26	\   	2.04	,	2.23	\   	2.02	,	2.20	\   	1.99	,	2.17	\   	1.97	,	2.14	\   	1.95	,	2.12	\   	1.93	,	2.10	\   	1.91	,	2.07	\   	1.90	,	2.05	\   	1.88	,	2.03	\   	1.87	,	2.01	\   	1.85	,	1.99	\   	1.84	,	1.98	\   	1.82	,	1.96	\   	1.81	,	1.94	\   	1.80	,	1.93	\   	1.79	,	1.91	\   	1.78	,	1.90	\   	1.76	,	1.88	\   	1.75	,	1.87	\   	1.74	,	1.86	\   	1.73	,	1.85	\   	1.73	,	1.83	\   	1.72	,	1.82	\   	1.71	,	1.81	\   	1.70	,	1.80	\   	1.69	,	1.79	\   	1.68	,	1.78	\   	1.68	,	1.77	\   	1.67	,	1.76	\   	1.66	,	1.75	\   	1.65	,	1.74	\   	1.65	,	1.74	\   	1.64	,	1.73	\   	1.63	,	1.72	\   	1.63	,	1.71	\   	1.62	,	1.70	\   	1.62	,	1.70	\   	1.61	,	1.69	\   	1.60	,	1.68	\   	1.60	,	1.68	\   	1.59	,	1.67	\   	1.59	,	1.66	\   	1.58	,	1.66	\   	1.58	,	1.65	\   	1.57	,	1.64	\   	1.57	,	1.64	\   	1.56	,	1.63	\   	1.56	,	1.63	\   	1.55	,	1.62	\   	1.55	,	1.62	\   	1.54	,	1.61	\   	1.54	,	1.61	\   	1.53	,	1.60	\   	1.53	,	1.60	\   	1.53	,	1.59	\   	1.52	,	1.59	\   	1.52	,	1.58	\   	1.51	,	1.58	\   	1.51	,	1.57	\   	1.50	,	1.57	\   	1.50	,	1.57	\   	1.50	,	1.56	\   	1.49	,	1.56	\   	1.49	,	1.55	\   	1.49	,	1.55	\   	1.48	,	1.55	\   	1.48	,	1.54	\   	1.47	,	1.54	\   	1.47	,	1.54	\   	1.47	,	1.53	\   	1.46	,	1.53	)"
fullmax5 = strtoreal(colshape(colshape(tokens(s_fullmax5), 2)[.,1], 2))

s_fullmax10 = 
"	19.40	,	.	\   	10.14	,	11.92	\   	7.18	,	8.39	\   	5.72	,	6.64	\   	4.85	,	5.60	\   	4.27	,	4.90	\   	3.86	,	4.40	\   	3.55	,	4.03	\   	3.31	,	3.73	\   	3.12	,	3.50	\   	2.96	,	3.31	\   	2.83	,	3.15	\   	2.71	,	3.01	\   	2.62	,	2.89	\   	2.53	,	2.79	\   	2.46	,	2.70	\   	2.39	,	2.62	\   	2.33	,	2.55	\   	2.28	,	2.49	\   	2.23	,	2.43	\   	2.19	,	2.38	\   	2.15	,	2.33	\   	2.11	,	2.29	\   	2.08	,	2.25	\   	2.05	,	2.21	\   	2.02	,	2.18	\   	1.99	,	2.14	\   	1.97	,	2.11	\   	1.94	,	2.08	\   	1.92	,	2.06	\   	1.90	,	2.03	\   	1.88	,	2.01	\   	1.86	,	1.99	\   	1.85	,	1.97	\   	1.83	,	1.95	\   	1.81	,	1.93	\   	1.80	,	1.91	\   	1.79	,	1.89	\   	1.77	,	1.88	\   	1.76	,	1.86	\   	1.75	,	1.85	\   	1.74	,	1.83	\   	1.72	,	1.82	\   	1.71	,	1.81	\   	1.70	,	1.80	\   	1.69	,	1.78	\   	1.68	,	1.77	\   	1.67	,	1.76	\   	1.66	,	1.75	\   	1.66	,	1.74	\   	1.65	,	1.73	\   	1.64	,	1.72	\   	1.63	,	1.71	\   	1.62	,	1.70	\   	1.62	,	1.69	\   	1.61	,	1.69	\   	1.60	,	1.68	\   	1.60	,	1.67	\   	1.59	,	1.66	\   	1.58	,	1.65	\   	1.58	,	1.65	\   	1.57	,	1.64	\   	1.57	,	1.63	\   	1.56	,	1.63	\   	1.55	,	1.62	\   	1.55	,	1.61	\   	1.54	,	1.61	\   	1.54	,	1.60	\   	1.53	,	1.60	\   	1.53	,	1.59	\   	1.52	,	1.59	\   	1.52	,	1.58	\   	1.52	,	1.58	\   	1.51	,	1.57	\   	1.51	,	1.57	\   	1.50	,	1.56	\   	1.50	,	1.56	\   	1.49	,	1.55	\   	1.49	,	1.55	\   	1.49	,	1.54	\   	1.48	,	1.54	\   	1.48	,	1.53	\   	1.48	,	1.53	\   	1.47	,	1.53	\   	1.47	,	1.52	\   	1.46	,	1.52	\   	1.46	,	1.51	\   	1.46	,	1.51	\   	1.45	,	1.51	\   	1.45	,	1.50	\   	1.45	,	1.50	\   	1.44	,	1.50	\   	1.44	,	1.49	\   	1.44	,	1.49	\   	1.44	,	1.49	\   	1.43	,	1.48	\   	1.43	,	1.48	\   	1.43	,	1.48	\   	1.42	,	1.48	\   	1.42	,	1.47	)"
fullmax10 = strtoreal(colshape(colshape(tokens(s_fullmax10), 2)[.,1], 2))

s_fullmax20 =
" 	15.39	,	.	\   	8.16	,	9.41	\   	5.87	,	6.79	\   	4.75	,	5.47	\   	4.08	,	4.66	\   	3.64	,	4.13	\   	3.32	,	3.74	\   	3.08	,	3.45	\   	2.89	,	3.22	\   	2.74	,	3.03	\   	2.62	,	2.88	\   	2.51	,	2.76	\   	2.42	,	2.65	\   	2.35	,	2.56	\   	2.28	,	2.48	\   	2.22	,	2.40	\   	2.17	,	2.34	\   	2.12	,	2.28	\   	2.08	,	2.23	\   	2.04	,	2.19	\   	2.01	,	2.15	\   	1.98	,	2.11	\   	1.95	,	2.07	\   	1.92	,	2.04	\   	1.89	,	2.01	\   	1.87	,	1.98	\   	1.85	,	1.96	\   	1.83	,	1.93	\   	1.81	,	1.91	\   	1.79	,	1.89	\   	1.77	,	1.87	\   	1.76	,	1.85	\   	1.74	,	1.83	\   	1.73	,	1.82	\   	1.72	,	1.80	\   	1.70	,	1.79	\   	1.69	,	1.77	\   	1.68	,	1.76	\   	1.67	,	1.74	\   	1.66	,	1.73	\   	1.65	,	1.72	\   	1.64	,	1.71	\   	1.63	,	1.70	\   	1.62	,	1.69	\   	1.61	,	1.68	\   	1.60	,	1.67	\   	1.59	,	1.66	\   	1.58	,	1.65	\   	1.58	,	1.64	\   	1.57	,	1.63	\   	1.56	,	1.62	\   	1.56	,	1.62	\   	1.55	,	1.61	\   	1.54	,	1.60	\   	1.54	,	1.59	\   	1.53	,	1.59	\   	1.52	,	1.58	\   	1.52	,	1.57	\   	1.51	,	1.57	\   	1.51	,	1.56	\   	1.50	,	1.56	\   	1.50	,	1.55	\   	1.49	,	1.54	\   	1.49	,	1.54	\   	1.48	,	1.53	\   	1.48	,	1.53	\   	1.47	,	1.52	\   	1.47	,	1.52	\   	1.47	,	1.51	\   	1.46	,	1.51	\   	1.46	,	1.51	\   	1.45	,	1.50	\   	1.45	,	1.50	\   	1.45	,	1.49	\   	1.44	,	1.49	\   	1.44	,	1.48	\   	1.44	,	1.48	\   	1.43	,	1.48	\   	1.43	,	1.47	\   	1.43	,	1.47	\   	1.42	,	1.46	\   	1.42	,	1.46	\   	1.42	,	1.46	\   	1.41	,	1.45	\   	1.41	,	1.45	\   	1.41	,	1.45	\   	1.40	,	1.44	\   	1.40	,	1.44	\   	1.40	,	1.44	\   	1.40	,	1.44	\   	1.39	,	1.43	\   	1.39	,	1.43	\   	1.39	,	1.43	\   	1.39	,	1.42	\   	1.38	,	1.42	\   	1.38	,	1.42	\   	1.38	,	1.42	\   	1.38	,	1.41	\   	1.37	,	1.41	\   	1.37	,	1.41	)"
fullmax20 = strtoreal(colshape(colshape(tokens(s_fullmax20), 2)[.,1], 2))

s_fullmax30 =
 "	12.76	,	.	\   	6.97	,	8.01	\   	5.11	,	5.88	\   	4.19	,	4.78	\   	3.64	,	4.12	\   	3.27	,	3.67	\   	3.00	,	3.35	\   	2.80	,	3.10	\   	2.64	,	2.91	\   	2.52	,	2.76	\   	2.41	,	2.63	\   	2.33	,	2.52	\   	2.25	,	2.43	\   	2.19	,	2.35	\   	2.13	,	2.29	\   	2.08	,	2.22	\   	2.04	,	2.17	\   	2.00	,	2.12	\   	1.96	,	2.08	\   	1.93	,	2.04	\   	1.90	,	2.01	\   	1.87	,	1.97	\   	1.84	,	1.94	\   	1.82	,	1.92	\   	1.80	,	1.89	\   	1.78	,	1.87	\   	1.76	,	1.84	\   	1.74	,	1.82	\   	1.73	,	1.80	\   	1.71	,	1.79	\   	1.70	,	1.77	\   	1.68	,	1.75	\   	1.67	,	1.74	\   	1.66	,	1.72	\   	1.64	,	1.71	\   	1.63	,	1.70	\   	1.62	,	1.68	\   	1.61	,	1.67	\   	1.60	,	1.66	\   	1.59	,	1.65	\   	1.58	,	1.64	\   	1.57	,	1.63	\   	1.57	,	1.62	\   	1.56	,	1.61	\   	1.55	,	1.60	\   	1.54	,	1.59	\   	1.54	,	1.59	\   	1.53	,	1.58	\   	1.52	,	1.57	\   	1.52	,	1.56	\   	1.51	,	1.56	\   	1.50	,	1.55	\   	1.50	,	1.54	\   	1.49	,	1.54	\   	1.49	,	1.53	\   	1.48	,	1.53	\   	1.48	,	1.52	\   	1.47	,	1.51	\   	1.47	,	1.51	\   	1.46	,	1.50	\   	1.46	,	1.50	\   	1.45	,	1.49	\   	1.45	,	1.49	\   	1.44	,	1.48	\   	1.44	,	1.48	\   	1.44	,	1.47	\   	1.43	,	1.47	\   	1.43	,	1.47	\   	1.42	,	1.46	\   	1.42	,	1.46	\   	1.42	,	1.45	\   	1.41	,	1.45	\   	1.41	,	1.45	\   	1.41	,	1.44	\   	1.40	,	1.44	\   	1.40	,	1.44	\   	1.40	,	1.43	\   	1.39	,	1.43	\   	1.39	,	1.43	\   	1.39	,	1.42	\   	1.39	,	1.42	\   	1.38	,	1.42	\   	1.38	,	1.41	\   	1.38	,	1.41	\   	1.37	,	1.41	\   	1.37	,	1.40	\   	1.37	,	1.40	\   	1.37	,	1.40	\   	1.36	,	1.40	\   	1.36	,	1.39	\   	1.36	,	1.39	\   	1.36	,	1.39	\   	1.36	,	1.38	\   	1.35	,	1.38	\   	1.35	,	1.38	\   	1.35	,	1.38	\   	1.35	,	1.37	\   	1.34	,	1.37	\   	1.34	,	1.37	\   	1.34	,	1.37	)"
fullmax30 = strtoreal(colshape(colshape(tokens(s_fullmax30), 2)[.,1], 2))


s_limlsize10 = 
 "	16.38	,	.	\   	8.68	,	7.03	\   	6.46	,	5.44	\   	5.44	,	4.72	\   	4.84	,	4.32	\   	4.45	,	4.06	\   	4.18	,	3.90	\   	3.97	,	3.78	\   	3.81	,	3.70	\   	3.68	,	3.64	\   	3.58	,	3.60	\   	3.50	,	3.58	\   	3.42	,	3.56	\   	3.36	,	3.55	\   	3.31	,	3.54	\   	3.27	,	3.55	\   	3.24	,	3.55	\   	3.20	,	3.56	\   	3.18	,	3.57	\   	3.21	,	3.58	\   	3.39	,	3.59	\   	3.57	,	3.60	\   	3.68	,	3.62	\   	3.75	,	3.64	\   	3.79	,	3.65	\   	3.82	,	3.67	\   	3.85	,	3.74	\   	3.86	,	3.87	\   	3.87	,	4.02	\   	3.88	,	4.12	\   	3.89	,	4.19	\   	3.89	,	4.24	\   	3.90	,	4.27	\   	3.90	,	4.31	\   	3.90	,	4.33	\   	3.90	,	4.36	\   	3.90	,	4.38	\   	3.90	,	4.39	\   	3.90	,	4.41	\   	3.90	,	4.43	\   	3.90	,	4.44	\   	3.90	,	4.45	\   	3.90	,	4.47	\   	3.90	,	4.48	\   	3.90	,	4.49	\   	3.90	,	4.50	\   	3.90	,	4.51	\   	3.90	,	4.52	\   	3.90	,	4.53	\   	3.90	,	4.54	\   	3.90	,	4.55	\   	3.90	,	4.56	\   	3.90	,	4.56	\   	3.90	,	4.57	\   	3.90	,	4.58	\   	3.90	,	4.59	\   	3.90	,	4.59	\   	3.90	,	4.60	\   	3.90	,	4.61	\   	3.90	,	4.61	\   	3.90	,	4.62	\   	3.90	,	4.62	\   	3.90	,	4.63	\   	3.90	,	4.63	\   	3.89	,	4.64	\   	3.89	,	4.64	\   	3.89	,	4.64	\   	3.89	,	4.65	\   	3.89	,	4.65	\   	3.89	,	4.65	\   	3.89	,	4.66	\   	3.89	,	4.66	\   	3.89	,	4.66	\   	3.89	,	4.66	\   	3.88	,	4.66	\   	3.88	,	4.66	\   	3.88	,	4.66	\   	3.88	,	4.66	\   	3.88	,	4.66	\   	3.88	,	4.66	\   	3.88	,	4.66	\   	3.87	,	4.66	\   	3.87	,	4.66	\   	3.87	,	4.66	\   	3.87	,	4.66	\   	3.87	,	4.66	\   	3.86	,	4.65	\   	3.86	,	4.65	\   	3.86	,	4.65	\   	3.86	,	4.64	\   	3.85	,	4.64	\   	3.85	,	4.64	\   	3.85	,	4.63	\   	3.85	,	4.63	\   	3.84	,	4.62	\   	3.84	,	4.62	\   	3.84	,	4.61	\   	3.84	,	4.60	\   	3.83	,	4.60	\   	3.83	,	4.59	)"
limlsize10 = strtoreal(colshape(colshape(tokens(s_limlsize10), 2)[.,1], 2))

s_limlsize15 = 
 "	8.96	,	.	\   	5.33	,	4.58	\   	4.36	,	3.81	\   	3.87	,	3.39	\   	3.56	,	3.13	\   	3.34	,	2.95	\   	3.18	,	2.83	\   	3.04	,	2.73	\   	2.93	,	2.66	\   	2.84	,	2.60	\   	2.76	,	2.55	\   	2.69	,	2.52	\   	2.63	,	2.48	\   	2.57	,	2.46	\   	2.52	,	2.44	\   	2.48	,	2.42	\   	2.44	,	2.41	\   	2.41	,	2.40	\   	2.37	,	2.39	\   	2.34	,	2.38	\   	2.32	,	2.38	\   	2.29	,	2.37	\   	2.27	,	2.37	\   	2.25	,	2.37	\   	2.24	,	2.37	\   	2.22	,	2.38	\   	2.21	,	2.38	\   	2.20	,	2.38	\   	2.19	,	2.39	\   	2.18	,	2.39	\   	2.19	,	2.40	\   	2.22	,	2.41	\   	2.33	,	2.42	\   	2.40	,	2.42	\   	2.45	,	2.43	\   	2.48	,	2.44	\   	2.50	,	2.45	\   	2.52	,	2.54	\   	2.53	,	2.55	\   	2.54	,	2.66	\   	2.55	,	2.73	\   	2.56	,	2.78	\   	2.57	,	2.82	\   	2.57	,	2.85	\   	2.58	,	2.87	\   	2.58	,	2.89	\   	2.58	,	2.91	\   	2.59	,	2.92	\   	2.59	,	2.93	\   	2.59	,	2.94	\   	2.59	,	2.95	\   	2.59	,	2.96	\   	2.60	,	2.97	\   	2.60	,	2.98	\   	2.60	,	2.98	\   	2.60	,	2.99	\   	2.60	,	2.99	\   	2.60	,	3.00	\   	2.60	,	3.00	\   	2.60	,	3.01	\   	2.60	,	3.01	\   	2.60	,	3.02	\   	2.61	,	3.02	\   	2.61	,	3.02	\   	2.61	,	3.03	\   	2.61	,	3.03	\   	2.61	,	3.03	\   	2.61	,	3.03	\   	2.61	,	3.04	\   	2.61	,	3.04	\   	2.61	,	3.04	\   	2.60	,	3.04	\   	2.60	,	3.04	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.60	,	3.05	\   	2.59	,	3.05	\   	2.59	,	3.05	\   	2.59	,	3.05	\   	2.59	,	3.05	\   	2.59	,	3.05	\   	2.59	,	3.04	\   	2.58	,	3.04	\   	2.58	,	3.04	\   	2.58	,	3.04	\   	2.58	,	3.04	\   	2.58	,	3.03	\   	2.57	,	3.03	\   	2.57	,	3.03	\   	2.57	,	3.03	\   	2.57	,	3.02	\   	2.56	,	3.02	\   	2.56	,	3.02	)"
limlsize15 = strtoreal(colshape(colshape(tokens(s_limlsize15), 2)[.,1], 2))

s_limlsize20 = 
 "	6.66	,	.	\   	4.42	,	3.95	\   	3.69	,	3.32	\   	3.30	,	2.99	\   	3.05	,	2.78	\   	2.87	,	2.63	\   	2.73	,	2.52	\   	2.63	,	2.43	\   	2.54	,	2.36	\   	2.46	,	2.30	\   	2.40	,	2.25	\   	2.34	,	2.21	\   	2.29	,	2.17	\   	2.25	,	2.14	\   	2.21	,	2.11	\   	2.18	,	2.09	\   	2.14	,	2.07	\   	2.11	,	2.05	\   	2.09	,	2.03	\   	2.06	,	2.02	\   	2.04	,	2.01	\   	2.02	,	1.99	\   	2.00	,	1.98	\   	1.98	,	1.98	\   	1.96	,	1.97	\   	1.95	,	1.96	\   	1.93	,	1.96	\   	1.92	,	1.95	\   	1.90	,	1.95	\   	1.89	,	1.95	\   	1.88	,	1.94	\   	1.87	,	1.94	\   	1.86	,	1.94	\   	1.85	,	1.94	\   	1.84	,	1.94	\   	1.83	,	1.94	\   	1.82	,	1.94	\   	1.81	,	1.95	\   	1.81	,	1.95	\   	1.80	,	1.95	\   	1.79	,	1.95	\   	1.79	,	1.96	\   	1.78	,	1.96	\   	1.78	,	1.97	\   	1.80	,	1.97	\   	1.87	,	1.98	\   	1.92	,	1.98	\   	1.95	,	1.99	\   	1.97	,	2.00	\   	1.99	,	2.00	\   	2.00	,	2.01	\   	2.01	,	2.09	\   	2.02	,	2.11	\   	2.03	,	2.18	\   	2.04	,	2.23	\   	2.04	,	2.27	\   	2.05	,	2.29	\   	2.05	,	2.31	\   	2.06	,	2.33	\   	2.06	,	2.34	\   	2.07	,	2.35	\   	2.07	,	2.36	\   	2.07	,	2.37	\   	2.08	,	2.38	\   	2.08	,	2.39	\   	2.08	,	2.39	\   	2.08	,	2.40	\   	2.09	,	2.40	\   	2.09	,	2.41	\   	2.09	,	2.41	\   	2.09	,	2.41	\   	2.09	,	2.42	\   	2.09	,	2.42	\   	2.09	,	2.42	\   	2.09	,	2.43	\   	2.10	,	2.43	\   	2.10	,	2.43	\   	2.10	,	2.43	\   	2.10	,	2.44	\   	2.10	,	2.44	\   	2.10	,	2.44	\   	2.10	,	2.44	\   	2.10	,	2.44	\   	2.09	,	2.44	\   	2.09	,	2.44	\   	2.09	,	2.45	\   	2.09	,	2.45	\   	2.09	,	2.45	\   	2.09	,	2.45	\   	2.09	,	2.45	\   	2.09	,	2.45	\   	2.09	,	2.45	\   	2.08	,	2.45	\   	2.08	,	2.45	\   	2.08	,	2.45	\   	2.08	,	2.45	\   	2.08	,	2.45	\   	2.07	,	2.44	\   	2.07	,	2.44	\   	2.07	,	2.44	)"
limlsize20 = strtoreal(colshape(colshape(tokens(s_limlsize20), 2)[.,1], 2))

s_limlsize25 = 
 "	5.53	,	.	\   	3.92	,	3.63	\   	3.32	,	3.09	\   	2.98	,	2.79	\   	2.77	,	2.60	\   	2.61	,	2.46	\   	2.49	,	2.35	\   	2.39	,	2.27	\   	2.32	,	2.20	\   	2.25	,	2.14	\   	2.19	,	2.09	\   	2.14	,	2.05	\   	2.10	,	2.02	\   	2.06	,	1.99	\   	2.03	,	1.96	\   	2.00	,	1.93	\   	1.97	,	1.91	\   	1.94	,	1.89	\   	1.92	,	1.87	\   	1.90	,	1.86	\   	1.88	,	1.84	\   	1.86	,	1.83	\   	1.84	,	1.81	\   	1.83	,	1.80	\   	1.81	,	1.79	\   	1.80	,	1.78	\   	1.78	,	1.77	\   	1.77	,	1.77	\   	1.76	,	1.76	\   	1.75	,	1.75	\   	1.74	,	1.75	\   	1.73	,	1.74	\   	1.72	,	1.73	\   	1.71	,	1.73	\   	1.70	,	1.73	\   	1.69	,	1.72	\   	1.68	,	1.72	\   	1.67	,	1.71	\   	1.67	,	1.71	\   	1.66	,	1.71	\   	1.65	,	1.71	\   	1.65	,	1.71	\   	1.64	,	1.70	\   	1.63	,	1.70	\   	1.63	,	1.70	\   	1.62	,	1.70	\   	1.62	,	1.70	\   	1.61	,	1.70	\   	1.61	,	1.70	\   	1.61	,	1.70	\   	1.60	,	1.70	\   	1.60	,	1.70	\   	1.59	,	1.70	\   	1.59	,	1.70	\   	1.59	,	1.70	\   	1.58	,	1.70	\   	1.58	,	1.71	\   	1.58	,	1.71	\   	1.57	,	1.71	\   	1.59	,	1.71	\   	1.60	,	1.71	\   	1.63	,	1.72	\   	1.65	,	1.72	\   	1.67	,	1.72	\   	1.69	,	1.72	\   	1.70	,	1.76	\   	1.71	,	1.81	\   	1.72	,	1.87	\   	1.73	,	1.91	\   	1.74	,	1.94	\   	1.74	,	1.96	\   	1.75	,	1.98	\   	1.75	,	1.99	\   	1.76	,	2.01	\   	1.76	,	2.02	\   	1.77	,	2.03	\   	1.77	,	2.04	\   	1.78	,	2.04	\   	1.78	,	2.05	\   	1.78	,	2.06	\   	1.79	,	2.06	\   	1.79	,	2.07	\   	1.79	,	2.07	\   	1.79	,	2.08	\   	1.80	,	2.08	\   	1.80	,	2.09	\   	1.80	,	2.09	\   	1.80	,	2.09	\   	1.80	,	2.09	\   	1.80	,	2.10	\   	1.80	,	2.10	\   	1.80	,	2.10	\   	1.80	,	2.10	\   	1.80	,	2.10	\   	1.80	,	2.11	\   	1.80	,	2.11	\   	1.80	,	2.11	\   	1.80	,	2.11	\   	1.80	,	2.11	\   	1.80	,	2.11	)"
limlsize25 = strtoreal(colshape(colshape(tokens(s_limlsize25), 2)[.,1], 2))

if        (choice == 1) {
	st_matrix(temp, ivbias5)
} else if (choice == 2) {
	st_matrix(temp, ivbias10)
} else if (choice == 3) {
	st_matrix(temp, ivbias20)
} else if (choice == 4) {
	st_matrix(temp, ivbias30)
} else if (choice == 5) {
	st_matrix(temp, ivsize10)
} else if (choice == 6) {
	st_matrix(temp, ivsize15)
} else if (choice == 7) {
	st_matrix(temp, ivsize20)
} else if (choice == 8) {
	st_matrix(temp, ivsize25)
} else if (choice == 9) {
	st_matrix(temp, fullrel5)
} else if (choice == 10) {
	st_matrix(temp, fullrel10)
} else if (choice == 11) {
	st_matrix(temp, fullrel20)
} else if (choice == 12) {
	st_matrix(temp, fullrel30)
} else if (choice == 13) {
	st_matrix(temp, fullmax5)
} else if (choice == 14) {
	st_matrix(temp, fullmax10)
} else if (choice == 15) {
	st_matrix(temp, fullmax20)
} else if (choice == 16) {
	st_matrix(temp, fullmax30)
} else if (choice == 17) {
	st_matrix(temp, limlsize10)
} else if (choice == 18) {
	st_matrix(temp, limlsize15)
} else if (choice == 19) {
	st_matrix(temp, limlsize20)
} else if (choice == 20) {
	st_matrix(temp, limlsize25)
}
} // end of program cdsy

mata mlib create livreg2, dir(PERSONAL) replace
mata mlib add livreg2 *()
mata mlib index
mata describe using livreg2
end
