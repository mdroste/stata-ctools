"""Executable differential cases: options are exercised on three data families.

Reference-only translations are explicit (e.g. qreg's nested VCE syntax).
Tuning options must preserve reference results; they are not acceptance smoke tests.
"""

import hashlib
import json
import re
import tempfile
from pathlib import Path

from validation import ROOT
from validation.stata import launch_driver

CASES = []


def case(command, name, options, reference, candidate, oracle="data", **kwargs):
    CASES.append(
        dict(
            command=command,
            name=name,
            options=options.split(),
            reference=reference,
            candidate=candidate,
            oracle=oracle,
            **kwargs,
        )
    )


def estimate(command, ref, lhs, options, reference_options=None, names=None, **kwargs):
    for i, opt in enumerate(options):
        ropt = reference_options[i] if reference_options is not None else opt
        case(
            command,
            f'{command}: {opt or "default"}',
            (names or [""] * len(options))[i],
            f"{ref} {lhs}, {ropt}",
            f"{command} {lhs}, {opt}",
            "estimate",
            **kwargs,
        )


for alg in ["auto", "lsd", "msd", "timsort", "sample", "counting", "merge", "ips4o"]:
    for tuning in ["", "stream(2)", "nostream", "nosortedby"]:
        case(
            "csort",
            f"sort {alg} {tuning}",
            "algorithm threads verbose " + tuning.split("(")[0],
            "sort g s, stable",
            f"csort g s, algorithm({alg}) {tuning} threads(2) verbose",
        )
for qualifier in ["", "if id<451", "in 20/480"]:
    for opt in ["", "noextend"]:
        case(
            "cencode",
            f"encode {qualifier} {opt}",
            "generate label noextend threads verbose",
            f"encode s {qualifier}, gen(out) label(encoded) {opt}",
            f"cencode s {qualifier}, gen(out) label(encoded) {opt} threads(2) verbose",
            setup="encode s, gen(__labelseed) label(encoded)\ndrop __labelseed",
        )
    for opt in ["", "maxlength(3)"]:
        case(
            "cdecode",
            f"decode {qualifier} {opt}",
            "generate maxlength threads verbose",
            f"decode group {qualifier}, gen(out) {opt}",
            f"cdecode group {qualifier}, gen(out) {opt} threads(2) verbose",
        )
    for opt in ["", "float", "force", "ignore(x) force", "percent force", "dpcomma force"]:
        case(
            "cdestring",
            f"destring {qualifier} {opt}",
            "generate float force ignore percent dpcomma threads verbose",
            f"coverage_destring_reference number {qualifier}, gen(out) {opt}",
            f"cdestring number {qualifier}, gen(out) {opt} threads(2) verbose",
        )
    for opt in ["", "epolate", "by(group)", "by(group) epolate"]:
        case(
            "cipolate",
            f"interpolation {qualifier} {opt}",
            "generate epolate by threads verbose",
            f"ipolate y time {qualifier}, gen(out) {opt}",
            f"cipolate y time {qualifier}, gen(out) {opt} threads(2) verbose",
            sortby="id",
        )
    for opt in [
        'parse(",")',
        'parse(",") notrim',
        'parse(",") limit(2)',
        'parse(",") destring force float',
        'parse(",") destring force percent',
        'parse(",") destring ignore("%x")',
    ]:
        case(
            "csplit",
            f"split {qualifier} {opt}",
            "generate parse notrim limit destring force float percent ignore threads verbose",
            f"split tokens {qualifier}, gen(out) {opt}",
            f"csplit tokens {qualifier}, gen(out) {opt} threads(2) verbose",
        )
    for opt in [
        "cuts(5 95) suffix(_w)",
        "cuts(1 99) trim suffix(_w)",
        "cuts(5 95) by(group) suffix(_w)",
        "cuts(5 95) replace",
    ]:
        case(
            "cwinsor",
            f"winsor {qualifier} {opt}",
            "cuts suffix trim by replace threads verbose",
            f"winsor2 y {qualifier}, {opt}",
            f"cwinsor y {qualifier}, {opt} threads(2) verbose",
        )
for opt in ["p(5) q(95) prefix(w_)", "p(5) q(95) trim prefix(w_)"]:
    case(
        "cwinsor",
        "winsor aliases " + opt,
        "p q prefix",
        "winsor2 y, cuts(5 95) suffix(_w)" + (" trim" if "trim" in opt else ""),
        "cwinsor y, " + opt,
        ref_after="rename y_w w_y",
    )

# Full VCE comparison for every density kernel and bandwidth. These are native
# qreg options inside vce(), not unsupported methods to omit from comparison.
for den in [
    "fitted",
    "residual",
    "kernel",
    "kernel(epanechnikov)",
    "kernel(epan2)",
    "kernel(biweight)",
    "kernel(cosine)",
    "kernel(gaussian)",
    "kernel(parzen)",
    "kernel(rectangle)",
    "kernel(triangle)",
]:
    for bw in ["hsheather", "bofinger", "chamberlain"]:
        case(
            "cqreg",
            f"qreg {den} {bw}",
            "denmethod bwmethod vce quantile",
            f"qreg y x z, quantile(.4) vce(iid, {den} {bw})",
            f"cqreg y x z, quantile(.4) denmethod({den}) bwmethod({bw})",
            "estimate",
            predict="xb stdp",
            scalars="df_r q sum_adev sum_rdev q_v f_r",
        )
estimate(
    "cqreg",
    "qreg",
    "y x z",
    [
        "quantile(.2)",
        "quantile(.8) vce(robust)",
        "nopreprocess(1) threads(2) verbose tolerance(1e-12) maxiter(500)",
    ],
    ["quantile(.2)", "quantile(.8) vce(robust)", ""],
    ["quantile", "quantile vce", "nopreprocess threads verbose tolerance maxiter"],
    predict="xb stdp",
)

olsopts = [
    "absorb(g)",
    "absorb(g h)",
    "absorb(g h) vce(robust)",
    "absorb(g h) vce(cluster g)",
    "absorb(g h) dofadjustments(none)",
    "absorb(g h) tolerance(1e-12) iterate(20000)",
    "absorb(g h) nostandardize",
    "absorb(g h) quad",
    "absorb(g h) threads(2) verbose",
]
olsrefs = [
    o.replace("nostandardize", "").replace("quad", "").replace("threads(2) verbose", "")
    for o in olsopts
]
estimate(
    "creghdfe",
    "reghdfe",
    "y x z",
    olsopts,
    olsrefs,
    [
        "absorb",
        "absorb",
        "vce",
        "vce",
        "dofadjustments",
        "tolerance iterate",
        "nostandardize",
        "quad",
        "threads verbose",
    ],
    predict="xb stdp",
    scalars="df_r df_m r2 rss",
)
for weight in ["aw=w", "fw=fw", "pw=w"]:
    estimate(
        "creghdfe",
        "reghdfe",
        f"y x z [{weight}]",
        ["absorb(g h) vce(robust)"],
        names=["absorb vce"],
        predict="xb",
    )
for residual in ["resid", "resid2(saved_resid)", "residuals(saved_resid)"]:
    estimate(
        "creghdfe",
        "reghdfe",
        "y x z",
        [f"absorb(g h) {residual}"],
        ["absorb(g h) residuals(saved_resid)"],
        [residual.split("(")[0]],
        predict="xb residuals d xbd dresiduals",
    )

iv_options = [
    "absorb(g)",
    "absorb(g) robust",
    "absorb(g) cluster(g)",
    "absorb(g) vce(cluster g)",
    "absorb(g) first",
    "absorb(g) ffirst",
    "absorb(g) rf",
    "absorb(g) small",
    "absorb(g) liml",
    "absorb(g) fuller(1)",
    "absorb(g) kclass(1)",
    "absorb(g) gmm2s robust",
    "absorb(g) cue robust",
    "absorb(g) coviv",
    "absorb(g) bw(3) kernel(bartlett)",
    "absorb(g) dkraay(3)",
    "absorb(g) kiefer",
    "absorb(g) nocons",
    "absorb(g) level(90)",
    "absorb(g) noheader nofooter nooutput",
    "absorb(g) title(Coverage) subtitle(Reference) depname(outcome)",
    "absorb(g) noid",
    "absorb(g) orthog(z2)",
    "absorb(g) endogtest(endog)",
    "absorb(g) redundant(z2)",
    "absorb(g) partial(x)",
    "absorb(g) fwl(x)",
    "absorb(g) dofminus(1)",
    "absorb(g) sdofminus(1)",
    "absorb(g) partial(x) nopartialsmall",
    "absorb(g) eform(Effect)",
    "absorb(g) plus noomitted vsquish noemptycells allbaselevels",
    "absorb(g) omitted baselevels",
    "absorb(g) savefirst savefprefix(cov_)",
    "absorb(g) saverf saverfprefix(covrf_)",
    "absorb(g) residuals(saved_resid)",
    "absorb(g) residuals2",
    "absorb(g) tolerance(1e-12) maxiter(2000)",
    "absorb(g) threads(2) verbose",
]
iv_names = [
    "absorb",
    "robust",
    "cluster",
    "vce",
    "first",
    "ffirst",
    "rf",
    "small",
    "liml",
    "fuller",
    "kclass",
    "gmm2s robust",
    "cue robust",
    "coviv",
    "bw kernel",
    "dkraay",
    "kiefer",
    "noconstant",
    "level",
    "noheader nofooter nooutput",
    "title subtitle depname",
    "noid",
    "orthog",
    "endogtest",
    "redundant",
    "partial",
    "fwl",
    "dofminus",
    "sdofminus",
    "nopartialsmall",
    "eform",
    "plus noomitted vsquish noemptycells allbaselevels",
    "omitted baselevels",
    "savefirst savefprefix",
    "saverf saverfprefix",
    "residuals",
    "residuals2",
    "tolerance maxiter",
    "threads verbose",
]
estimate(
    "civreghdfe",
    "ivreghdfe",
    "y x (endog=z z2)",
    iv_options,
    [o.replace("threads(2) verbose", "") for o in iv_options],
    iv_names,
    predict="xb stdp",
)
# PPML option interactions already have native comparisons in the advanced suite;
# repeat numerical/display/tuning aliases across independent distributions here.
ppmlopts = [
    "absorb(g h)",
    "absorb(g h) vce(robust)",
    "absorb(g h) cluster(g)",
    "absorb(g h) keepsingletons",
    "noabsorb",
    "absorb(g h) exposure(w)",
    "absorb(g h) offset(z2)",
    "absorb(g h) dofadjustments(none)",
    "absorb(g h) eform irr nolog timeit",
    "absorb(g h) verbose(1) threads(2)",
    "absorb(g h) tolerance(1e-10) iterate(2000) irlstolerance(1e-10) irlsmaxiter(2000)",
    "absorb(g h) septolerance(1e-10) maxiterations(2000) itolerance(1e-10)",
]
ppmlrefs = [
    o.replace("threads(2)", "")
    .replace("irlstolerance(1e-10)", "")
    .replace("irlsmaxiter(2000)", "")
    .replace("septolerance(1e-10)", "")
    .replace("maxiterations(2000)", "maxiter(2000)")
    .replace("iterate(2000)", "maxiter(2000)")
    .replace("itolerance(1e-10)", "itol(1e-10)")
    for o in ppmlopts
]
estimate(
    "cppmlhdfe",
    "ppmlhdfe",
    "count_y x z",
    ppmlopts,
    ppmlrefs,
    [
        "absorb",
        "vce",
        "cluster",
        "keepsingletons",
        "noabsorb",
        "exposure",
        "offset",
        "dofadjustments",
        "eform irr nolog timeit",
        "verbose threads",
        "tolerance iterate irlstolerance irlsmaxiter",
        "septolerance maxiterations itolerance",
    ],
    predict="xb",
    scalars="df_m ll deviance",
)

# Range statistics compare every output variable, including counts/missingness.
for stat in [
    "count",
    "mean",
    "sum",
    "sd",
    "variance",
    "min",
    "max",
    "median",
    "first",
    "last",
    "skewness",
    "kurtosis",
]:
    for opt in [
        "interval(time -2 2)",
        "interval(time -2 2) by(group)",
        "interval(time -2 2) by(group) excludeself",
    ]:
        case(
            "crangestat",
            f"range {stat} {opt}",
            "interval by excludeself threads verbose",
            f"rangestat ({stat}) y, {opt}",
            f"crangestat ({stat}) y, {opt} threads(2) verbose",
        )
for options in [
    "neighbor(1)",
    "neighbor(3) ties",
    "noreplacement",
    "noreplacement descending",
    "caliper(.05) common",
    "radius caliper(.1)",
    *[
        f"kernel kerneltype({k}) bwidth(.15)"
        for k in ["epan", "normal", "biweight", "uniform", "tricube"]
    ],
]:
    case(
        "cpsmatch",
        "matching " + options,
        "outcome pscore neighbor ties noreplacement descending caliper common radius kernel kerneltype bwidth threads verbose",
        f"psmatch2 treat, pscore(ps) outcome(y) {options}",
        f"cpsmatch treat, pscore(ps) outcome(y) {options} threads(2) verbose",
        "matching",
    )
for method in ["logit", "probit"]:
    case(
        "cpsmatch",
        method + " propensity",
        "outcome logit probit",
        "psmatch2 treat x z, outcome(y) " + ("logit" if method == "logit" else ""),
        f"cpsmatch treat x z, outcome(y) {method}",
        "matching",
    )

# Sampling uses different RNGs; compare stratified counts here, with a separate
# repeated-seed distribution test for randomness, inclusion, and replacement.
for selection in ["", "if id<451", "in 20/480"]:
    for percent in [0, 25, 100]:
        case(
            "csample",
            f"sample {percent} {selection}",
            "by threads verbose",
            f'sample {percent} {selection.replace("in 20/480","if inrange(id,20,480)")}, by(group)',
            f"csample {percent} {selection}, by(group) threads(2) verbose",
            "sampling",
        )
    case(
        "csample",
        "sample count " + selection,
        "count by",
        f'sample 25 {selection.replace("in 20/480","if inrange(id,20,480)")}, count by(group)',
        f"csample {selection}, count(25) by(group)",
        "sampling",
    )
    for n in ["", "25"]:
        case(
            "cbsample",
            f"bootstrap {n} {selection}",
            "strata threads verbose",
            f"bsample {n} {selection}, strata(group)",
            f"cbsample {n} {selection}, strata(group) threads(2) verbose",
            "sampling",
        )

merge_setup = """tempfile usingdata
preserve
keep id g s y
rename y using_y
replace id=id+100
sort id
save `usingdata'
restore
sort id"""
for kind in ["1:1", "m:1", "1:m", "m:m", "1:1 _n"]:
    keys = "id" if kind != "1:1 _n" else ""
    for opt in [
        "",
        "keep(match)",
        "keep(master match) generate(status)",
        "nogenerate",
        "keepusing(using_y)",
        "sorted",
        "nolabel nonotes noreport",
        "update",
        "update replace",
        "force",
    ]:
        case(
            "cmerge",
            f"merge {kind} {opt}",
            "keep generate nogenerate keepusing sorted nolabel nonotes noreport update replace force threads verbose",
            f'merge {kind} {keys} using "`usingdata\'", {opt}',
            f'cmerge {kind} {keys} using "`usingdata\'", {opt} threads(2) verbose',
            setup=merge_setup,
            sortby="id",
        )

join_setup = """tempfile usingdata
preserve
keep id time g s y
rename id uid
rename y using_y
keep in 1/100
save `usingdata'
restore
keep in 101/160"""
for by in ["", "by(g)", "by(s)", "by(g s)"]:
    for opt in ["", "keepusing(uid using_y)", "all prefix(u_)", "all suffix(_other)"]:
        case(
            "crangejoin",
            f"rangejoin {by} {opt}",
            "by keepusing all prefix suffix threads verbose",
            f'rangejoin time lo hi using "`usingdata\'", {by} {opt}',
            f'crangejoin time lo hi using "`usingdata\'", {by} {opt} threads(2) verbose',
            setup=join_setup,
        )

for opt in [
    "nquantiles(5)",
    "nquantiles(25)",
    "nquantiles(10) controls(z)",
    "nquantiles(10) absorb(g)",
    "nquantiles(10) controls(z) absorb(g)",
    "discrete",
    "linetype(none)",
    "linetype(lfit)",
    "linetype(qfit)",
    "linetype(connect)",
    "reportreg",
    "title(Coverage) ytitle(Outcome) xtitle(Regressor)",
    "colors(red blue)",
    "msymbols(Oh)",
    "mlabels(group)",
    "legend(off)",
]:
    case(
        "cbinscatter",
        "bins " + opt,
        "nquantiles controls absorb discrete linetype reportreg title ytitle xtitle colors msymbols mlabels legend nograph threads verbose",
        "binscatter y x, " + opt + " nodraw",
        "cbinscatter y x, " + opt + " nograph threads(2) verbose",
        "bins",
    )
for w in ["aw=w", "fw=fw", "pw=w", "iw=w"]:
    case(
        "cbinscatter",
        "bins weight " + w,
        "nquantiles",
        f"binscatter y x [{w}], nquantiles(10) nodraw",
        f"cbinscatter y x [{w}], nquantiles(10) nograph",
        "bins",
    )
for opt in [
    "",
    "novarnames",
    "quote",
    "nolabel",
    "datafmt",
    "delimiter(;)",
    "noquoteif",
    "mmap",
    "nofsync",
    "direct",
    "prefault",
    "crlf",
    "noparallel",
    "threads(2) verbose",
]:
    native = opt if opt in ["", "novarnames", "quote", "nolabel", "datafmt", "delimiter(;)"] else ""
    raw = " rawfields" if opt == "noquoteif" else ""
    case(
        "cexport",
        "CSV writer " + opt,
        "delimiter novarnames quote nolabel datafmt noquoteif mmap nofsync direct prefault crlf noparallel threads verbose replace",
        f"coverage_export_readback, ctools(0) kind(delimited) opts({native}){raw}",
        f"coverage_export_readback, ctools(1) kind(delimited) opts({opt})",
    )
for opt in [
    "",
    "firstrow(variables)",
    "firstrow(varlabels)",
    "sheet(Coverage)",
    "cell(B3)",
    "missing(NA)",
    "nolabel",
    "locale(en_US)",
    "threads(2) verbose",
]:
    native = "" if opt == "threads(2) verbose" else opt
    case(
        "cexport",
        "Excel writer " + opt,
        "firstrow sheet cell missing nolabel locale threads verbose replace",
        f"coverage_export_readback, ctools(0) kind(excel) opts({native})",
        f"coverage_export_readback, ctools(1) kind(excel) opts({opt})",
    )
import_setup = """tempfile incoming
export delimited id y s group using "`incoming'", replace"""
for opt in [
    "",
    "varnames(1)",
    "varnames(nonames)",
    "case(upper)",
    "case(lower)",
    "encoding(utf8)",
    "bindquotes(strict)",
    "bindquotes(nobind)",
    "stripquotes(yes)",
    "stripquotes(no)",
    "rowrange(3:100)",
    "colrange(2:3)",
    "asfloat",
    "asdouble",
    "stringcols(_all)",
    "numericcols(2)",
    "decimalseparator(.)",
    "groupseparator(,)",
    "maxquotedrows(5)",
    "emptylines(include)",
    "parselocale(en_US)",
    "favorstrfixed",
    "threads(2) verbose",
]:
    native = "" if opt == "threads(2) verbose" else opt
    case(
        "cimport",
        "CSV reader " + opt,
        "varnames case encoding bindquotes stripquotes rowrange colrange asfloat asdouble stringcols numericcols decimalseparator groupseparator maxquotedrows emptylines parselocale favorstrfixed threads verbose clear",
        f'import delimited using "`incoming\'", clear {native}',
        f'cimport delimited using "`incoming\'", clear {opt}',
        setup=import_setup,
    )
for opt in ["cluster(g)", "cluster(g) strata(group)"]:
    case(
        "cbsample",
        "cluster bootstrap " + opt,
        "cluster strata",
        f"bsample 10, {opt}",
        f"cbsample 10, {opt}",
        "sampling",
        setup="replace g=ceil(id/3)\nreplace group=mod(g,4)",
    )
for opt in ["", "strata(group)", "cluster(g)", "cluster(g) strata(group)"]:
    case(
        "cbsample",
        "bootstrap frequencies " + opt,
        "weight strata cluster",
        f"coverage_weighted_bootstrap, ctools(0) {opt}",
        f"coverage_weighted_bootstrap, ctools(1) {opt}",
        "sampling",
        setup="replace g=ceil(id/3)\nreplace group=mod(g,4)",
    )

for residual in ["residuals(saved_resid)", "residuals2"]:
    estimate(
        "civreghdfe",
        "ivreghdfe",
        "y x (endog=z z2)",
        ["absorb(g) " + residual],
        ["absorb(g) residuals(saved_resid)"],
        [residual.split("(")[0]],
        predict="residuals xbd",
    )
for dopt in ["d(saved_d)", "d2"]:
    estimate(
        "cppmlhdfe",
        "ppmlhdfe",
        "count_y x z",
        ["absorb(g h) " + dopt],
        ["absorb(g h) d(saved_d)"],
        [dopt.split("(")[0]],
        predict="xb xbd d mu eta stdp anscombe deviance pearson response scores working",
    )
# Solver controls must converge to the same independent native solution.
ppml_controls = [
    "guess(simple)",
    "guess(ols)",
    "guess(variable w)",
    "standardize_data(0)",
    "min_ok(3)",
    "use_exact_solver(1)",
    "use_exact_partial(0)",
    "use_heuristic_tol(0)",
    "start_inner_tol(.001)",
    "remove_collinear_variables(0)",
    "use_step_halving(1) step_halving_memory(.8) max_step_halving(4)",
    "mu_tol(1e-8)",
    "relu_tol(1e-6) relu_zero_tol(1e-10) relu_maxiter(200) relu_strict(1) relu_accelerate(1)",
    "simplex_tol(1e-12) simplex_maxiter(2000)",
    "acceleration(cg)",
    "acceleration(none)",
    "acceleration(sd)",
    "acceleration(aitken)",
    "acceleration(hybrid)",
    "acceleration(lsmr) btol(1e-12) conlim(1e9)",
    "transform(kaczmarz)",
    "transform(cimmino)",
    "poolsize(1)",
    "separation(fe)",
    "separation(simplex)",
    "separation(relu)",
    "separation(mu)",
    "separation(fe simplex relu mu)",
]
for opt in ppml_controls:
    names = " ".join(__import__("re").findall(r"([a-z_]+)\(", opt))
    case(
        "cppmlhdfe",
        "PPML control " + opt,
        names,
        "ppmlhdfe count_y x z, absorb(g h) tol(1e-12) use_exact_solver(1) use_exact_partial(1) min_ok(3)",
        "cppmlhdfe count_y x z, absorb(g h) tol(1e-12) " + opt,
        "estimate",
        predict="xb",
        scalars="df_m ll deviance",
    )
# Initialization must not change CUE's converged result.
case(
    "civreghdfe",
    "CUE starting values",
    "b0 cue",
    "ivreghdfe y x (endog=z z2), absorb(g) cue robust",
    "civreghdfe y x (endog=z z2), absorb(g) cue robust b0(initial_b)",
    "estimate",
    setup="quietly ivreghdfe y x (endog=z z2), absorb(g) robust\nmatrix initial_b=e(b)",
)
case(
    "cdestring",
    "destring replace",
    "replace",
    "destring number, replace",
    "cdestring number, replace",
)

for cmd, ref in [("creghdfe", "reghdfe"), ("civreghdfe", "ivreghdfe")]:
    lhs = "y x z" if cmd == "creghdfe" else "y x (endog=z z2)"
    estimate(
        cmd, ref, lhs, ["absorb(g) residuals(saved_resid)"], names=["residuals"], predict="score"
    )

# Inputs for native time-series covariance options must have declared time.
for c in CASES:
    if c["command"] == "civreghdfe" and any(
        o in c["candidate"] for o in ["bw(3)", "dkraay(3)", "kiefer"]
    ):
        c["setup"] = "bysort g (id): gen time_iv=_n\nxtset g time_iv"
    if c["command"] == "cbinscatter" and "discrete" in c["candidate"]:
        c["setup"] = "replace x=round(x,.5)"
    if c["command"] == "cbinscatter" and ("[pw=w]" in c["reference"] or "[iw=w]" in c["reference"]):
        # Bin weighted means use the same normalized weights; binscatter only
        # accepts aw/fw, so translate ctools' pw/iw extension to aw.
        c["reference"] = c["reference"].replace("[pw=w]", "[aw=w]").replace("[iw=w]", "[aw=w]")
    if c["command"] == "cbinscatter" and "mlabels(group)" in c["reference"]:
        # mlabels is a ctools presentation extension, not native binscatter syntax.
        c["reference"] = c["reference"].replace("mlabels(group)", "")

for opt in [
    "by(group)",
    "by(group) controls(z)",
    "by(group) absorb(g)",
    "by(group) controls(z) absorb(g)",
]:
    case(
        "cbinscatter",
        "grouped bins " + opt,
        "by controls absorb",
        "binscatter y x, nquantiles(10) " + opt + " nodraw",
        "cbinscatter y x, nquantiles(10) " + opt + " nograph",
        "bins",
    )
case(
    "cexport",
    "date format override",
    "datestring",
    "coverage_export_readback, ctools(0) kind(delimited) opts(datafmt)",
    "coverage_export_readback, ctools(1) kind(delimited) opts(datestring(%tdCCYY-NN-DD))",
    setup="replace y=id-300 if !missing(y)\nformat y %tdCCYY-NN-DD",
)
# Exported formats and workbook edits have specialized exact metadata checks
# in validate_io_excel_options.do and validate_io_formats.do.
for opt in ["charset(utf8)", 'delimiters(",")', 'delimiters(",") collapsedelimiters']:
    case(
        "cimport",
        "CSV parser alias " + opt,
        "charset delimiters collapsedelimiters",
        'import delimited using "`incoming\'", clear '
        + opt.replace("charset(utf8)", "encoding(utf8)"),
        f'cimport delimited using "`incoming\'", clear {opt}',
        setup=import_setup,
    )

# Diagnostics must exist; two missing scalars are not evidence of equivalence.
for c in CASES:
    if c["command"] == "civreghdfe":
        if "orthog(z2)" in c["candidate"]:
            c["required"] = "cstat cstatp:cstat_p cstatdf:cstat_df"
        if "endogtest(endog)" in c["candidate"]:
            c["required"] = "estat:endogtest estatp:endogtest_p estatdf:endogtest_df"
        if "redundant(z2)" in c["candidate"]:
            c["required"] = "redstat:redund redp:redund_p reddf:redund_df"
        if " cue " in c["reference"] + " ":
            c["reference"] = "coverage_unitscale " + c["reference"]
        if "first" in c["options"] or "ffirst" in c["options"]:
            c["reference"] = c["reference"].replace("ivreghdfe ", "coverage_iv_firstref ", 1)
            c["required"] = "F_first1" + (" partial_r2_1" if "ffirst" in c["options"] else "")
for stage in ["first", "rf"]:
    case(
        "civreghdfe",
        f"IV saved {stage} estimates",
        f"save{stage} save" + ("fprefix" if stage == "first" else "rfprefix"),
        f"coverage_saved_iv, ctools(0) stage({stage})",
        f"coverage_saved_iv, ctools(1) stage({stage})",
        "estimate",
        required="df_r",
        predict="xb",
    )

# Activate conversion branches with actual percent signs, decimal commas,
# ignored characters and malformed input, rather than clean numeric strings.
for c in CASES:
    if c["command"] == "cdestring":
        opt = c["candidate"]
        if "percent" in opt:
            c["setup"] = 'replace number=number+"%" if !missing(y)'
        elif "dpcomma" in opt:
            c["setup"] = 'replace number=subinstr(number,".",",",.) if !missing(y)'
        elif "ignore(x)" in opt:
            c["setup"] = 'replace number="x"+number+"x" if !missing(y)'
        elif "force" in opt:
            c["setup"] = 'replace number="not_numeric" if mod(id,11)==0'
# Conflicts, missing updates, imported labels and notes make the merge switches
# observable. Duplicate-key multiplicity is separately exercised by legacy tests.
merge_overlap_setup = (
    merge_setup.replace(
        "rename y using_y",
        'clonevar using_y=y\nreplace y=y+10\nlabel define using_groups 0 "from using" 1 "using one" 2 "using two"\nlabel values using_y using_groups\nnotes: using dataset note\nnotes y: using variable note',
    )
    + "\nreplace y=. if mod(id,7)==0"
)
for c in CASES:
    if c["command"] == "cmerge":
        c["setup"] = merge_overlap_setup
# reghdfe's results lose accuracy for outcomes far from unit scale; see
# coverage_unitscale, which estimates at unit scale and rescales.
for c in CASES:
    if c["oracle"] == "estimate" and c["reference"].startswith("reghdfe "):
        c["reference"] = "coverage_unitscale " + c["reference"]
# Equal propensity scores force the ties branch to allocate fractional weights.
case(
    "cpsmatch",
    "exact score ties",
    "pscore outcome neighbor ties",
    "psmatch2 treat, pscore(ps) outcome(y) neighbor(3) ties",
    "cpsmatch treat, pscore(ps) outcome(y) neighbor(3) ties",
    "matching",
    setup="replace ps=round(ps,.1)\nreplace ps=.05 if ps==0\nreplace ps=.95 if ps==1",
)

for cmd, native, var in [("cencode", "encode", "s"), ("cdecode", "decode", "group")]:
    case(
        cmd,
        native + " replace adapter",
        "replace",
        f"coverage_codec_replace {var}, command({native})",
        f"{cmd} {var}, replace",
    )
case(
    "cmerge",
    "merge assert success",
    "assert",
    'merge 1:1 id using "`usingdata\'", assert(match master using)',
    'cmerge 1:1 id using "`usingdata\'", assert(match master using)',
    setup=merge_overlap_setup,
)
case(
    "cmerge",
    "preserve master observation order",
    "preserve_order",
    'coverage_merge_preserved using "`usingdata\'"',
    'cmerge 1:1 id using "`usingdata\'", preserve_order(1)',
    setup=merge_overlap_setup + "\ngsort -id",
)
case(
    "creghdfe",
    "mobility groups across disconnected components",
    "groupvar",
    "reghdfe y x z, absorb(g h) groupvar(mobility)",
    "creghdfe y x z, absorb(g h) groupvar(mobility)",
    setup="replace g=g+100*group\nreplace h=h+100*group",
)

"""Audit public option declarations against executable differential cases.

This is an inventory gate, not a claim that a passing call proves equivalence.
The runner separately checks all case/fixture results. Existing specialized
suites remain authoritative for format metadata, random draws, and diagnostics.
"""


def inventory(root=ROOT):
    found = {}
    for path in sorted((root / "build").glob("c*.ado")):
        text = re.sub(r"/\*.*?\*/", "", path.read_text(), flags=re.S)
        text = re.sub(r"///[^\n]*\n\s*", " ", text)
        program = None
        opts = {}
        for line in text.splitlines():
            m = re.match(r"\s*program define (\w+)", line)
            if m:
                program = m[1]
            if program not in (path.stem, path.stem + "_excel", "_" + path.stem + "_impl"):
                continue
            m = re.match(r"\s*(?:capture noisily |capture )?syntax\b(.*)", line)
            if not m or "," not in m[1]:
                continue
            args = m[1].split(",", 1)[1]
            # Remove balanced option argument declarations, preserving names.
            while "(" in args:
                new = re.sub(r"\([^()]*\)", "", args)
                if new == args:
                    raise ValueError("unbalanced syntax: " + str(path))
                args = new
            for token in re.findall(r"[A-Za-z_][A-Za-z_0-9:]*", args):
                token = token.replace(":", "")
                if token == "options":
                    continue
                opts[token.lower()] = token
        if path.stem == "cbsample":
            opts["weight"] = "weight"
        if not opts:
            raise ValueError("public command has no audited syntax: " + str(path))
        found[path.stem] = opts
    # Statistical-format dispatch options are composed in private adapters.
    for cmd, helper in [("cimport", "_cio_import"), ("cexport", "_cio_export")]:
        for line in (root / f"build/{helper}.ado").read_text().splitlines():
            if "local options" not in line:
                continue
            for token in re.findall(r"\b([A-Z][A-Za-z]*(?:\([^)]*\))?)", line):
                word = token.split("(")[0]
                if word != "CLEAR":
                    found[cmd][word.lower()] = word
    shp = (root / "build/_cio_shp_export.ado").read_text()
    for line in shp.splitlines():
        if not re.match(r"\s*syntax\b", line):
            continue
        args = re.sub(r"\([^()]*\)", "", line.split(",", 1)[1])
        for token in re.findall(r"[A-Za-z_][A-Za-z_0-9]*", args):
            found["cexport"][token.lower()] = token
    return found


def exercised(case, option, declaration):
    if option not in case["options"]:
        return False
    minimum = "".join(c for c in declaration if c.isupper()).lower() or option
    # Stata abbreviations must be prefixes, never substring mentions.
    tokens = re.findall(r"\b[a-zA-Z_][a-zA-Z_0-9]*\b", case["candidate"].lower())
    return any(option.startswith(t) and len(t) >= len(minimum) for t in tokens)


def records():
    legacy = json.loads((ROOT / "validation/option_evidence.json").read_text())
    rows = []
    for cmd, options in inventory(ROOT).items():
        for option, declaration in options.items():
            cases = [
                i
                for i, c in enumerate(CASES)
                if c["command"] == cmd and exercised(c, option, declaration)
            ]
            if cmd.endswith("_p"):
                parent = cmd[:-2]
                cases = [
                    i
                    for i, c in enumerate(CASES)
                    if c["command"] == parent and option in c.get("predict", "").split()
                ]
            evidence = legacy.get(cmd, {}).get(option)
            if evidence:
                if (
                    evidence["kind"] not in ("differential", "contract", "rejection", "unsupported")
                    or not evidence["reason"]
                    or not evidence["evidence"]
                ):
                    raise ValueError(f"invalid evidence: {cmd}.{option}")
                for e in evidence["evidence"]:
                    path = ROOT / e["path"]
                    if not path.exists() or e["anchor"] not in path.read_text():
                        raise ValueError(f"stale evidence: {cmd}.{option}: {e}")
            rows.append(
                dict(
                    command=cmd,
                    option=option,
                    case_ids=cases,
                    status="differential" if cases else evidence["kind"] if evidence else "gap",
                    evidence=evidence,
                )
            )
    extras = {(cmd, opt) for cmd, opts in legacy.items() for opt in opts} - {
        (r["command"], r["option"]) for r in rows
    }
    if extras:
        raise ValueError(f"stale option entries: {extras}")
    return rows


def report(rows):
    lines = [
        "# Check command options",
        "",
        "Complete the [test setup](../validation/README.md#before-you-start), then run:",
        "",
        "```sh",
        "python3 validation/run.py --stata option_coverage",
        "```",
        "",
        "To check one command and its options, replace `cencode` below:",
        "",
        "```sh",
        "python3 validation/run.py cencode",
        "```",
        "",
        "| Command | Option | Check |",
        "|---|---|---|",
    ]
    labels = {
        "differential": "Compare command results",
        "contract": "Check option behavior",
        "rejection": "Check error handling",
        "unsupported": "Check rejection of an unsupported option",
        "gap": "**Missing test**",
    }
    for row in rows:
        lines.append(f'| `{row["command"]}` | `{row["option"]}` | {labels[row["status"]]} |')
    return "\n".join(lines) + "\n"


"""Run every explicit option comparison, retaining per-case failures and coverage."""


def quote(s):
    return '`"' + s + "\"'"


def source(cases, log):
    lines = [
        f'cd "{ROOT}"',
        "clear all",
        "set more off",
        "set linesize 255",
        "capture log close _all",
        f'log using "{log}", text replace',
        "capture noisily {",
        "quietly do validation/support/validate_setup.do",
        "quietly do validation/support/coverage_helpers.do",
        "local failed_cases=0",
        "di \"REFERENCE_STATA_VERSION=`c(stata_version)' OS=`c(os)'\"",
    ]
    for cmd in [
        "qreg",
        "reghdfe",
        "ivreghdfe",
        "ppmlhdfe",
        "psmatch2",
        "binscatter",
        "rangestat",
        "rangejoin",
        "winsor2",
    ]:
        lines.append("which " + cmd)
    for index, c in cases:
        lines += ["forvalues variant=1/3 {", "quietly coverage_fixture `variant'"]
        for line in c.get("setup", "").splitlines():
            lines.append("quietly " + line)
        opts = ""
        for key in ("sortby", "predict", "scalars", "required", "ref_after"):
            if c.get(key):
                opts += " " + key.replace("ref_after", "refafter") + "(" + quote(c[key]) + ")"
        lines += [
            "local before=$TESTS_FAILED",
            "global FAILURE_COUNT=0",
            f'capture noisily coverage_{c["oracle"]}, name("option-{index} variant-`variant\'") '
            + "reference("
            + c["reference"]
            + ") candidate("
            + c["candidate"]
            + ")"
            + opts,
            "local case_rc=_rc",
            f'if `case_rc\' test_fail "option-{index}" "harness error r(`case_rc\')"',
            "local bad=($TESTS_FAILED>`before')",
            "if `bad' print_failures",
            "local failed_cases=`failed_cases'+`bad'",
            f"di \"OPTION_CASE={index} VARIANT=`variant' FAILED=`bad' RC=`case_rc'\"",
            "}",
        ]
    lines += [
        'di "OPTION_COVERAGE_PASSED=$TESTS_PASSED FAILED=$TESTS_FAILED SKIPPED=$TESTS_SKIPPED"',
        "local rc=cond(`failed_cases',9,0)",
        'di "OPTION_COVERAGE_COMPLETE"',
        "}",
        "local harness_rc=_rc",
        "if `harness_rc' local rc=`harness_rc'",
        "log close",
        "exit, clear",
    ]
    return "\n".join(lines) + "\n"


def validate_cases(text, cases):
    results = {}
    for i, v, f, rc in re.findall(
        r"^OPTION_CASE=(\d+) VARIANT=(\d+) FAILED=(\d+) RC=(\d+)$", text, re.M
    ):
        if f not in ("0", "1"):
            raise ValueError("invalid failure flag")
        key = (int(i), int(v))
        if key in results:
            raise ValueError("duplicate case result")
        results[key] = dict(failed=int(f), rc=int(rc))
    expected = {(i, v) for i, _ in cases for v in (1, 2, 3)}
    if set(results) != expected:
        raise ValueError(f"incomplete option cases: missing {expected-set(results)}")
    return results


def run(output, commands=None):
    from types import SimpleNamespace

    a = SimpleNamespace(output=Path(output), commands=commands)
    cases = [(i, c) for i, c in enumerate(CASES) if not a.commands or c["command"] in a.commands]
    if not cases:
        raise SystemExit("No cases selected")
    a.output.mkdir(parents=True, exist_ok=True)
    log = a.output / "results.log"
    inputs = list((ROOT / "build").glob("*.ado")) + list((ROOT / "build").glob("*mac_arm.plugin"))
    hashes = lambda: {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs
    }
    before_hashes = hashes()
    launch_error = None
    with tempfile.TemporaryDirectory(prefix="ctools-options-") as t:
        driver = Path(t) / "options.do"
        driver.write_text(source(cases, log))
        try:
            launch_driver(driver, log)
        except ValueError as error:
            launch_error = str(
                error
            )  # Retain failed case details, but never turn launcher failure into success.
    text = log.read_text(errors="replace")
    results = validate_cases(text, cases)
    report = []
    for i, c in cases:
        for variant in (1, 2, 3):
            report.append(dict(index=i, variant=variant, **c, **results[(i, variant)]))
    after_hashes = hashes()
    (a.output / "provenance.json").write_text(
        json.dumps(dict(before=before_hashes, after=after_hashes), indent=2) + "\n"
    )
    (a.output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    failures = [r for r in report if r["failed"] or r["rc"]]
    for r in failures:
        print(f'FAIL {r["name"]} variant={r["variant"]} rc={r["rc"]}')
    print(
        f"{len(report)-len(failures)}/{len(report)} option/data cases passed; {len(failures)} failed; log: {log}"
    )
    if before_hashes != after_hashes:
        raise SystemExit("Build inputs changed during validation; rerun on a stable build")
    if failures or launch_error:
        raise SystemExit(1)
