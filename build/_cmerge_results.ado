*! version 1.2.0 20261029 github.com/mdroste/stata-ctools
*! Private helper for cmerge.
program define _cmerge_results
    gettoken macname 0 : 0
    gettoken colon 0 : 0
    local codes ""
    foreach w of local 0 {
        local l = strlen(`"`w'"')
        local c ""
        if inlist(`"`w'"', "1", "2", "3", "4", "5") local c `w'
        else if substr("masters", 1, max(3, `l')) == `"`w'"' local c 1
        else if substr("usings", 1, max(2, `l')) == `"`w'"' local c 2
        else if substr("matches", 1, max(3, `l')) == `"`w'"' local c 3
        else if substr("matched", 1, max(3, `l')) == `"`w'"' local c 3
        else if substr("match_updates", 1, max(8, `l')) == `"`w'"' local c 4
        else if substr("match_conflicts", 1, max(8, `l')) == `"`w'"' local c 5
        else {
            di as err `"`w':  invalid {it:resulttype}"'
            di as err "{p 4 4 2}"
            di as err "results, specified in options {bf:assert()} and {bf:keep()},"
            di as err "are the integers 1 through 5, or"
            di as err "{bf:master} (equivalent to 1),"
            di as err "{bf:using} (2),"
            di as err "{bf:match} (3),"
            di as err "{bf:match_update} (4), or"
            di as err "{bf:match_conflict} (5).  The last two arise only when"
            di as err "option {cmd:update} is specified."
            di as err "{p_end}"
            exit 198
        }
        local codes `codes' `c'
    }
    local codes : list uniq codes
    local codes : list sort codes
    c_local `macname' "`codes'"
end
