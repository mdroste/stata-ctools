*! version 1.2.6 20261111 github.com/mdroste/stata-ctools
*! ctools: a file name with a leading ~ read as the home directory
* _ctools_home lname name sets the caller's local lname to name, with a
* leading ~ (alone, or before / or \) replaced by the home directory: HOME,
* or USERPROFILE on Windows, as Stata's own file commands read it. The
* plugin opens names as they are, so ~/data.csv would name a folder ~ in the
* working directory. Other names, URLs included, are left as they are.
program define _ctools_home
    version 14.1
    args lname name
    local home : environment HOME
    if c(os) == "Windows" local home : environment USERPROFILE
    mata: st_local("name", _ctools_home_name(st_local("name"), st_local("home")))
    c_local `lname' `"`macval(name)'"'
end

version 14.1
mata:
string scalar _ctools_home_name(string scalar name, string scalar home)
{
    if (home == "" | substr(name, 1, 1) != "~") return(name)
    if (strlen(name) > 1 & !anyof(("/", "\"), substr(name, 2, 1))) return(name)
    return(home + substr(name, 2, .))
}
end
