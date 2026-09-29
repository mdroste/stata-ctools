clear all
set more off
log using "temp/io_formats/xpf.log", text replace
import sasxport5 "temp/io_formats/formats.xpf", novallabels
describe
list, noobs
log close
exit, clear
