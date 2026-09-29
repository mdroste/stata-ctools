
libname datapath '/Users/Mike/Documents/GitHub/stata-ctools/' ;
libname xptfile xport '/Users/Mike/Documents/GitHub/stata-ctools/temp/io_formats/rename5.xpt' ;

proc copy in = xptfile out = datapath ;

proc format library = work ;
    value MYLONGVA
           1 = 'Yes'
           2 = 'No'
           3 = 'Other' ;

quit ;
