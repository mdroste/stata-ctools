/*******************************************************************************
 * validate_cimport_helper.do
 *
 * Helper file for validate_cimport.do - creates test data files
 * This file is called by validate_cimport.do to create test CSV files
 ******************************************************************************/

/*******************************************************************************
 * Create Test Data Files
 ******************************************************************************/
print_section "Creating Test Data Files"

capture file close _all

* Basic CSV
file open fh using "temp/basic.csv", write replace
file write fh "id,name,value,category" _n
file write fh "1,Alpha,100.5,A" _n
file write fh "2,Beta,200.25,B" _n
file write fh "3,Gamma,300.75,A" _n
file write fh "4,Delta,400.1,C" _n
file write fh "5,Epsilon,500.9,B" _n
file close fh

* Tab-delimited (use _tab for actual tab characters)
file open fh using "temp/tab_delim.tsv", write replace
file write fh "id" _tab "name" _tab "value" _n
file write fh "1" _tab "Alpha" _tab "100" _n
file write fh "2" _tab "Beta" _tab "200" _n
file write fh "3" _tab "Gamma" _tab "300" _n
file close fh

* Semicolon-delimited
file open fh using "temp/semicolon.csv", write replace
file write fh "id;name;value" _n
file write fh "1;Alpha;100" _n
file write fh "2;Beta;200" _n
file close fh

* No header
file open fh using "temp/noheader.csv", write replace
file write fh "1,Alpha,100" _n
file write fh "2,Beta,200" _n
file write fh "3,Gamma,300" _n
file close fh

* Quoted fields
file open fh using "temp/quoted.csv", write replace
file write fh `"id,name,value"' _n
file write fh `"1,"Smith, John",100"' _n
file write fh `"2,"Doe, Jane",200"' _n
file write fh `"3,"O'Brien, Pat",300"' _n
file close fh

* Mixed case headers
file open fh using "temp/mixedcase.csv", write replace
file write fh "ID,FirstName,LastValue" _n
file write fh "1,John,100" _n
file write fh "2,Jane,200" _n
file close fh

* Numeric types
file open fh using "temp/numerics.csv", write replace
file write fh "int_col,float_col,big_col" _n
file write fh "1,1.5,1000000000" _n
file write fh "2,2.5,2000000000" _n
file write fh "3,3.5,3000000000" _n
file close fh

* Missing values
file open fh using "temp/missing.csv", write replace
file write fh "id,value,name" _n
file write fh "1,100,Alpha" _n
file write fh "2,,Beta" _n
file write fh "3,300," _n
file write fh "4,," _n
file close fh

* European format (semicolon + comma decimal)
file open fh using "temp/european.csv", write replace
file write fh "id;amount;price" _n
file write fh "1;1.234,56;99,99" _n
file write fh "2;2.345,67;199,99" _n
file write fh "3;3.456,78;299,99" _n
file close fh

* UTF-8 international
file open fh using "temp/utf8_intl.csv", write replace
file write fh "name,city,country" _n
file write fh "Jose Garcia,Sao Paulo,Brasil" _n
file write fh "Francois Muller,Zurich,Schweiz" _n
file write fh "Soren Orsted,Kobenhavn,Danmark" _n
file close fh

* Large file (10K rows)
file open fh using "temp/large_10k.csv", write replace
file write fh "id,value,category" _n
forvalues i = 1/10000 {
    local cat = mod(`i', 10) + 1
    file write fh "`i',`=runiform()',cat`cat'" _n
}
file close fh

* Zipcode-style data (leading zeros)
file open fh using "temp/zipcodes.csv", write replace
file write fh "zipcode,city,value" _n
file write fh "01234,Boston,100" _n
file write fh "00501,Holtsville,200" _n
file write fh "90210,Beverly Hills,300" _n
file close fh

* Row/col range test file
file open fh using "temp/range_test.csv", write replace
file write fh "a,b,c,d,e" _n
forvalues i = 1/100 {
    file write fh "`i',`=`i'*2',`=`i'*3',`=`i'*4',`=`i'*5'" _n
}
file close fh

/*******************************************************************************
 * PATHOLOGICAL TEST DATA FILES
 ******************************************************************************/
* (silent) "  Creating pathological test data files..."

* Empty file (0 bytes)
file open fh using "temp/empty_file.csv", write replace
file close fh

* Header only (no data rows)
file open fh using "temp/header_only.csv", write replace
file write fh "id,name,value" _n
file close fh

* Single empty row after header
file open fh using "temp/single_empty_row.csv", write replace
file write fh "id,name,value" _n
file write fh ",," _n
file close fh

* Empty rows in the middle of data
file open fh using "temp/empty_rows_middle.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100" _n
file write fh ",," _n
file write fh "3,Gamma,300" _n
file write fh ",," _n
file write fh "5,Epsilon,500" _n
file close fh

* Blank lines (no commas) in the middle
file open fh using "temp/blank_lines_middle.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100" _n
file write fh "" _n
file write fh "3,Gamma,300" _n
file write fh "" _n
file write fh "5,Epsilon,500" _n
file close fh

* All columns empty (just commas)
file open fh using "temp/all_columns_empty.csv", write replace
file write fh "a,b,c,d,e" _n
forvalues i = 1/10 {
    file write fh ",,,," _n
}
file close fh

* One column entirely missing
file open fh using "temp/one_col_all_missing.csv", write replace
file write fh "id,empty_col,value" _n
file write fh "1,,100" _n
file write fh "2,,200" _n
file write fh "3,,300" _n
file write fh "4,,400" _n
file write fh "5,,500" _n
file close fh

* One row entirely missing (all empty values)
file open fh using "temp/one_row_all_missing.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100" _n
file write fh "2,Beta,200" _n
file write fh ",," _n
file write fh "4,Delta,400" _n
file write fh "5,Epsilon,500" _n
file close fh

* Inconsistent column counts (fewer columns in some rows)
file open fh using "temp/fewer_cols.csv", write replace
file write fh "a,b,c,d,e" _n
file write fh "1,2,3,4,5" _n
file write fh "1,2,3" _n
file write fh "1,2,3,4,5" _n
file write fh "1,2" _n
file close fh

* Inconsistent column counts (more columns in some rows)
file open fh using "temp/extra_cols.csv", write replace
file write fh "a,b,c" _n
file write fh "1,2,3" _n
file write fh "1,2,3,4,5" _n
file write fh "1,2,3" _n
file write fh "1,2,3,4,5,6,7" _n
file close fh

* Trailing commas on all rows
file open fh using "temp/trailing_commas.csv", write replace
file write fh "id,name,value," _n
file write fh "1,Alpha,100," _n
file write fh "2,Beta,200," _n
file write fh "3,Gamma,300," _n
file close fh

* Trailing commas on some rows
file open fh using "temp/trailing_commas_some.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100," _n
file write fh "2,Beta,200" _n
file write fh "3,Gamma,300," _n
file close fh

* Leading whitespace in fields
file open fh using "temp/leading_whitespace.csv", write replace
file write fh "id,name,value" _n
file write fh "1, Alpha,100" _n
file write fh "2,  Beta, 200" _n
file write fh "3,   Gamma,  300" _n
file close fh

* Trailing whitespace in fields
file open fh using "temp/trailing_whitespace.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha ,100 " _n
file write fh "2,Beta  ,200  " _n
file write fh "3,Gamma   ,300" _n
file close fh

* Mixed leading and trailing whitespace
file open fh using "temp/mixed_whitespace.csv", write replace
file write fh " id , name , value " _n
file write fh " 1 , Alpha , 100 " _n
file write fh " 2 ,  Beta  , 200 " _n
file close fh

* Duplicate header names
file open fh using "temp/duplicate_headers.csv", write replace
file write fh "id,name,value,name,value" _n
file write fh "1,Alpha,100,Beta,200" _n
file write fh "2,Gamma,300,Delta,400" _n
file close fh

* Special characters in headers
file open fh using "temp/special_header_chars.csv", write replace
file write fh "id#1,name@2,value$3" _n
file write fh "1,Alpha,100" _n
file write fh "2,Beta,200" _n
file close fh

* Headers with spaces
file open fh using "temp/headers_with_spaces.csv", write replace
file write fh "my id,first name,total value" _n
file write fh "1,John,100" _n
file write fh "2,Jane,200" _n
file close fh

* Numeric headers
file open fh using "temp/numeric_headers.csv", write replace
file write fh "2020,2021,2022,2023" _n
file write fh "100,200,300,400" _n
file write fh "110,220,330,440" _n
file close fh

* Very long header names
file open fh using "temp/long_headers.csv", write replace
local long_name = "a" * 100
file write fh "id,`long_name',value" _n
file write fh "1,test,100" _n
file write fh "2,test2,200" _n
file close fh

* File ending without newline
file open fh using "temp/no_final_newline.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100" _n
file write fh "2,Beta,200"
file close fh

* File with Windows line endings (CRLF)
file open fh using "temp/crlf_endings.csv", write replace
file write fh "id,name,value" _char(13) _n
file write fh "1,Alpha,100" _char(13) _n
file write fh "2,Beta,200" _char(13) _n
file close fh

* File with mixed line endings
file open fh using "temp/mixed_line_endings.csv", write replace
file write fh "id,name,value" _n
file write fh "1,Alpha,100" _char(13) _n
file write fh "2,Beta,200" _n
file write fh "3,Gamma,300" _char(13) _n
file close fh

* Extremely long lines (many columns)
file open fh using "temp/very_wide.csv", write replace
local header ""
forvalues i = 1/100 {
    if `i' > 1 local header "`header',"
    local header "`header'v`i'"
}
file write fh "`header'" _n
forvalues j = 1/5 {
    local line ""
    forvalues i = 1/100 {
        if `i' > 1 local line "`line',"
        local line "`line'`i'"
    }
    file write fh "`line'" _n
}
file close fh

* Very long string values
file open fh using "temp/very_long_strings.csv", write replace
file write fh "id,text" _n
forvalues i = 1/5 {
    local long_text = "x" * 500
    file write fh "`i',`long_text'" _n
}
file close fh

* Quoted fields with embedded commas
file open fh using "temp/embedded_commas.csv", write replace
file write fh "id,name,address" _n
file write fh `"1,John,"123 Main St, Apt 4""' _n
file write fh `"2,Jane,"456 Oak Ave, Suite 100""' _n
file write fh `"3,Bob,"789 Pine Rd, Building A, Floor 2""' _n
file close fh

* Quoted fields with embedded quotes (escaped)
file open fh using "temp/embedded_quotes.csv", write replace
file write fh "id,name,quote" _n
file write fh `"1,John,"He said ""hello""!"' _n
file write fh `"2,Jane,"She replied ""goodbye""."' _n
file close fh

* Quoted fields with embedded newlines
file open fh using "temp/embedded_newlines.csv", write replace
file write fh "id,name,notes" _n
file write fh `"1,John,"Line 1"' _n
file write fh `"Line 2""' _n
file write fh `"2,Jane,"Single line""' _n
file close fh

* Mix of quoted and unquoted fields
file open fh using "temp/mixed_quoting.csv", write replace
file write fh "id,name,value" _n
file write fh `"1,"John Smith",100"' _n
file write fh "2,Jane,200" _n
file write fh `"3,"Bob ""The Builder""",300"' _n
file close fh

* Only whitespace in some fields
file open fh using "temp/whitespace_only_fields.csv", write replace
file write fh "id,name,value" _n
file write fh "1,   ,100" _n
file write fh "2,	,200" _n
file write fh "3, 	 ,300" _n
file close fh

* Numeric-looking strings (phone numbers, IDs with leading zeros)
file open fh using "temp/numeric_strings.csv", write replace
file write fh "id,phone,zipcode,ssn" _n
file write fh "1,555-123-4567,01234,123-45-6789" _n
file write fh "2,800-555-0000,00501,987-65-4321" _n
file write fh "3,123-456-7890,90210,111-22-3333" _n
file close fh

* Scientific notation
file open fh using "temp/scientific_notation.csv", write replace
file write fh "id,small,large,mixed" _n
file write fh "1,1e-10,1e10,1.23e5" _n
file write fh "2,2.5e-8,2.5e8,-3.14e-3" _n
file write fh "3,0.0,1E10,1.0E-5" _n
file close fh

* Negative numbers
file open fh using "temp/negative_numbers.csv", write replace
file write fh "id,value,balance" _n
file write fh "1,-100,-50.25" _n
file write fh "2,-0,-0.0" _n
file write fh "3,-999999,-0.00001" _n
file close fh

* Mixed positive and negative
file open fh using "temp/mixed_signs.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "2,-100" _n
file write fh "3,0" _n
file write fh "4,-0" _n
file write fh "5,+100" _n
file close fh

* Boolean-like values
file open fh using "temp/boolean_like.csv", write replace
file write fh "id,flag1,flag2,flag3,flag4" _n
file write fh "1,true,TRUE,True,1" _n
file write fh "2,false,FALSE,False,0" _n
file write fh "3,yes,YES,Yes,Y" _n
file write fh "4,no,NO,No,N" _n
file close fh

* Date-like values
file open fh using "temp/date_formats.csv", write replace
file write fh "id,date1,date2,date3,date4" _n
file write fh "1,2023-01-15,01/15/2023,15-Jan-2023,20230115" _n
file write fh "2,2023-12-31,12/31/2023,31-Dec-2023,20231231" _n
file write fh "3,1999-06-01,06/01/1999,01-Jun-1999,19990601" _n
file close fh

* Time-like values
file open fh using "temp/time_formats.csv", write replace
file write fh "id,time1,time2,time3" _n
file write fh "1,12:30:45,12:30,12:30:45.123" _n
file write fh "2,00:00:00,0:00,00:00:00.000" _n
file write fh "3,23:59:59,23:59,23:59:59.999" _n
file close fh

* Currency-like values
file open fh using "temp/currency_values.csv", write replace
file write fh "id,price1,price2,price3" _n
file write fh `"1,$100.00,"$1,000.00","$10,000.00""' _n
file write fh `"2,EUR 50.00,"1,234.56 EUR","EUR9,999.99""' _n
file close fh

* Percentage values
file open fh using "temp/percentage_values.csv", write replace
file write fh "id,pct1,pct2,pct3" _n
file write fh "1,50%,0.50,50.0%" _n
file write fh "2,100%,1.00,100.00%" _n
file write fh "3,0%,0.00,0.0%" _n
file close fh

* Very small dataset (1 row, 1 col)
file open fh using "temp/tiny_1x1.csv", write replace
file write fh "x" _n
file write fh "1" _n
file close fh

* Very small dataset (1 row, many cols)
file open fh using "temp/tiny_1xN.csv", write replace
file write fh "a,b,c,d,e,f,g,h,i,j" _n
file write fh "1,2,3,4,5,6,7,8,9,10" _n
file close fh

* Very small dataset (many rows, 1 col)
file open fh using "temp/tiny_Nx1.csv", write replace
file write fh "x" _n
forvalues i = 1/100 {
    file write fh "`i'" _n
}
file close fh

* Single row CSV
file open fh using "temp/single_row.csv", write replace
file write fh "a,b,c" _n
file write fh "1,2,3" _n
file close fh

* Single column CSV (20 rows)
file open fh using "temp/single_col.csv", write replace
file write fh "value" _n
forvalues i = 1/20 {
    file write fh "`i'" _n
}
file close fh

* All same value
file open fh using "temp/all_same_value.csv", write replace
file write fh "a,b,c" _n
forvalues i = 1/50 {
    file write fh "1,1,1" _n
}
file close fh

* Alternating missing pattern
file open fh using "temp/alternating_missing.csv", write replace
file write fh "a,b,c" _n
forvalues i = 1/20 {
    if mod(`i', 2) == 1 {
        file write fh "`i',,`i'" _n
    }
    else {
        file write fh ",`i'," _n
    }
}
file close fh

* Checkerboard missing pattern
file open fh using "temp/checkerboard_missing.csv", write replace
file write fh "a,b,c,d" _n
forvalues i = 1/10 {
    if mod(`i', 2) == 1 {
        file write fh "`i',,`i'," _n
    }
    else {
        file write fh ",`i',,`i'" _n
    }
}
file close fh

* First column all missing
file open fh using "temp/first_col_missing.csv", write replace
file write fh "empty,b,c" _n
forvalues i = 1/10 {
    file write fh ",`i',`=`i'*2'" _n
}
file close fh

* Last column all missing
file open fh using "temp/last_col_missing.csv", write replace
file write fh "a,b,empty" _n
forvalues i = 1/10 {
    file write fh "`i',`=`i'*2'," _n
}
file close fh

* First row all missing (after header)
file open fh using "temp/first_row_missing.csv", write replace
file write fh "a,b,c" _n
file write fh ",," _n
forvalues i = 2/10 {
    file write fh "`i',`=`i'*2',`=`i'*3'" _n
}
file close fh

* Last row all missing
file open fh using "temp/last_row_missing.csv", write replace
file write fh "a,b,c" _n
forvalues i = 1/9 {
    file write fh "`i',`=`i'*2',`=`i'*3'" _n
}
file write fh ",," _n
file close fh

* Only string data
file open fh using "temp/strings_only.csv", write replace
file write fh "name,city,country" _n
file write fh "John,New York,USA" _n
file write fh "Jane,London,UK" _n
file write fh "Bob,Paris,France" _n
file write fh "Alice,Tokyo,Japan" _n
file close fh

* Mixed numeric and string with ambiguous types
file open fh using "temp/ambiguous_types.csv", write replace
file write fh "id,mixed" _n
file write fh "1,100" _n
file write fh "2,abc" _n
file write fh "3,200" _n
file write fh "4,def" _n
file close fh

* Column that starts numeric then becomes string
file open fh using "temp/type_change.csv", write replace
file write fh "id,value" _n
file write fh "1,100" _n
file write fh "2,200" _n
file write fh "3,three hundred" _n
file write fh "4,400" _n
file close fh

* Very deep file (many rows, few columns)
file open fh using "temp/very_deep.csv", write replace
file write fh "id,value" _n
forvalues i = 1/5000 {
    file write fh "`i',`=runiform()'" _n
}
file close fh

* Square file (equal rows and columns)
file open fh using "temp/square_10x10.csv", write replace
local header ""
forvalues i = 1/10 {
    if `i' > 1 local header "`header',"
    local header "`header'v`i'"
}
file write fh "`header'" _n
forvalues j = 1/10 {
    local line ""
    forvalues i = 1/10 {
        if `i' > 1 local line "`line',"
        local line "`line'`=`j'*`i''"
    }
    file write fh "`line'" _n
}
file close fh

* Sparse data (mostly missing)
file open fh using "temp/sparse_data.csv", write replace
file write fh "a,b,c,d,e" _n
forvalues i = 1/50 {
    local line ""
    forvalues j = 1/5 {
        if `j' > 1 local line "`line',"
        if runiform() < 0.2 {
            local line "`line'`=runiformint(1,100)'"
        }
    }
    file write fh "`line'" _n
}
file close fh

* Dense data (no missing values)
file open fh using "temp/dense_data.csv", write replace
file write fh "a,b,c,d,e" _n
forvalues i = 1/100 {
    file write fh "`i',`=`i'+1',`=`i'+2',`=`i'+3',`=`i'+4'" _n
}
file close fh

* (silent) "  Pathological test data files created."
