/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * Copyright 2004 Komarov Valery
 * Copyright 2006 Christophe Leitienne
 * Copyright 2008-2017 David Hoerl
 * Copyright 2013 Bob Colbert
 * Copyright 2013-2018 Evan Miller
 *
 * This file is part of libxls -- A multiplatform, C/C++ library for parsing
 * Excel(TM) files.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS
 * IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef XLS_STRUCT_INC
#define XLS_STRUCT_INC

#include "../libxls/ole.h"

#define XLS_RECORD_EOF          0x000A
#define XLS_RECORD_DEFINEDNAME  0x0018
#define XLS_RECORD_NOTE         0x001C
#define XLS_RECORD_1904         0x0022
#define XLS_RECORD_FILEPASS     0x002F
#define XLS_RECORD_CONTINUE     0x003C
#define XLS_RECORD_WINDOW1      0x003D
#define XLS_RECORD_CODEPAGE     0x0042
#define XLS_RECORD_OBJ          0x005D
#define XLS_RECORD_MERGEDCELLS  0x00E5
#define XLS_RECORD_DEFCOLWIDTH  0x0055
#define XLS_RECORD_COLINFO      0x007D
#define XLS_RECORD_BOUNDSHEET   0x0085
#define XLS_RECORD_PALETTE      0x0092
#define XLS_RECORD_MULRK        0x00BD
#define XLS_RECORD_MULBLANK     0x00BE
#define XLS_RECORD_RSTRING      0x00D6
#define XLS_RECORD_DBCELL       0x00D7
#define XLS_RECORD_XF           0x00E0
#define XLS_RECORD_MSODRAWINGGROUP   0x00EB
#define XLS_RECORD_MSODRAWING   0x00EC
#define XLS_RECORD_SST          0x00FC
#define XLS_RECORD_LABELSST     0x00FD
#define XLS_RECORD_EXTSST       0x00FF
#define XLS_RECORD_TXO          0x01B6
#define XLS_RECORD_HYPERREF     0x01B8
#define XLS_RECORD_BLANK        0x0201
#define XLS_RECORD_NUMBER       0x0203
#define XLS_RECORD_LABEL        0x0204
#define XLS_RECORD_BOOLERR      0x0205
#define XLS_RECORD_STRING       0x0207 // only follows a formula
#define XLS_RECORD_ROW          0x0208
#define XLS_RECORD_INDEX        0x020B
#define XLS_RECORD_ARRAY        0x0221 // Array-entered formula
#define XLS_RECORD_DEFAULTROWHEIGHT    0x0225
#define XLS_RECORD_FONT         0x0031 // spec says 0x0231 but Excel expects 0x0031
#define XLS_RECORD_FONT_ALT     0x0231
#define XLS_RECORD_WINDOW2      0x023E
#define XLS_RECORD_RK           0x027E
#define XLS_RECORD_STYLE        0x0293
#define XLS_RECORD_FORMULA      0x0006
#define XLS_RECORD_FORMULA_ALT  0x0406 // Apple Numbers bug
#define XLS_RECORD_FORMAT       0x041E
#define XLS_RECORD_BOF          0x0809

#define BLANK_CELL  XLS_RECORD_BLANK  // compat

#if defined(_AIX) || defined(__sun)
#pragma pack(1)
#else
#pragma pack(push, 1)
#endif

typedef struct BOF
{
    XLSWORD id;
    XLSWORD size;
}
BOF;

typedef struct BIFF
{
    XLSWORD ver;
    XLSWORD type;
    XLSWORD id_make;
    XLSWORD year;
    XLSDWORD flags;
    XLSDWORD min_ver;
}
BIFF;

typedef struct WIND1
{
    XLSWORD xWn;
    XLSWORD yWn;
    XLSWORD dxWn;
    XLSWORD dyWn;
    XLSWORD grbit;
    XLSWORD itabCur;
    XLSWORD itabFirst;
    XLSWORD ctabSel;
    XLSWORD wTabRatio;
}
WIND1;

typedef struct BOUNDSHEET
{
    XLSDWORD	filepos;
    XLSBYTE	type;
    XLSBYTE	visible;
    char	name[1];
}
BOUNDSHEET;

typedef struct ROW
{
    XLSWORD	index;
    XLSWORD	fcell; // first cell, 0-indexed
    XLSWORD	lcell; // last cell, 1-indexed
    XLSWORD	height;
    XLSWORD	notused;
    XLSWORD	notused2; //used only for BIFF3-4
    XLSWORD	flags;
    XLSWORD	xf;
}
ROW;

typedef struct COL
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf;
}
COL;


typedef struct FORMULA // BIFF8
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf;
	// next 8 bytes either a IEEE double, or encoded on a byte basis
    XLSBYTE	resid;
    XLSBYTE	resdata[5];
    XLSWORD	res;
    XLSWORD	flags;
    XLSBYTE	chn[4]; // BIFF8
    XLSWORD	len;
    XLSBYTE	value[1]; //var
}
FORMULA;

typedef struct FARRAY // BIFF8
{
    XLSWORD	row1;
    XLSWORD	row2;
    XLSBYTE	col1;
    XLSBYTE	col2;
    XLSWORD	flags;
    XLSBYTE	chn[4]; // BIFF8
    XLSWORD	len;
    XLSBYTE	value[1]; //var
}
FARRAY;

typedef struct RK
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf;
    XLSDWORD   value;
}
RK;

typedef struct MULRK
{
    XLSWORD	row;
    XLSWORD	col;
	struct {
		XLSWORD	xf;
		XLSDWORD   value;
	}		rk[1];
	//XLSWORD	last_col;
}
MULRK;

typedef struct MULBLANK
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf[1];
	//XLSWORD	last_col;
}
MULBLANK;

typedef struct BLANK
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf;
}
BLANK;

typedef struct LABEL
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf;
    XLSBYTE	value[1]; // var
}
LABEL;

typedef struct BOOLERR
{
    XLSWORD    row;
    XLSWORD    col;
    XLSWORD    xf;
    XLSBYTE    value;
    XLSBYTE    iserror;
}
BOOLERR;

typedef struct SST
{
    XLSDWORD	num;
    XLSDWORD	numofstr;
    XLSBYTE	strings[1];
}
SST;

typedef struct XF5
{
    XLSWORD	font;
    XLSWORD	format;
    XLSWORD	type;
    XLSWORD	align;
    XLSWORD	color;
    XLSWORD	fill;
    XLSWORD	border;
    XLSWORD	linestyle;
}
XF5;

typedef struct XF8
{
    XLSWORD	font;
    XLSWORD	format;
    XLSWORD	type;
    XLSBYTE	align;
    XLSBYTE	rotation;
    XLSBYTE	ident;
    XLSBYTE	usedattr;
    XLSDWORD	linestyle;
    XLSDWORD	linecolor;
    XLSWORD	groundcolor;
}
XF8;

typedef struct BR_NUMBER
{
    XLSWORD	row;
    XLSWORD	col;
    XLSWORD	xf;
    double value;
}
BR_NUMBER;

typedef struct COLINFO
{
    XLSWORD	first;
    XLSWORD	last;
    XLSWORD	width;
    XLSWORD	xf;
    XLSWORD	flags;
/* There should be an unused XLSWORD field at the end here. However, some files in
 * the wild report it as a XLSBYTE, which results in a boundary-check parse error.
 * Since the value is ignored anyway, we'll just pretend it was never there.
 *
 * See issue https://github.com/evanmiller/libxls/issues/27
 */
}
COLINFO;

typedef struct MERGEDCELLS
{
    XLSWORD	rowf;
    XLSWORD	rowl;
    XLSWORD	colf;
    XLSWORD	coll;
}
MERGEDCELLS;

typedef struct FONT
{
    XLSWORD	height;
    XLSWORD	flag;
    XLSWORD	color;
    XLSWORD	bold;
    XLSWORD	escapement;
    XLSBYTE	underline;
    XLSBYTE	family;
    XLSBYTE	charset;
    XLSBYTE	notused;
    char    name[1];
}
FONT;

typedef struct FORMAT
{
    XLSWORD	index;
    char	value[1];
}
FORMAT;

#pragma pack(pop)

//---------------------------------------------------------

typedef	struct st_sheet
{
    XLSDWORD count;        // Count of sheets
    struct st_sheet_data
    {
        XLSDWORD filepos;
        XLSBYTE visibility;
        XLSBYTE type;
        char * name;
    }
    * sheet;
}
st_sheet;

typedef	struct st_font
{
    XLSDWORD count;		// Count of FONT's
    struct st_font_data
    {
        XLSWORD	height;
        XLSWORD	flag;
        XLSWORD	color;
        XLSWORD	bold;
        XLSWORD	escapement;
        XLSBYTE	underline;
        XLSBYTE	family;
        XLSBYTE	charset;
        char *	name;
    }
    * font;
}
st_font;

typedef struct st_format
{
    XLSDWORD count;		// Count of FORMAT's
    struct st_format_data
    {
         XLSWORD index;
         char *value;
    }
    * format;
}
st_format;

typedef	struct st_xf
{
    XLSDWORD count;	// Count of XF
    //	XF** xf;
    struct st_xf_data
    {
        XLSWORD	font;
        XLSWORD	format;
        XLSWORD	type;
        XLSBYTE	align;
        XLSBYTE	rotation;
        XLSBYTE	ident;
        XLSBYTE	usedattr;
        XLSDWORD	linestyle;
        XLSDWORD	linecolor;
        XLSWORD	groundcolor;
    }
    * xf;
}
st_xf;


typedef	struct st_sst
{
    XLSDWORD count;
    XLSDWORD lastid;
    XLSDWORD continued;
    XLSDWORD lastln;
    XLSDWORD lastrt;
    XLSDWORD lastsz;
    struct str_sst_string
    {
        char * str;
    }
    * string;
}
st_sst;


typedef	struct st_cell
{
    XLSDWORD count;
    struct st_cell_data
    {
        XLSWORD	id;
        XLSWORD	row;
        XLSWORD	col;
        XLSWORD	xf;
        char *	str;		// String value;
        double	d;
        int32_t	l;
        XLSWORD	width;		// Width of col
        XLSWORD	colspan;
        XLSWORD	rowspan;
        XLSBYTE	isHidden;	// Is cell hidden
    }
    * cell;
}
st_cell;


typedef	struct st_row
{
    //	XLSDWORD count;
    XLSWORD lastcol;	// numCols - 1
    XLSWORD lastrow;	// numRows - 1
    struct st_row_data
    {
        XLSWORD index;
        XLSWORD fcell;
        XLSWORD lcell;
        XLSWORD height;
        XLSWORD flags;
        XLSWORD xf;
        XLSBYTE xfflags;
        st_cell cells;
    }
    * row;
}
st_row;


typedef	struct st_colinfo
{
    XLSDWORD count;				// Count of COLINFO
    struct st_colinfo_data
    {
        XLSWORD	first;
        XLSWORD	last;
        XLSWORD	width;
        XLSWORD	xf;
        XLSWORD	flags;
    }
    * col;
}
st_colinfo;

typedef struct xlsWorkBook
{
    //FILE*		file;
    OLE2Stream*	olestr;
    int32_t		filepos;		// position in file

    //From Header (BIFF)
    XLSBYTE		is5ver;
    XLSBYTE		is1904;
    XLSWORD		type;
    XLSWORD		activeSheetIdx;	// index of the active sheet

    //Other data
    XLSWORD		codepage;		// Charset codepage
    char*		charset;
    st_sheet	sheets;
    st_sst		sst;			// SST table
    st_xf		xfs;			// XF table
    st_font		fonts;
    st_format	formats;		// FORMAT table

	char		*summary;		// ole file
	char		*docSummary;	// ole file

    void        *converter;
    void        *utf16_converter;
    void        *utf8_locale;
}
xlsWorkBook;

typedef struct xlsWorkSheet
{
    XLSDWORD		filepos;
    XLSWORD		defcolwidth;
    st_row		rows;
    xlsWorkBook *workbook;
    st_colinfo	colinfo;
}
xlsWorkSheet;

#ifdef __cplusplus
typedef struct st_cell::st_cell_data xlsCell;
typedef	struct st_row::st_row_data xlsRow;
#else
typedef struct st_cell_data xlsCell;
typedef	struct st_row_data xlsRow;
#endif

typedef struct xls_summaryInfo
{
	XLSBYTE		*title;
	XLSBYTE		*subject;
	XLSBYTE		*author;
	XLSBYTE		*keywords;
	XLSBYTE		*comment;
	XLSBYTE		*lastAuthor;
	XLSBYTE		*appName;
	XLSBYTE		*category;
	XLSBYTE		*manager;
	XLSBYTE		*company;
}
xlsSummaryInfo;

typedef void (*xls_formula_handler)(XLSWORD bof, XLSWORD len, XLSBYTE *formula);

#endif
