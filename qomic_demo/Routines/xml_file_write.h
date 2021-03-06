/*
 *  This file is part of FREGENE program.
 *  Copyright (c) Clive J Hoggart   <c.hoggart@imperial.ac.uk>
 *                Taane G Clark     <tgc@well.ox.ac.uk>
 *                Marc Chadeau      <m.chadeau@imperial.ac.uk>
 *
 *  FREGENE is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 2, as
 *  published by the Free Software Foundation.
 *
 *  FREGENE is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */ 



#ifndef XML_FILE_WRITE_H

#define XML_FILE_WRITE_H

#define MA_NO_ERROR			0 
#define MA_NO_ACTION			2  
#define MA_CLOSE_FILE_ERROR	   15
#define MA_XML_FILE_ERROR			-13

#include <stdio.h>

class MaXmlTagWrite
{
protected :
	FILE *XF;

public :
	//MaXmlTagWrite();
	MaXmlTagWrite();
	MaXmlTagWrite(FILE *f, char version[10]="1.0");
	MaXmlTagWrite(char *s, char version[10]="1.0");

	~MaXmlTagWrite();

	FILE *GetXmlFile() {return XF;}
	void SetXmlFile(FILE *f); 
	FILE *OpenXmlFile(char *s);
	int  CloseXmlFile(FILE *f);
	//writing xml declaration
	int WriteDeclaration(char version[10]="1.0");
	

	//writing xml DocType
	int WriteDocType();

	
	//writing xml comment
	int WriteComment(char *s);

	//writing string : redefinition of fwrite
	int WriteString(char *s);


	//writing Space
	int WriteSpace();

	//writing Tabulation
	int WriteTab();

	//writing line
	int WriteLine(char *s);

	//writing Empty Element
	int WriteEmptyElement(char *s);


	//writing begin tag
	int StartTag(char *s);
	//wrinting end tag
	int EndTag(char *s);

	//writing begin tag CDATA
	int StartCDATA();
	//wrinting end tag CDATA
	int EndCDATA();



	//writing begin tag head
	int StartHead();
	//wrinting end tag head
	int EndHead();

	//writing begin tag manip name
	int StartManip();
	//wrinting end tag manip name
	int EndManip();


	//writing begin tag Data in case of binary
	int StartBinData();
	//wrinting end tag Data in case of binary
	int EndBinData();

	//writing begin tag Data in case of ascii
	int StartAData();
	//wrinting end tag Data in case of ascii
	int EndAData();


	//write fixed tag
	int WriteAuthor(char *s);
	int WriteDate();
	int WriteDataFormat();           //Binaries Data or ascii data
	int WriteType(char *c);			   //Type de fichier
	int WriteDescription(char *des); //Description de ce que contient le champs DATA
	
};

//========================================================================
//========================================================================>>>



#endif 

