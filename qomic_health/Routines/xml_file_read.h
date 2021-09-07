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

#ifndef XML_FILE_READ_H

#define XML_FILE_READ_H

#define MA_NO_ERROR 0 
#define MA_NO_ACTION 2  
#define MA_CLOSE_FILE_ERROR 15
#define MA_XML_FILE_ERROR -13

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifndef Mymax
#define Mymax(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef Mymin
#define Mymin(a,b)            (((a) < (b)) ? (a) : (b))
#endif

class MaXmlTagRead
{
   protected :
      FILE *XF;
   
   
   public :
      MaXmlTagRead(FILE *f);
   MaXmlTagRead(char *s);
   ~MaXmlTagRead();
   
   void  SetXmlFile(FILE *f); 
   FILE* OpenXmlFile(char *s);
   int  CloseXmlFile(FILE *f);
   int  CloseXmlFile();
                
   bool GetTagPos(const char *tag, long start, long, long &pos_found);
   bool GetTagFirstPos(const char *tag, long start, long, long &first_pos);
   bool ReadTag(const char *tag, long start, long end, char *content, unsigned long content_size);
   bool ReadTag(const char *tag, long start, long end, char *content, unsigned long content_size, long &pos_found);
   bool ReadSizeTag(char *tag, long start, long end, unsigned long &content_size);
   char *GetTagEndOf(char *tag);
   
};


#endif
