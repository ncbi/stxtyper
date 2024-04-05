// stxtyper.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   stx typing using protein reference sequences
*   
* Dependencies: NCBI BLAST, gunzip (optional)
*
* Release changes:
*  1.0.19 03/26/2024          BlastAlignment::targetAlign is removed
*  1.0.18 03/19/2024 PD-4910  Element symbol is <stx type>_operon, Element name contains operon quality attribute"

          Sequence name -> Element name in header
          Sequence name, now Element name, should be type/subtype and include info when not complete e.g.,:
              stx1a operon
              stx2c operon
              Partial stx2 operon
              stx2 operon with frameshift
              Novel stx2 operon
              stx2 operon with internal stop
          Element subtype should be STX_TYPE
          Subclass should be to type only where element symbol is to type only
              E.g. a partial stx2a should have Element symbol of stx2 and a subclass of STX2
          Target length should be in nucleotide sequence coordinates
          Reference sequence length and % coverage of reference sequence should be blank

*   1.0.17 03/18/2024 PD-4910  
*   1.0.16 03/13/2024 PD-4926  --amrfinder: <stx type>_operon, "Gene symbol" -> "Element symbol"
*   1.0.15 03/11/2024 PD-4924  dead stxA2j EFK5907329.1 is replaced by EMA1832120.1
*   1.0.14 03/05/2024 PD-4918  --print_node: print AMRFinderPlus hierarchy node
*   1.0.13 03/05/2024 PD-4910  --amrfinder prints output in the AMRFinderPlus format
*   1.0.12 02/29/2024          TsvOut.live() is used
*   1.0.11 02/28/2024 PD-4911  wrong QC for log file output
*                     PD-4897  single-subunit operons are de-redundified
*   1.0.10 02/22/2024 PD-4901  multiple types for the same protein sequence are allowed; only Flemming's reference proteins are used
*   1.0.9  02/16/2024 PD-4901  new database (includes Flemming's data)
*   1.0.8  02/16/2024 PD-4892, PD-4898  steps: (1) find operons where A subtype = B subtype and operon identity >= threshold
*                                              (2) find other operons where operon identity >= threshold for each subunit
*                                              (3) find other operons with relaxed intergenic region size
*   1.0.7  02.15/2024 PD-4897  extend intergenic region for partial operons
*   1.0.6  02/13/2024 PD-4874  --translation_table is removed
*          02/13/2024 PD-4894  EXTENDED operon type
*          02/13/2024 PD-4892  Subunits A and B are not preferred to be of the same stx class
*   1.0.5  02/12/2024 PD-4891  -v == --version
*                              stable choice of a blast hit among equivalent ones
*   1.0.4  02/08/2024 PD-4891  stable choice of a blast hit among equivalent ones
*          02/08/2024 PD-4887  a stand-alone blast hit suppresses a covered blast hit if the reported gene symbol is the same
*                              "STANDARD" -> "COMPLETE"
*                              PARTIAL_CONTIG_END is detected for separate A and B subunits
*   1.0.3  02/07/2024 PD-4874  strand is '+' or '-'
*   1.0.2  02/06/2024          blastx -> tblastn
*   1.0.1  02/05/2024 PD-4874  github.com/vbrover/stxtyper
*   1.0.0  11/21/2023 PD-4798
*
*/
   
   
#ifdef _MSC_VER
  #error "UNIX is required"
#endif

#undef NDEBUG 

#include "common.hpp"
#include "tsv.hpp"
using namespace Common_sp;

#include "common.inc"



namespace 
{


string input_name;
bool amrfinder = false;
bool print_node = false;
map<string,double> stxClass2identity;

// PAR
constexpr size_t intergenic_max {36};  // Max. intergenic region in the reference set + 2
constexpr size_t slack = 30;  

const string stxS ("stx");
const string na ("na");



string stxType_reported_operon2elementName (const string &stxType_reported,
                                            const string &operon)
{
  string elementName (stxType_reported  + " operon");
       if (operon == "FRAMESHIFT")
    elementName += " with frameshift";
  else if (operon == "INTERNAL_STOP")
    elementName += " with internal stop";
  else if (contains (operon, "PARTIAL"))
    elementName = "Partial " + elementName;
  else if (operon == "EXTENDED")
    elementName = "Extended " + elementName;
  else if (contains (operon, "NOVEL"))
    elementName = "Novel " + elementName;
    
  return elementName;
}



struct BlastAlignment 
{
  size_t length {0}, nident {0}  // aa
       ,    refStart {0},    refEnd {0},    refLen {0}
       , targetStart {0}, targetEnd {0}, targetLen {0};
    // Positions are 0-based
    // targetStart < targetEnd
  bool stopCodon {false};
  bool frameshift {false};

  // target    
  string targetName; 
  string targetSeq;  
  bool targetStrand {true}; 
    // false <=> negative
//size_t targetAlign {0};
    // bp
  
  // Reference
  // Whole sequence ends with '*'
  string refAccession;
  string refSeq;
  // Function of refAccession
  string stxType;
  string stxClass;
    // Function of stxType
  string stxSuperClass;
    // Function of stxClass
  char subunit {'\0'};
    // 'A' or 'B'
  
  bool reported {false};


  BlastAlignment (const string &line)
    {
      {
        string sseqid;
        {
    	    istringstream iss (line);
    	    iss >> targetName >> sseqid >> targetStart >> targetEnd >> targetLen >> refStart >> refEnd >> refLen >> targetSeq >> refSeq;
  	  // format:  qseqid       sseqid    qstart         qend         qlen         sstart      send      slen      qseq         sseq
      // blast:                          62285          63017        88215        105         837       837          
        }
  	    QC_ASSERT (! targetSeq. empty ());	
  	    
  	    string famId;
        try
        {	
  		    famId        = rfindSplit (sseqid, '|');  
  		    refAccession = rfindSplit (sseqid, '|');
  		  }
  		  catch (const exception &e)
  		  {
  		  	throw runtime_error (string ("Bad StxTyper database\n") + e. what () + "\n" + line);
  		  }
        QC_ASSERT (famId. size () == 6);
        QC_ASSERT (isLeft (famId, stxS));
        subunit = famId [3];
        stxType = famId. substr (4);
      }
      ASSERT (stxType. size () == 2);
      	      
	    length = targetSeq. size ();
	    nident = 0;
	    QC_ASSERT (targetSeq. size () == refSeq. size ());
	    FFOR (size_t, i, targetSeq. size ())
	      if (targetSeq [i] == refSeq [i])
	        nident++;

      stxClass = stxType;
      if (   stxType == "2a"
          || stxType == "2c"
          || stxType == "2d"
         )
        stxClass = "2";
      
      stxSuperClass = stxClass. substr (0, 1);

	    QC_ASSERT (refStart < refEnd);

	    QC_ASSERT (targetStart != targetEnd);
	    targetStrand = targetStart < targetEnd;  
	    if (! targetStrand)
	      swap (targetStart, targetEnd);
	    
	    QC_ASSERT (refStart >= 1);
	    QC_ASSERT (targetStart >= 1);
	    refStart--;
	    targetStart--;
	    
    //targetAlign = targetEnd - targetStart;
    //QC_ASSERT (targetAlign_aa % 3 == 0);
    //targetAlign_aa /= 3;
	    
	    const size_t stopCodonPos = targetSeq. find ('*');
	    if (stopCodonPos != string::npos && stopCodonPos + 1 < targetSeq. size ())
	      stopCodon = true;
    }
  void qc () const
    {
      if (! qc_on)
        return;
      QC_ASSERT (length);
      QC_ASSERT (nident);
      QC_ASSERT (nident <= length);
	    QC_ASSERT (targetStart < targetEnd);
	    QC_ASSERT (targetEnd <= targetLen);
      QC_ASSERT (refStart < refEnd);
	    QC_ASSERT (refEnd <= refLen);
	    if (! frameshift)
	    {
  	    QC_ASSERT (nident <= refEnd - refStart);
  	    QC_ASSERT (refEnd - refStart <= length);	    
  	  }
      QC_ASSERT (! targetName. empty ());
      QC_ASSERT (contains (stxClass2identity, stxClass));
      QC_ASSERT (isLeft (stxType, stxClass));
      QC_ASSERT (subunit == 'A' || subunit == 'B');
      QC_ASSERT (! refAccession. empty ());
      QC_ASSERT (! targetSeq. empty ());
      QC_ASSERT (! refSeq. empty ());
      QC_ASSERT (targetSeq. size () == refSeq. size ());
      QC_IMPLY (! frameshift, length == targetSeq. size ());
      QC_ASSERT (stxType. size () == 2);
    }
  void saveTsvOut (TsvOut& td,
                   bool verboseP) const 
    { if (! td. live ())
        return;
      const string stxType_reported (verboseP ? getGenesymbol () : (stxS + stxType. substr (0, 1)));
      const string operon (frameshift 
                             ? "FRAMESHIFT"
                             : stopCodon 
                               ? "INTERNAL_STOP"
                               : truncated () || otherTruncated ()
                                 ? "PARTIAL_CONTIG_END"
                                 : verboseP && getRelCoverage () == 1.0
                                   ? "COMPLETE_SUBUNIT"
                                   : getExtended ()
                                     ? "EXTENDED"
                                     : "PARTIAL"
                          );
      const char strand (targetStrand ? '+' : '-');
      const double refCoverage = getRelCoverage () * 100.0;
      const double refIdentity = getIdentity ()    * 100.0; 
      // td     
      if (! input_name. empty ())
  	    td << input_name;
      if (amrfinder)
      {
        const string subunitS (1, subunit);
        string subclass (stxType_reported /*stxS + stxType*/);
        strUpper (subclass);
        td << na               // 1 "Protein identifier"  
           << targetName       // 2 "Contig id"
           << targetStart + 1  // 3 "Start"
           << targetEnd        // 4 "Stop"
           << strand           // 5 "Strand"
           << stxType_reported + "_operon" // 6 "Element symbol"
           << stxType_reported_operon2elementName (stxType_reported, operon)    // 7 "Element name"
           << "plus"           // 8 "Scope"
           << "VIRULENCE"      // 9 "Element type"
           << "STX_TYPE"       //10 "Element subtype"
           << subclass. substr (0, 4)   //11 "Class"
           << subclass         //12 "Subclass"
           << operon           //13 "Method"  
           << targetEnd - targetStart /*targetAlign*/      //14 "Target length" 
           << noString /*refLen*/  //15 "Reference sequence length"
           << noString /*refCoverage*/      //16 "% Coverage of reference sequence"
           << refIdentity      //17 "% Identity to reference sequence"
           << length           //18 "Alignment length"
           << refAccession     //19 "Accession of closest sequence"
           << "Shiga toxin " + stxS + stxType + " subunit " + subunitS //20 "Name of closest sequence"
           << na               //21 "HMM id"
           << na               //22 "HMM description"
           ;
        if (print_node)
          td << getGenesymbol ();
      }
      else
      {
        td << targetName
           << stxType_reported
           << operon
           << noString
           << targetStart + 1
           << targetEnd
           << strand;
        if (subunit == 'B')
          td << noString
             << noString
             << noString;
        td << refAccession
           << refIdentity
           << refCoverage;
        if (subunit == 'A')
          td << noString
             << noString
             << noString;
      }
      td. newLn ();
    }
    

  string getGenesymbol () const
    { return stxS + subunit + stxType; }
  void merge (const BlastAlignment &prev)
    { ASSERT (targetName   == prev. targetName);
      ASSERT (refAccession == prev. refAccession);
      ASSERT (targetStrand == prev. targetStrand);
      ASSERT (targetLen    == prev. targetLen);
      ASSERT (refLen       == prev. refLen);
      ASSERT (targetStart > prev. targetStart);
      targetStart = prev. targetStart;
      if (targetStrand)
        refStart = prev. refStart;
      else
        refEnd = prev. refEnd;
      length += prev. length;  // Approximately
      nident += prev. nident;  // Approximately
    //targetAlign += prev. targetAlign;
      if (prev. stopCodon)
        stopCodon = true;
      frameshift = true;
    }
  size_t getFrame () const
    { return (targetStart % 3) + 1; }
  double getIdentity () const 
    { return (double) nident / (double) (length); }
  size_t getAbsCoverage () const 
    { return refEnd - refStart; }
  double getRelCoverage () const 
    { return (double) getAbsCoverage () / (double) refLen; }
  size_t getDiff () const
    { return refStart + (refLen - refEnd) + (length - nident); }
  bool truncated () const
    { return    (targetStart           <= 3 /*Locus::end_delta*/ && ((targetStrand && refStart)            || (! targetStrand && refEnd + 1 < refLen)))
             || (targetLen - targetEnd <= 3 /*Locus::end_delta*/ && ((targetStrand && refEnd + 1 < refLen) || (! targetStrand && refStart)));
    }
  bool otherTruncated () const
    { constexpr size_t missed_max = intergenic_max + 3 * 20 /*min. domain length*/;  // PAR
      return    (targetStrand == (subunit == 'B') && targetStart           <= missed_max)
             || (targetStrand == (subunit == 'A') && targetLen - targetEnd <= missed_max);
    }
  bool getExtended () const
    { ASSERT (! truncated ());
      return ! refStart && refEnd + 1 == refLen; 
    }
  bool insideEq (const BlastAlignment &other) const
    { return    targetStart >= other. targetStart 
             && targetEnd   <= other. targetEnd;
    }  
  string refMap (size_t len) const
    { QC_ASSERT (refLen <= len);
      string s = string (refStart, '-'); 
      FFOR (size_t, i, refSeq. size ())
        if (refSeq [i] != '-')
          s += targetSeq [i];
      s += string (len - refEnd, '-'); 
      return s;
    }
  static bool frameshiftLess (const BlastAlignment* a,
                              const BlastAlignment* b)
    { ASSERT (a);
      ASSERT (b);
      ASSERT (! a->reported);
      ASSERT (! b->reported);
      LESS_PART (*a, *b, targetName);
      LESS_PART (*a, *b, targetStrand);
      LESS_PART (*a, *b, refAccession);
      LESS_PART (*a, *b, targetStart);
      LESS_PART (*a, *b, targetEnd);
      return false;
    }
  static bool sameTypeLess (const BlastAlignment* a,
                            const BlastAlignment* b) 
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, targetName);
      LESS_PART (*a, *b, targetStrand);
      LESS_PART (*a, *b, stxClass);
      LESS_PART (*a, *b, subunit);
      LESS_PART (*a, *b, targetStart);
      LESS_PART (*a, *b, getDiff ());
      LESS_PART (*a, *b, refAccession);
      return false;
    }
  static bool less (const BlastAlignment* a,
                    const BlastAlignment* b) 
    // = sameTypeLess(), but without stxClass
    { ASSERT (a);
      ASSERT (b);
    //LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, targetName);
      LESS_PART (*a, *b, targetStrand);
      LESS_PART (*a, *b, subunit);
      LESS_PART (*a, *b, targetStart);
      LESS_PART (*a, *b, getDiff ());
      LESS_PART (*a, *b, refAccession);
      return false;
    }
  static bool reportLess (const BlastAlignment* a,
                          const BlastAlignment* b)  
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, targetName);
      LESS_PART (*a, *b, targetStrand);
      LESS_PART (*b, *a, getAbsCoverage ());
      LESS_PART (*a, *b, getDiff ());
      LESS_PART (*a, *b, targetStart);
      LESS_PART (*a, *b, refAccession);
      return false;
    }
};



struct Operon
{
  const BlastAlignment* al1 {nullptr};
    // !nullptr
  const BlastAlignment* al2 {nullptr};
  // al1->targetEnd < al2->targetStart
  

  Operon () = default;
  Operon (const BlastAlignment& al1_arg,
          const BlastAlignment& al2_arg)
    : al1 (& al1_arg)
    , al2 (& al2_arg)
    {}
  explicit Operon (const BlastAlignment& al1_arg)
    : al1 (& al1_arg)
    {}
  void qc () const
    { if (! qc_on)
        return;
      QC_ASSERT (al1);
      al1->qc ();
      QC_ASSERT (al1->reported);
      if (! al2)
        return;
      al2->qc ();
      QC_ASSERT (al1->targetName   == al2->targetName);
      QC_ASSERT (al1->targetStrand == al2->targetStrand);
      QC_ASSERT (al1->targetEnd    <  al2->targetStart);
      QC_ASSERT (al1->subunit      != al2->subunit);
      QC_ASSERT (al2->reported);
    }
  void saveTsvOut (TsvOut& td,
                   bool verboseP) const 
    { ASSERT (al1);
      if (! td. live ())
        return;
      if (al2)
      {
        string stxType (getStxType (verboseP));
        const string standard ("COMPLETE");
        const bool novel =    al1->stxClass != al2->stxClass 
                           || getIdentity () < stxClass2identity [al1->stxClass]
                           || stxType. size () <= 1;
        const string operonType =    getA () -> frameshift
                                  || getB () -> frameshift
                                    ? "FRAMESHIFT"
                                    :    getA () -> stopCodon 
                                      || getB () -> stopCodon 
                                      ? "INTERNAL_STOP"
                                      :    getA () -> truncated () 
                                        || getB () -> truncated ()
                                        ? "PARTIAL_CONTIG_END"
                                        : partial ()
                                          ? "PARTIAL"  
                                          :    getA () -> getExtended ()
                                            || getB () -> getExtended ()
                                            ? "EXTENDED"
                                            : novel
                                              ? standard + "_NOVEL"
                                              : standard;
        if (! verboseP)
        {
          ASSERT (stxType. size () <= 2);
          if (operonType == standard)
            { ASSERT (stxType. size () == 2); }
          else
            if (stxType. size () == 2)
              stxType. erase (1);  
        }
        const string targetName (al1->targetName);
        const size_t start = al1->targetStart + 1;
        const size_t stop  = al2->targetEnd;
  	    const char strand (al1->targetStrand ? '+' : '-');
  	    const string stxType_reported (stxS + stxType);
  	    const double refIdentity = getIdentity () * 100.0;
  	    // td
        if (! input_name. empty ())
  	      td << input_name;
        if (amrfinder)
        {
          const string genesymbol (al1->stxType == al2->stxType ? stxS + al1->stxType : stxType_reported);
          string subclass (stxType_reported /*genesymbol*/);
          strUpper (subclass);
          const size_t targetAlign = al2->targetEnd - al1->targetStart;
        //const size_t refLen = al1->refLen + al2->refLen;
        //const double refCoverage = double (al1->getAbsCoverage () + al2->getAbsCoverage ()) / double (refLen) * 100.0;
          const size_t alignmentLen = al1->length + al2->length;
          const string refAccessions (al1->refAccession + ", " + al2->refAccession);
          const string fam (al1->getGenesymbol () + ", " + al2->getGenesymbol ());
          td << na                // 1 "Protein identifier"  
             << targetName        // 2 "Contig id"
             << start             // 3 "Start"
             << stop              // 4 "Stop"
             << strand            // 5 "Strand"
             << stxType_reported + "_operon"  // 6 "Element symbol"
             << stxType_reported_operon2elementName (stxType_reported, operonType)    // 7 "Element name"
             << "plus"            // 8 "Scope"
             << "VIRULENCE"       // 9 "Element type"
             << "STX_TYPE"        //10 "Element subtype"
             << subclass. substr (0, 4)   //11 "Class"
             << subclass          //12 "Subclass"
             << operonType        //13 "Method"  
             << targetAlign       //14 "Target length" 
             << noString /*refLen*/  //15 "Reference sequence length"
             << noString /*refCoverage*/  //16 "% Coverage of reference sequence"
             << refIdentity       //17 "% Identity to reference sequence"
             << alignmentLen      //18 "Alignment length"
             << refAccessions     //19 "Accession of closest sequence"
             << "Shiga toxin " + genesymbol //20 "Name of closest sequence"
             << na                //21 "HMM id"
             << na                //22 "HMM description"
             ;
          if (print_node)
            td << fam;
        }
        else
    	    td << targetName
             << stxType_reported
             << operonType
             << refIdentity
             << start
             << stop
             << strand
             // Approximately if frameshift
             << getA () -> refAccession
             << getA () -> getIdentity () * 100.0
             << getA () -> getRelCoverage () * 100.0
             << getB () -> refAccession
             << getB () -> getIdentity () * 100.0
             << getB () -> getRelCoverage () * 100.0
             ;
        td. newLn ();
      }
      else
        al1->saveTsvOut (td, verboseP);
    }


private:
  const BlastAlignment* getA () const
    { return al1->targetStrand ? al1 : al2; }
  const BlastAlignment* getB () const
    { return al1->targetStrand ? al2 : al1; }
  bool hasAl2 () const
    { return al2; }
  string getRefAccession2 () const
    { if (al2)
        return al2->refAccession;
      return noString;
    }
  string getStxType (bool verboseP) const
    { if (! al2)
        return al1->stxType;
      if (al1->stxClass != al2->stxClass)
      {
      //return al1->stxClass + "/" + al2->stxClass;  ??  // order alphabetically
        if (al1->stxSuperClass == al2->stxSuperClass)
          return al1->stxSuperClass;  
        return noString;
      }
      if (al1->stxClass != "2")
        return al1->stxType; 
      const string a (getA () -> refMap (319 + 1));
      const string b (getB () -> refMap ( 89 + 1));
      if (   (a [312] == 'F' || a [312] == 'S')
          && (a [318] == 'K' || a [318] == 'E')
          && b [34] == 'D'
         )
        return "2a";
      if (    a [312] == 'F'
          && (a [318] == 'K' || a [318] == 'E')
          && b [34] == 'N'
         )
        return "2c";
      if (   a [312] == 'S'
          && a [318] == 'E'
          && b [34] == 'N'
         )
        return "2d";
      if (verboseP)
        return string ("2 ") + a [312] + a [318] + b [34];
      return "2";
    }
  bool partial () const
    { return    (getA () -> getRelCoverage () < 1.0 && ! getA () -> getExtended ())
             || (getB () -> getRelCoverage () < 1.0 && ! getB () -> getExtended ());
    }
public:
  double getIdentity () const
    { return double (al1->nident + al2->nident) / double (al1->length/*refLen*/ + al2->length/*refLen*/); }
  bool insideEq (const Operon &other) const
    { return    al1->targetStrand        == other. al1->targetStrand
    	       && al1->targetStart + slack >= other. al1->targetStart 
             && al2->targetEnd           <= other. al2->targetEnd + slack;
    }
  bool operator< (const Operon &other) const
    { LESS_PART (*this, other, al1->targetName);
      LESS_PART (other, *this, getIdentity ());
      LESS_PART (*this, other, hasAl2 ());
      LESS_PART (*this, other, al1->refAccession);
      LESS_PART (*this, other, hasAl2 ());
      LESS_PART (*this, other, getRefAccession2 ());
      return false;
    }
  static bool reportLess (const Operon &a,
                          const Operon &b)
    { LESS_PART (a, b, al1->targetName);
      LESS_PART (a, b, al1->targetStart);
      LESS_PART (a, b, al1->targetEnd);
      LESS_PART (b, a, al1->targetStrand);
      LESS_PART (a, b, al1->refAccession);
      LESS_PART (a, b, hasAl2 ());
      LESS_PART (a, b, getRefAccession2 ());
      return false;
    }
};



void goodBlasts2operons (const VectorPtr<BlastAlignment> &goodBlastAls, 
                         Vector<Operon> &operons, 
                         bool sameType,
                         bool strong,
                         TsvOut &logTd)
{
  IMPLY (sameType, strong);
  
  LOG ("\nGood blasts:");

  size_t start = 0;
  FFOR (size_t, i, goodBlastAls. size ())
  {
    const BlastAlignment* alB = goodBlastAls [i];
    ASSERT (alB);
    if (alB->reported)
      continue;
    alB->saveTsvOut (logTd, true);
    if (alB->subunit != 'B')
      continue;
    while (   start < i 
           && ! (   goodBlastAls [start] -> targetName   == alB->targetName
                 && goodBlastAls [start] -> targetStrand == alB->targetStrand
                 && (   ! sameType 
                     || goodBlastAls [start] -> stxClass == alB->stxClass
                    )
                )
          )
      start++;
    FOR_START (size_t, j, start, i)
    {
      const BlastAlignment* alA = goodBlastAls [j];
      ASSERT (alA);
      if (alA->reported)
        continue;
      ASSERT (alA->targetName   == alB->targetName);
      ASSERT (alA->targetStrand == alB->targetStrand);
      IMPLY (sameType, alA->stxClass == alB->stxClass);
      ASSERT (alA->subunit <= alB->subunit);
      if (alA->subunit == alB->subunit)
        break;
      ASSERT (alA->subunit == 'A');
      const BlastAlignment* al1 = alA;
      const BlastAlignment* al2 = alB;
      if (! al1->targetStrand)
        swap (al1, al2);
      if (   al1->targetEnd <= al2->targetStart  
          && al2->targetStart - al1->targetEnd <= intergenic_max * (strong ? 1 : 2)  // PAR  // PD-4897
         )
      {
        Operon op (*al1, *al2);
        LOG ("Operon:\t" + to_string (op. getIdentity ()) + "\t" + to_string (stxClass2identity [op. al1->stxClass]));
        op. saveTsvOut (logTd, true);  
        if (   ! strong 
            || (   op. getIdentity () >= stxClass2identity [op. al1->stxClass]
                && op. getIdentity () >= stxClass2identity [op. al2->stxClass]
               )
           )
        {
          operons << std::move (op);
          var_cast (al1) -> reported = true;
          var_cast (al2) -> reported = true;
        }
      }
    }
  }
  
  LOG ("# Operons: " + to_string (operons. size ()));
  LOG ("\nSuppress goodBlastAls by operons");

  for (const BlastAlignment* al : goodBlastAls)
  {
    ASSERT (al);
    if (! al->reported)
      for (const Operon& op : operons)
      {
        ASSERT (op. al2);
        if (   al->targetName          == op. al1->targetName
            && al->targetStart + slack >= op. al1->targetStart 
            && al->targetEnd           <= op. al2->targetEnd + slack
            && al->targetStrand        == op. al1->targetStrand
           ) 
        {
          var_cast (al) -> reported = true;
          break;
        }
      }
  }
}



// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Determine stx type(s) of a genome, print .tsv-file", true, false, true, true)
    {
    	addKey ("nucleotide", "Input nucleotide FASTA file (can be gzipped)", "", 'n', "NUC_FASTA");
    //addKey ("translation_table", "NCBI genetic code for translated BLAST", "11", 't', "TRANSLATION_TABLE");
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report, for example it can be an assembly name", "", '\0', "NAME");
      addKey ("output", "Write output to OUTPUT_FILE instead of STDOUT", "", 'o', "OUTPUT_FILE");
    	addKey ("blast_bin", "Directory for BLAST. Deafult: $BLAST_BIN", "", '\0', "BLAST_DIR");
    	addFlag ("amrfinder", "Print output in the nucleotide AMRFinderPlus format");
    	addFlag ("print_node", "Print AMRFinderPlus hierarchy node");

      version = SVN_REV;
    }



  void shellBody () const final
  {
    const string fName      = shellQuote (getArg ("nucleotide"));
    const uint   gencode    =             /*arg2uint ("translation_table")*/ 11; 
                 input_name =             getArg ("name");
    const string output     =             getArg ("output");
          string blast_bin  =             getArg ("blast_bin");
                 amrfinder  =             getFlag ("amrfinder");
                 print_node =             getFlag ("print_node");
    
    if (contains (input_name, '\t'))
      throw runtime_error ("NAME cannot contain a tab character");
    if (print_node && ! amrfinder)
      throw runtime_error ("--print_node requires --amrfinder");


    stderr << "Software directory: " << shellQuote (execDir) << '\n';
    stderr << "Version: " << version << '\n'; 
    
		const string logFName (tmp + "/log"); 
    const string qcS (qc_on ? " -qc" : "");


    #define BLASTX 0

    // blast_bin
    if (blast_bin. empty ())
    	if (const char* s = getenv ("BLAST_BIN"))
    		blast_bin = string (s);
    if (! blast_bin. empty ())
    {
	    addDirSlash (blast_bin);
	  #if BLASTX
	    prog2dir ["blastx"]      = blast_bin;
	  #else
	    prog2dir ["tblastn"]     = blast_bin;  
	    prog2dir ["makeblastdb"] = blast_bin;  
	  #endif
	  }

    const string dna_flat = uncompress (fName,  "dna_flat");
    
  #if BLASTX
    size_t nDna = 0;
    size_t dnaLen_max = 0;
    size_t dnaLen_total = 0;
  #endif
    {
      prog2dir ["fasta_check"] = execDir;
      exec (fullProg ("fasta_check") + dna_flat + "  -hyphen  -ambig  " + qcS + "  -log " + logFName + " > " + tmp + "/nseq", logFName); 
    	const StringVector vec (tmp + "/nseq", (size_t) 10, true); 
    	if (vec. size () != 3)
        throw runtime_error ("fasta_check failed: " + vec. toString ("\n"));
    #if BLASTX
      nDna         = str2<size_t> (vec [0]);
      dnaLen_max   = str2<size_t> (vec [1]);
      dnaLen_total = str2<size_t> (vec [2]);
    #endif
    }
  #if BLASTX
    QC_ASSERT (nDna);
    QC_ASSERT (dnaLen_max);
    QC_ASSERT (dnaLen_total);
  #endif

	//stderr. section ("Running blast");
	  const string blastOut (tmp + "/blast");
		{
			const Chronometer_OnePass_cerr cop ("blast");
 			// Database: created by ~brovervv/code/database/stx.prot.sh
    #if BLASTX
 			findProg ("blastx");
  		const string blast_fmt ("-outfmt '6 qseqid sseqid qstart qend qlen sstart send slen qseq sseq'");
			exec (fullProg ("blastx") + " -query " + dna_flat + " -db " + execDir + "stx.prot  " 
			      + "-comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000  -word_size 5  -query_gencode " + to_string (gencode) + " "
			      + getBlastThreadsParam ("blastx", min (nDna, dnaLen_total / 10002)) 
			      + " " + blast_fmt + " -out " + blastOut + " > /dev/null 2> " + tmp + "/blast-err", tmp + "/blast-err");
 		#else
 			findProg ("makeblastdb");
 			exec (fullProg ("makeblastdb") + "-in " + dna_flat + "  -dbtype nucl  -out " + tmp + "/db  -logfile " + tmp + "/db.log  > /dev/null", tmp + "db.log");
 			findProg ("tblastn");
  		const string blast_fmt ("-outfmt '6 sseqid qseqid sstart send slen qstart qend qlen sseq qseq'");
			exec (fullProg ("tblastn") + " -query " + execDir + "stx.prot  -db " + tmp + "/db  "
			      + "-comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000  -word_size 5  -db_gencode " + to_string (gencode) 
			    //+ "  -task tblastn-fast  -threshold 100  -window_size 15"  // from amrfinder.cpp: Reduces time by 9% 
			    //+ "   -num_threads 10  -mt_mode 1"  // Reduces time by 30%
			      + " " + blast_fmt + " -out " + blastOut + " > /dev/null 2> " + tmp + "/blast-err", tmp + "/blast-err");
		#endif
		}


    // stxClass2identity[]
    stxClass2identity ["1a"] = 0.983;
    stxClass2identity ["1c"] = 0.983;
    stxClass2identity ["1d"] = 0.983;
    stxClass2identity ["1e"] = 0.983;
    stxClass2identity ["2"]  = 0.98;  
    stxClass2identity ["2b"] = 0.98;  
    stxClass2identity ["2e"] = 0.98;
    stxClass2identity ["2f"] = 0.98;
    stxClass2identity ["2g"] = 0.98;
    stxClass2identity ["2h"] = 0.98;
    stxClass2identity ["2i"] = 0.98;
    stxClass2identity ["2j"] = 0.98;
    stxClass2identity ["2k"] = 0.985;
    stxClass2identity ["2l"] = 0.985;
    stxClass2identity ["2m"] = 0.98;
    stxClass2identity ["2n"] = 0.98;
    stxClass2identity ["2o"] = 0.98;
    
    
    Cout out (output);
    TsvOut td (& *out, 2, false);
    TsvOut logTd (logPtr, 2, false);

    
    if (! input_name. empty ())
      td << "name";
    if (amrfinder)
    {
      td << /* 1*/ "Protein identifier"  
         << /* 2*/ "Contig id"
         << /* 3*/ "Start"
         << /* 4*/ "Stop"
         << /* 5*/ "Strand"
         << /* 6*/ "Element symbol"  // PD-4924
         << /* 7*/ "Element name"  // PD-4910
         << /* 8*/ "Scope"
         << /* 9*/ "Element type"
         << /*10*/ "Element subtype"
         << /*11*/ "Class"
         << /*12*/ "Subclass"
         << /*13*/ "Method"
         << /*14*/ "Target length"
         << /*15*/ "Reference sequence length"
         << /*16*/ "% Coverage of reference sequence"
         << /*17*/ "% Identity to reference sequence"
         << /*18*/ "Alignment length"
         << /*19*/ "Accession of closest sequence"
         << /*20*/ "Name of closest sequence"
         << /*21*/ "HMM id"
         << /*22*/ "HMM description"
         ;
      if (print_node)
        td << "Hierarchy node";  
    }
    else
      td << "target_contig"
         << "stx_type"
         << "operon" 
         << "identity"
         << "target_start"
         << "target_stop"
         << "target_strand"
         << "A_reference"
         << "A_identity"
         << "A_coverage"
         << "B_reference"
         << "B_identity"
         << "B_coverage"
         ; 
    td. newLn ();

	  VectorOwn<BlastAlignment> blastAls;   
    {
      LineInput f (blastOut);
  	  while (f. nextLine ())
  	  {
  	    const Unverbose unv;
  	    LOG (f. line);
  	    auto al = new BlastAlignment (f. line);
  	    al->qc ();  
 	      blastAls << al;
  	  }
  	}
  	
  	LOG ("# All stx blasts: " + to_string (blastAls. size ()));
    LOG ("Finding frame shifts:");
	  {
      // Multiple frame shifts are possible
      blastAls. sort (BlastAlignment::frameshiftLess); 
      const BlastAlignment* prev = nullptr;
      for (const BlastAlignment* al : blastAls)
      {        
        ASSERT (al);
        if (   prev
            && al->targetName   == prev->targetName
            && al->targetStrand == prev->targetStrand
            && al->refAccession == prev->refAccession
            && al->targetStart  >  prev->targetStart
            && (int) al->targetStart - (int) prev->targetEnd < 10  // PAR
            && al->getFrame () != prev->getFrame ()
           )
        {
          var_cast (al) -> merge (*prev);
          al->qc ();
          var_cast (prev) -> reported = true;
        }
        al->saveTsvOut (logTd, true);
        prev = al;
      }
    }
    
    LOG ("All blasts:");
	  VectorPtr<BlastAlignment> goodBlastAls;   
	  {
      blastAls. sort (BlastAlignment::sameTypeLess);
	    size_t start = 0;
      FFOR (size_t, i, blastAls. size ())
      {
        const BlastAlignment* al = blastAls [i];
        ASSERT (al);
        while (   start < i 
               && ! (   blastAls [start] -> targetName   == al->targetName
                     && blastAls [start] -> targetStrand == al->targetStrand
                     && blastAls [start] -> stxClass     == al->stxClass
                     && blastAls [start] -> subunit      == al->subunit
                     && blastAls [start] -> targetEnd    >  al->targetStart
                    )
              )
          start++;
        if (al->reported)
          break; 
        al->saveTsvOut (logTd, true);
        bool suppress = false;
        FOR_START (size_t, j, start, i)
        {
          const BlastAlignment* prev = blastAls [j];
          ASSERT (prev);
          ASSERT (! prev->reported);
          ASSERT (al->targetName   == prev->targetName);
          ASSERT (al->targetStrand == prev->targetStrand);
          ASSERT (al->stxClass     == prev->stxClass);
          ASSERT (al->subunit      == prev->subunit);
          if (   al->insideEq (*prev)
              && al->getDiff () >= prev->getDiff ()
            //&& al->getIdentity () <= prev-> getIdentity ()
             )
          {
            suppress = true;
            break;
          }
        }
        if (! suppress)
          goodBlastAls << al;
      }
    }
    
    Vector<Operon> operons;

    LOG ("\nSame type operons:");
    goodBlasts2operons (goodBlastAls, operons, true, true, logTd);
    
    goodBlastAls. sort (BlastAlignment::less);

    LOG ("\nStrong operons:");
    goodBlasts2operons (goodBlastAls, operons, false, true, logTd);

    LOG ("\nWeak operons:");
    goodBlasts2operons (goodBlastAls, operons, false, false, logTd);
   	  
  	LOG ("\ngoodOperons");
    Vector<Operon> goodOperons;
    {    
      operons. sort ();
      for (const Operon& op : operons)
      {
     	  op. saveTsvOut (logTd, true); 
       	op. qc ();     
        bool found = false;
        for (const Operon& goodOp : goodOperons)
          if (   op. al1->targetName == goodOp. al1->targetName
              && op. insideEq (goodOp)
              && goodOp. getIdentity () >= op. getIdentity ()
             )
          {
            found = true;
            break;
          }
        if (! found)
          goodOperons << op;          
      }      
    }

    // De-redundify single-subunit operons
	  {
	    size_t start = 0;
      FFOR (size_t, i, goodBlastAls. size ())
      {
        const BlastAlignment* al = goodBlastAls [i];
        ASSERT (al);
        if (al->reported)
          continue; 
        while (   start < i 
               && ! (   goodBlastAls [start] -> targetName   == al->targetName
                     && goodBlastAls [start] -> targetStrand == al->targetStrand
                     && goodBlastAls [start] -> subunit      == al->subunit
                     && goodBlastAls [start] -> targetEnd    >  al->targetStart
                    )
              )
          start++;
        al->saveTsvOut (logTd, true);
        FOR_START (size_t, j, start, i)
        {
          const BlastAlignment* prev = goodBlastAls [j];
          ASSERT (prev);
          ASSERT (al->targetName   == prev->targetName);
          ASSERT (al->targetStrand == prev->targetStrand);
          ASSERT (al->subunit      == prev->subunit);
          if (   al->insideEq (*prev)
              && al->getDiff () >= prev->getDiff ()
            //&& al->getIdentity () <= prev-> getIdentity ()
             )
          {
            var_cast (al) -> reported = true;
            break;
          }
        }
      }
    }

  	LOG ("\ngoodBlastAls -> goodOperons (single-subunit)");
    goodBlastAls. sort (BlastAlignment::reportLess); 
    FFOR (size_t, i, goodBlastAls. size ())
    {
      const BlastAlignment* al1 = goodBlastAls [i];
      ASSERT (al1);
      if (al1->reported)
        continue;
      Operon op (*al1);
      goodOperons << std::move (op);
      FFOR_START (size_t, j, i + 1, goodBlastAls. size ())
      {
        const BlastAlignment* al2 = goodBlastAls [j];
        ASSERT (al2);          
        if (! (   al2->targetName   == al1->targetName
               && al2->targetStrand == al1->targetStrand
              )
           )
          break;
        if (   ! al2->reported
            && al2->insideEq (*al1)
            && (   al2->stxType [0] == al1->stxType [0] 
                || al2->getDiff () >= al1->getDiff ()
              //|| al1->getIdentity () >= al2->getIdentity ()
               )
           )
          var_cast (al2) -> reported = true;
      }
    }

    // Report
    goodOperons. sort (Operon::reportLess);     
  	for (const Operon& op : goodOperons)
   	  op. saveTsvOut (td, false);
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



