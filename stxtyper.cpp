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
*  1.0.44 05/14/2025 PD-5329  New struct Hsp
*         04/25/2025 PD-5301  C++ refactoring: struct Hsp; trailing stop codons are not counted as nident
*  1.0.43 03/25/2025          -help or -version with other parameters is an error
*  1.0.42 02/20/2025 PD-5246  bioinformatically functional operons are preferred over weak operons disregarding the percent of identity
*  1.0.41 02/18/2025          Single subunit operon: PARTIAL_CONTIG_END/PARTIAL bug
*  1.0.40 02/04/2025 PD-5231  PARTIAL_CONTIG_END < EXTENDED
*  1.0.39 01/31/2025 PD-5231  suppress single-subunit operons overlappng with two-subunit operons
*  1.0.38 01/30/2025 PD-5231  select matches to reference proteins as if the matches had been at nucleotide level
*                             BlastAlignment::{positives --> sx}; Xs are not counted as identities
*  1.0.37 01/29/2025 PD-5198  operon type priorities: ... < PARTIAL_CONTIG_END < PARTIAL < AMBIGUOUS [complete operon] < COMPLETE_NOVEL < COMPLETE
*  1.0.36 01/14/2025 PD-5215  re-enable reporting "Name of closest sequence" for two-subunit operons
*  1.0.35 01/14/2025 PD-5215  "Closest reference accession" field has two accessions separated by "," for two-subunit operons
*  1.0.34 12/23/2024 PD-5209  make the AMRFinderPlus columns "Closest reference accession" and "Closest reference name" to be NA (per PD-4910)
*  1.0.33 12/23/2024 PD-5205  replace ", " by "," in the AMRFinderPlus column "Closest reference accession"
*  1.0.32 12/20/2024 PD-5201  change empty fields to NA
*  1.0.31 12/17/2024 PD-5181  COMPLETE and COMPLETE_NOVEL is preferred over the other operon types
*                             operons with higher identity and higher coverage are preferred
*  1.0.30 12/14/2024          bug in --debug
*                    PD-5191  strong operons should not be partial
*                             --threads: reduces time by 30% for large input DNA
*                             bug in frame shifts detection (for multiple frame shifts in the same protein)
*  1.0.29 12/13/2024 PD-5192  only operons with >= 80% identity are reported   // Only BLAST HSPs with identity >= 80% are considered
*  1.0.28 12/03/2024          tblastn -gapextend 2
*         10/30/2024          colorizeDir()
*  1.0.27 10/23/2024 PD-5155  "Hierarchy node" with mixed types is <stx1>::<stx2>
*  1.0.26 10/22/2024 PD-5085  change column "Element length" to "Target length"
*  1.0.25 08/16/2024 PD-5085  AMRFinderPlus column names to match MicroBIGG-E
*  1.0.24 08/05/2024 PD-5076  "na" -> "NA"
*  1.0.23 07/29/2024 PD-5064  AMBIGUOUS operon type
*  1.0.22 07/25/2024          first codon L|I|V -> M
*  1.0.21 07/15/2024 PD-5038  --nucleotide_output 
*  1.0.20 05/21/2024 PD-5002  {A|B}_reference_subtype
*  1.0.19 03/26/2024          BlastAlignment::targetAlign is removed
*  1.0.18 03/19/2024 PD-4910  element symbol is <stx type>_operon, Element name contains operon quality attribute"

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
*          02/13/2024 PD-4892  subunits A and B are not preferred to be of the same stx class
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
#include "seq.hpp"
using namespace Seq_sp;
#include "amrfinder_columns.hpp"

#include "common.inc"


#undef PROT_MATCH  // 0 <=> nucleotide level matching to protein reference sequences  // PD-5231



namespace 
{


string input_name;
bool amrfinder = false;
bool print_node = false;
map<string,double> stxClass2identity;

// PAR
// bp
constexpr size_t intergenic_max {36};  // = max. intergenic region in the reference set + 2
constexpr size_t slack = 30;  
//
constexpr double identity_min = 0.8;  // PD-5192

const string stxS ("stx");



string stxType_reported_operon2elementName (const string &stxType_reported,
                                            const string &operon)
{
  string elementName (stxType_reported  + " operon");
       if (operon == frameshift_Name)
    elementName += " with frameshift";
  else if (operon == internalStop_Name)
    elementName += " with internal stop";
  else if (contains (operon, "PARTIAL"))
    elementName = "Partial " + elementName;
  else if (operon == "EXTENDED")
    elementName = "Extended " + elementName;
  else if (contains (operon, "NOVEL"))
    elementName = "Novel " + elementName;
    
  return elementName;
}



struct BlastAlignment final : Hsp
// qseq: whole sequence ends with '*'
{
  // Function of qseqid
  string stxType;
  string stxClass;
    // Function of stxType
  string stxSuperClass;
    // Function of stxClass
  char subunit {'\0'};
    // 'A' or 'B'
  string subClass;  // = as in AMRFinderPlus report
  
  bool reported {false};


  BlastAlignment (const string &line)
    : Hsp (line, true/*qProt*/, false/*sProt*/, true/*aProt*/, true/*qStopCodon*/, true/*bacterialStartCodon*/)
    {
      {
  	    string famId;
        try
        {	
  		    subClass = rfindSplit (qseqid, '|');  
  		    famId    = rfindSplit (qseqid, '|');  
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
      	      
      stxClass = stxType;
      if (   stxType == "2a"
          || stxType == "2c"
          || stxType == "2d"
         )
        stxClass = "2";
      
      stxSuperClass = stxClass. substr (0, 1);
    }
  void qc () const final
    {
      if (! qc_on)
        return;
      Hsp::qc ();
      QC_ASSERT (contains (stxClass2identity, stxClass));
      QC_ASSERT (isLeft (stxType, stxClass));
      QC_ASSERT (subunit == 'A' || subunit == 'B');
      QC_ASSERT (subClass. size () > stxS. size ());
      QC_ASSERT (stxType. size () == 2);
    }
  void saveTsvOut (TsvOut& td,
                   bool verboseP) const 
    { if (! td. live ())
        return;
      const string stxType_reported (verboseP ? getGenesymbol () : (stxS + stxType. substr (0, 1)));
      const string quality (findDisruption (Disruption::eFrameshift)
                              ? frameshift_Name
                              : sInternalStop 
                                ? internalStop_Name
                                : sTruncated () || otherTruncated ()
                                  ? partialContigEnd_Name
                                  : debugP && verboseP && qComplete ()
                                    ? "COMPLETE_SUBUNIT"
                                    : c_extended ()
                                      ? "EXTENDED"
                                      : partial_Name
                           );
      const char strand (strand2char (sInt. strand));
      const double refCoverage = qRelCoverage () * 100.0;
      const double refIdentity = relIdentity ()  * 100.0; 
      // td     
      if (! input_name. empty ())
  	    td << input_name;
      if (amrfinder)
      {
        const string subunitS (1, subunit);
        string subclass (stxType_reported /*stxS + stxType*/);
        strUpper (subclass);
        td << na               // 1 "Protein identifier"  
           << sseqid           // 2 "Contig id"
           << sInt. start + 1  // 3 "Start"
           << sInt. stop       // 4 "Stop"
           << strand           // 5 "Strand"
           << stxType_reported + "_operon" // 6 "Element symbol"
           << stxType_reported_operon2elementName (stxType_reported, quality)    // 7 "Element name"
           << "plus"           // 8 "Scope"
           << "VIRULENCE"      // 9 "Element type"
           << "STX_TYPE"       //10 "Element subtype"
           << subclass. substr (0, 4)   //11 "Class"
           << subclass         //12 "Subclass"
           << quality          //13 "Method"  
           << sInt. len ()     /*targetAlign*/      //14 "Target length" 
           << na /*qlen*/      //15 "Reference sequence length"
           << na /*refCoverage*/   //16 "% Coverage of reference sequence"
           << refIdentity      //17 "% Identity to reference sequence"
           << length           //18 "Alignment length"
           << qseqid           //19 "Accession of closest sequence"
           << "Shiga toxin " + stxS + stxType + " subunit " + subunitS //20 "Name of closest sequence"
           << na               //21 "HMM id"
           << na               //22 "HMM description"
           ;
        if (print_node)
          td << getGenesymbol ();
      }
      else
      {
        td << sseqid
           << stxType_reported
           << quality
           << na
           << sInt. start + 1
           << sInt. stop
           << strand;
        if (subunit == 'B')
          td << na
             << na
             << na
             << na;
        td << qseqid
           << subClass
           << refIdentity
           << refCoverage;
        if (subunit == 'A')
          td << na
             << na
             << na
             << na;
      }
      td. newLn ();
    }
    

  bool c_extended () const 
    { return    ! qInt. start          // N-terminus is complete
             && c_complete == efalse;  // Only "*" (stop codon) is missing
    }
  string getGenesymbol () const
    { return stxS + subunit + stxType; }
  bool otherTruncated () const
    { constexpr size_t missed_max = intergenic_max + 3 * 20 /*min. domain length*/;  // PAR
      return    ((sInt. strand == 1) == (subunit == 'B') && sInt. start       <= missed_max)
             || ((sInt. strand == 1) == (subunit == 'A') && slen - sInt. stop <= missed_max);
    }
  static bool sameClassLess (const BlastAlignment* a,
                             const BlastAlignment* b) 
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, sseqid);
      LESS_PART (*a, *b, sInt. strand);
      LESS_PART (*a, *b, stxClass);
      LESS_PART (*a, *b, subunit);
      LESS_PART (*a, *b, sInt. start);
      return false;
    }
  static bool less (const BlastAlignment* a,
                    const BlastAlignment* b) 
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, sseqid);
      LESS_PART (*a, *b, sInt. strand);
      LESS_PART (*a, *b, subunit);
      LESS_PART (*a, *b, sInt. start);
      LESS_PART (*b, *a, relIdentity ());  
      LESS_PART (*b, *a, qRelCoverage ());  
      // Tie resolution
      LESS_PART (*a, *b, qseqid);  
      return false;
    }
  static bool reportLess (const BlastAlignment* a,
                          const BlastAlignment* b)  
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, sseqid);
      LESS_PART (*a, *b, sInt. strand);
    //LESS_PART (*b, *a, getAbsCoverage ()); 
      LESS_PART (*b, *a, relIdentity ());  
      LESS_PART (*b, *a, qRelCoverage ());  
      LESS_PART (*a, *b, sInt. start);
      // Tie resolution
      LESS_PART (*a, *b, qseqid);  
      return false;
    }
};



struct Operon
{
  const BlastAlignment* al1 {nullptr};
    // !nullptr
  const BlastAlignment* al2 {nullptr};
  // al1->send <= al2->sstart
  

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
      QC_ASSERT (al1->sseqid       == al2->sseqid);
      QC_ASSERT (al1->sInt. strand == al2->sInt. strand);
      QC_ASSERT (al1->sInt. stop   <= al2->sInt. start);
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
                           || relIdentity () < stxClass2identity [al1->stxClass]
                           || stxType. size () <= 1;
        const string quality =    al1->findDisruption (Disruption::eFrameshift)  
                               || al2->findDisruption (Disruption::eFrameshift)
                                 ? frameshift_Name
                                 :    al1->sInternalStop
                                   || al2->sInternalStop
                                     ? internalStop_Name
                                     :    al1->sTruncated () 
                                       || al2->sTruncated ()
                                         ? partialContigEnd_Name
                                         :    al1->c_extended ()
                                           || al2->c_extended ()
                                             ? "EXTENDED" 
                                             :    ! al1->qComplete ()
                                               || ! al2->qComplete ()
                                               ? partial_Name
                                               // complete operon types
                                               : novel
                                                 ? xs ()
                                                   ? "AMBIGUOUS"
                                                   : standard + "_NOVEL"
                                                 : standard;
        if (! verboseP)
        {
          ASSERT (stxType. size () <= 2);
          if (quality == standard)
            { ASSERT (stxType. size () == 2); }
          else
            if (stxType. size () == 2)
              stxType. erase (1);  
        }
        const string sseqid (al1->sseqid);
        const size_t start = al1->sInt. start + 1;
        const size_t stop  = al2->sInt. stop;
  	    const char strand = strand2char (al1->sInt. strand);
  	    const string stxType_reported (stxS + stxType);
  	    const double refIdentity = relIdentity () * 100.0;
  	    // td
        if (! input_name. empty ())
  	      td << input_name;
        if (amrfinder)
        {
          const string genesymbol (al1->stxType == al2->stxType ? stxS + al1->stxType : stxType_reported);
          string subclass (stxType_reported /*genesymbol*/);
          strUpper (subclass);
          const size_t targetAlign = al2->sInt. stop - al1->sInt. start;
        //const size_t qlen = al1->qlen + al2->qlen;
        //const double refCoverage = double (al1->getAbsCoverage () + al2->getAbsCoverage ()) / double (qlen) * 100.0;
          const size_t alignmentLen = al1->length + al2->length;
          const string refAccessions (al1->qseqid + "," + al2->qseqid);  // No space: PD-5205
          const string fam (al1->getGenesymbol () + fusion_infix + al2->getGenesymbol ());
          td << na                // 1 "Protein identifier"  
             << sseqid            // 2 "Contig id"
             << start             // 3 "Start"
             << stop              // 4 "Stop"
             << strand            // 5 "Strand"
             << stxType_reported + "_operon"  // 6 "Element symbol"
             << stxType_reported_operon2elementName (stxType_reported, quality)    // 7 "Element name"
             << "plus"            // 8 "Scope"
             << "VIRULENCE"       // 9 "Element type"
             << "STX_TYPE"        //10 "Element subtype"
             << subclass. substr (0, 4)   //11 "Class"
             << subclass          //12 "Subclass"
             << quality           //13 "Method"  
             << targetAlign       //14 "Target length" 
             << na /*qlen*/       //15 "Reference sequence length"
             << na /*refCoverage*/ //16 "% Coverage of reference sequence"
             << refIdentity       //17 "% Identity to reference sequence"
             << alignmentLen      //18 "Alignment length"
             << refAccessions     //19 "Accession of closest sequence"        // PD-5209
             << "Shiga toxin " + genesymbol //20 "Name of closest sequence"   // PD-5209
             << na                //21 "HMM id"
             << na                //22 "HMM description"
             ;
          if (print_node)
            td << fam;
        }
        else
    	    td << sseqid
             << stxType_reported
             << quality
             << refIdentity
             << start
             << stop
             << strand
             << getA () -> qseqid
             << getA () -> subClass
             << getA () -> relIdentity () * 100.0
             << getA () -> qRelCoverage () * 100.0
             << getB () -> qseqid
             << getB () -> subClass
             << getB () -> relIdentity () * 100.0
             << getB () -> qRelCoverage () * 100.0
             ;
        td. newLn ();
      }
      else
        al1->saveTsvOut (td, verboseP);
    }


  const BlastAlignment* getA () const
    { return al1->sInt. strand == 1 ? al1 : al2; }
  const BlastAlignment* getB () const
    { return al1->sInt. strand == 1 ? al2 : al1; }
private:
  bool hasAl2 () const
    { return al2; }
  string getRefAccession2 () const
    { if (al2)
        return al2->qseqid;
      return noString;
    }
  string getStxType (bool verboseP) const
    { if (! al2)
        return al1->stxType;
      if (al1->stxClass != al2->stxClass)
      {
      //return al1->stxClass + fusion_infix + al2->stxClass;  // order alphabetically
        if (al1->stxSuperClass == al2->stxSuperClass)
          return al1->stxSuperClass;  
        return noString;
      }
      if (al1->stxClass != "2")
        return al1->stxType; 
      const string a (getA () -> qMap (319 + 1));
      const string b (getB () -> qMap ( 89 + 1));
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
  size_t xs () const
    { ASSERT (al2);
      return al1->sx + al2->sx;
    }
  size_t getTargetEnd () const
    { return al2 ? al2->sInt. stop : al1->sInt. stop; }
  size_t getNident () const
    { return al2 ? al1->nident + al2->nident : al1->nident; }
  size_t getLength () const
    { return al2 ? al1->length + al2->length : al1->length; }
public:
  double relIdentity () const
    { return (double) getNident () / (double) getLength (); }
  double getRelCoverage () const 
    { ASSERT (al2);
      return double (al1->qAbsCoverage () + al2->qAbsCoverage ()) / double (al1->qlen + al2->qlen); 
    }
  bool perfect () const
    { ASSERT (al2);
      return    al1->perfect ()
             && al2->perfect ();
    }
  bool insideEq (const Operon &other,
                 size_t slack_arg) const
    { return    al1->sInt. strand            == other. al1->sInt. strand
    	       && al1->sInt. start + slack_arg >= other. al1->sInt. start 
             && getTargetEnd ()              <= other. getTargetEnd () + slack_arg;
    }
  bool betterEq (const Operon &other) const
    {  if (al1->sseqid != other. al1->sseqid)
         return false;
       if (! other. insideEq (*this, 3 * slack))  // PAR
         return false;
     #ifdef PROT_MATCH
       if (perfect () > other. perfect ())
         return true;
       if (perfect () < other. perfect ())
         return false;
     #endif
       return relIdentity () >= other. relIdentity ();
    }
  bool operator< (const Operon &other) const
    // Ordering by quality
    { ASSERT (al2);
      LESS_PART (*this, other, al1->sseqid);
    #ifdef PROT_MATCH
      LESS_PART (other, *this, perfect ());
    #endif
      LESS_PART (other, *this, relIdentity ());
      LESS_PART (other, *this, getRelCoverage ());
      // Tie resolution
      LESS_PART (*this, other, al1->qseqid);
      LESS_PART (*this, other, getRefAccession2 ());
      return false;
    }
  static bool reportLess (const Operon &a,
                          const Operon &b)
    { LESS_PART (a, b, al1->sseqid);
      LESS_PART (a, b, al1->sInt. start);
      LESS_PART (a, b, al1->sInt. stop);
      LESS_PART (b, a, al1->sInt. strand);
      LESS_PART (a, b, hasAl2 ());
      // Tie resolution
      LESS_PART (a, b, al1->qseqid);
      LESS_PART (a, b, getRefAccession2 ());
      return false;
    }
};



VectorPtr<BlastAlignment> processDisruptions (VectorPtr<Hsp> &hsps)
{
  VectorPtr<BlastAlignment> newAls;
  
  if (hsps. empty ())
    return newAls;
  
  VectorPtr<Hsp> firstOrigHsps;   // Subset of hsps
  const Vector<Hsp> mergedHsps (Hsp::merge (hsps, firstOrigHsps, nullptr/*sm*/, 20/*intronScore*/, true/*bacteria*/));  // PAR
//ASSERT (mergedHsps. size () <= hsps. size ());
  ASSERT (mergedHsps. size () == firstOrigHsps. size ());
  ASSERT (hsps. containsAll (firstOrigHsps));      
        
  VectorOwn<BlastAlignment> newAls_;
  bool hasDisruptions = false;
  FFOR (size_t, i, mergedHsps. size ())
  {
    const Hsp* origHsp = firstOrigHsps [i];  // Matches mergedHsps[i]
    ASSERT (origHsp);
    ASSERT (! origHsp->merged);
    auto al = new BlastAlignment (* static_cast <const BlastAlignment*> (origHsp));
    * static_cast <Hsp*> (al) = std::move (mergedHsps [i]);
    al->qc ();
    newAls_ << al;
    if (! al->disrs. empty ())
      hasDisruptions = true;
  }
  
  if (hasDisruptions)
  {
    for (const Hsp* hsp : hsps)
      const_static_cast <BlastAlignment*> (hsp) -> reported = true;
    newAls = newAls_;
    newAls_. clear ();
  }
      
  return newAls;
}



void paretoBest (const VectorPtr<BlastAlignment> &blastAls, 
                 size_t start, 
                 size_t end)
{
  ASSERT (start <= end);
  ASSERT (end <= blastAls. size ());

  FOR_START (size_t, i, start, end)
    if (! blastAls [i] -> reported)
      FOR_START (size_t, j, start, end)
        if (   i != j
            &&   blastAls [j] -> qBetterEq (* blastAls [i])
            && ! blastAls [i] -> qBetterEq (* blastAls [j])
           )
        {
          var_cast (blastAls [i]) -> reported = true;
          break;
        }
}



void goodBlasts2operons (const VectorPtr<BlastAlignment> &goodBlastAls, 
                         Vector<Operon> &operons, 
                         bool sameClass,
                         ebool strong,
                         TsvOut &logTd)
{
  IMPLY (sameClass, strong == etrue);
  
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
           && ! (   goodBlastAls [start] -> sseqid       == alB->sseqid
                 && goodBlastAls [start] -> sInt. strand == alB->sInt. strand
                 && (   ! sameClass 
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
      ASSERT (alA->sseqid  == alB->sseqid);
      ASSERT (alA->sInt. strand == alB->sInt. strand);
      IMPLY (sameClass, alA->stxClass == alB->stxClass);
      ASSERT (alA->subunit <= alB->subunit);
      if (alA->subunit == alB->subunit)
        break;
      ASSERT (alA->subunit == 'A');
      const BlastAlignment* al1 = alA;
      const BlastAlignment* al2 = alB;
      if (al1->sInt. strand == -1)
        swap (al1, al2);
      if (   al1->sInt. stop <= al2->sInt. start  
          && al2->sInt. start - al1->sInt. stop <= intergenic_max * (strong == efalse ? 2 : 1)  // PAR  // PD-4897
         )
      {
        Operon op (*al1, *al2);
        LOG ("Operon: " + to_string (op. relIdentity ()) + " " + to_string (stxClass2identity [op. al1->stxClass]) + ":");
        op. saveTsvOut (logTd, true);  
      #if 1
        bool good = false;
        switch (strong)
        {
          case etrue: good =    op. relIdentity () >= stxClass2identity [op. al1->stxClass]
                             && op. relIdentity () >= stxClass2identity [op. al2->stxClass];
            break;
          case enull: good = op. perfect ();
            break;
          default: good = true;
        }
        if (good)
      #else
        if (   strong != etrue
            || (   op. relIdentity () >= stxClass2identity [op. al1->stxClass]
                && op. relIdentity () >= stxClass2identity [op. al2->stxClass]
              #ifdef PROT_MATCH
                && op. perfect ()
              #endif
               )
           )
      #endif
        {
          operons << std::move (op);
          LOG ("Added");
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
        if (   al->sseqid              == op. al1->sseqid
            && al->sInt. start + slack >= op. al1->sInt. start 
            && al->sInt. stop          <= op. al2->sInt. stop + slack
            && al->sInt. strand        == op. al1->sInt. strand
           ) 
        {
          var_cast (al) -> reported = true;
          break;
        }
      }
  }
  
  if (qc_on)
    for (const Operon& op : operons)
    {
      ASSERT (op. al1);
      ASSERT (op. al2);
      ASSERT (op. al1->reported);
      ASSERT (op. al2->reported);
    }    
}



// ThisApplication

struct ThisApplication final : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Determine stx type(s) of a genome, print .tsv-file", true, true/*threadsUsed*/, true, true)
    {
    	addKey ("nucleotide", "Input nucleotide FASTA file (can be gzipped)", "", 'n', "NUC_FASTA");
    //addKey ("translation_table", "NCBI genetic code for translated BLAST", "11", 't', "TRANSLATION_TABLE");
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report, for example it can be an assembly name", "", '\0', "NAME");
      addKey ("output", "Write output to OUTPUT_FILE instead of STDOUT", "", 'o', "OUTPUT_FILE");
    	addKey ("blast_bin", "Directory for BLAST. Deafult: $BLAST_BIN", "", '\0', "BLAST_DIR");
    	addFlag ("amrfinder", "Print output in the nucleotide AMRFinderPlus format");
    	addFlag ("print_node", "Print AMRFinderPlus hierarchy node");
      addKey ("nucleotide_output", "Output nucleotide FASTA file of reported nucleotide sequences", "", '\0', "NUC_FASTA_OUT");
        // Flag: only for COMPLETE_NOVEL ??

      version = SVN_REV;
    }



  void shellBody () const final
  {
    const string fName      = shellQuote (getArg ("nucleotide"));
    const uint   gencode    =           /*arg2uint ("translation_table")*/ 11; 
                 input_name =             getArg ("name");
    const string output     =             getArg ("output");
          string blast_bin  =             getArg ("blast_bin");
                 amrfinder  =             getFlag ("amrfinder");
                 print_node =             getFlag ("print_node");
    const string dna_out    = shellQuote (getArg ("nucleotide_output"));
    
    if (contains (input_name, '\t'))
      throw runtime_error ("NAME cannot contain a tab character");
    if (print_node && ! amrfinder)
      throw runtime_error ("--print_node requires --amrfinder");


    const bool screen = ! isRedirected (cerr);

    stderr << "Software directory: " << colorizeDir (execDir, screen) << '\n';
    stderr << "Version: " << version << '\n'; 
    
		const string logFName (tmp + "/log"); 
    const string qcS (qc_on ? " -qc" : noString);


    // blast_bin
    if (blast_bin. empty ())
    	if (const char* s = getenv ("BLAST_BIN"))
    		blast_bin = string (s);
    if (! blast_bin. empty ())
    {
	    addDirSlash (blast_bin);
	    prog2dir ["tblastn"]     = blast_bin;  
	    prog2dir ["makeblastdb"] = blast_bin;  
	  }

    const string dna_flat = uncompress (fName,  "dna_flat");
    
    {
      prog2dir ["fasta_check"] = execDir;
      exec (fullProg ("fasta_check") + dna_flat + "  -hyphen  -ambig  " + qcS + "  -log " + logFName + " > " + tmp + "/nseq", logFName); 
    	const StringVector vec (tmp + "/nseq", (size_t) 10, true); 
    	if (vec. size () != 3)
        throw runtime_error ("fasta_check failed: " + vec. toString ("\n"));
    }

	//stderr. section ("Running blast");
	  const string blastOut (tmp + "/blast");
		{
			const Chronometer_OnePass_cerr cop ("blast");
 			// Database: created by ~brovervv/code/database/stx.prot.sh
 			findProg ("makeblastdb");
 			exec (fullProg ("makeblastdb") + "-in " + dna_flat + "  -dbtype nucl  -out " + tmp + "/db  -logfile " + tmp + "/db.log  > /dev/null", tmp + "db.log");
 			findProg ("tblastn");
			exec (fullProg ("tblastn") + " -query " + execDir + "stx.prot  -db " + tmp + "/db"
			      + Hsp::blastp_par_fast + "  -gapextend 2  -db_gencode " + to_string (gencode) 
			      + "  -mt_mode 1  -num_threads " + to_string (threads_max)  // Reduces time by 30% for large DNA
			      + " " + Hsp::format_par (true) + " -out " + blastOut + " > /dev/null 2> " + tmp + "/blast-err", tmp + "/blast-err");
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
    

    const string tmpOut (tmp + "/out");
    OFStream fOut (tmpOut);
    TsvOut td (& fOut, 2, false);
    TsvOut logTd (logPtr, 2, false);
    logTd. usePound = false;

    
    if (! input_name. empty ())
      td << "name";
    if (amrfinder)
    {
      td << /* 1*/ prot_colName 
         << /* 2*/ contig_colName
         << /* 3*/ start_colName
         << /* 4*/ stop_colName
         << /* 5*/ strand_colName
         << /* 6*/ genesymbol_colName
         << /* 7*/ elemName_colName
         << /* 8*/ scope_colName
         << /* 9*/ type_colName
         << /*10*/ subtype_colName
         << /*11*/ class_colName
         << /*12*/ subclass_colName
         << /*13*/ method_colName
         << /*14*/ targetLen_colName
         << /*15*/ refLen_colName
         << /*16*/ refCov_colName
         << /*17*/ refIdent_colName
         << /*18*/ alignLen_colName
         << /*19*/ closestRefAccession_colName
         << /*20*/ closestRefName_colName
         << /*21*/ hmmAccession_colName
         << /*22*/ hmmDescr_colName
         ;
      if (print_node)
        td << hierarchyNode_colName;  
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
         << "A_reference_subtype"
         << "A_identity"
         << "A_coverage"
         << "B_reference"
         << "B_reference_subtype"
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
      blastAls. sort (Hsp::less); 
      const BlastAlignment* prev = nullptr;
      VectorPtr<Hsp> hsps;  // Subset of blastAls
      FFOR (size_t, i, blastAls. size ())  // Fixed blastAls.size()
      {
        const BlastAlignment* al = blastAls [i];
        ASSERT (al);
        if (   prev
            && ! (   al->sseqid       == prev->sseqid
                  && al->sInt. strand == prev->sInt. strand
                  && al->qseqid       == prev->qseqid
                 )
           )
        {
          blastAls << processDisruptions (hsps);
          hsps. clear ();
        } 
        hsps << al;  
        prev = al;     
      }
      blastAls << processDisruptions (hsps);
    }
    
    // Pareto-best
    LOG ("All blasts:");
	  {
      blastAls. sort (BlastAlignment::sameClassLess);
	    size_t start = 0;
      FFOR (size_t, i, blastAls. size ())
      {
        const BlastAlignment* al = blastAls [i];
        ASSERT (al);
        al->saveTsvOut (logTd, true);
        if (! (   blastAls [start] -> sseqid        == al->sseqid
               && blastAls [start] -> sInt. strand  == al->sInt. strand
               && blastAls [start] -> stxClass      == al->stxClass
               && blastAls [start] -> subunit       == al->subunit
              )
           )
        {
          paretoBest (blastAls, start, i);
          start = i;
        }
      }
      paretoBest (blastAls, start, blastAls. size ());
    }
    
	  VectorPtr<BlastAlignment> goodBlastAls;   
	  goodBlastAls. reserve (blastAls. size ());
	  for (const BlastAlignment* al : blastAls)
	    if (! al->reported)
	      goodBlastAls << al;
    
    Vector<Operon> operons;

    LOG ("\nSame type operons:");
    goodBlasts2operons (goodBlastAls, operons, true, etrue, logTd);
    
    goodBlastAls. sort (BlastAlignment::less);

    LOG ("\nStrong operons:");
    goodBlasts2operons (goodBlastAls, operons, false, etrue, logTd);

    LOG ("\nWeak complete operons:");
    goodBlasts2operons (goodBlastAls, operons, false, enull, logTd);  // PD-5246

    LOG ("\nWeak non-complete operons:");
    goodBlasts2operons (goodBlastAls, operons, false, efalse, logTd);
   	  
  	LOG ("\ngoodOperons");
    Vector<Operon> goodOperons;
    {    
      operons. sort ();
      for (const Operon& op : operons)
      {
     	  op. saveTsvOut (logTd, true); 
       	op. qc ();     
       	if (op. relIdentity () < identity_min)
       	  continue;
        bool betterFound = false;
        for (const Operon& goodOp : goodOperons)
          if (goodOp. betterEq (op))
          {
            betterFound = true;
            break;
          }
        if (! betterFound)
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
        al->saveTsvOut (logTd, true);
        if (! (   goodBlastAls [start] -> sseqid        == al->sseqid
               && goodBlastAls [start] -> sInt. strand  == al->sInt. strand
               && goodBlastAls [start] -> subunit       == al->subunit
              )
           )
        {
          paretoBest (goodBlastAls, start, i);
          start = i;
        }
      }
      paretoBest (goodBlastAls, start, goodBlastAls. size ());
    }

  	LOG ("\ngoodBlastAls -> goodOperons (single-subunit)");
    goodBlastAls. sort (BlastAlignment::reportLess); 
    FFOR (size_t, i, goodBlastAls. size ())
    {
      const BlastAlignment* al1 = goodBlastAls [i];
      ASSERT (al1);
      if (al1->reported)
        continue;
     	if (al1->relIdentity () < identity_min)
     	  continue;
      Operon op (*al1);
      bool good = true;
      for (const Operon& op_good : goodOperons)
        if (op_good. betterEq (op))  
        {
          good = false;
          break;
        }
      if (good)
        goodOperons << std::move (op);
    }


    // Report
    goodOperons. sort (Operon::reportLess);     
  	for (const Operon& op : goodOperons)
   	  op. saveTsvOut (td, false);

    // Output
    {
      TextTable tt (tmpOut);
      tt. qc ();
      {
        Cout out (output);
   		  tt. saveText (*out);
   		}
      if (! emptyArg (dna_out))
      {
        const StringVector columns {"target_contig", "target_start", "target_stop", "target_strand", "stx_type", "operon"};
        tt. filterColumns (columns);
        tt. saveHeader = false;
        tt. qc ();
        const string extract (tmp + "/extract");
        tt. saveFile (extract);
        prog2dir ["fasta_extract"] = execDir;
        exec (fullProg ("fasta_extract") + dna_flat + " " + extract + qcS + " -log " + logFName + " > " + dna_out, logFName);  
      }
    }
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



