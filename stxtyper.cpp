// stxtyper.cpp
// --> amr_report.cpp ??

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
//#include "../cpp/genetics/alignment.hpp"  ??

#include "common.inc"



// PAR!
//#define DATA_VER_MIN "2024-01-31.1"  // ??



namespace 
{


string input_name;
map<string,double> stxClass2identity;

constexpr size_t slack = 30;  // PAR
const string stxS ("stx");



struct BlastAlignment 
{
  // BLASTX alignment
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
  
  bool reported {false};


  BlastAlignment (const string &line,
                  bool &stx)
    {
      stx = true;

	    // Cf. amr_report.cpp
      {
        string sseqid;
        {
    	    istringstream iss (line);
    	    iss >> targetName >> sseqid >> length >> nident >> targetStart >> targetEnd >> targetLen >> refStart >> refEnd >> refLen >> targetSeq >> refSeq;
  	  // format:  qseqid       sseqid    length    nident         qstart         qend         qlen      sstart      send      slen         qseq    sseq
      // blastx:                            733       733          62285        63017        88215         105       837       837          
        }
  	    QC_ASSERT (! targetSeq. empty ());	
  	    
  	    // sseqid
  	    // 0|EED0303793.1|1|1|stxB2j|stxB2j||1|STX2J|STX2|Shiga_toxin_Stx2j_subunit_B
  	    string classS;
  	    string famId;
        try
        {	
  	//  /*product      =*/                    rfindSplit (sseqid, '|'); 
  		    classS       =                      rfindSplit (sseqid, '|'); 
  	//  /*subclass     =*/                    rfindSplit (sseqid, '|'); 
  	//  /*reportable   =(uchar)*/ str2<int>  (rfindSplit (sseqid, '|')); 
  	//  /*resistance   =*/                    rfindSplit (sseqid, '|'); 
  	//  /*gene         =*/                    rfindSplit (sseqid, '|');  
  		    famId        =                      rfindSplit (sseqid, '|');  
  	//  /*parts        = (size_t)*/str2<int> (rfindSplit (sseqid, '|'));
  	//  /*part         = (size_t)*/str2<int> (rfindSplit (sseqid, '|'));
  		    refAccession =                      rfindSplit (sseqid, '|');
  	//  /*gi           =*/                    str2<long> (sseqid);  // dummy
  		  }
  		  catch (const exception &e)
  		  {
  		  	throw runtime_error (string ("Bad AMRFinder database\n") + e. what () + "\n" + line);
  		  }
  		  if (! (   classS == "STX1" 
  		         || classS == "STX2" 
  		        )
  		     )
  		  {
  		    stx = false;
  		    return;
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

	    QC_ASSERT (refStart < refEnd);

	    QC_ASSERT (targetStart != targetEnd);
	    targetStrand = targetStart < targetEnd;  
	    if (! targetStrand)
	      swap (targetStart, targetEnd);
	    
	    QC_ASSERT (refStart >= 1);
	    QC_ASSERT (targetStart >= 1);
	    refStart--;
	    targetStart--;
	    
	    const size_t stopCodonPos = targetSeq. find ('*');
	    if (stopCodonPos != string::npos && stopCodonPos + 1 < targetSeq. size ())
	      stopCodon = true;
    }
  void qc () const
    {
      if (! qc_on)
        return;
      QC_ASSERT (length);
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
  void saveTsvOut (TsvOut& td) const 
    { if (! input_name. empty ())
  	    td << input_name;
      td << targetName
         << (stxS + stxType [0])
         << (frameshift 
              ? "FRAMESHIFT"
              : stopCodon 
                ? "INTERNAL_STOP"
                : truncated ()
                  ? "PARTIAL_CONTIG_END"
                  : "PARTIAL"
            )
         << noString
         << targetStart + 1
         << targetEnd
         << targetStrand;
      if (subunit == 'B')
        td << noString
           << noString
           << noString;
      td << refAccession
         << getIdentity () * 100.0
         << getCoverage () * 100.0;
      if (subunit == 'A')
        td << noString
           << noString
           << noString;
      td. newLn ();
    }
    

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
      if (prev. stopCodon)
        stopCodon = true;
      frameshift = true;
    }
  size_t getFrame () const
    { return (targetStart % 3) + 1; }
  double getIdentity () const 
    { return (double) nident / (double) (length /*refEnd - refStart*/); }
  double getCoverage () const 
    { return (double) (refEnd - refStart) / (double) refLen; }
  bool truncated () const
    { return    targetStart           <= 3 /*Locus::end_delta*/
             || targetLen - targetEnd <= 3 /*Locus::end_delta*/;
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
  static bool frameshiftLess (const BlastAlignment& a,
                              const BlastAlignment& b)
    { LESS_PART (a, b, targetName);
      LESS_PART (a, b, targetStrand);
      LESS_PART (a, b, refAccession);
      LESS_PART (a, b, targetStart);
      LESS_PART (a, b, targetEnd);
      return false;
    }
  bool operator< (const BlastAlignment &other) const
    { LESS_PART (*this, other, targetName);
      LESS_PART (*this, other, targetStrand);
      LESS_PART (*this, other, stxClass);
      LESS_PART (*this, other, subunit);
      LESS_PART (*this, other, targetStart);
      LESS_PART (other, *this, nident);
      return false;
    }
  // For VectorPtr
  static bool stxLess (const BlastAlignment* a,
                       const BlastAlignment* b)
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, targetName);
      LESS_PART (*a, *b, targetStrand);
      LESS_PART (*a, *b, subunit);
      LESS_PART (*a, *b, targetStart);
      LESS_PART (*b, *a, nident);
      return false;
    }
  static bool reportLess (const BlastAlignment* a,
                          const BlastAlignment* b)  
    { ASSERT (a);
      ASSERT (b);
      LESS_PART (*a, *b, reported);
      LESS_PART (*a, *b, targetName);
      LESS_PART (*a, *b, targetStrand);
      LESS_PART (*b, *a, nident);
      LESS_PART (*a, *b, targetStart);
      return false;
    }
};



struct Operon
{
  static constexpr size_t intergenic_max {36};  // PAR ??
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
  void saveTsvOut (TsvOut& td) const 
    { ASSERT (al1);
      if (al2)
      {
        string stxType (getStxType ());
        const bool novel =    al1->stxClass != al2->stxClass 
                           || getIdentity () < stxClass2identity [al1->stxClass]
                           || stxType. size () <= 1;
        const string operonType =    getA () -> frameshift
                                  || getB () -> frameshift
                                    ? "FRAMESHIFT"
                                    :    getA () -> stopCodon 
                                      || getB () -> stopCodon 
                                      ? "INTERNAL_STOP"
                                      : partial ()
                                        ? truncated ()
                                          ? "PARTIAL_CONTIG_END"
                                          : "PARTIAL"  // Use Blast Rules ??
                                        : novel
                                          ? "NOVEL"
                                          : "STANDARD";
        if (   operonType != "STANDARD"
            && stxType. size () >= 2
           )
          stxType. erase (1);  
        if (! input_name. empty ())
  	      td << input_name;
  	    td << al1->targetName
           << (stxS + stxType)
           << operonType
           << getIdentity () * 100.0
           << al1->targetStart + 1
           << al2->targetEnd
           << al1->targetStrand
           // Approximately if frameshift
           << getA () -> refAccession
           << getA () -> getIdentity () * 100.0
           << getA () -> getCoverage () * 100.0
           << getB () -> refAccession
           << getB () -> getIdentity () * 100.0
           << getB () -> getCoverage () * 100.0
           ;
        td. newLn ();
      }
      else
        al1->saveTsvOut (td);
    }


private:
  const BlastAlignment* getA () const
    { return al1->targetStrand ? al1 : al2; }
  const BlastAlignment* getB () const
    { return al1->targetStrand ? al2 : al1; }
  string getStxType () const
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
      if (verbose ())
        return string ("2 ") + a [312] + a [318] + b [34];
      return "2";
    }
  bool partial () const
    { return    getA () -> getCoverage () < 1.0
             || getB () -> getCoverage () < 1.0;
    }
  bool truncated () const
    { return    getA () -> truncated ()
             || getB () -> truncated ();
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
      return false;
    }
  static bool reportLess (const Operon &a,
                          const Operon &b)
    { LESS_PART (a, b, al1->targetName);
      LESS_PART (a, b, al1->targetStart);
      LESS_PART (a, b, al1->targetEnd);
      LESS_PART (b, a, al1->targetStrand);
      return false;
    }
};




// ThisApplication

struct ThisApplication : ShellApplication
{
  ThisApplication ()
    : ShellApplication ("Determine stx type(s) of a genome, print .tsv-file", true, true, true)
    {
    	addKey ("nucleotide", "Input nucleotide FASTA file (can be gzipped)", "", 'n', "NUC_FASTA");
    //addKey ("database", "Alternative directory with AMRFinder database. Default: $AMRFINDER_DB", "", 'd', "DATABASE_DIR");
    //addFlag ("database_version", "Print database version", 'V');
    	addKey ("translation_table", "NCBI genetic code for translated BLAST", "11", 't', "TRANSLATION_TABLE");
      addKey ("name", "Text to be added as the first column \"name\" to all rows of the report, for example it can be an assembly name", "", '\0', "NAME");
      addKey ("output", "Write output to OUTPUT_FILE instead of STDOUT", "", 'o', "OUTPUT_FILE");
    	addKey ("blast_bin", "Directory for BLAST. Deafult: $BLAST_BIN", "", '\0', "BLAST_DIR");

      version = SVN_REV;
    }



  void shellBody () const final
  {
    const string  fName            = shellQuote (getArg ("nucleotide"));
  //      string  db               =             getArg ("database");
  //const bool    database_version =             getFlag ("database_version");
    const uint    gencode          =             arg2uint ("translation_table"); 
                  input_name       =             getArg ("name");
    const string  output           =             getArg ("output");
          string  blast_bin        =             getArg ("blast_bin");
    
    if (contains (input_name, '\t'))
      throw runtime_error ("NAME cannot contain a tab character");


  #if 0
    if (database_version)
      cout   << "Software directory: " << shellQuote (execDir) << endl;
    else
  #endif
      stderr << "Software directory: " << shellQuote (execDir) << '\n';
  #if 0
    if (database_version)
      cout   << "Software version: " << version << endl; 
    else
	    stderr << "Software version: " << version << '\n'; 
	#endif
    
    OFStream::prepare (output);
        
		const string logFName (tmp + "/log");  // Command-local log file


  #if 0
    string defaultDb;
    #ifdef CONDA_DB_DIR
    // we're in condaland
      if (const char* s = getenv("CONDA_PREFIX")) {
        defaultDb = string (s) + "/share/amrfinderplus/data/latest";
      } else if (const char* s = getenv("PREFIX")) {
        const Warning warning (stderr);
        stderr << "This was compiled for running under bioconda, but $CONDA_PREFIX was not found" << '\n';
        defaultDb = string (s) + "/share/amrfinderplus/data/latest";
        stderr << "Reverting to $PREFIX: " << defaultDb;
      } else {
        const Warning warning (stderr);
        stderr << "This was compiled for running under bioconda, but $CONDA_PREFIX was not found" << '\n';
        stderr << "Reverting to hard coded directory: " << CONDA_DB_DIR "/latest";
        defaultDb = CONDA_DB_DIR "/latest";
      }
    #else
    // not in condaland
      defaultDb = execDir + "data/latest";
    #endif
    ASSERT (isRight (defaultDb, "/latest"));
        
		// db
		if (db. empty ())
		{
    	if (const char* s = getenv ("AMRFINDER_DB"))
    		db = s;
    	else
			  db = defaultDb;
		}
		ASSERT (! db. empty ());		  
  #endif


    // blast_bin
    if (blast_bin. empty ())
    	if (const char* s = getenv ("BLAST_BIN"))
    		blast_bin = string (s);
    if (! blast_bin. empty ())
    {
	    addDirSlash (blast_bin);
	    prog2dir ["blastx"] = blast_bin;
	  }


  #if 0
    const string downloadLatestInstr ("\nTo download the latest version to the default directory run: amrfinder -u");
    
		if (! directoryExists (db))
		  throw runtime_error ("No valid AMRFinder database is found.\nThis directory (or symbolic link to directory) is not found: " + db + downloadLatestInstr);
    if (database_version)
      cout   << "Database directory: " << shellQuote (path2canonical (db)) << endl;
    else
		  stderr << "Database directory: " << shellQuote (path2canonical (db)) << '\n';
    setSymlink (db, tmp + "/db", true);

		if (! fileExists (db + "/AMRProt.phr"))
			throw runtime_error ("The BLAST database for AMRProt was not found. Use amrfinder -u to download and prepare database for AMRFinderPlus");


		// PD-3051
		{
  	  istringstream versionIss (version);
  		const SoftwareVersion softwareVersion (versionIss);
  		const SoftwareVersion softwareVersion_min (db + "/database_format_version.txt"); 
  		const DataVersion dataVersion (db + "/version.txt");
  		istringstream dataVersionIss (DATA_VER_MIN); 
  		const DataVersion dataVersion_min (dataVersionIss);  
      if (database_version)
        cout   << "Database version: " << dataVersion. str () << endl;
      else
        stderr << "Database version: " << dataVersion. str () << '\n';
      if (softwareVersion < softwareVersion_min)
        throw runtime_error ("Database requires software version at least " + softwareVersion_min. str ());
      if (dataVersion < dataVersion_min)
        throw runtime_error ("Software requires database version at least " + dataVersion_min. str () + downloadLatestInstr);
      if (database_version)
        return;
    }
  #endif


    const string qcS (qc_on ? " -qc" : "");

    const string dna_flat = uncompress (fName,  "dna_flat");
    
    size_t nDna = 0;
    size_t dnaLen_max = 0;
    size_t dnaLen_total = 0;
    {
      prog2dir ["fasta_check"] = execDir;
      exec (fullProg ("fasta_check") + dna_flat + "  -hyphen  -ambig  " + qcS + "  -log " + logFName + " > " + tmp + "/nseq", logFName); 
    	const StringVector vec (tmp + "/nseq", (size_t) 10, true); 
    	if (vec. size () != 3)
        throw runtime_error ("fasta_check failed: " + vec. toString ("\n"));
      nDna         = str2<size_t> (vec [0]);
      dnaLen_max   = str2<size_t> (vec [1]);
      dnaLen_total = str2<size_t> (vec [2]);
    }
    QC_ASSERT (nDna);
    QC_ASSERT (dnaLen_max);
    QC_ASSERT (dnaLen_total);

		stderr. section ("Running blastx");
		{
			const Chronometer_OnePass_cerr cop ("blastx");
  		#define BLAST_FMT  "-outfmt '6 qseqid sseqid length nident qstart qend qlen sstart send slen qseq sseq'"
 			findProg ("blastx");
 			// Optmize ??
			exec (fullProg ("blastx") + " -query " + dna_flat + " -db " + execDir + "stx.prot  "  // tmp + "/db/AMRProt" /* /db/stx ??*/  + "  " 
			      + "-comp_based_stats 0  -evalue 1e-10  -seg no  -max_target_seqs 10000  -word_size 5  -query_gencode " + to_string (gencode) + " "
			      + getBlastThreadsParam ("blastx", min (nDna, dnaLen_total / 10002)) 
			      + " " BLAST_FMT " -out " + tmp + "/blastx > /dev/null 2> " + tmp + "/blastx-err", tmp + "/blastx-err");
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
    
    
    ostream* tsvOut = & cout;
    unique_ptr<ostream> tsvOutS;
    if (! output. empty ())
    {
      tsvOutS. reset (new OFStream (output));
      tsvOut = tsvOutS. get ();
    }
    ASSERT (tsvOut);
    TsvOut td (*tsvOut, 2, false);
    
    if (! input_name. empty ())
      td << "name";
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

	  Vector<BlastAlignment> blastAls;   
    {
      LineInput f (tmp + "/blastx");
  	  while (f. nextLine ())
  	  {
  	    const Unverbose unv;
  	    if (verbose ())
  	      cout << f. line << endl;
  	    bool stx = true;
  	    BlastAlignment al (f. line, stx);
  	    if (! stx)
  	      continue;
  	    al. qc ();  
 	      blastAls << std::move (al);
  	  }
  	}
  	if (verbose ())
  	  cout << "# All stx blasts: " << blastAls. size () << endl;

    if (verbose ())
      cout << "Finding frame shifts:" << endl;
	  {
      // Multiple frame shifts are possible
      const BlastAlignment* prev = nullptr;
      blastAls. sort (BlastAlignment::frameshiftLess); 
      for (const BlastAlignment& al : blastAls)
      {        
        if (verbose ())
          al. saveTsvOut (td);
        if (   prev
            && al. targetName   == prev->targetName
            && al. targetStrand == prev->targetStrand
            && al. refAccession == prev->refAccession
            && al. targetStart  >  prev->targetStart
            && (int) al. targetStart - (int) prev->targetEnd < 10  // PAR
            && al. getFrame () != prev->getFrame ()
           )
        {
          var_cast (al). merge (*prev);
          al. qc ();
          var_cast (prev) -> reported = true;
        }
        prev = & al;
      }
    }
    
    if (verbose ())
      cout << "All blasts:" << endl;
	  VectorPtr<BlastAlignment> goodBlastAls;   
	  {
	    size_t start = 0;
      blastAls. sort ();
      FFOR (size_t, i, blastAls. size ())
      {
        const BlastAlignment& al = blastAls [i];
        if (verbose ())
          al. saveTsvOut (td);
        while (   start < i 
               && ! (   blastAls [start]. targetName   == al. targetName
                     && blastAls [start]. targetStrand == al. targetStrand
                     && blastAls [start]. stxClass     == al. stxClass
                     && blastAls [start]. subunit      == al. subunit
                     && blastAls [start]. targetEnd    >  al. targetStart
                    )
              )
          start++;
        if (al. reported)
          continue;
        bool suppress = false;
        FOR_START (size_t, j, start, i)
        {
          const BlastAlignment& prev = blastAls [j];
          ASSERT (al. targetName   == prev. targetName);
          ASSERT (al. targetStrand == prev. targetStrand);
          ASSERT (al. stxClass     == prev. stxClass);
          ASSERT (al. subunit      == prev. subunit);
          if (   al. insideEq (prev)
              && al. nident <= prev. nident
             )
          {
            suppress = true;
            break;
          }
        }
        if (! suppress)
          goodBlastAls << & al;
      }
    }
    
    if (verbose ())
      cout << endl << "Good blasts:" << endl;
    Vector<Operon> operons;
	  {
	    size_t start = 0;
      FFOR (size_t, i, goodBlastAls. size ())
      {
        const BlastAlignment* alB = goodBlastAls [i];
        ASSERT (alB);
        if (verbose ())
          alB->saveTsvOut (td);
        if (alB->subunit != 'B')
          continue;
        while (   start < i 
               && ! (   goodBlastAls [start] -> targetName   == alB->targetName
                     && goodBlastAls [start] -> targetStrand == alB->targetStrand
                     && goodBlastAls [start] -> stxClass     == alB->stxClass
                    )
              )
          start++;
        FOR_START (size_t, j, start, i)
        {
          const BlastAlignment* alA = goodBlastAls [j];
          ASSERT (alA);
          ASSERT (alA->targetName   == alB->targetName);
          ASSERT (alA->targetStrand == alB->targetStrand);
          ASSERT (alA->stxClass     == alB->stxClass);
          ASSERT (alA->subunit      <= alB->subunit);
          if (alA->subunit == alB->subunit)
            break;
          ASSERT (alA->subunit == 'A');
          const BlastAlignment* al1 = alA;
          const BlastAlignment* al2 = alB;
          if (! al1->targetStrand)
            swap (al1, al2);
          if (   al1->targetEnd < al2->targetStart
              && al2->targetStart - al1->targetEnd <= Operon::intergenic_max
             )
          {
            Operon op (*al1, *al2);
            if (verbose ())
            {
              cout << "Operon:" << '\t' << op. getIdentity () << '\t' << stxClass2identity [op. al1->stxClass] << endl;  
              op. saveTsvOut (td);  
            }
            if (op. getIdentity () >= stxClass2identity [op. al1->stxClass])
            {
              operons << std::move (op);
              var_cast (al1) -> reported = true;
              var_cast (al2) -> reported = true;
            }
          }
        }
      }
    }
  	if (verbose ())
  	  cout << "# Operons: " << operons. size () << endl;

    if (verbose ())
      cout << endl << "stx1/stx2 operons:" << endl;
	  {
	    size_t start = 0;
	    goodBlastAls. sort (BlastAlignment::stxLess);
      FFOR (size_t, i, goodBlastAls. size ())
      {
        const BlastAlignment* alB = goodBlastAls [i];
        ASSERT (alB);
        if (alB->reported)
          break;
        if (verbose ())
          alB->saveTsvOut (td);
        if (alB->subunit != 'B')
          continue;
        while (   start < i 
               && ! (   goodBlastAls [start] -> targetName   == alB->targetName
                     && goodBlastAls [start] -> targetStrand == alB->targetStrand
                    )
              )
          start++;
        FOR_START (size_t, j, start, i)
        {
          const BlastAlignment* alA = goodBlastAls [j];
          ASSERT (alA);
        //ASSERT (! alA->reported);
          ASSERT (alA->targetName   == alB->targetName);
          ASSERT (alA->targetStrand == alB->targetStrand);
          ASSERT (alA->subunit      <= alB->subunit);
          if (alA->subunit == alB->subunit)
            break;
          ASSERT (alA->subunit == 'A');
          const BlastAlignment* al1 = alA;
          const BlastAlignment* al2 = alB;
          if (! al1->targetStrand)
            swap (al1, al2);
          if (   al1->targetEnd < al2->targetStart
              && al2->targetStart - al1->targetEnd <= Operon::intergenic_max
             )
          {
            Operon op (*al1, *al2);
            if (verbose ())
              op. saveTsvOut (td);  
            operons << std::move (op);
            var_cast (al1) -> reported = true;
            var_cast (al2) -> reported = true;
          }
        }
      }
    }
  	if (verbose ())
  	  cout << "# Operons: " << operons. size () << endl;
   	  
  	if (verbose ())
  	  cout << endl << "goodOperons" << endl;
    Vector<Operon> goodOperons;
    {    
      operons. sort ();
      for (const Operon& op : operons)
      {
        if (verbose ())
       	  op. saveTsvOut (td); 
       	op. qc ();     
        bool found = false;
        // Order by targetName, targetStart for speed ??
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

  	if (verbose ())
  	  cout << endl << "Suppress goodBlastAls by goodOperons" << endl;
    for (const BlastAlignment* al : goodBlastAls)
    {
      ASSERT (al);
      if (! al->reported)
        for (const Operon& op : goodOperons)
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

  	if (verbose ())
  	  cout << endl << "goodBlastAls -> goodOperons" << endl;
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
            && al2->nident <= al1->nident
           )
          var_cast (al2) -> reported = true;
      }
    }

    // Report
    goodOperons. sort (Operon::reportLess);     
  	for (const Operon& op : goodOperons)
   	  op. saveTsvOut (td);
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



