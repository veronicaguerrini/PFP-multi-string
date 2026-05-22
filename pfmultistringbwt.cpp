/* **************************************************************************
 * pfmultistringbwt.cpp
 * Output the multi-string BWT, the SA (option -S) or the DA (option -D)
 * computed using the prefix-free parsing technique
 *  
 * Usage:
 *   pfbwt[NT][64].x -h
 * for usage info
 *
 * See newscan.cpp for a description of what it does  
 * 
 **************************************************************************** */
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <semaphore.h>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
extern "C" {
#include "gsa/gsacak.h"
#include "utils.h"
}

using namespace std;
// using namespace __gnu_cxx;

//for SAP values
unsigned char M[]={0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01};
#define get_sap(i) ( (sap[(i)/8]&M[(i)%8]) ? 1 : 0 )

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
   char *basename;
   string parseExt =  EXTPARSE;    // extension final parse file  
   string occExt =    EXTOCC;      // extension occurrences file  
   string dictExt =   EXTDICT;     // extension dictionary file  
   //04_2026 string lastExt =   EXTLST;      // extension file containing last chars   
   string saExt =     EXTSAI;      // extension file containing sa info   
   int w = 10;            // sliding window size and its default 
   int th = 0;            // number of helper threads, default none 
   bool SA = false;       // output all SA values
   bool DA = false;       // output all DA values 04_2026
   int sampledSA = 0;// output sampled SA values
};

// mask for sampled SA: start of a BWT run, end of a BWT run, or both 
#define START_RUN 1
#define END_RUN 2


static long get_num_words(uint8_t *d, long n);
static long binsearch(uint_t x, uint_t a[], long n);
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid);
//04_2026
static void compute_dict_bwt_sap(uint8_t *d, long dsize,long dwords, int w, uint_t **sat, uint8_t **sapt);
static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, uint32_t *ilist, uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts);
static void fwrite_chars_same_suffix_sa(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, uint32_t *ilist, uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
                                     int_t suffixLen, FILE *safile, uint8_t *bwsainfo,long);
static void fwrite_chars_same_suffix_da(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, uint32_t *ilist, uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
                                     FILE *dafile, uint32_t *bwdainfo,long); //04_2026
static void fwrite_chars_same_suffix_ssa(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, uint32_t *ilist, uint32_t *istart, FILE *fbwt, long &easy_bwts, long &hard_bwts,
                                     int_t suffixLen, FILE *ssafile, FILE *esafile, uint8_t *bwsainfo,long, int &, uint64_t &, int);
static uint8_t *load_bwsa_info(Args &arg, long n);
static uint32_t *load_bwda_info(Args &arg, long n); //04_2026

// class representing the suffix of a dictionary word
// instances of this class are stored to a heap to handle the hard bwts
struct SeqId {
  uint32_t id;       // lex. id of the dictionary word to which the suffix belongs
  int remaining;     // remaining copies of the suffix to be considered  
  uint32_t *bwtpos;  // list of bwt positions of this dictionary word
  uint8_t char2write;// char to be written (is the one preceeding the suffix)

  // constructor
  SeqId(uint32_t i, int r, uint32_t *b, int8_t c) : id(i), remaining(r), bwtpos(b) {
    char2write = c;
  }

  // advance to the next bwt position, return false if there are no more positions 
  bool next() {
    remaining--;
    bwtpos += 1;
    return remaining>0;
  }
  bool operator<(const SeqId& a) const;
};

bool SeqId::operator<(const SeqId& a) const {
    return *bwtpos > *(a.bwtpos);
}



/* *******************************************************************
 * Computation of the final BWT
 * 
 * istart[] and islist[] are used together. For each dictionary word i 
 * (considered in lexicographic order) for k=istart[i]...istart[i+1]-1
 * ilist[k] contains the ordered positions in BWT(P) containing word i 
 * ******************************************************************* */
void bwt(Args &arg, uint8_t *d, uint8_t *bitdict, long dsize, // dictionary and its size 
         uint32_t *ilist, long psize, // ilist and their size 
         uint32_t *istart, long dwords) // starting point in ilist for each word and # words
{  
  // possibly read bwsa info file and open sa output file
  uint8_t *bwsainfo = load_bwsa_info(arg,psize);
  uint32_t *bwdainfo = load_bwda_info(arg,psize);
  FILE *safile=NULL;
  FILE *dafile=NULL;
  // open the necessary SA/DA files (possibly none)
  if(arg.SA) safile = open_aux_file(arg.basename,EXTSA,"wb");
  if(arg.DA) dafile = open_aux_file(arg.basename,"da","wb"); 

  // compute sa and bwt of d and do some checking on them 
  uint_t *sa; uint8_t *sap;
  compute_dict_bwt_sap(d,dsize,dwords,arg.w,&sa,&sap); 

  // derive eos from sa. for i=0...dwords-1, eos[i] is the eos position of string i in d
  uint_t *eos = sa+1;
  // note that sa[0] = dsize-1 the position of the final 0x0 in d[]
  assert((long) sa[0]==dsize-1);
  // the bwt algorithm orders the 0x1 symbols in position order 
  // (ie they are considered different symbols)
  for(int i=0;i<dwords-1;i++)
    assert(eos[i]<eos[i+1]); // important since we'll use binary search on eos[]

  // open output file
  FILE *fbwt = open_aux_file(arg.basename,"bwt","wb");
    
  //We need to consider the dictionary word suffixes differently according to bitdict:
  //if bitdict[i]=0, the dict word IS NOT ending with $, then all the proper suffixes of length >=w MUST be considered
  //otherwise, the dict word IS ending with $, then all the proper suffixes of length >=1 MUST be considered
  
  time_t start = time(NULL);

  long easy_bwts = 0;
  long hard_bwts = 0;
  long next;
  uint32_t seqid;

  long i=dwords+1;
  // vector to save seqid and the corresponding BWT char 
  vector<uint32_t> id2merge; 
  vector<uint8_t> char2write;
  int_t suffixLen=0;
  //loop for $0x1... suffixes (i.e. $-suffixes)
  assert(d[sa[i]+1]==EndOfWord);
  //loop for the $-suffixes
  while(i<dsize && d[sa[i]+1]==EndOfWord){ 
    suffixLen = getlen(sa[i],eos,dwords,&seqid);
    assert(bitdict[seqid]);
    id2merge.push_back(seqid);           // sequence to consider
    char2write.push_back(d[sa[i]-1]);    // corresponding char
    i++;
  }
  printf("Processed %ld suffixes corresponding to 0-suffixes.\n",i-dwords-1);
  // output to fbwt the bwt chars corresponding to the current dictionary suffix, and, if requested, some SA values 
  if(arg.SA)
    fwrite_chars_same_suffix_sa(id2merge,char2write,ilist,istart,fbwt,easy_bwts,hard_bwts,suffixLen,safile,bwsainfo,psize);
  else if(arg.DA!=0)
    fwrite_chars_same_suffix_da(id2merge,char2write,ilist,istart,fbwt,easy_bwts,hard_bwts,dafile,bwdainfo,psize);
  else 
    fwrite_chars_same_suffix(id2merge,char2write,ilist,istart,fbwt,easy_bwts,hard_bwts);
  //Clear vectors
  id2merge.clear();
  char2write.clear();

  for(; i<dsize; i=next ) {
    // we are considering d[sa[i]....]
    next = i+1;  // prepare for next iteration  
    // compute length of this suffix and sequence it belongs
    //seqid identifies a word in the lex sorted dictionary and thus it can be used as index of bitdict
    suffixLen = getlen(sa[i],eos,dwords,&seqid);

    // ignore suffixes of lenght < w if bitdict=0
    if(bitdict[seqid]==0 && suffixLen<arg.w) continue; 
    // ignore suffixes that are full words
    if(sa[i]==0 || d[sa[i]-1]==EndOfWord) continue;

    id2merge.push_back(seqid);           // sequence to consider
    char2write.push_back(d[sa[i]-1]);    // corresponding char

    //we need to check for equal suffixes starting at i
    if(next<dsize){
      assert(sa[next]>0);
      while(d[sa[next]-1]!=EndOfWord && get_sap(next)==1) {
        int_t nextsuffixLen = getlen(sa[next],eos,dwords,&seqid);
        assert(nextsuffixLen==suffixLen);
        id2merge.push_back(seqid);           // sequence to consider
        char2write.push_back(d[sa[next]-1]);  // corresponding char
        next++;
      }
    }

    // output to fbwt the bwt chars corresponding to the current dictionary suffix, and, if requested, some SA values 
    if(arg.SA)
      fwrite_chars_same_suffix_sa(id2merge,char2write,ilist,istart,fbwt,easy_bwts,hard_bwts,suffixLen,safile,bwsainfo,psize);
    else if(arg.DA!=0)
      fwrite_chars_same_suffix_da(id2merge,char2write,ilist,istart,fbwt,easy_bwts,hard_bwts,dafile,bwdainfo,psize);
    else  
      fwrite_chars_same_suffix(id2merge,char2write,ilist,istart,fbwt,easy_bwts,hard_bwts);
    
    //Clear vectors
    id2merge.clear();
    char2write.clear();
  }

  cout << "Easy bwt chars: " << easy_bwts << endl;
  cout << "Hard bwt chars: " << hard_bwts << endl;
  cout << "Generating the final BWT took " << difftime(time(NULL),start) << " wall clock seconds\n";    
  fclose(fbwt);
  delete[] sap;
  delete[] sa;
  if(arg.SA) free(bwsainfo);
  if(arg.DA) free(bwdainfo);
  if(arg.SA) fclose(safile);
  if(arg.DA) fclose(dafile);
}  

void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-h  \tshow help and exit" << endl
        << "\t-S  \tcompute full suffix array" << endl
        << "\t-D  \tcompute full document array" << endl;
  exit(1);
}

void parseArgs( int argc, char** argv, Args& arg ) {
  int c;
  extern char *optarg;
  extern int optind;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

   string sarg;
   while ((c = getopt( argc, argv, "t:w:sehSD") ) != -1) {
      switch(c) {
        case 'S':
        arg.SA = true; break;
        case 'D':
        arg.DA = true; break; //04_2026
        case 'w':
        sarg.assign( optarg );
        arg.w = stoi( sarg ); break;
        case 'h':
           print_help(argv, arg); exit(1);
        case '?':
        cout << "Unknown option. Use -h for help." << endl;
        exit(1);
      }
   }
   // the only input parameter is the file name
   arg.basename = NULL; 
   if (argc == optind+1) 
     arg.basename = argv[optind];
   else {
      cout << "Invalid number of arguments" << endl;
      print_help(argv,arg);
   }
   // check check sa and da
   if(arg.SA && arg.DA!=0) {
     cout << "Computing both the full SA and the DA is not implemented yet";
     exit(1);
   }
   /////
   if(arg.w <4) {
     cout << "Windows size must be at least 4\n";
     exit(1);
   }
}


int main(int argc, char** argv)
{
  time_t start = time(NULL);  

  // translate command line parameters
  Args arg;
  parseArgs(argc, argv, arg);
  // read dictionary file 
  FILE *g = open_aux_file(arg.basename,EXTDICT,"rb");
  fseek(g,0,SEEK_END);
  long dsize = ftell(g);
  if(dsize<0) die("ftell dictionary");
  if(dsize<=1+arg.w) die("invalid dictionary file");
  cout  << "Dictionary file size: " << dsize << endl;
  #if !M64
  if(dsize > 0x7FFFFFFE) {
    printf("Dictionary size greater than  2^31-2!\n");
    printf("Please use 64 bit version\n");
    exit(1);
  }
  #endif

  uint8_t *d = new uint8_t[dsize];  
  rewind(g);
  long e = fread(d,1,dsize,g);
  if(e!=dsize) die("fread");
  fclose(g);
  
  // read occ file
  g = open_aux_file(arg.basename,EXTOCC,"rb");
  fseek(g,0,SEEK_END);
  e = ftell(g);
  if(e<0) die("ftell occ file");
  if(e%4!=0) die("invalid occ file");
  int dwords = e/4;
  cout  << "Dictionary words: " << dwords << endl;
  uint32_t *occ = new uint32_t[dwords+1];  // dwords+1 since istart overwrites occ
  rewind(g);
  e = fread(occ,4,dwords,g);
  if(e!=dwords) die("fread 2");
  fclose(g);
  assert(dwords==get_num_words(d,dsize));

  // read bitdict file 
  g = open_aux_file(arg.basename,"bitdict","rb");
  uint8_t *bitdict = new uint8_t[dwords];
  e = fread(bitdict,1,dwords,g);
  if(e!=dwords) die("fread bitdict");
  fclose(g);
  //---

  // read ilist file 
  g = open_aux_file(arg.basename,EXTILIST,"rb");
  fseek(g,0,SEEK_END);
  e = ftell(g);
  if(e<0) die("ftell ilist file");
  if(e%4!=0) die("invalid ilist file");
  long psize = e/4;
  cout  << "Parsing size: " << psize << endl;
  if(psize>0xFFFFFFFEL) die("More than 2^32 -2 words in the parsing");
  uint32_t *ilist = new uint32_t[psize];  
  rewind(g);
  e = fread(ilist,4,psize,g);
  if(e!=psize) die("fread 3");
  fclose(g);
  assert(ilist[0]==0); // EOF is in PBWT[0] 

  // convert occ entries into starting positions inside ilist
  // ilist also contains the position of EOF but we don't care about it since it is not in dict 
  uint32_t last=1; // starting position in ilist of the smallest dictionary word  
  for(long i=0;i<dwords;i++) {
    uint32_t tmp = occ[i];
    occ[i] = last;
    last += tmp;
  }
  assert(last==psize);
  occ[dwords]=psize;
  
  // compute and write the final bwt 
  bwt(arg,d,bitdict,dsize,ilist,psize,occ,dwords); // version not using threads 04_2024

  delete[] ilist;
  delete[] occ;
  delete[] d;  
  cout << "==== Elapsed time: " << difftime(time(NULL),start) << " wall clock seconds\n";      
  return 0;
}

// --------------------- aux functions ----------------------------------

static uint8_t *load_bwsa_info(Args &arg, long n)
{  
  // maybe sa info is not really needed 
  if(arg.SA==false and arg.sampledSA==0) return NULL;
  // open .bwsa file for reading and .bwlast for writing
  FILE *fin = open_aux_file(arg.basename,EXTBWSAI,"rb");
  // allocate and load the bwsa array
  uint8_t *sai = (uint8_t *) malloc(n*IBYTES);
  if(sai==NULL) die("malloc failed (BWSA INFO)"); 
  long s = fread(sai,IBYTES,n,fin);
  if(s!=n) die("bwsa info read");
  if(fclose(fin)!=0) die("bwsa info file close");
  return sai;
}

static uint32_t *load_bwda_info(Args &arg, long n)
{  
  // maybe da info is not really needed 
  if(arg.DA==false) return NULL;
  // open .bwda file for reading
  FILE *fin = open_aux_file(arg.basename,"bwdai","rb");
  // allocate and load the bwsa array
  uint32_t *dai = (uint32_t *) malloc(n*sizeof(uint32_t));
  if(dai==NULL) die("malloc failed (BWDA INFO)"); 
  long s = fread(dai,sizeof(uint32_t),n,fin);
  if(s!=n) die("bwda info read");
  if(fclose(fin)!=0) die("bwda info file close");
  return dai;
}

// compute the number of words in a dictionary
static long get_num_words(uint8_t *d, long n)
{
  long i,num=0;
  for(i=0;i<n;i++)
    if(d[i]==EndOfWord) num++;
  assert(d[n-1]==EndOfDict);
  return num;
}

// binary search for x in an array a[0..n-1] that doesn't contain x
// return the lowest position that is larger than x
static long binsearch(uint_t x, uint_t a[], long n)
{
  long lo=0; long hi = n-1;
  while(hi>lo) {
    assert( ((lo==0) || x>a[lo-1]) && x< a[hi]);
    int mid = (lo+hi)/2;
    assert(x!=a[mid]);  // x is not in a[]
    if(x<a[mid]) hi = mid;
    else lo = mid+1;
  }
  assert(((hi==0) || x>a[hi-1]) && x< a[hi]);
  return hi; 
}


// return the length of the suffix starting in position p.
// also write to seqid the id of the sequence containing that suffix 
// n is the # of distinct words in the dictionary, hence the length of eos[]
static int_t getlen(uint_t p, uint_t eos[], long n, uint32_t *seqid)
{
  assert(p<eos[n-1]);
  *seqid = binsearch(p,eos,n);
  assert(eos[*seqid]> p); // distance between position p and the next $
  return eos[*seqid] - p;
}

// compute the SA and SAP array for the set of (unique) dictionary words using gSACA-K
static void compute_dict_bwt_sap(uint8_t *d, long dsize,long dwords, int w, uint_t **sat, uint8_t **sapt) // output parameters
{
  uint_t *sa = new uint_t[dsize];
  uint8_t *sap = (uint8_t *)malloc(dsize/8+1);
  (void) dwords; (void) w;

  cout  << "Each SA entry: " << sizeof(*sa) << " bytes\n";
  cout  << "Each SAP entry: " << sizeof(*sap) << " bytes\n";

  cout << "Computing SA and SAP of dictionary" << endl;

  cout << "dwords=" << dwords<< endl;
  time_t  start = time(NULL);
  gsacak_sap(d,sa,sap,dsize);
  cout << "Computing SA/SAP took " << difftime(time(NULL),start) << " wall clock seconds\n";  

  // ------ do some checking on the sa
  assert(d[dsize-1]==EndOfDict);
  assert(sa[0]==(unsigned long)dsize-1);// sa[0] is the EndOfDict symbol 
  for(long i=0;i<dwords;i++) 
    assert(d[sa[i+1]]==EndOfWord); // there are dwords EndOfWord symbols 
 
  assert(sa[dwords]==(unsigned long)dsize-2);  
  assert(d[sa[dwords+1]]>Dollar);  // end of Dollar chars in the first column

  // copy sa and lcp address
  *sat = sa;  *sapt = sap;  
}


// write to the bwt all the characters preceding a given suffix
// doing a merge operation if necessary
static void fwrite_chars_same_suffix(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, 
                                    uint32_t *ilist, uint32_t *istart,
                                    FILE *fbwt, long &easy_bwts, long &hard_bwts)
{
  size_t numwords = id2merge.size(); // numwords dictionary words contain the same suffix
  bool samechar=true;

  for(size_t i=1;(i<numwords)&&samechar;i++)
    samechar = (char2write[i-1]==char2write[i]); 
  if(samechar) {
    for(size_t i=0; i<numwords; i++) {
      uint32_t s = id2merge[i];
      for(long j=istart[s];j<istart[s+1];j++){ //write occ[s] times the BWT char 
        if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
      }
      easy_bwts +=  istart[s+1]- istart[s]; 
    }
  }
  else {  // many words, many chars...    
    vector<SeqId> heap; // create heap
    for(size_t i=0; i<numwords; i++) {
      uint32_t s = id2merge[i]; //s is the dict word ID
      heap.push_back(SeqId(s,istart[s+1]-istart[s], ilist+istart[s], char2write[i]));
      //s                        lex. id of the dictionary word to which the suffix belongs
      //istart[s+1]-istart[s]    remaining copies of the suffix to be considered  (occ)
      //ilist+istart[s]          list of bwt positions of this dictionary word
      //char2write               char to be written (is the one preceeding the suffix)
    }
    std::make_heap(heap.begin(),heap.end()); //sort according to bwt positions of the dictionary words in the parse
    while(heap.size()>0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      if(fputc(s.char2write,fbwt)==EOF) die("BWT write error 2");
      hard_bwts += 1;
      // remove top 
      pop_heap(heap.begin(),heap.end());
      heap.pop_back();
      // if remaining positions, reinsert to heap
      if(s.next()) {
        heap.push_back(s);
        push_heap(heap.begin(),heap.end());
      }
    }
  }
}

static void fwrite_chars_same_suffix_da(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, 
                                    uint32_t *ilist, uint32_t *istart,
                                    FILE *fbwt, long &easy_bwts, long &hard_bwts, FILE *dafile, uint32_t *bwdainfo, long n)
{
  size_t numwords = id2merge.size(); // numwords dictionary words contain the same suffix

  if(numwords==1) {
    uint32_t s = id2merge[0];
    for(long j=istart[s];j<istart[s+1];j++) {
      if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
      uint64_t da = bwdainfo[ilist[j]]-1; //0-based
      if(fwrite(&da,sizeof(uint32_t),1,dafile)!=1) die("DA write error 1");
    }
    easy_bwts +=  istart[s+1]- istart[s];
  }
  else {  // many words, many chars...     
    vector<SeqId> heap; // create heap
    for(size_t i=0; i<numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s,istart[s+1]-istart[s], ilist+istart[s], char2write[i]));
    }
    std::make_heap(heap.begin(),heap.end());
    while(heap.size()>0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      if(fputc(s.char2write,fbwt)==EOF) die("BWT write error 2");
      uint64_t da = bwdainfo[*(s.bwtpos)]-1; //0-based
      if(fwrite(&da,sizeof(uint32_t),1,dafile)!=1) die("DA write error 2");     
      hard_bwts += 1;
      // remove top 
      pop_heap(heap.begin(),heap.end());
      heap.pop_back();
      // if remaining positions, reinsert to heap
      if(s.next()) {
        heap.push_back(s);
        push_heap(heap.begin(),heap.end());
      }
    }
  }
}



// write to the bwt all the characters preceding a given suffix
// and the corresponding SA entries doing a merge operation
static void fwrite_chars_same_suffix_sa(vector<uint32_t> &id2merge,  vector<uint8_t> &char2write, 
                                    uint32_t *ilist, uint32_t *istart,
                                    FILE *fbwt, long &easy_bwts, long &hard_bwts,
                                    int_t suffixLen, FILE *safile, uint8_t *bwsainfo, long n)
{
  size_t numwords = id2merge.size(); // numwords dictionary words contain the same suffix

  if(numwords==1) {
    uint32_t s = id2merge[0];
    for(long j=istart[s];j<istart[s+1];j++) {
      if(fputc(char2write[0],fbwt)==EOF) die("BWT write error 1");
      uint64_t sa = get_myint(bwsainfo,n,ilist[j]) - suffixLen;
      if(fwrite(&sa,SABYTES,1,safile)!=1) die("SA write error 1");
    }
    easy_bwts +=  istart[s+1]- istart[s];
  }
  else {  // many words, many chars...     
    vector<SeqId> heap; // create heap
    for(size_t i=0; i<numwords; i++) {
      uint32_t s = id2merge[i];
      heap.push_back(SeqId(s,istart[s+1]-istart[s], ilist+istart[s], char2write[i]));
    }
    std::make_heap(heap.begin(),heap.end());
    while(heap.size()>0) {
      // output char for the top of the heap
      SeqId s = heap.front();
      if(fputc(s.char2write,fbwt)==EOF) die("BWT write error 2");
      uint64_t sa = get_myint(bwsainfo,n,*(s.bwtpos)) - suffixLen;
      if(fwrite(&sa,SABYTES,1,safile)!=1) die("SA write error 2");    
      hard_bwts += 1;
      // remove top 
      pop_heap(heap.begin(),heap.end());
      heap.pop_back();
      // if remaining positions, reinsert to heap
      if(s.next()) {
        heap.push_back(s);
        push_heap(heap.begin(),heap.end());
      }
    }
  }
}
