#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
}

using namespace std;
// using namespace __gnu_cxx;

// =============== algorithm limits =================== 
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;

// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT32_MAX)
typedef uint32_t occ_int_t;

#define EoString '$'
int pos_rel=0;

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  string str;
  occ_int_t occ;
  word_int_t rank=0;
};

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   int w = 10;            // sliding window size and its default 
   int p = 100;           // modulus for establishing stopping w-tuples 
   bool SAinfo = false;   // compute SA information
   bool DAinfo = false;   // compute SA information
   int th=0;              // number of helper threads
   int verbose=0;         // verbosity level 
};


// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
  int *window;
  int asize;
  const uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime 
  
  KR_window(int w): wsize(w) {
    asize = 256;
    asize_pot = 1;
    for(int i=1;i<wsize;i++) 
      asize_pot = (asize_pot*asize)% prime; // ugly linear-time power algorithm  
    // alloc and clear window
    window = new int[wsize];
    reset();
  }
  
  // init window, hash, and tot_char 
  void reset() {
    for(int i=0;i<wsize;i++) window[i]=0;
    // init hash value and related values
    hash=tot_char=0;  
  }
  
  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    // complex expression to avoid negative numbers 
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution  
    hash = (asize*hash + c) % prime;      //  add char i 
    window[k]=c;
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash; 
  }
  // debug only 
  string get_window() {
    string w = "";
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }
  
  ~KR_window() {
    delete[] window;
  } 

};
// -----------------------------------------------------------

// compute 64-bit KR hash of a string 
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    const uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k=0;k<s.size();k++) {
      int c = (unsigned char) s[k];
      assert(c>=0 && c< 256);
      hash = (256*hash + c) % prime;    //  add char k
    } 
    return hash; 
}



// save current word in the freq map and update it leaving only the 
// last minsize chars which is the overlap with next word  
static void save_update_word(Args& arg, string& w, map<uint64_t,word_stats>& freq, FILE *tmp_parse_file, FILE *sa, FILE *da, uint64_t &pos, uint32_t& num_startingWithDollar, uint32_t s) //2026
{
  size_t minsize = arg.w; 

  assert(pos_rel==0 || w.size() > minsize); 
  if(w.size() <= minsize) return; 
  // save overlap 
  string overlap(w.substr(w.size() - minsize)); // keep last minsize chars
  
  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");

  // update frequency table for current hash
  if(freq.find(hash)==freq.end()) {
     if(w[0]==EoString)
        num_startingWithDollar++; 
      freq[hash].occ = 1; // new hash
      freq[hash].str = w; 
  }
  else {
      freq[hash].occ += 1; // known hash
      if(freq[hash].occ <=0) {
        cerr << "Emergency exit! Maximum # of occurence of dictionary word (";
        cerr<< MAX_WORD_OCC << ") exceeded\n";
        exit(1);
      }
      if(freq[hash].str != w) {
        cerr << "Emergency exit! Hash collision for strings:\n";
        cerr << freq[hash].str << "\n  vs\n" <<  w << endl;
        exit(1);
      }
  }
  if(pos_rel==0) {
      pos += w.size()-1; // -1 is for the initial $ of the first word
      pos_rel=1;
  } else pos += w.size() -minsize;
  if(sa) if(fwrite(&pos,IBYTES,1,sa)!=1) die("Error writing to sa info file");
  if(da) if(fwrite(&s,sizeof(num_startingWithDollar),1,da)!=1) die("Error writing to da info file"); //04_2026
  
  // keep only the overlapping part of the window
  w.assign(overlap);
}

static void save_update_word_2(Args& arg, string& w, map<uint64_t,word_stats>& freq, FILE *tmp_parse_file, FILE *sa, FILE *da, uint64_t &pos, uint32_t& num_startingWithDollar, uint32_t s) //2026
{
  size_t minsize = arg.w; 
  if (w.size() <= minsize)
    cerr << "w = " << w << endl;
  assert(w.size() > minsize);
  // save overlap 
  string overlap(w.substr(w.size() - 1)); // keep last char
  
  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");

  // update frequency table for current hash
  if(freq.find(hash)==freq.end()) {
      if(w[0]==EoString) //2026
        num_startingWithDollar++;
      freq[hash].occ = 1; // new hash
      freq[hash].str = w; 
  }
  else {
      freq[hash].occ += 1; // known hash
      if(freq[hash].occ <=0) {
        cerr << "Emergency exit! Maximum # of occurence of dictionary word (";
        cerr<< MAX_WORD_OCC << ") exceeded\n";
        exit(1);
      }
      if(freq[hash].str != w) {
        cerr << "Emergency exit! Hash collision for strings:\n";
        cerr << freq[hash].str << "\n  vs\n" <<  w << endl;
        exit(1);
      }
  }

  if(pos_rel==0)  pos += w.size()-1; // -1 is for the initial $ of the first word
  else pos += w.size() -minsize;
  if(sa) if(fwrite(&pos,IBYTES,1,sa)!=1) die("Error writing to sa info file");
  if(da) if(fwrite(&s,sizeof(num_startingWithDollar),1,da)!=1) die("Error writing to da info file"); //04_2026
  
  // keep only the overlapping part of the window
  w.assign(overlap);
}


// prefix free parse of file fnam. w is the window size, p is the modulus 
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args& arg, map<uint64_t,word_stats>& wordFreq, uint32_t& num_tot_seqs,uint32_t& num_starting)
{
  uint64_t tot_char_read=0; //num char read
  //open a, possibly compressed, input file
  string fnam = arg.inputFileName;
  #ifdef GZSTREAM 
  igzstream f(fnam.c_str());
  #else
  ifstream f(fnam);
  #endif    
  if(!f.rdbuf()->is_open()) {// is_open does not work on igzstreams 
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + fnam);
  }

  // open the 1st pass parsing file 
  FILE *g = open_aux_file(arg.inputFileName.c_str(),EXTPARS0,"wb");
  FILE *sa_file = NULL, *da_file=NULL;
  if(arg.SAinfo) 
    sa_file = open_aux_file(arg.inputFileName.c_str(),EXTSAI,"wb");
  if(arg.DAinfo) 
    da_file = open_aux_file(arg.inputFileName.c_str(),"dai","wb");
  
  // main loop on the chars of the input file
  int c;
  uint64_t pos = 0; // ending position +1 of previous word in the original text, used for computing sa_info 
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file 
  
  string word("");
  // init empty KR window: constructor only needs window size
  KR_window krw(arg.w);
  while( (c = f.get()) != EOF ) {
    word.append(1,EoString);
    krw.reset();
    pos_rel=0;
    word.append(1,c);
    uint64_t hash = krw.addchar(c);
    assert(krw.tot_char==1);
    while( (c = f.get()) != EoString) { 
      if(c<=Dollar) {
        cerr << "Invalid char found in input file. Exiting...\n"; exit(1);
      }
      word.append(1,c);
      hash = krw.addchar(c);
      if(hash%arg.p==0) {
        save_update_word(arg,word,wordFreq,g,sa_file,da_file,pos,num_starting,num_tot_seqs+1); 
      }
    }
    word.append(1,EoString);
    save_update_word_2(arg,word,wordFreq,g,sa_file,da_file,pos,num_starting,num_tot_seqs+1);
    num_tot_seqs++; //2026
    word="";
    tot_char_read+=krw.tot_char+1;
    if(pos!=tot_char_read){
      cerr << "Pos: " << pos << " tot " << tot_char_read << endl;
    }
    assert(pos==tot_char_read);
  }
  // close input and output files 
  if(sa_file) if(fclose(sa_file)!=0) die("Error closing SA file");
  if(da_file) if(fclose(da_file)!=0) die("Error closing DA file");
  
  if(fclose(g)!=0) die("Error closing parse file");
  if(pos!=tot_char_read) cerr << "Pos: " << pos << " tot " << tot_char_read << endl;
  
  f.close();
  return tot_char_read;
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a <= *b;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, map<uint64_t,word_stats> &wfreq, vector<const string *> &sortedDict)
{
  assert(sortedDict.size() == wfreq.size());
  FILE *fdict, *fwlen=NULL, *focc=NULL;
  FILE *fbit=NULL; //04_2026
  // open dictionary and occ files
  fdict = open_aux_file(arg.inputFileName.c_str(),EXTDICT,"wb");
  focc = open_aux_file(arg.inputFileName.c_str(),EXTOCC,"wb");
  fbit = open_aux_file(arg.inputFileName.c_str(),"bitdict","wb"); 
  
  word_int_t wrank = 1, num_bitdict=0; // current word rank (1 based)
  for(auto x: sortedDict) {          // *x is the string representing the dictionary word
    const char *word = (*x).data();       // current dictionary word
    size_t len = (*x).size();  // offset and length of word
    uint8_t bit=0;
    if(word[len-1]==EoString){ //dict word ending with EoString
      bit=1;
      num_bitdict++;
    }
    //---
    //assert(len>(size_t)arg.w);
    uint64_t hash = kr_hash(*x);
    auto& wf = wfreq.at(hash);
    assert(wf.occ>0);
    size_t s = fwrite(word,1,len, fdict);
    if(s!=len) die("Error writing to DICT file");
    
    s = fwrite(&bit,1,1, fbit);
    if(s!=1) die("Error writing to BITDICT file");  
    //---
    if(fputc(EndOfWord,fdict)==EOF) die("Error writing EndOfWord to DICT file");
    s = fwrite(&wf.occ,sizeof(wf.occ),1, focc);
    if(s!=1) die("Error writing to OCC file");
    assert(wf.rank==0); // word should have no rank at this time
    wf.rank = wrank++;  // save the rank of the current word
  }
  cout << "Number of 1-bit in .bitdict: " << num_bitdict << endl;
  if(fclose(fbit)!=0) die("Error closing BITDICT file"); 
  if(fputc(EndOfDict,fdict)==EOF) die("Error writing EndOfDict to DICT file");
  if(fclose(focc)!=0) die("Error closing OCC file");
  if(fclose(fdict)!=0) die("Error closing DICT file");
}

void remapParse(Args &arg, map<uint64_t,word_stats> &wfreq, word_int_t offset) 
{
  // open parse files. the old parse can be stored in a single file or in multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), EXTPARS0, arg.th);
  FILE *newp = open_aux_file(arg.inputFileName.c_str(), EXTPARSE, "wb");
  
  // recompute occ as an extra check 
  vector<occ_int_t> occ(wfreq.size()+offset+1,0); 
  uint64_t hash;
  word_int_t rank, rank_tmp,r=1; 
  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    rank_tmp=wfreq.at(hash).rank;
    if (wfreq.at(hash).str[0]==EoString){      
      r++;
      rank=rank_tmp;
    }
    else{
      rank = rank_tmp+offset;
    }
    occ[rank]++;
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("Error writing to new parse file");
  }
  size_t s = fwrite(&r,sizeof(rank),1,newp);

  if(s!=1) die("Error writing to new parse file");
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");

  // check old and recomputed occ's coincide 
  for(auto& x : wfreq){
    if (x.second.str[0]==EoString)
      assert(x.second.occ==occ[x.second.rank]);
    else
      assert(x.second.occ==occ[x.second.rank+offset]);
  }
}
 



void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl       
        << "\t-h  \tshow help and exit" << endl
        << "\t-s  \tcompute suffix array info" << endl
        << "\t-d  \tcompute document array info" << endl; //04_2026
  #ifdef GZSTREAM
  cout << "If the input file is gzipped it is automatically extracted\n";
  #endif
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
   while ((c = getopt( argc, argv, "p:w:sdht:vc") ) != -1) {
      switch(c) {
        case 's':
        arg.SAinfo = true; break;
        case 'd':
        arg.DAinfo = true; break; 
        case 'w':
        sarg.assign( optarg );
        arg.w = stoi( sarg ); break;
        case 'p':
        sarg.assign( optarg );
        arg.p = stoi( sarg ); break;
        case 'v':
           arg.verbose++; break;
        case 'h':
           print_help(argv, arg); exit(1);
        case '?':
        cout << "Unknown option. Use -h for help." << endl;
        exit(1);
      }
   }
   // the only input parameter is the file name 
   if (argc == optind+1) {
     arg.inputFileName.assign( argv[optind] );
   }
   else {
      cout << "Invalid number of arguments" << endl;
      print_help(argv,arg);
   }
   // check algorithm parameters 
   if(arg.w <4) {
     cout << "Windows size must be at least 4\n";
     exit(1);
   }
   if(arg.p<10) {
     cout << "Modulus must be at leas 10\n";
     exit(1);
   } 
}



int main(int argc, char** argv)
{
  // translate command line parameters
  Args arg;
  parseArgs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;  

  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;  
  // init sorted map counting the number of occurrences of each word
  map<uint64_t,word_stats> wordFreq;  
  uint64_t totChar;
  word_int_t s=0,g=0;
  // ------------ parsing input file 
  try {
      totChar = process_file(arg,wordFreq,s,g); 
  }
  catch(const std::bad_alloc&) {
      cout << "Out of memory (parsing phase)... emergency exit\n";
      die("bad alloc exception");
  }
  // first report 
  uint64_t totDWord = wordFreq.size();
  cout << "Total input symbols: " << totChar << endl;
  cout << "Found " << totDWord << " distinct words" <<endl;
  cout << "Parsing took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  // check # distinct words
  if(totDWord>MAX_DISTINCT_WORDS) {
    cerr << "Emergency exit! The number of distinc words (" << totDWord << ")\n";
    cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
    exit(1);
  }
  cout << "Total number of sequences: " << s << "\nDictionary words starting with separators: " << g << endl; //2026
  // -------------- second pass  
  start_wc = time(NULL);
  // create array of dictionary words
  vector<const string *> dictArray;
  dictArray.reserve(totDWord);
  // fill array
  uint64_t sumLen = 0;
  uint64_t totWord = 0;
  for (auto& x: wordFreq) {
    sumLen += x.second.str.size();
    totWord += x.second.occ;
    dictArray.push_back(&x.second.str);
  }
  assert(dictArray.size()==totDWord);
  cout << "Sum of lenghts of dictionary words: " << sumLen << endl; 
  cout << "Total number of words: " << totWord << endl; 
  // sort dictionary
  sort(dictArray.begin(), dictArray.end(),pstringCompare);
  // write plain dictionary and occ file, also compute rank for each hash 
  cout << "Writing plain dictionary and occ file\n";
  writeDictOcc(arg, wordFreq, dictArray);
  dictArray.clear(); // reclaim memory
  cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
    
  FILE *info = open_aux_file(arg.inputFileName.c_str(), "info", "wb");
  assert(s>=g); 
  //write # words in the parse starting with EoString plus a dummy word
  s+=1;
  if(fwrite(&s,sizeof(word_int_t),1,info)!=1) die("info write error"); 
  s-=g;
  //write offset s-g 
  if(fwrite(&s,sizeof(word_int_t),1,info)!=1) die("info write error"); 
  if(fclose(info)!=0) die("Error closing info file"); 

  // remap parse file
  start_wc = time(NULL);
  cout << "Generating remapped parse file together with the map\n";
  remapParse(arg, wordFreq,s);
  cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";        
  return 0;
}

