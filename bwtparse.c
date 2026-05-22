#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>
#include "gsa/gsacak.h"
#include "utils.h" 

// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT32_MAX -1)
typedef uint32_t word_int_t;

// Input: file.parse, [file.sai], [file.dai]
// Output: file.ilist, [file.bwsai], [file.bwdai]

// Compute the SA of a parsing and its associated BWT

// If the -s option is used, also remap the sa values in the .sai 
// file (that uses IBYTES bytes per entry) producing the .bwsai file

// If the -d option is used, also remap the da values in the .dai 
// file (that uses 4 bytes per entry) producing the .bwdai file

typedef uint_t sa_index_t;

// -------------------------------------------------------------
// struct containing command line parameters and other globals
typedef struct {
   char *basename;
   bool SAinfo;
   bool DAinfo;
   int th;    // number of segments for the last and sa files
} Args;


// read the parse file, add a 0 EOS symbol and return a pointer 
// to a new allocate uint32_t array containing it. Store size in *tsize
// note *tsize is the number of elements in the parsing, but we add a 0
// symbol at the end so the returned array has *tsize+1 elements 
static uint32_t *read_parse(char *basename, long *tsize) 
{  
  FILE *parse = open_aux_file(basename,EXTPARSE,"rb");
  // get file size
  if(fseek(parse,0,SEEK_END)!=0) die("parse fseek");
  long nn = ftell(parse);
  // check input file is OK
  if(nn%4!=0) {
    printf("Invalid input file: size not multiple of 4\n");
    exit(1);
  }
  #if !M64
  // if in 32 bit mode, the number of words is at most 2^31-2
  if(nn/4 > 0x7FFFFFFE) {
    printf("Input containing more than 2^31-2 phrases!\n");
    printf("Please use 64 bit version\n");
    exit(1);
  }
  #else
  // if in 64 bit mode, the number of words is at most 2^32-2 (for now)
  if(nn/4 > 0xFFFFFFFEu) {
    printf("Input containing more than 2^32-2 phrases!\n");
    printf("This is currently a hard limit\n");
    exit(1);
  }
  #endif
  printf("Parse file contains %ld words\n",nn/4);
  long n = nn/4;
  // ------ allocate and read text file, len is n+1 for the EOS 
  uint32_t *Text = malloc((n+1)*sizeof(*Text));
  if(Text==NULL) die("malloc failed (Text)");
  rewind(parse);

  // read the array in one shot
  assert(sizeof(*Text)==4);
  size_t s = fread(Text,sizeof(*Text),n,parse);
  if(s!=n) {
    char *msg=NULL;
    int e= asprintf(&msg,"read parse error: %zu vs %ld\n", s,n); 
    (void) e; die(msg);
  }
  if(fclose(parse)!=0) die("parse file close");
  Text[n]=0; // sacak needs a 0 eos 
  *tsize= n;
  return Text;
}  

static uint32_t *read_map(char *basename, long *tsize) 
{  
  FILE *map_file = open_aux_file(basename,"map","rb");
  // get file size
  if(fseek(map_file,0,SEEK_END)!=0) die("map fseek");
  long nn = ftell(map_file);
  // check input file is OK
  if(nn%4!=0) {
    printf("Invalid input file: size not multiple of 4\n");
    exit(1);
  }
  #if !M64
  // if in 32 bit mode, the number of words is at most 2^31-2
  if(nn/4 > 0x7FFFFFFE) {
    printf("Input containing more than 2^31-2 phrases!\n");
    printf("Please use 64 bit version\n");
    exit(1);
  }
  #else
  // if in 64 bit mode, the number of words is at most 2^32-2 (for now)
  if(nn/4 > 0xFFFFFFFEu) {
    printf("Input containing more than 2^32-2 phrases!\n");
    printf("This is currently a hard limit\n");
    exit(1);
  }
  #endif
  printf("map file contains %ld words\n",nn/4);
  long n = nn/4;
  // ------ allocate and read map file
  uint32_t *MAP = malloc(n*sizeof(*MAP));
  if(MAP==NULL) die("malloc failed (MAP)");
  rewind(map_file);
  // read the array in one shot
  assert(sizeof(*MAP)==4);
  size_t s = fread(MAP,sizeof(*MAP),n,map_file);
  if(s!=n) {
    char *msg=NULL;
    int e= asprintf(&msg,"read map error: %zu vs %ld\n", s,n); 
    (void) e; die(msg);
  }
  if(fclose(map_file)!=0) die("map file close"); 
  *tsize= n;
  return MAP;
} 
static void print_help(char *name)
{
  printf("Usage: %s <basename> [options]\n\n", name);
  puts("Compute the BWT of basename.parse and store its inverted list occurrence");
  puts("  Options:");
  puts("\t-h  \tshow help and exit");
  puts("\t-s  \tpermute also basename."EXTSAI " file");
  puts("\t-d  \tpermute also basename.dai file");
  exit(1);
}

static void parseArgs(int argc, char** argv, Args *arg ) {
  extern int optind, opterr, optopt;
  extern char *optarg;  
  int c;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

  arg->SAinfo = false;
  arg->DAinfo = false; //04_2026
  arg->th = 0;
  while ((c = getopt( argc, argv, "sdht:") ) != -1) {
    switch(c) {
      case 's':
      arg->SAinfo = true; break;
      case 'd':
      arg->DAinfo = true; break;
      case 'h':
         print_help(argv[0]); exit(1);
      case '?':
      puts("Unknown option. Use -h for help.");
      exit(1);
    }
  }
  // read base name as the only non-option parameter 
  if (argc!=optind+1)
    print_help(argv[0]);
  arg->basename = argv[optind];
}

static sa_index_t *compute_SA(uint32_t *Text, long n, long k) 
{
  sa_index_t *SA = malloc(n*sizeof(*SA));
  if(SA==NULL) die("malloc failed  (SA)");
  printf("Computing SA of size %ld over an alphabet of size %ld\n",n,k);
  int depth = sacak_int(Text, SA, n, k);
  if(depth>=0)
    printf("SA computed with depth: %d\n", depth);
  else
    die("Error computing the SA"); 
  return SA;
}

  
static uint8_t *load_sa_info(Args *arg, long n)
{  
  // maybe sa info was not required 
  if(arg->SAinfo==false) return NULL;
  // open .sai file for reading
  mFile *fin = mopen_aux_file(arg->basename,EXTSAI,arg->th);
  // allocate and load the sa info array
  uint8_t *sai = malloc(n*IBYTES);
  if(sai==NULL) die("malloc failed (SA INFO)"); 
  size_t s = mfread(sai,IBYTES,n,fin);
  if(s!=n) die("sa info read");
  if(mfclose(fin)!=0) die("sa info file close");
  return sai;
}

static uint32_t *load_da_info(Args *arg, long n)
{  
  // maybe da info was not required 
  if(arg->DAinfo==false) return NULL;
  // open .dai file for reading
  mFile *fin = mopen_aux_file(arg->basename,"dai",arg->th);
  // allocate and load the sa info array
  uint32_t *dai = malloc(n*sizeof(uint32_t));
  if(dai==NULL) die("malloc failed (DA INFO)"); 
  size_t s = mfread(dai,sizeof(uint32_t),n,fin);
  if(s!=n) die("da info read");
  if(mfclose(fin)!=0) die("da info file close");
  return dai;
}

static FILE *open_sa_out(Args *arg)
{
  if(arg->SAinfo==false) return NULL;
  return open_aux_file(arg->basename,EXTBWSAI,"wb"); 
}

//04_2026
static FILE *open_da_out(Args *arg)
{
  if(arg->DAinfo==false) return NULL;
  return open_aux_file(arg->basename,"bwdai","wb"); 
}


int main(int argc, char *argv[])
{
  uint32_t *Text; // array of parsing symbols
  uint32_t *Map; // array of "initial" words  
  long n;         // length of Text[] (not including final 0 symbol)
  long check;    // length of MAP[]
  size_t s;
  Args arg;

  // read arguments 
  parseArgs(argc,argv,&arg);
  // start measuring wall clock time 
  time_t start_wc = time(NULL);
  // read parse file
  Text = read_parse(arg.basename,&n);
  
  // ------- compute largest input symbol (ie alphabet size-1)
  long k=0;
  for(long i=0;i<n;i++) {
    if(Text[i]>k) k = Text[i];
  }
  // -------- alloc and compute SA of the parse
  sa_index_t *SA = compute_SA(Text,n+1,k+1);
  
  // load sa/da info file, if requested
  uint8_t *sa_info = load_sa_info(&arg,n-1);
  uint32_t *da_info = load_da_info(&arg,n-1); //04_2026
 
  FILE *sa_out = open_sa_out(&arg);
  FILE *da_out = open_da_out(&arg);

  // transform SA->BWT inplace and write remapped last array, and possibly sainfo
  sa_index_t *BWTsa = SA; // BWT overlapping SA
  assert(n>1);
  // first BWT symbol
  assert(SA[0]==n);

  for(long i=1;i<=n;i++) {
    if(SA[i]==0) {  
      assert(i==1);  // Text[0]=$abc... is the second lex word 
      BWTsa[i-1] = 0;   // eos in BWT, there is no phrase in D corresponding to this symbol so we write dummy values
      if(arg.SAinfo) write_myint(0,sa_out); // dummy end of word position, it is never used an 0 does not appear elsewhere in sa_out
      if(arg.DAinfo) fwrite(&SA[i],sizeof(uint32_t),1,da_out); // dummy end of word position 
    }
    else {
      if(arg.SAinfo) get_and_write_myint(sa_info,n-1,SA[i]-1,sa_out); // ending position of BWT symbol in original text
      if(arg.DAinfo) fwrite(da_info+SA[i]-1,sizeof(uint32_t),1,da_out);
      BWTsa[i-1] = Text[SA[i]-1];
    }
  }
  if(arg.SAinfo) {
    if(fclose(sa_out)!=0) die("sa_out close");
    free(sa_info);
  } 
  if(arg.DAinfo){
    if(fclose(da_out)!=0) die("da_out close");
    free(da_info);
  } 
  //Upload s and offset value
  word_int_t s_offset[2];
  FILE *info = open_aux_file(arg.basename,"info","rb");
  s = fread(s_offset,sizeof(word_int_t),2,info);
  if(s!=2) die("not enough info data!");
  printf("info file contains %d and %d\n",s_offset[0],s_offset[1]);
  fclose(info);
  //Upload the .map file
  Map = read_map(arg.basename,&check);

  assert(check+1==s_offset[0]);
  //Update BWT values
  uint32_t *BWT = Text;
  for(long i=0;i<n;i++){
    if(BWTsa[i]>s_offset[0])
      BWT[i] = BWTsa[i]-s_offset[1];
    else if(BWTsa[i]>0)
      BWT[i] = Map[BWTsa[i]-1];
    else
      BWT[i] = BWTsa[i];
  }
  free(Map);

  //Reduce k accordingly
  k-=s_offset[1]; 
  //------------/ 04_2026

  // read # of occ of each char from file .occ
  uint32_t *occ = malloc((k+1)*sizeof(*occ)); // extra space for the only occ of 0
  if(occ==NULL) die("malloc failed (OCC)");
  FILE *occin = open_aux_file(arg.basename,"occ","rb");
  s = fread(occ+1,sizeof(*occ), k,occin);
  if(s!=k) die("not enough occ data!");
  occ[0] = 1; // we know there is somewhere a 0 BWT entry 
  fclose(occin);
  // create F vector
  uint32_t *F = malloc((k+1)*sizeof(*F));
  if(F==NULL) die("malloc failed (F)");
  // init F[] using occ[]
  F[0] = 0;
  for(int i=1;i<=k;i++)
    F[i]=F[i-1]+occ[i-1];
  assert(F[k]+occ[k]==n);
  puts("---- computing inverted list ----");
  // ----- compute inverse list overwriting SA
  uint32_t *IList = (uint32_t *) SA+1;
  for(long i=0;i<n;i++) {
    IList[F[BWT[i]]++] = i;
    occ[BWT[i]]--;
  }
  // ---check
  assert(IList[0]==0); // EOF is in BWT[0] since P[0] = $xxx is the smallest word and appears once
  assert(BWT[IList[0]]==0);
  for(long i=0;i<=k;i++) 
    assert(occ[i]==0);
  // ---save Ilist   
  FILE *ilist = open_aux_file(arg.basename,EXTILIST,"wb");
  s = fwrite(IList,sizeof(*IList),n,ilist);
  if(s!=n) die("Ilist write");
  fclose(ilist);
  printf("---- %ld ilist positions written (%ld bytes) ----\n",n+1,(n+1)*4l);
  // deallocate
  free(F);
  free(occ);
  free(SA);
  free(Text);
  printf("==== Elapsed time: %.0lf wall clock seconds\n", difftime(time(NULL),start_wc));  
  return 0;
}
