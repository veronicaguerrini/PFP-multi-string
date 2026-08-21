# PFP-multi-string
Tool to build the multi-string BWT and optionally the Document Array of string collections using the Prefix-Free Parsing preprocessing (PFP).

### Download and Compile

```sh
git clone https://github.com/veronicaguerrini/PFP-multi-string
cd PFP-multi-string
make
```
### Usage

To build the multi-string BWT for file *test.fasta* just type
```sh
./PFP-multi-string test.fasta
```

The input file cannot contain the characters 0, 1 or 2, which are used internally by the algorithm, and the character '$', which is internally used as end-of-string marker.
The character '$' is the end-of-string symbol in the output multi-string BWT.
The PFP parameters can be modified using the command line options `-w` for the sliding window size and `-p` for the hash modulus.

Using option `-D` it is possible to compute the Document Array at the same time of the multi-string BWT.
```sh
./PFP-multi-string -D test.txt
```
### Authors

* Veronica Guerrini
* Felipe A. Louza
* Giovanni Manzini
* Giovanna Rosone
