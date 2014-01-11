/****************************************************************************************
 * This program retrives a sequence from a fasta file.
 * Version 1.00
 * 1.01 (13/Jan/2010): added the  [-s <seq-id>] [-i <seq-list>] -f <seq-file> option for 
 * reading the required sequence ids from file 
 * 
 * Written by Itai Sharon (itai.sharon@gmail.com), 12/Nov/2010
 ****************************************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <set>
#include <map>
#include <string>

using namespace std;
typedef pair<size_t, size_t> Segment;


// Just in case the compiler is a little bit old fashioned...
//typedef enum {false, true} bool;

/****************************************************************************************/
void usage(const char* prog_name)
{
	fprintf(stderr, "Usage: %s <seq-id>[,start,end] <seq-file>\n", prog_name);
	fprintf(stderr, "       or\n");
	fprintf(stderr, "       %s [-s <seq-id>] [-i <seq-list>] -f <seq-file>\n", prog_name);
}

/****************************************************************************************/
int main(int argc, const char* argv[]) 
{
	set<string>	seq_ids;
	FILE*		fp = NULL;
	const char*	infile = NULL;
	bool		include_match = true;
	map<string, Segment>	segments;


	/*** 1. Read command line parameters ***/
	if(argc == 3) {
		infile = argv[2];
		int start, end;
		char str[1024];
		if(sscanf(argv[1], "%s.%d,%d", str, &start, &end) == 3) {
			segments[string(str)] = Segment(start, end);
		}
		seq_ids.insert(argv[1]);
	}
	else if(argc >= 5) {
		// 3rd parameter must be -f, 4th is the input fasta file
		if(strcmp(argv[3], "-f")) {
			usage(argv[0]);
			return -1;
		}
		infile = argv[4];

		// If this is a single sequence (-s) - just read it
		if(!strcmp(argv[1], "-s")) {
			seq_ids.insert(argv[2]);
		}
		// Otherwise - if the ids are read from a file - read them
		else if(!strcmp(argv[1], "-i")) {
			FILE* fp = fopen(argv[2], "r");
			if(!fp) {
				fprintf(stderr, "Failed to open %s for reading\n\n", argv[2]);
				return -1;
			}
			char	line[4096];
			while(fscanf(fp, "%s", line) == 1) {
				seq_ids.insert(line);
			}
			fclose(fp);
		}
		// Otherwise - this is an error
		else {
			usage(argv[0]);
			return(-1);
		}
		if(argc == 6) {
			if(!strcmp(argv[5], "-v")) {
				include_match = false;
			}
			else {
				usage(argv[0]);
				return -1;
			}
		}
		else if(argc != 5) {
			usage(argv[0]);
			return -1;
		}
	}
	else {
		usage(argv[0]);
		return(-1);
	}


	fp = fopen(infile, "r");
	if(!fp) {
		fprintf(stderr, "Failed to open %s for reading\n\n", infile);
		return -1;
	}

	/*** 2. Read sequences, collect whichever is relevant ***/
	char	line[4096];
	bool	out = false;
	size_t pos = 1;
	while(fgets(line, 4096, fp)) {
		if(*line == '>') {
			char *p=line+1;
			while(!isspace(*p))
				p++;
			char ch = *p;
			*p = 0;
			set<string>::iterator it = seq_ids.find(line+1); 
			if(it != seq_ids.end()) {
//				fprintf(stderr, ".");
				out = include_match;
				seq_ids.erase(it);
			}
			else if((out == include_match) && (seq_ids.size() == 0) && include_match)
				break;
			else
				out = !include_match;
			*p = ch;
			if(out)
				printf("%s", line);
		}
		else if(out) {
			for(int i=0; line[i]!=0; i++) {
				if((line[i] >= 'a') && (line[i] <= 'z'))
					line[i] -= ('a'-'A');
			}
			printf("%s", line);
		}
	}
	fclose(fp);
	fprintf(stderr, "\n");

	if(seq_ids.size() > 0) {
		fprintf(stderr, "\nUnable to find the following sequences:\n");
		for(set<string>::iterator it = seq_ids.begin(); it != seq_ids.end(); it++)
			fprintf(stderr, "%s\n", it->c_str()); 
		
	}
	return 0;
}
