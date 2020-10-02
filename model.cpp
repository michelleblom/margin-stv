#include<fstream>
#include<algorithm>
#include<iostream>
#include<sstream>
#include<time.h>
#include "model.h"

using namespace std;
typedef boost::char_separator<char> boostcharsep;

const double GAP = 0;

void GetTime(struct mytimespec *t)
{
	#ifndef _WIN32
	struct timeval tv;
	gettimeofday(&tv, NULL);
	t->seconds = tv.tv_sec + (tv.tv_usec/1E6);
	#else
	t->seconds = clock()/CLOCKS_PER_SEC;
	#endif
}

template <typename T>
T ToType(const std::string &s) 
{ 
	try
	{
		return boost::lexical_cast<T>(s); 
	}
	catch(...)
	{
		throw STVException("Lexical cast problem");
	}
}

void Split(const string &line, const boostcharsep &sep, vector<string> &r)
{
	try
	{
		boost::tokenizer<boostcharsep> tokens(line, sep);
		r.assign(tokens.begin(), tokens.end());

		for(int i = 0; i < r.size(); ++i)
		{
			boost::algorithm::trim(r[i]);
		}
	}
	catch(exception &e)
	{
		stringstream ss;
		ss << e.what();
		throw STVException(ss.str());
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		throw STVException("Error in Split Function.");
	}
}

double ReadMargin(const char *path){
	double margin = numeric_limits<int>::max();
	try{
		ifstream infile(path);
		string line;
		getline(infile, line);
		margin = ToType<double>(line);
	}
	catch(exception &e)
	{
		throw e;
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		cout << "Unexpected error reading in margin." << endl;
		return false;
	}

	return margin;
}

bool ReadConfig(const char *path, Config &config, Candidates &cand)
{
	try
	{
		boostcharsep sp(",");
	
		ifstream infile(path);
		string line;
		getline(infile, line);

		vector<string> columns;
		Split(line, sp, columns);

		config.ncandidates = ToType<int>(columns[0]);
		config.nseats = ToType<int>(columns[1]);
		//config.glb = ToType<double>(columns[3]);
		//config.gub = ToType<double>(columns[4]);

		if(columns.size() == 6){
			config.quota = ToType<double>(columns[5]);
		}

		getline(infile,line);
		columns.clear();
		Split(line, sp, columns);

		if(columns.size() != config.ncandidates)
		{
			infile.close();
			cout << "Incorrect specification of candidates in config."
				<< endl;
			return false;
		}

		infile.close();

		for(int i = 0; i < columns.size(); ++i)
		{
			int id = ToType<int>(columns[i]);
			Candidate c;
			c.index = i;
			c.id = id;

			cand.push_back(c);

			config.id2index.insert(pair<int,int>(id,i));
		}
	}
	catch(exception &e)
	{
		throw e;
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		cout << "Unexpected error reading in config file." << endl;
		return false;
	}

	return true;

}


// Assume input format:
// (first_id, second_id, third_id, ...) : #appears
bool ReadBallots(const char *path, Ballots &ballots, Candidates &candidates,
	Config &config)
{
	try
	{
		cout << "Reading Ballots" << endl;
		boostcharsep sp(",():");
		
		int cntr = ballots.size();

		ifstream infile(path);
		string line;
		while(getline(infile, line))
		{
			vector<string> columns;
			Split(line, sp, columns);

			Ballot b;
			b.tag = cntr;
			b.votes = ToType<double>(columns.back());

			for(int i = 0; i < columns.size()-1; ++i)
			{
				if(columns[i] == "") continue;
				int ccode = ToType<int>(columns[i]);
				int index = config.id2index.find(ccode)->second;
					
				if(find(b.prefs.begin(),
					b.prefs.end(), index) != b.prefs.end())
				{
					continue;
				}

				b.prefs.push_back(index);
			}

			if(b.prefs.empty()) continue;

			Candidate &cand = candidates[b.prefs.front()];
			bool found_ballot = false;
			int bid = -1;
			for(int j = 0; j < cand.ballots.size(); ++j){
				const Ballot &bf = ballots[cand.ballots[j]];
				if(bf.prefs == b.prefs){
					bid = bf.tag;
					found_ballot = true;
					break;
				}
			}
			if(found_ballot){
				ballots[bid].votes += b.votes;
				cand.sum_votes += b.votes;
				config.totalvotes += b.votes;
				continue;
			}
			else{
				cand.sum_votes += b.votes;
				cand.ballots.push_back(b.tag);

				config.totalvotes += b.votes;
	
				ballots.push_back(b);
				++cntr;
			}
		}

		cout << "Finished reading ballots" << endl;
		infile.close();
	}
	catch(exception &e)
	{
		throw e;
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		cout << "Unexpected error reading in ballots." << endl;
		return false;
	}

	return true;
}
