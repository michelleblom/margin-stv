/*
    Copyright (C) 2017-2020  Michelle Blom

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


#ifndef _MODEL_H
#define _MODEL_H

#include<vector>
#include<set>
#include<exception>
#include<boost/lexical_cast.hpp>
#include<boost/filesystem.hpp>
#include<boost/tokenizer.hpp>
#include<boost/algorithm/string.hpp>
#include<map>

#ifndef _WIN32
#include<sys/time.h>
#endif


typedef std::vector<int> Ints;
typedef std::set<int> SInts;
typedef std::vector<Ints> Ints2d;
typedef std::vector<SInts> SInts2d;
typedef std::vector<Ints2d> Ints3d;
typedef std::vector<SInts2d> SInts3d;

typedef std::vector<double> Doubles;
typedef std::vector<Doubles> Doubles2d;

extern const double GAP;

struct mytimespec
{
	double seconds;
};

void GetTime(struct mytimespec* t);


struct Config
{
	int ncandidates;
	int nseats;
	double totalvotes;
	double quota;
	double glb;
	double gub;

	double tlimit_wsol;
	double tlimit_wosol;
	double tlimit_leaf;

	int threads;
	int subthreads;

	int DIV_BILIN;
	int FIX_UNTIL;

	bool deleteonly;
	bool addonly;

	std::map<int,int> id2index;

	Config() : ncandidates(0), nseats(0), totalvotes(0),
        	quota(0), glb(0), gub(-1), tlimit_wsol(-1),
		tlimit_wosol(-1), tlimit_leaf(-1), threads(1), subthreads(1),
		DIV_BILIN(1), FIX_UNTIL(0), deleteonly(false), addonly(false) {}
};


class STVException
{
	private:
		const std::string message;

	public:
		STVException(const STVException &me) : message(me.what()){}
		STVException(const std::string &str) : message(str) {}
		const std::string& what() const { return message; }
};

struct IntDouble
{
	int id;
	double weight;

	IntDouble(int _id, double _weight) : 
		id(_id), weight(_weight) {}
};

typedef std::vector<IntDouble> IDS;
typedef std::vector<IDS> IDS2d;
typedef std::list<IntDouble> L_IDS;

struct Candidate
{
	int id;
	int index;
	double sum_votes;

	Ints ballots;

	// For simulation
	double sim_votes;
	IDS bweights;

	int standing;
	int seat;
	int surplus;

	double max_votes;

	Candidate() : id(0), index(0), sum_votes(0), sim_votes(0), 
		standing(1), seat(-1), surplus(0), max_votes(0) {}
};

typedef std::vector<Candidate> Candidates;

struct Ballot
{
	int tag;
	double votes;
	Ints prefs;
};

typedef std::vector<Ballot> Ballots;

template <typename T>
T ToType(const std::string &s);

typedef std::map<std::vector<int>,int> I2Map;

struct Node{
	double dist;
	double dist1;

	Ints order_c;
	Ints order_a;
	Ints remcand;
	std::set<int> elected;

	int seatsleft;

	Ballots rev_ballots;
	I2Map ballotmap;
	Ints bid2newid;

	Ints cand_proc;
	Ints cand_act;
	Ints cand_equota;

	SInts3d poss_tally;
	Ints3d d_ij_slist;

	int last_elim;
	int last_q;

	Doubles2d max_votes;

	SInts2d c_clumped_order;
	Ints a_clumped_order; 

	Node(int ncand, int nballots) {
		bid2newid.resize(nballots, -1);
		max_votes.resize(ncand);
		poss_tally.resize(ncand);
		d_ij_slist.resize(ncand);
		last_elim = -1;
		last_q = -1;
		dist = -1;
		seatsleft = 0;
		dist1 = 0;

		cand_proc.resize(ncand, -1);
		cand_act.resize(ncand, -1);
		cand_equota.resize(ncand, -1);
	}

	void Reset(){
		last_elim = -1;
		last_q = -1;
		dist = -1;
		seatsleft = 0;
		dist1 = 0;

		Ballots().swap(rev_ballots);
		I2Map().swap(ballotmap);
		SInts().swap(elected);
		Ints().swap(remcand);

		rev_ballots.clear();
		ballotmap.clear();
		elected.clear();
		remcand.clear();

		for(int i = 0; i < cand_proc.size(); ++i){
			cand_proc[i] = -1;
			cand_act[i] = -1;
			cand_equota[i] = -1;

			Doubles().swap(max_votes[i]);
			SInts2d().swap(poss_tally[i]);
			Ints2d().swap(d_ij_slist[i]);

			max_votes[i].clear();
			poss_tally[i].clear();
			d_ij_slist[i].clear();
		}

		for(int i = 0; i < bid2newid.size(); ++i){
			bid2newid[i] = -1;
		}
	}

	void ClearEqClassData(){
		Ballots().swap(rev_ballots);
		I2Map().swap(ballotmap);

		rev_ballots.clear();
		ballotmap.clear();
		for(int i = 0; i < cand_proc.size(); ++i){
			cand_proc[i] = -1;
			cand_act[i] = -1;
			cand_equota[i] = -1;

			Doubles().swap(max_votes[i]);
			SInts2d().swap(poss_tally[i]);
			Ints2d().swap(d_ij_slist[i]);

			max_votes[i].clear();
			poss_tally[i].clear();
			d_ij_slist[i].clear();
		}


	}
};


bool ReadBallots(const char *path, Ballots &ballots,
	Candidates &candidates, Config &config);

bool ReadConfig(const char *path, Config &config,
	Candidates &cand);

double ReadMargin(const char *path);

#endif
