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


#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include<iostream>
#include<fstream>
#include<stdio.h>
#include "cplex_utils.h"

#include<stdio.h>

#include "stv_distance.h"
#include "sim_stv.h"

#define EPSILON 0.000001

using namespace std;

typedef map<int,double> IDMAP;
typedef vector<IDMAP> IDMAPs;
typedef vector<IDMAPs> IDMAP2d;

struct TimerRecord{
	bool aborted;
	mytimespec start;
	mytimespec end;

	double bestobj;
	double best_dual;
	mytimespec besttime;
};	


ILOMIPINFOCALLBACK6(TimerCallback, TimerRecord&, timer, double, lb, double, ub, 
	const Config&, config, bool, isleaf, ostream&,log)
{
	if(timer.aborted) return;
	GetTime(&timer.end);

	if(hasIncumbent()){
		double inc_obj = getIncumbentObjValue();
		double dual_obj = getBestObjValue();
	
    	if(dual_obj >= ub){
			timer.aborted = true;
			abort();
    	}

		if(inc_obj - dual_obj <= 1.0){
			timer.aborted = true;
			abort();
		}

		if(inc_obj < timer.bestobj-0.99){
			timer.bestobj = inc_obj;
			GetTime(&timer.besttime);
		}
		else if(dual_obj > timer.best_dual + 0.99){
			timer.best_dual = dual_obj;
			GetTime(&timer.besttime);
		}
		else if(!isleaf){
			if(config.tlimit_wsol > 0 && timer.end.seconds - 
				timer.besttime.seconds > config.tlimit_wsol){
				timer.aborted = true;
				abort();
			}
		}
		else{
			if(config.tlimit_leaf > 0 && timer.end.seconds - 
				timer.besttime.seconds > config.tlimit_leaf){
				timer.aborted = true;
				abort();
			}
		}
	}
	else if(!isleaf && config.tlimit_wosol > 0 && timer.end.seconds -
		timer.start.seconds > config.tlimit_wosol){
		abort();
	}
	else if(isleaf && config.tlimit_leaf > 0 && timer.end.seconds -
		timer.start.seconds > config.tlimit_leaf){
		abort();
	}
}



int next_mod(vector<int> &mask, int n){
	int i;
	for(i = 0; i < n && mask[i]; ++i){
		mask[i] = 0;
	}

	if(i < n){
		mask[i] = 1;
		return 1;
	}
	return 0;
}

void DefineOrders(const Node &node, int j, const Ints &mask, 
	Ints2d &possible_orders, const Ints2d &rorder_c,
	const Ints &rorder_a){

	for(int i = 0; i < rorder_c[j].size(); ++i){
		Ints order;
		order.push_back(rorder_c[j][i]);
		possible_orders.push_back(order);
	}

	int size = 1;

	for(int k = j+1; k < rorder_c.size(); ++k){
		if(mask[k] == 1){
			Ints2d additional_orders;
			for(int m = 0; m < possible_orders.size(); ++m){
				Ints &order = possible_orders[m];

				for(int l = 1; l < rorder_c[k].size(); ++l){
					Ints temp = order;
				   	temp.push_back(rorder_c[k][l]);
		 			additional_orders.push_back(temp);
				}
		
				order.push_back(rorder_c[k][0]);
			}
			for(int m = 0; m < additional_orders.size(); ++m){
				possible_orders.push_back(additional_orders[m]);
			}
			++size;
		}	
	}

	if(size < rorder_c.size()){
		Ints2d additional_orders;
		for(int i = 0; i < node.remcand.size(); ++i){
			const int cand = node.remcand[i];
			for(int k = 0; k < possible_orders.size(); ++k){
				Ints order = possible_orders[k];

				if(node.cand_proc[order[order.size()-1]] ==
					node.order_c.size()-1){
					continue;
				}
				order.push_back(cand);
				additional_orders.push_back(order);
			}
		}
		for(int m = 0; m < additional_orders.size(); ++m){
			possible_orders.push_back(additional_orders[m]);
		}
	}
}

int GetRevBallotID(const Node &node, const Ints &prefs, const Ints2d &rorder_c,
	const Ints &rorder_a){
	// First, we have to recreate the equivalence class pref order
	Ints equiv;
	equiv.push_back(prefs[0]);
	if(node.cand_proc[prefs[0]] == -1 ||
		node.cand_proc[prefs[0]] == rorder_c.size()-1){
		return node.ballotmap.find(equiv)->second;
	}

	for(int i = 1; i < prefs.size(); ++i){
		if(equiv.size() == rorder_c.size()){
			break;
		}
		const int cand = prefs[i];
		const int cp = node.cand_proc[cand];

		if(cp == -1){
			equiv.push_back(cand);
			return node.ballotmap.find(equiv)->second;
		}

		bool include = true;
		for(int j = 0; j < equiv.size(); ++j){
			const int ecand = equiv[j];
			const int ep = node.cand_proc[ecand];
			if(ep >= cp){
				include = false;
				break;
			}
		}

		if(include){
			equiv.push_back(cand);
		}

		if(cp ==rorder_c.size()-1){
			break;
		}
	}

	return node.ballotmap.find(equiv)->second;
}




void Preliminaries(const Ballots &ballots, 
	const Candidates &cand, const Config &config, Node &node,
	bool compbounds, const set<int> &elected_o, std::ofstream &log){
	
	//double comp_lower_bound = 0;
	const int ncand = node.order_c.size();

	node.last_elim = ncand - 1;
	for(int j = ncand-1; j >= 0; --j){
		if(node.order_a[j] == 0) break;

		--node.last_elim;
	}

	node.last_q = ncand - 1;
	for(int j = ncand-1; j >= 0; --j){
		if(node.order_a[j] == 1) break;

		--node.last_q;
	}


	for(int i = 0; i < config.ncandidates; ++i){
		node.max_votes[i].resize(ncand, 0);
		int lelim = -1;
		for(int j = 0; j < ncand; ++j){
			if(!node.order_a[j]) lelim = j;

			if(node.order_c[j] == i){
				node.cand_proc[i] = j;
				node.cand_act[i] = node.order_a[j];
				node.cand_equota[i] = lelim+1;
				break;
			}
		}
	}

	return;
}


void CreateEquivalenceClasses_RL(const Ballots &ballots, 
	const Candidates &cand, const Config &config, Node &node,
	const set<int> &elected_o, std::ofstream &log, const Ints2d &rorder_c,
	const Ints &rorder_a){
	try{
		const int ncand = rorder_c.size();
	
		int last_elim = ncand - 1;
		for(int j = ncand-1; j >= 0; --j){
			if(node.order_a[j] == 0) break;

			--last_elim;
		}

		int last_q = ncand - 1;
		for(int j = ncand-1; j >= 0; --j){
			if(rorder_a[j] == 1) break;

			--last_q;
		}

		for(int i = 0; i < config.ncandidates; ++i){
			node.poss_tally[i].clear();
			node.poss_tally[i].resize(ncand);
			node.d_ij_slist[i].clear();
			node.d_ij_slist[i].resize(ncand);
		}

		for(int i = 0; i < config.ncandidates; ++i){
			int lelim = -1;
			for(int j = 0; j < ncand; ++j){
				if(!rorder_a[j]) lelim = j;

				if(find(rorder_c[j].begin(), rorder_c[j].end(), i) !=
					rorder_c[j].end()){
					node.cand_proc[i] = j;
					node.cand_act[i] = node.order_a[j];
					node.cand_equota[i] = lelim+1;
					break;
				}
			}
		}

		if(config.deleteonly){
			for(int i = 0; i < ballots.size(); ++i){
				const Ballot &blt = ballots[i];
				node.rev_ballots.push_back(blt);
				int first = blt.prefs[0];
				for(int l = 0; l < ncand; ++l){
					node.poss_tally[first][l].insert(i);
					if(node.cand_proc[first] == l)
						break;
				}
			}
		}
		else{
			Ints mask(config.ncandidates, 0);

			int cntr = 0;
			while(next_mod(mask, ncand)){
				int j = -1;
				for(int i = 0; i < ncand; ++i){
					if(mask[i]){
						j = i;
						break;
					}
				}

				if(j < 0) continue;
	
				Ints2d possible_orders;
				DefineOrders(node, j, mask, possible_orders, rorder_c, rorder_a);

				for(int k = 0; k < possible_orders.size(); ++k){
					Ballot b;
					b.tag = cntr++;
					b.votes = 0;
					b.prefs = possible_orders[k];

					node.rev_ballots.push_back(b);
					node.ballotmap.insert(pair<vector<int>,int>(b.prefs, b.tag));	
				
					int first = possible_orders[k][0];
					for(int l = 0; l < ncand; ++l){
						node.poss_tally[first][l].insert(b.tag);
						if(node.cand_proc[first] == l)
							break;
					}
				}
			}

			for(int i = 0; i < node.remcand.size(); ++i){
				const int cand = node.remcand[i];
				Ballot b;
				b.tag = cntr++;
				b.votes = 0;
				b.prefs.push_back(cand);

				node.rev_ballots.push_back(b);
				node.ballotmap.insert(pair<vector<int>,int>(b.prefs, b.tag));	

				for(int l = 0; l < ncand; ++l){
					node.poss_tally[cand][l].insert(b.tag);
				}
			}

			for(int b = 0; b < ballots.size(); ++b){
				const Ballot &bt = ballots[b];
				const int rid = GetRevBallotID(node, bt.prefs,rorder_c, rorder_a);

				node.rev_ballots[rid].votes += bt.votes;
				node.bid2newid[b] = rid;
			}
		}
		
		// iterate over each round
		for(int j = 0; j < ncand; ++j){
			// look at who is eliminated/elected in round j
			const Ints cs = rorder_c[j];
			const int a = rorder_a[j];

			for(int m = 0; m < cs.size(); ++m){
				const int c = cs[m];

				// the votes they could distribute to others are
				// those in poss_tally[c][j].
				const SInts &tally = node.poss_tally[c][j];
			 
				for(SInts::const_iterator it = tally.begin();
					it != tally.end(); ++it){	
					// Determine who this ballot may be distributed to
					const int bid = *it;
					const Ints &prefs = node.rev_ballots[bid].prefs;
					for(int k = 0; k < prefs.size(); ++k){
						const int kc = prefs[k];
						// if next preferred candidate has already been
						// eliminated/elected -> skip
						if(node.cand_proc[kc] >= 0 && 
							node.cand_proc[kc] <= j) continue;

						// otherwise, candidate 'kc' could get ballots
						// with id 'bid' in round 'j'
						node.d_ij_slist[kc][j].push_back(bid);
						if(j < ncand - 1){
							// so, candidate 'kc's tally in round 
							// j + 1 may have signature 'bid' in it.

							for(int l = j+1; l < ncand; ++l){
								node.poss_tally[kc][l].insert(bid);
								if(node.cand_proc[kc] == l)
									break;
							}
						}
						// it's possible this ballots could skip 'kc'
						// because 'kc' already has a quota, so we
						// should keep working our way through the 
						// preference list. But, if either candidate 'c'
						// or 'kc' is being eliminated, then we know
						// we can stop here.
						if(!a || node.cand_act[kc] == 0) break;
					}
				}		
			}
		}

		node.dist = max(node.dist, -1.0);
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		throw STVException("Unexpected error in creating eq classes.");
	}
}


// Function for simulating STV election
double distance(const Ballots &ballots, const Candidates &cand, 
	const Config &config, Node &node, double upperbound,
	double tleft, bool lprelax, const set<int> &elected_o, ofstream &log,
	const Ints &origc, const Ints &origa, double &def_ub)
{
	if(origc.size() == node.order_c.size()){
		bool skip = true;
		for(int i = 0; i < origc.size(); ++i){
			if(origc[i] != node.order_c[i] ||
				origa[i] != node.order_a[i]){
				skip = false;
				break;
			}
		}
		if(skip){
			return 0;
		}
	}

	double dist = -1;

	// if order_a[0] = 0 it means candidate 0 was eliminated
	// if order_a[0] = 1 if means candidate 1 was elected.
	try{
		Ints2d rorder_c;
		Ints rorder_a; 
		const int ncand = node.order_c.size();	
		Ints elim;
		for(int i = 0; i < node.order_c.size(); ++i){
			if(node.order_a[i]){
				if(!elim.empty()){
					if(elim.size() > 1){
						Ints elim1;
						for(int j = 0; j < elim.size()-1; ++j)
							elim1.push_back(elim[j]);
						Ints elim2;
						elim2.push_back(elim.back());
						rorder_c.push_back(elim1);	
						rorder_a.push_back(0);	
						rorder_c.push_back(elim2);	
						rorder_a.push_back(0);
					}
					else{
						rorder_c.push_back(elim);
						rorder_a.push_back(0);
					}
					elim.clear();
				}

				Ints st;
				st.push_back(node.order_c[i]);
				rorder_c.push_back(st);
				rorder_a.push_back(1);
			}
			else{
				elim.push_back(node.order_c[i]);
			}
		}

		if(!elim.empty()){
			if(elim.size() > 1){
				Ints elim1;
				for(int j = 0; j < elim.size()-1; ++j)
					elim1.push_back(elim[j]);
				Ints elim2;
				elim2.push_back(elim.back());
				rorder_c.push_back(elim1);	
				rorder_a.push_back(0);	
				rorder_c.push_back(elim2);	
				rorder_a.push_back(0);
			}
			else{
				rorder_c.push_back(elim);
				rorder_a.push_back(0);
			}
			elim.clear();
		}

		CreateEquivalenceClasses_RL(ballots, cand, config, node,
			elected_o, log, rorder_c, rorder_a);

		int rlast_q = 0;
		int rlast_e = 0; 

		int cntr = 0;
		for(int j = 0; j < rorder_c.size(); ++j){
			cntr += rorder_c[j].size();

			if(rorder_a[j]){
				rlast_q = j;
			}
			else{
				rlast_e = j;
			}
		}	


		IloEnv env;
		IloModel cmodel(env);

		// Assume ballots contains all possible rankings, even 
		// if the number of times that ranking was voted for is '0'
		const int sigs = node.rev_ballots.size();

		IloNumVarArray ps(env, sigs);
		IloNumVarArray ms(env, sigs);
		IloNumVarArray ys(env, sigs);

		char varname[500];

		IloExpr balance(env);

		double lb = max(0.0, max(config.glb, node.dist));
		double ub = max(lb, upperbound);
		double quota_lb = config.quota;
		double quota_ub = config.quota;

		if(config.deleteonly){
			quota_lb = 1+ ((config.totalvotes - ub)/(config.nseats+1));
		}

		if(config.addonly){
			quota_ub = 1+ ((config.totalvotes + ub)/(config.nseats+1));
		}

		IloExpr obj(env);
		IloExpr totvotes(env);
		for(int i = 0; i < sigs; ++i){
			const double ns = node.rev_ballots[i].votes;
			// p_s variable: number of ballots modified so that their
			// new signature is 's' 
			sprintf(varname, "vps_%d", i);
			if(config.deleteonly){
			    ps[i] = IloNumVar(env, 0, 0, ILOFLOAT, varname);
			}
			else if(lprelax){
			    ps[i] = IloNumVar(env, 0, ub, ILOFLOAT, varname);
			}
			else{
			    ps[i] = IloNumVar(env, 0, ub, ILOINT, varname);
			}

			// m_s variable: number of ballots whose signature in the
			// original profile is 's', but are modified to something 
			// other than 's' in the new profile
			sprintf(varname, "vms_%d", i);
			if(config.addonly){
			    ms[i] = IloNumVar(env, 0, 0, ILOFLOAT, varname);
			}
			if(lprelax){
			    ms[i] = IloNumVar(env, 0, min(ns, ub), ILOFLOAT, varname);
			}
			else{
			    ms[i] = IloNumVar(env, 0, min(ns, ub), ILOINT, varname);
			}

			// y_s variable: total number of ballots with signature 's'
			// in the new election profile
			sprintf(varname, "vys_%d", i);

			ys[i] = IloNumVar(env, 0, min(ns+ub,config.totalvotes), 
			    ILOFLOAT, varname);

			// n_s = total number of ballots with signature 's' in
			// the original profile.
			// constraint: n_s + p_s - m_s = y_s
			//     rewrite: n_s = y_s + m_s - p_s
			if(!(config.deleteonly || config.addonly)){
				cmodel.add(ns == ys[i] + ms[i] - ps[i]);
				balance += (ps[i] - ms[i]);
				obj += ps[i];
			}

			if(config.deleteonly){
				cmodel.add(ys[i] == ns - ms[i]);
				obj += ms[i];
			}
			else if(config.addonly){
				cmodel.add(ys[i] == ns + ps[i]);
				obj += ps[i];
			}
			totvotes += ys[i];
		}
		cmodel.add(obj <= ub);
		cmodel.add(obj >= lb);

		cmodel.add(IloMinimize(env, obj));
		if(!(config.deleteonly || config.addonly)){
			cmodel.add(balance == 0);
		}

		// Variable to represent quota
		IloNumVar vquota(env, quota_lb, quota_ub, ILOFLOAT, "quota");
		if(config.deleteonly || config.addonly){
			cmodel.add(vquota == 1 + (totvotes/(config.nseats+1)));
		}

		const int ncandr = rorder_c.size();

		NumVarArray2D s_ij(env, config.ncandidates);
		NumVarArray2D nots_ij(env, config.ncandidates);
		NumVarArray2D v_ij(env, config.ncandidates);

		NumVarArray3D d_ij_s(env, config.ncandidates);
		Ints3d set_d_ij_s(config.ncandidates);

		NumVarArray3D y_ij_s(env, config.ncandidates);
		Ints3d set_y_ij_s(config.ncandidates);

		IloNumVarArray w_k(env, ncandr);
		IloNumVarArray tw_k(env, ncandr);

		NumVarArray2D w_k_bins(env, ncandr);

		ExprArray1D twcons(env, ncandr);

		for(int i = 0; i < ncandr; ++i){
			twcons[i] = IloExpr(env);

			if(ncand == config.ncandidates && (rlast_q == i-1 || 
				rlast_e == i-1)) break;

			if(rorder_a[i]){
				sprintf(varname, "vwk_%d", i);
				w_k[i] = IloNumVar(env, 0, 1, ILOFLOAT, varname);

				const int MTV = config.totalvotes; 
					
				sprintf(varname, "vtwk_%d", i);
				tw_k[i] = IloNumVar(env, 0, MTV, ILOFLOAT, varname);

				twcons[i] += tw_k[i];

				if(config.DIV_BILIN > 0){
                	w_k_bins[i] = IloNumVarArray(env, config.DIV_BILIN);
					double left = 0;
                	double gap = 1.0/(double)config.DIV_BILIN;

					IloExpr sum_over_wbins(env);
					IloExpr cless(env);
					IloExpr cmore(env);

                	for(int j = 0; j < config.DIV_BILIN; ++j){
                    	sprintf(varname, "vwk_%d_bin_%d", i, j);
						w_k_bins[i][j] = IloNumVar(env, 0, 1, ILOINT, varname);

						sum_over_wbins += w_k_bins[i][j];
						cmore += left * w_k_bins[i][j];
						cless += (left+gap) * w_k_bins[i][j];

						left += gap;
                	}

					cmodel.add(100*w_k[i] >= 100*cmore);
					cmodel.add(100*w_k[i] <= 100*cless);
					cmodel.add(sum_over_wbins == 1);
				}
			}
		}

		Ints2d s_setvalue(config.ncandidates);

		for(int i = 0; i < config.ncandidates; ++i){
			s_ij[i] = IloNumVarArray(env, ncandr);
			nots_ij[i] = IloNumVarArray(env, ncandr);
			v_ij[i] = IloNumVarArray(env, ncandr);

			y_ij_s[i] = NumVarArray2D(env, ncandr);
			set_y_ij_s[i].resize(ncandr);

			s_setvalue[i].resize(ncandr, -1);

			for(int j = 0; j < ncandr; ++j){
				if(ncand == config.ncandidates && (rlast_q == j-1 || 
					rlast_e == j-1)) break;

				// variable s_{i,j}: binary, whether or not candidate 'i'
				// has a quota at the start of round 'j'

				// co: which candidate is being eliminated/elected in round j
				const int co = rorder_c[j][0]; 
				// whether co is being eliminated/elected.
				const int ca = rorder_a[j];

				if(i == rorder_c[j][0] && ca){
					s_setvalue[i][j] = 1;	
				}
				else if(!ca){
					// candidate being eliminated this round: no candidate
					// can have a quota
					s_setvalue[i][j] = 0;
				}
				else if(node.cand_proc[i] >= 0 && 
					node.cand_proc[i] < j && node.cand_act[i]){
					s_setvalue[i][j] = 1;
				}
				else if(node.cand_act[i] == 0){
					// if candidate 'i' will be eliminated at some point
					// then they will never have a quota: s_{i,j} = 0
					s_setvalue[i][j] = 0;
				}
				else if(node.cand_equota[i] > j){
					// if the round in which candidate 'i' can have a 
					// quota is greater than j, then s_{i,j} = 0
					s_setvalue[i][j] = 0;
				}
				else{
					// otherwise, candidate 'i' might be able to have
					// a quota at the start of round j.
					sprintf(varname, "vs_ij_%d_%d",i,j);
					s_ij[i][j] = IloNumVar(env, 0, 1, ILOINT, varname);
					sprintf(varname, "vns_ij_%d_%d",i,j);
					nots_ij[i][j] = IloNumVar(env, 0, 1, ILOINT, varname);

					cmodel.add(s_ij[i][j] + nots_ij[i][j] == 1);
				}
		
				// variable v_{i,j} to represent the number of votes candidate
				// i has in their tally AT THE START OF round j
				sprintf(varname, "v_ij_%d_%d",i,j);
				
				const double MV = config.totalvotes;

				if(!ca){
					// if the candidate is being eliminated in round 'j' we 
					// know they will have somewhere between 0 and the quota
					// (minus a little bit) number of votes.
					if(config.deleteonly || config.addonly){
						v_ij[i][j] = IloNumVar(env, 0, min(MV,quota_ub), 
							ILOFLOAT, varname);
						cmodel.add(v_ij[i][j] <= vquota - 0.001);
					}
					else{
						v_ij[i][j] = IloNumVar(env, 0, min(MV, 
							(config.quota-0.001)), ILOFLOAT, varname);
					}
				}
				else if(co == i && ca){
					// if we know candidate i is being elected in round j then
					// we know they have somewhere between the quota and the 
					// total number of avail votes in their tally
					if(config.deleteonly || config.addonly){
						v_ij[i][j] = IloNumVar(env, 0, min(MV, 
							config.totalvotes),ILOFLOAT, varname);
						cmodel.add(v_ij[i][j] >= vquota);
					}
					else{
						v_ij[i][j] = IloNumVar(env, config.quota, min(MV, 
							config.totalvotes),ILOFLOAT, varname);
					}
				}
				else{
					// otherwise, we don't know much: could have anywhere
					// between 0 and total number of available votes.
					v_ij[i][j] = IloNumVar(env, 0, MV, ILOFLOAT, varname);
				}

				IloExpr setv(env);

				if(s_setvalue[i][j] == 0){
					if(config.deleteonly || config.addonly){
						cmodel.add(v_ij[i][j] <= vquota - 0.001);
					}
					else{
						cmodel.add(v_ij[i][j] <= min(MV,config.quota-0.001));
					}
				}
				else if(s_setvalue[i][j] == 1){
					if(config.deleteonly || config.addonly){
						cmodel.add(v_ij[i][j] >= vquota);
					}
					else{
						cmodel.add(v_ij[i][j] >= config.quota);
					}
				}
				else{
					if(j > 0 && s_setvalue[i][j-1] == -1){
						cmodel.add(s_ij[i][j] >= s_ij[i][j-1]);
					}
					if(config.deleteonly || config.addonly){
						// Creating and adding constraint: v[i][j] >= s[i][j] Q
						// z = s_ij * vquota
						IloNumVar zsijvq(env, 0, quota_ub, ILOFLOAT);
						cmodel.add(zsijvq >= quota_lb * s_ij[i][j]);
						cmodel.add(zsijvq <= quota_ub * s_ij[i][j]);
						cmodel.add(zsijvq >= vquota-quota_ub*(1 - s_ij[i][j]));
						cmodel.add(zsijvq <= vquota-quota_lb*(1 - s_ij[i][j]));
						cmodel.add(zsijvq <= vquota+(1 - s_ij[i][j])*quota_ub);

						cmodel.add(v_ij[i][j] >= zsijvq); 

						// Creating and adding constraint: v[i][j] < 
						// nots[i][j] Q + s[i][j] * max
						// z = nots[i][j] * Q
						IloNumVar znsij(env, 0, quota_ub, ILOFLOAT);
						cmodel.add(znsij >= quota_lb * nots_ij[i][j]);
						cmodel.add(znsij <= quota_ub * nots_ij[i][j]);
						cmodel.add(znsij >= vquota-quota_ub*(1-nots_ij[i][j]));
						cmodel.add(znsij <= vquota-quota_lb*(1-nots_ij[i][j]));
						cmodel.add(znsij <= vquota+(1-nots_ij[i][j])*quota_ub);

						cmodel.add(v_ij[i][j] <= znsij - 0.001*nots_ij[i][j]
							+ s_ij[i][j] * min(MV, config.totalvotes));
						// Creating and adding constraint: v[i][j] >= s[i][j] Q
						cmodel.add(v_ij[i][j] >= config.quota * s_ij[i][j]);

						// Creating and adding constraint: v[i][j]<nots[i][j]Q+
						//     s[i][j] * max
						cmodel.add(v_ij[i][j] <= nots_ij[i][j] * (config.quota-
							0.001) + s_ij[i][j]*min(MV,config.totalvotes));
					}
					else{
						// Creating and adding constraint: v[i][j] >= s[i][j] Q
						cmodel.add(v_ij[i][j] >= config.quota * s_ij[i][j]);

						// Creating and adding constraint: v[i][j] < nots[i][j] Q +
						//     s[i][j] * max
						cmodel.add(v_ij[i][j] <= nots_ij[i][j] * (config.quota-
							0.001) + s_ij[i][j] * min(MV, config.totalvotes));
					}
				}

				// deal with ballots that could possible be in candidates
				// tally in this round.
				const SInts &ts = node.poss_tally[i][j];
				y_ij_s[i][j] = IloNumVarArray(env, sigs);
				set_y_ij_s[i][j].resize(sigs, 0);

				for(SInts::const_iterator it = ts.begin(); it!=ts.end();++it){
					const int bid = *it;
					const Ballot &bl = node.rev_ballots[bid];
					if(j > 0 && bl.prefs[0] != i){
						// variable y_{i,j,s}: number of votes of signature 's'
						// in candidate 'i's tally at the start of round 'j'
						sprintf(varname, "yijs_%d_%d_%d", i, j, bid);
						y_ij_s[i][j][bid] = IloNumVar(env, 0, min(bl.votes+ub,
							config.totalvotes), ILOFLOAT, varname);

						set_y_ij_s[i][j][bid] = 1;

						// eqv => v_{i,j} - sum_{sigs that could be in i's tally
						//   at the start of round j} y_{i,j,s} = 0
						setv += y_ij_s[i][j][bid];

						// we'll define the constraints 'defining' this
						// variable later.
					}
					else{
						y_ij_s[i][j][bid] = ys[bid];
						set_y_ij_s[i][j][bid] = 1;
						setv += ys[bid];
					}
				}

				cmodel.add(v_ij[i][j] == setv);

				if(co == i) break;
			}
		}
			
		for(int i = 0; i < config.ncandidates; ++i){
			d_ij_s[i] = NumVarArray2D(env, ncandr);
			set_d_ij_s[i].resize(ncandr);

			for(int j = 0; j < ncandr; ++j){
				if(ncand == config.ncandidates && (rlast_q == j-1 || 
					rlast_e == j-1)) break;

				const int co = rorder_c[j][0];
				const int ca = rorder_a[j];

				// let's work out which signatures could possibly be
				// distributed to candidate i in round j.
				if(find(rorder_c[j].begin(),rorder_c[j].end(), i) !=
					rorder_c[j].end()) break;

				SInts ds;
				ds.insert(node.d_ij_slist[i][j].begin(), 
					node.d_ij_slist[i][j].end());

				d_ij_s[i][j] = IloNumVarArray(env, sigs);
				set_d_ij_s[i][j].resize(sigs, 0);
				for(SInts::const_iterator kt = ds.begin(); kt != ds.end();++kt){
					const int k = *kt;
					const int bid = *kt;
					sprintf(varname, "d_%d_%d_%d", i, j, bid);

					// need to work out list of candidates that the ballots
					// could potentially go to over candidate 'i'
					Ints posscands;
					const Ballot &blt = node.rev_ballots[bid];
					for(int l = 0; l < blt.prefs.size(); ++l){
						const int cand = blt.prefs[l];
						if(node.cand_proc[cand] >= 0 && 
							node.cand_proc[cand] <= j){
							continue;
						}

						posscands.push_back(cand);
						if(cand == i){
							break;
						}
					}	
			
					const double M = config.totalvotes;

					// if cand was eliminated in j, then we know that their
					// votes have been distributed to only one other candidate
					if(!ca){
						d_ij_s[i][j][bid] = IloNumVar(env,0,M,ILOFLOAT,varname);
						set_d_ij_s[i][j][bid] = 1;
						IloExpr deqy(env);

						// so d_ij_s[i][j][bid] == y_ij_s[co][j][bid]
						for(Ints::const_iterator et = rorder_c[j].begin();
							et != rorder_c[j].end(); ++et){
							if(set_y_ij_s[*et][j][bid] == 1)
								deqy += y_ij_s[*et][j][bid];
						}
						cmodel.add(deqy == d_ij_s[i][j][bid]);
					}
					else if(posscands.size() > 0){	
						// d_ij_s[i][j][bid] == w_c * y_ij_s[co][j][bid]
						//      * Prod(s_vars for cands in posscands) * (s_var
						//      -1)
						// prod = min(s1, s2, ..., sn);
					
						IloNumVar z;	
						bool isz0 = false;

						// add constraint to define z=y_ij_s[co][j][bid]*prod
						if(posscands.size() == 1 && node.cand_act[i] == 0){
							// prod == 1
							z = y_ij_s[co][j][bid];
						}
						else{
							bool isprod1 = false;
							bool isprod0 = false;
	
							if(posscands.size() == 1){
								if(s_setvalue[i][j] == 0)
									isprod1 = true;
								else if(s_setvalue[i][j] == 1)
									isprod0 = true;
							}

							if(isprod1){
								z = y_ij_s[co][j][bid];
							}
							else if(isprod0){
								isz0 = true;
							}
							else{
								sprintf(varname,"y_ij_sBYprod_%d_%d_%d",i,j,k);
								z = IloNumVar(env, 0, M, ILOFLOAT, varname);

								IloNumVar prod;
								
								if(posscands.size() == 1){
									prod = nots_ij[i][j];
								}
								else{
				   					sprintf(varname, "prod_%d_%d_%d", i, j, k);
									prod = IloNumVar(env,0,1,ILOINT,varname);

									// prod >= sum(posscands) s_ij[cl][j]-(n-1)
									int coef = 0;
									for(int l = 0; l < posscands.size(); ++l){
										const int cl = posscands[l];
										if(l == posscands.size() - 1){
											if(s_setvalue[cl][j] == 0)
												++coef;
										}
										else if(s_setvalue[cl][j] == 1)
											++coef;
									}

									IloExpr setprod(env);
									setprod -= prod;

									for(int l = 0; l < posscands.size(); ++l){
										const int cl = posscands[l];
										if(l == posscands.size() - 1){
											if(s_setvalue[cl][j] == -1)
												setprod += nots_ij[cl][j];
										}
										else if(s_setvalue[cl][j] == -1)
											setprod += s_ij[cl][j];
									}
	
									cmodel.add(setprod <= -coef + 
										(int)posscands.size() - 1);
								}	

								// ****
								// z <= M * prod
								cmodel.add(z <= M * prod);

								// z <= y_ij_s[co][j][bid] 
								cmodel.add(z <= y_ij_s[co][j][bid]);

								// z >= y_ij_s[co][j][bid] - (1 - prod) * M
								cmodel.add(z >= y_ij_s[co][j][bid]-(1-prod)*M);
							}
						}					

						if(isz0){
							set_d_ij_s[i][j][bid] = 0;
						}
						else{
							sprintf(varname, "d_%d_%d_%d", i, j, bid);
							d_ij_s[i][j][bid] = IloNumVar(env,0,M,
								ILOFLOAT,varname);
							set_d_ij_s[i][j][bid] = 1;
							

							// Now deal with (prod of two continuous vars): 
							//     d_ij_s[i][j][bid] == w_c * z
							// w_c = x, z = y
							// x -> 0,1
							// y -> 0, M
                            
							if(config.DIV_BILIN == 0){
								cmodel.add(d_ij_s[i][j][bid] >= 0); 
								cmodel.add(d_ij_s[i][j][bid] >= M*w_k[j] + z - M); 
								cmodel.add(d_ij_s[i][j][bid] <= M*w_k[j]); 
								cmodel.add(d_ij_s[i][j][bid] <= z);
							} 
							else{
                            	// z = z^L + sum_{1..N} \deltay_n
								IloExpr setz(env);

								const double gap = 1.0/(double)config.DIV_BILIN;
                            	double wbval = 0;

								IloExpr bilin_1(env);
								bilin_1 += M*w_k[j];

								IloExpr bilin_2(env);
								bilin_2 += M*w_k[j];

								IloExpr bilin_3(env);
								IloExpr bilin_4(env);

                            	for(int l = 0; l < config.DIV_BILIN; ++l){
                                	sprintf(varname, "dyn_%d_%d_%d_%d", i, j, k, l);
									IloNumVar dyn(env, 0, M, ILOFLOAT, varname);

									cmodel.add(0.1 * dyn <= 0.1 * M * w_k_bins[j][l]);

									double xl = wbval;
									double xr = wbval + gap;

									bilin_1 += xr*dyn - M*xr * w_k_bins[j][l];
									bilin_2 += xl*dyn - M*xl * w_k_bins[j][l];

									bilin_3 += xr * dyn;
									bilin_4 += xl * dyn;
			
									setz += dyn;
									wbval += gap;
								}

								cmodel.add(setz == z);
								cmodel.add(d_ij_s[i][j][bid] >= bilin_1);
								cmodel.add(d_ij_s[i][j][bid] <= bilin_2);
								cmodel.add(100*d_ij_s[i][j][bid] <= 100*bilin_3);
								cmodel.add(100*d_ij_s[i][j][bid] >= 100*bilin_4);
							}
							twcons[j] -= z;
						}
					}
				}
			}
		}

		for(int i = 0; i < config.ncandidates; ++i){
			for(int j = 0; j < ncandr; ++j){ 
				if(ncand == config.ncandidates && (rlast_q == j-1 || 
					rlast_e == j-1)) break;
				const SInts &ts = node.poss_tally[i][j];
				for(SInts::const_iterator it = ts.begin();
					it != ts.end(); ++it){
					const int bid = *it;
					if(j > 0){
						// y_ij_s[i][j][bid] == y_ij_s[i][j-1][bid] +
						//     d_ij_s[i][j-1][bid] 
						
						IloExpr yeq(env);
						yeq += y_ij_s[i][j][bid];

						if(set_y_ij_s[i][j-1][bid] == 1)
							yeq -= y_ij_s[i][j-1][bid];

						if(set_d_ij_s[i][j-1][bid] == 1){
							yeq -= d_ij_s[i][j-1][bid];
						}

						cmodel.add(yeq == 0);
					}
				}

				if(find(rorder_c[j].begin(), rorder_c[j].end(), i) !=
					rorder_c[j].end()) break;
			}
		}

		for(int j = 0; j < ncandr; ++j){
			if(ncand == config.ncandidates && 
				rlast_q == j-1) break;
			if(ncand == config.ncandidates && 
				rlast_e == j-1) break;

			if(rorder_a[j]){
				const int c = rorder_c[j][0];
				cmodel.add(twcons[j] == 0);
				
				// in round j, candidate is elected with quota
				for(int k = 0; k < config.ncandidates; ++k){
					if(node.cand_proc[k] >= 0 &&  
						node.cand_proc[k] <= j){
						continue;
					}

					//const int ck = node.order_c[k];
					// add constraint v[c][j] >= v[k][j]
					cmodel.add(v_ij[c][j] >= v_ij[k][j]);
				}

				// correct expression is w[cand] * transferrablev[cand][j] =
				//      v[cand] - Q
				if(j != ncandr - 1){
					const int MTV = config.totalvotes;

                    // piecewise linear relaxation for w[cand] * tv
                    // tv ranges from 0 to MTV
                    const double gap = 1.0/(double)config.DIV_BILIN;
                    double wbval = gap;

					sprintf(varname,"xlin_wt_%d_%d",c,j);
					IloNumVar xlin(env, 0, MTV, ILOFLOAT, varname);
					// x = w[cand] [0 1]
					// y = tv [0 MTV]	

					if(config.DIV_BILIN == 0){
						cmodel.add(xlin >= 0);
						cmodel.add(xlin >= MTV * w_k[j] + tw_k[j] - MTV);
						cmodel.add(xlin <= MTV * w_k[j]);
						cmodel.add(xlin <= tw_k[j]);
					}
					else{
					/*	sprintf(varname,"xlin_wt_%d_%d",c,j);
						IloNumVar xlin(env, 0, MTV, ILOFLOAT, varname);*/

            	        // tv = tv^L + sum_{1..N} dtn
						IloExpr set_tv(env);

						IloExpr bilin_1(env);
						bilin_1 += MTV * w_k[j];

						IloExpr bilin_2(env);
						bilin_2 += MTV * w_k[j];

						IloExpr bilin_3(env);
						IloExpr bilin_4(env);

						for(int l = 0; l < config.DIV_BILIN; ++l){
                    	    sprintf(varname, "dtn_%d_%d_%d", c, j, l);
							IloNumVar dtn(env, 0, MTV, ILOFLOAT, varname);
							cmodel.add(0.1*dtn <= 0.1*MTV * w_k_bins[j][l]);

							set_tv += dtn;

							bilin_1 += wbval*dtn - MTV*wbval*w_k_bins[j][l];
							bilin_2  += (wbval-gap)*dtn - MTV*(wbval-gap)*w_k_bins[j][l];
							bilin_3 += wbval * dtn;
							bilin_4 += (wbval-gap) * dtn;

							wbval += gap;
						}
						cmodel.add(tw_k[j] == set_tv);
						cmodel.add(xlin >= bilin_1);
						cmodel.add(xlin <= bilin_2);
						cmodel.add(100*xlin <= 100*bilin_3);
						cmodel.add(100*xlin >= 100*bilin_4);
					}
				
					sprintf(varname, "tvgs_%d", c);
					IloNumVar tvgs(env, 0, 1, ILOINT, varname);

					// if transferablevotes<=(total votes-quota)+M*tvgs 
					if(config.deleteonly || config.addonly){
						cmodel.add(tw_k[j] <= (v_ij[c][j] - vquota) + 
							MTV * tvgs);

						// v[cand] - Q -> tvgs*(v[cand]-Q) + (1-tvgs)*tv

                 	   // linearise tvgs * v[cand]
						// A -> v[cand] [0 config.totalvotes]
						// x > tvgs  [0 1]
						sprintf(varname, "tvgs_z1_%d", c);
						IloNumVar z1(env,0,config.totalvotes,ILOFLOAT,varname);

						cmodel.add(z1 <= config.totalvotes * tvgs);
						cmodel.add(z1 <= v_ij[c][j]);
						cmodel.add(z1 >= v_ij[c][j]-(1-tvgs)*config.totalvotes);
						cmodel.add(z1 >= 0);

						// tvgs*tv
						// A -> tv [0 config.totalvotes]
						// x -> tvgs
						sprintf(varname, "tvgs_z2_%d", c);
						IloNumVar z2(env,0,config.totalvotes, ILOFLOAT,varname);
						cmodel.add(z2 <= config.totalvotes * tvgs);
						cmodel.add(z2 <= tw_k[j]);
						cmodel.add(z2 >= tw_k[j]-(1-tvgs) * config.totalvotes);

						// tvgs * vquota
						IloNumVar ztv(env, 0, quota_ub, ILOFLOAT);
						cmodel.add(ztv >= quota_lb * tvgs);
						cmodel.add(ztv <= quota_ub * tvgs);
						cmodel.add(ztv >= vquota - quota_ub*(1 - tvgs));
						cmodel.add(ztv <= vquota - quota_lb*(1 - tvgs));
						cmodel.add(ztv <= vquota + (1 - tvgs)*quota_ub);

						cmodel.add(xlin == z1 - ztv + tw_k[j] - z2);
					}
					else{
						cmodel.add(tw_k[j] <= (v_ij[c][j] - config.quota) + 
							MTV * tvgs);

						// linearise tvgs * v[cand]
						// A -> v[cand] [0 config.totalvotes]
						// x > tvgs  [0 1]
						sprintf(varname, "tvgs_z1_%d", c);
						IloNumVar z1(env,0,config.totalvotes,ILOFLOAT,varname);

						cmodel.add(z1 <= config.totalvotes * tvgs);
						cmodel.add(z1 <= v_ij[c][j]);
						cmodel.add(z1 >= v_ij[c][j]-(1-tvgs)*config.totalvotes);
						cmodel.add(z1 >= 0);

						// tvgs*tv
						// A -> tv [0 config.totalvotes]
						// x -> tvgs
						sprintf(varname, "tvgs_z2_%d", c);
						IloNumVar z2(env,0,config.totalvotes,ILOFLOAT,varname);
						cmodel.add(z2 <= config.totalvotes * tvgs);
						cmodel.add(z2 <= tw_k[j]);
						cmodel.add(z2 >= tw_k[j]-(1-tvgs)*config.totalvotes);

						cmodel.add(xlin==z1-config.quota*tvgs+tw_k[j]-z2);
					}
				}
			}
			else{
				// candidate(s) is/are eliminated
				if(rorder_c[j].size() == 1){
					for(int k = 0; k < config.ncandidates; ++k){
						if(node.cand_proc[k] >= 0 &&  
							node.cand_proc[k] <= j){
							continue;
						}

						const int co = rorder_c[j][0];

					// add constraint v[c][j] <= v[k][j]
						cmodel.add(v_ij[co][j] <= v_ij[k][j]);
					}
				}
				else{
					if(j < rorder_c.size() - 1){
						if(!(node.remcand.empty()  && (rlast_q == j ||
							rlast_e == j))){
							for(int k = 0; k < rorder_c[j].size(); ++k){
								const int co = rorder_c[j][k];
								for(int l = 0; l < config.ncandidates; ++l){
									if(node.cand_proc[l] >= 0 &&  
										node.cand_proc[l] <= j){
										continue;
									}
									// add constraint v[co][j] <= v[cl][j+1]
									cmodel.add(v_ij[co][j] <= v_ij[l][j+1]);
								}
							}
						}
					}
				}
			}
		}

		TimerRecord timer;
		timer.aborted = false;
		GetTime(&timer.start);
		GetTime(&timer.besttime);
		timer.bestobj = numeric_limits<double>::max();
		timer.best_dual = 0;
		
		IloCplex cplex(cmodel);
		if(node.remcand.size() <= 8){
			cplex.setOut(log);
		}
		else{
			cplex.setOut(env.getNullStream());
		}
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, tleft);
		cplex.setParam(IloCplex::MIPEmphasis, 1);
		cplex.setParam(IloCplex::Threads, 1);

		cplex.use(TimerCallback(env, timer, lb, ub, config, node.remcand.size() == 0, log));
		bool result = cplex.solve();

		if(cplex.getCplexStatus() == IloCplex::Infeasible){
			dist = -1;
		}
		else if(result){
 			dist = cplex.getObjValue();
			if(cplex.getBestObjValue() >= 0 &&
				cplex.getBestObjValue() <= ub){
				dist = max(lb, cplex.getBestObjValue());
			}
		} 
		else{
			if(cplex.getBestObjValue() >= 0 &&
				cplex.getBestObjValue() <= ub){
				dist = max(lb, cplex.getBestObjValue());
			}
			else{
				dist = lb;
			}
		}

/*		if(result && node.remcand.size() == 0 && dist < def_ub){
			Ballots newballots(ballots);
			
			for(int i = 0; i < sigs; ++i){
				const Ballot &bl = node.rev_ballots[i];
				int p = (!config.deleteonly) ? cplex.getValue(ps[i]) : 0;
				int m = (!config.addonly) ? cplex.getValue(ms[i]) : 0;

				if(p > 0){
					log << "(";
					for(int j = 0; j < bl.prefs.size(); ++j){
						log << bl.prefs[j] << " ";
					}
					log << ") + " << p << endl;
					newballots.push_back(bl);
				}
				else if(m > 0){
					log << "(";
					for(int j = 0; j < bl.prefs.size(); ++j){
						log << bl.prefs[j] << " ";
					}
					log << ") - " << m << endl;

					// remove m first preference votes from bl.prefs[0]
					int cand = bl.prefs[0];
					int rem = m;
					for(int j = 0; j < newballots.size(); ++j){
						Ballot &nb = newballots[j];
						if(nb.prefs[0] == cand){
							if(nb.votes >= rem){
								nb.votes -= rem;
								rem = 0;
							}
							else{
								rem -= nb.votes;
								nb.votes = 0;
							}
						}
						if(rem == 0)
							break;
					}
				}
			}
			Doubles votecounts;
			for(int i = 0; i < newballots.size(); ++i){
				newballots[i].tag = i;
				votecounts.push_back(newballots[i].votes);
			}

			Candidates candcopy(cand);
			Ints simorderc;
			Ints simordera;
			double mindiff = 0;
			bool alternate = false;

			SimSTV(newballots, votecounts, candcopy, config, simorderc,
				simordera, false, mindiff);  

			for(int i = 0; i < candcopy.size(); ++i){
				const Candidate &c = candcopy[i];
				if(c.seat != -1 && elected_o.find(i) == elected_o.end()){
					// We have found an alternate outcome
					alternate = true;
					break;
				}	
			}

			if(!alternate){
				if(!(config.deleteonly || config.addonly)){
					double extra = WEUB(newballots, votecounts, candcopy,
						config, elected_o, origc, origa);
					if(dist + extra <= def_ub){
						// Definite upper bound of dist + extra
						def_ub = dist + extra;
						log<<"Found definite upper bound of "<<def_ub << endl;
					}
					else{
						log << "Simulation finds candidate upper bound of: "
							<< dist+extra << " (not used)" << endl;
					}
				}
				else{
					log<<"Manipulation does not actually change result."<<endl;
				}
			}
			else{
				def_ub = dist;
				log << "Found definite upper bound of " << def_ub << endl;
			}
		}*/

		node.ClearEqClassData();
		cplex.end();
		cmodel.end();
		env.end();

	}
	catch(IloCplex::Exception e){
		log << e.getMessage() << endl;
		throw STVException("CPLEX error in STV distance calc");
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		throw STVException("Unexpected error in STV distance calculation.");
	}


	return dist;
}


