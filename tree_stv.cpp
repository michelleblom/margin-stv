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


#include<set>
#include<vector>
#include<list>
#include<iostream>
#include<fstream>
#include<cmath>
#include<boost/thread/thread.hpp>
#include<sstream>

#include "tree_stv.h"
#include "stv_distance.h"

using namespace std;

typedef list<Node> Fringe;
typedef vector<Node> Nodes;
typedef vector<Fringe> Fringes;

void do_join(boost::thread& t){
   	t.join();
}
	
void join_all(std::vector<boost::thread>& v){
   	std::for_each(v.begin(),v.end(),do_join);
}


void InsertIntoFringe(const Node &n, Fringe &fringe){
	Fringe::iterator it = fringe.begin();
	for( ; it != fringe.end(); ++it){
		if(n.dist < it->dist ||(n.dist <= it->dist && 
			n.dist1 < it->dist1)){
			break;
		}
	}
	
	fringe.insert(it, n);
}

Node& PopNode(Fringe &fringe){
	return *(fringe.begin());
}

void PrintNode(const Node &n, ostream &log){
	for(int i = 0; i < n.order_c.size(); ++i){
		log << n.order_c[i] << "--";
		if(n.order_a[i]){
			log << "q ";
		}
		else {
			log << "e ";
		}
	}

	log << " with distance " << n.dist << " ";
}

void PrintFringe(const Fringe &fringe, ostream &log){
	log << "--------------------------------" << endl;
	log << "Current state of priority queue: " << endl;
	log << "    Node ";
	PrintNode(fringe.front(), log);
	log << endl;
	log << "    .... " << endl;
	log << "    Node ";
	PrintNode(fringe.back(), log);
	log << endl;
	
/*	for(Fringe::const_iterator it = fringe.begin();
		it != fringe.end(); ++it){
		log << "    Node ";
		PrintNode(*it, log);
		log << endl;	
	}*/
	log << "--------------------------------" << endl;
}

// check if node2 is subsumed by node1
bool SameNode(const Node &node1, const Node &node2){
	if(node1.order_c.size() != node2.order_c.size())
		return false;

	for(int i = 0; i < node1.order_c.size(); ++i){
		if(node1.order_c[i] != node2.order_c[i]){
			return false;
		}

		if(node1.order_a[i] != node2.order_a[i]){
			return false;
		}	
	}
	return true;
}

bool Subsumes(const Config &config, const Node &node1, const Node &node2){
	if(node1.dist != node2.dist)
		return false;

	if(node2.c_clumped_order.size() > node1.c_clumped_order.size()){
		return false;
	}

	for(int i = 0; i < node2.c_clumped_order.size(); ++i){
		if(node1.a_clumped_order[i] != node2.a_clumped_order[i])
			return false;

		if(node1.c_clumped_order[i] != node2.c_clumped_order[i])
			return false;
	}

	return true;
}

void ComputeClumpedOrder(Node &node){
	node.c_clumped_order.clear();
	node.a_clumped_order.clear();

	Ints elim;
	for(int i = 0; i < node.order_c.size(); ++i){
		if(node.order_a[i]){
			if(!elim.empty()){
				SInts elim1;
				elim1.insert(elim.begin(), elim.end());
				node.c_clumped_order.push_back(elim1);	
				node.a_clumped_order.push_back(0);	
				elim.clear();
			}

			SInts st;
			st.insert(node.order_c[i]);
			node.c_clumped_order.push_back(st);
			node.a_clumped_order.push_back(1);
		}
		else{
			elim.push_back(node.order_c[i]);
		}
	}

	if(!elim.empty()){
		SInts elim1;
		elim1.insert(elim.begin(), elim.end());
		node.c_clumped_order.push_back(elim1);	
		node.a_clumped_order.push_back(0);	
		elim.clear();
	}
}

void GetChildren(const Node &n, Nodes &children, const set<int> &elected_o){

	if(n.seatsleft == n.remcand.size()){
		Node newn(n.cand_proc.size(), n.bid2newid.size());
		newn.order_c = n.order_c;
		newn.order_a = n.order_a;
		newn.elected = n.elected;

		for(int i = 0; i < n.seatsleft; ++i){
			newn.order_c.push_back(n.remcand[i]);
			newn.order_a.push_back(1);
			newn.elected.insert(n.remcand[i]);
		}

		newn.seatsleft = 0;
		newn.dist = n.dist; 
		newn.dist1 = n.dist1; 
		if(newn.elected != elected_o){
			children.push_back(newn);
		}

		return;
	}

	if(n.seatsleft == 1 && n.remcand.size() == 2){
		const int fc = n.remcand[0];
		const int sc = n.remcand[1];

		for(int o = 0; o <= 1; ++o){
			Node newn(n.cand_proc.size(), n.bid2newid.size());
			newn.order_c = n.order_c;
			newn.order_a = n.order_a;
			newn.elected = n.elected;
			newn.seatsleft = 0;
			newn.dist = n.dist;
			newn.dist1 = n.dist1;

			newn.order_c.push_back(fc);
			newn.order_a.push_back(o);

			newn.order_c.push_back(sc);
			newn.order_a.push_back(1-o);	
			
			if(o == 1){
				newn.elected.insert(fc);
			}
			else{
				newn.elected.insert(sc);
			}

			if(newn.elected != elected_o){
				children.push_back(newn);			
			}
		}

		for(int o = 0; o <= 1; ++o){
			Node newn(n.cand_proc.size(), n.bid2newid.size());
			newn.order_c = n.order_c;
			newn.order_a = n.order_a;
			newn.elected = n.elected;
			newn.seatsleft = 0;
			newn.dist = n.dist;
			newn.dist1 = n.dist1;

			newn.order_c.push_back(sc);
			newn.order_a.push_back(o);

			newn.order_c.push_back(fc);
			newn.order_a.push_back(1-o);	
			
			if(o == 1){
				newn.elected.insert(sc);
			}
			else{
				newn.elected.insert(fc);
			}

			if(newn.elected != elected_o){
				children.push_back(newn);			
			}
		}

	}
	else if(n.remcand.size() == 2 && n.seatsleft == 0){
		const int fc = n.remcand[0];
		const int sc = n.remcand[1];

		Node newn(n.cand_proc.size(), n.bid2newid.size());
		newn.order_c = n.order_c;
		newn.order_a = n.order_a;
		newn.elected = n.elected;
		newn.seatsleft = n.seatsleft;
		newn.dist = n.dist;
		newn.dist1 = n.dist1;

		newn.order_c.push_back(fc);
		newn.order_a.push_back(0);
		newn.order_c.push_back(sc);
		newn.order_a.push_back(0);	
		children.push_back(newn);		

		newn.order_c[newn.order_c.size()-2] = sc;	
		newn.order_c[newn.order_c.size()-1] = fc;	
		children.push_back(newn);		
	}
	else{
		for(int j = 0; j < n.remcand.size(); ++j){
			const int cand = n.remcand[j];

			// candidate can be elected or eliminated
			for(int o = 0; o <= 1; ++o){
				Node newn(n.cand_proc.size(), n.bid2newid.size());
				newn.order_c = n.order_c;
				newn.order_a = n.order_a;
				newn.elected = n.elected;
				newn.seatsleft = n.seatsleft;
				newn.dist = n.dist;
				newn.dist1 = n.dist1;

				for(int k = 0; k < n.remcand.size(); ++k){
					if(j != k){
						newn.remcand.push_back(n.remcand[k]);
					}
				}

				newn.order_c.push_back(cand);
				newn.order_a.push_back(o);	
			
				if(o == 1){
					newn.elected.insert(cand);
					--newn.seatsleft;

					if(newn.seatsleft == 0){
						for(int k = 0; k < newn.remcand.size(); ++k){
							newn.order_c.push_back(newn.remcand[k]);
							newn.order_a.push_back(0);
						}
						newn.remcand.clear();
					}
				}

				if(!(newn.seatsleft == 0 && newn.elected == elected_o)){
					children.push_back(newn);
				}
			}
		}
	}
}

void PruneFringe(Fringe &fringe, double ubound, ostream &log){
	for(Fringe::iterator it = fringe.begin(); it != fringe.end(); ){
		if(it->dist >= ubound){
			log << "Pruning node ";
			PrintNode(*it, log);
			log << endl;

			if(fringe.size() == 1){
				fringe.clear();
				return;
			}

			fringe.erase(it++);
		}
		else{
			++it;		
		}
	}
}


struct ThreadInputData {
	Ints order_c, order_a;
	set<int> elected_o;
	const Fringe *fringe;
};

struct ThreadOutputData {
	double upperbound;
	double bestdfupperbound;
	double timelimit;
	bool compbounds;
	Fringe fringe;
	Ints best_order_c;
	Ints best_order_a;
	int dtcounter;
};

struct SubThreadInputData{
	set<int> elected_o;
	double upperbound;
	double timelimit;
	Ints order_c;
	Ints order_a;
};

void SolveDistanceTo(Node &n, const Ballots &ballots,
	const Candidates &cands, const Config &config, 
	const SubThreadInputData &sidata,
	const string &logf, double &dtvalue, double &defub){

	const Ints &order_c = sidata.order_c;
	const Ints &order_a = sidata.order_a;

	ofstream log;
	log.open(logf.c_str(), ofstream::app);

	log << "Solving DT for ";
	PrintNode(n, log);
	log << endl;

	dtvalue = distance(ballots, cands, config,
		n, sidata.upperbound,  sidata.timelimit, 
		false, sidata.elected_o, log, order_c, order_a, defub);
	
	log << "Solved DT for ";
	PrintNode(n, log);
	log << " value " << dtvalue << endl;

	log.close();
}



void ExpandNode(const Node &node, const Ballots &ballots, 
	const Candidates &cands, const Config &config, 
	const ThreadInputData &idata, ThreadOutputData &tdata, 
	const string &logf){

	const set<int> &elected_o = idata.elected_o;

	double &upperbound = tdata.upperbound;
	double timelimit = tdata.timelimit;
	bool compbounds = tdata.compbounds;
    
	Fringe &fringe = tdata.fringe;
	Ints &best_order_c = tdata.best_order_c; 
	Ints &best_order_a = tdata.best_order_a;

	ofstream log;
	log.open(logf.c_str(), ofstream::app);

	log << "Expanding node: ";
	PrintNode(node, log);
	log << endl; 

	Nodes children;
	GetChildren(node, children, elected_o);

	int ca = 0;
	int cprocessed = 0;
	while(cprocessed < children.size()){
		boost::thread **tasks = new boost::thread*[config.subthreads];
		Doubles dtvalues(config.subthreads, -1);
		Doubles dfvalues(config.subthreads, tdata.bestdfupperbound);
		Ints corr_child(config.subthreads, -1);
		vector<string> logfilenames(config.subthreads, "");

		SubThreadInputData sidata;
		sidata.elected_o = idata.elected_o;
		sidata.timelimit = timelimit;
		sidata.upperbound = upperbound;
		sidata.order_c = idata.order_c;
		sidata.order_a = idata.order_a;

		int tcntr = 0;
		for(int j = cprocessed; j < children.size(); ++j){
			Node &n = children[j];

			Preliminaries(ballots, cands, config, n,
				compbounds, elected_o, log);

			if(n.dist >= 0 && n.dist >= upperbound){
				++cprocessed;
				continue;
			}

			stringstream ss;
			ss << tcntr << "_" << logf;
			
			logfilenames[tcntr] = ss.str();
			tasks[tcntr] = new boost::thread(
				SolveDistanceTo,
				boost::ref(n),
				boost::ref(ballots),
				boost::ref(cands),
				boost::ref(config),
				boost::ref(sidata),
				boost::ref(logfilenames[tcntr]),
				boost::ref(dtvalues[tcntr]),
				boost::ref(dfvalues[tcntr])
			);

			corr_child[tcntr] = j;
			++tcntr;
			++cprocessed;

			if(tcntr == config.subthreads)
				break;
		}

		for(int j = 0; j < tcntr; ++j)
			tasks[j]->join();

		
		for(int j = 0; j < tcntr; ++j){
			Node &n = children[corr_child[j]];
			tdata.dtcounter += 1;
		
			n.dist1 = dtvalues[j];
			if(n.dist1 >= 0) 
				n.dist = max(n.dist1, n.dist);
			else
				n.dist = n.dist1;

			tdata.bestdfupperbound=min(tdata.bestdfupperbound,dfvalues[j]);

			if(dfvalues[j] < upperbound){
				log << "New definite upper bound " <<
					dfvalues[j] << endl;
				upperbound = dfvalues[j];
			}

			if(n.dist >= 0 && n.dist < upperbound &&
				n.remcand.empty()){
				upperbound = n.dist;

				log << "Leaf node ";
				PrintNode(n, log);
				log << " found " << endl;
				log << "New upper bound: "<<upperbound << endl;

				best_order_c = n.order_c; 
				best_order_a = n.order_a;

				// Remove all nodes from queue with lower bound 
				// greater than or equal to upperbound
				PruneFringe(fringe, upperbound, log);	
				continue;
			}

			if(n.dist >= 0 && n.dist < upperbound-0.01){
				fringe.insert(fringe.end(), n);
				++ca;
			}
		}
	}
		
	if(ca == 0){
		log << "No children, node pruned" << endl;	
	}

	log.close();
}

bool RunTreeSTV(const Ballots &ballots, const Candidates &cands,
	const Config &config, const Ints &order_c, const Ints &order_a, 
	double upperbound, double timelimit, bool compbounds, const char *logf)
{
	mytimespec start;
	GetTime(&start);

	Fringe fringe;

	int dtcounter = 0;
	
	set<int> elected_o;
	for(int j = 0; j < order_c.size(); ++j){
		if(order_a[j]){
			elected_o.insert(order_c[j]);
		}
	}

	ofstream log(logf);

	double bestdfupperbound = upperbound;

	if(config.FIX_UNTIL > 0){
		Node newn(config.ncandidates, ballots.size());
		newn.seatsleft = config.nseats;
		Ints temp_order_c(config.FIX_UNTIL);
		Ints temp_order_a(config.FIX_UNTIL);
		Ints remcand;

		newn.dist = 0;
		newn.dist1 = 0;
		Preliminaries(ballots, cands, config, newn,
			compbounds, elected_o, log);

		for(int i = 0; i < config.FIX_UNTIL; ++i){
			temp_order_c[i] = order_c[i];	
			temp_order_a[i] = order_a[i];
			if(order_a[i] == 1){
				--newn.seatsleft;
				newn.elected.insert(order_c[i]);
			}	
		}

		for(int j = config.FIX_UNTIL; j < config.ncandidates; ++j){
			remcand.push_back(order_c[j]);
		}

		newn.order_c = temp_order_c;
		newn.order_a = temp_order_a;
		newn.remcand = remcand;

		InsertIntoFringe(newn, fringe);
	}
	else{
		// BUILD FRINGE
		for(int i = 0; i < config.ncandidates; ++i){
			// candidate can either be elected or eliminated
			Ints temp_order_c(1);
			temp_order_c[0] = i;
			Ints remcand;
			for(int j = 0; j < config.ncandidates; ++j){
				if(i != j) remcand.push_back(j);	
			}

			for(int o = 0; o <= 1; ++o){
				Ints temp_order_a(1);
				temp_order_a[0] = o;

				mytimespec tnow;
				GetTime(&tnow);
				double tleft = timelimit - (tnow.seconds - start.seconds);
				if(tleft < 10){
					log << "RunTreeSTV timed out" << endl;
					log.close();
					return false;
				}

				Node newn(config.ncandidates, ballots.size());
				newn.order_c = temp_order_c;
				newn.order_a = temp_order_a;
				newn.remcand = remcand;
				if(o == 1) newn.elected.insert(i);
				newn.seatsleft = config.nseats - o;
			
				Preliminaries(ballots, cands, config, newn,
					compbounds, elected_o, log);

				if(newn.dist >= upperbound)
				{
					continue;
				}	

				
				double dfub = bestdfupperbound;
				newn.dist1 = distance(ballots,  cands, config,
					newn, upperbound, tleft, false, 
					elected_o, log, order_c, order_a, dfub);

				if(newn.remcand.empty()){
					if(dfub < bestdfupperbound){
						log << "New definite upper bound" << dfub << endl;
						bestdfupperbound = dfub;
					}
				}

				++dtcounter;

				if(newn.dist1 >= 0) 
					newn.dist = max(newn.dist1, newn.dist);
				else
					newn.dist = newn.dist1;

				if(newn.dist >= 0 && newn.dist < upperbound){
					InsertIntoFringe(newn, fringe);
				}
			}	 
		}
	}

	Ints best_order_c;
	Ints best_order_a;

	double curr_ubound = upperbound;
	double lowerbound = 0;
	while(fringe.size() > 0){
		PrintFringe(fringe, log);

		log << "CURRENT UPPER BOUND = " << curr_ubound << endl;
		bool all_leaves = true;	
		double blower = curr_ubound;
		for(Fringe::const_iterator it = fringe.begin();
			it != fringe.end(); ++it){
			blower = min(blower, it->dist);
			if(!it->remcand.empty()){
				all_leaves = false;
			}
		}
		log << "BEST LOWER BOUND = " << blower << endl;
		lowerbound = blower;

		mytimespec tnow;
		GetTime(&tnow);
		log << "TOTAL TIME USED SO FAR: " << tnow.seconds -
				start.seconds << endl;
		log << "DistanceTo's solved: " << dtcounter << endl;

		if(all_leaves){
			log << "Fringe contains only leaves we couldn't handle" << endl;
			log << "LB " << blower << " UB " << curr_ubound << endl;
			break;
		}

		if(ceil(blower) == ceil(curr_ubound)){
			log << "Lower and Upper bound are close enough." << endl;
			log << "LB " << blower << " UB " << curr_ubound <<endl;
			break;
		}

		// Expand first config.threads nodes
		double tleft = timelimit - (tnow.seconds - start.seconds);
		if(tleft < 10){
			log << "RunTreeSTV timed out" << endl;
			log.close();
			return false;
		}
		const int threads = min((int)fringe.size(), config.threads);
		vector<ThreadOutputData> tdatas(threads);
		vector<string> logfiles;

		ThreadInputData idata;
		idata.order_a = order_a;
		idata.order_c = order_c;
		idata.elected_o = elected_o;
		idata.fringe = &fringe;

		for(int i = 0; i < threads; ++i){
			stringstream ss;
			ss << i << "_" << logf;

			logfiles.push_back(ss.str());
			tdatas[i].upperbound = curr_ubound;
			tdatas[i].bestdfupperbound = bestdfupperbound;
			tdatas[i].dtcounter = 0;
			tdatas[i].timelimit = tleft;
			tdatas[i].compbounds = compbounds;
			tdatas[i].best_order_c = best_order_c;
			tdatas[i].best_order_a = best_order_a;
		}	

		boost::thread **tasks = new boost::thread*[threads];

		Fringe::iterator bg = fringe.begin();
		for(int j = 0; j < threads; ++j){
			tasks[j] = new boost::thread(
				ExpandNode,
				boost::ref(*bg),
				boost::ref(ballots),
				boost::ref(cands),
				boost::ref(config),
				boost::ref(idata),
				boost::ref(tdatas[j]),
				boost::ref(logfiles[j])
			);
			++bg;
		}

		for(int j = 0; j < threads; ++j)
			tasks[j]->join();

		for(int j = 0; j < threads; ++j){
			const ThreadOutputData &td = tdatas[j];
			dtcounter += td.dtcounter;
			if(td.upperbound < curr_ubound){
				curr_ubound = td.upperbound;
				best_order_c = td.best_order_c;
				best_order_a = td.best_order_a;
			}

			if(td.bestdfupperbound < bestdfupperbound){
				bestdfupperbound = td.bestdfupperbound;
			}

			fringe.erase(fringe.begin());
		}

		PruneFringe(fringe, curr_ubound, log);

		for(int j = 0; j < threads; ++j){
			ThreadOutputData &td = tdatas[j];
			Fringe &newnodes = td.fringe;
			for(Fringe::iterator it = newnodes.begin();
				it != newnodes.end(); ++it){
				bool ignore = false;
				ComputeClumpedOrder(*it);
				for(Fringe::const_iterator itt = fringe.begin();
					itt != fringe.end(); ++itt){
					if(Subsumes(config, *itt, *it)){
						ignore = true;
						break;
					}
				}

				if(ignore){
					log << "Node subsumed: ";
					PrintNode(*it, log);
					log << endl;
					continue;	
				}
				if(it->dist < curr_ubound){
					log << "Adding node to fringe: ";
					PrintNode(*it, log);
					log << endl;

					InsertIntoFringe(*it, fringe);
				}
			}
		}

		log << "Size of fringe: " << fringe.size() << endl;
		for(int i = 0; i < threads; ++i){
			delete tasks[i];
		}
		delete[] tasks;
	}

	lowerbound = curr_ubound;
	for(Fringe::const_iterator it = fringe.begin();
		it != fringe.end(); ++it){
		lowerbound = min(lowerbound, it->dist);
	}

	mytimespec tnow;
	GetTime(&tnow);
	log << "TOTAL TIME USED SO FAR: " << tnow.seconds -
		start.seconds << endl;

	if(!best_order_c.empty()){
		log << "====================================" <<endl;
		log << "Minimal manipulation: " << curr_ubound << endl;
		log << "Manipulated order: ";
		for(int i = 0; i < best_order_c.size(); ++i){
			log << cands[best_order_c[i]].id << "--";
			if(best_order_a[i])
				log << "q ";
			else
				log << "e "; 
		}
		log << endl;
	}
	else{
		log << "All nodes pruned " << endl;
	}

	log << "Distance calls: " << dtcounter << endl;
	log << "CURRENT UPPER BOUND: " << curr_ubound << endl;
	log << "CURRENT LOWER BOUND: " << lowerbound << endl;
	log << "Best definite upper bound: " << bestdfupperbound << endl;
	log << "====================================" <<endl;
	log.close();
	return true;
}
