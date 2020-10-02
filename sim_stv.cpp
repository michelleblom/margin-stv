#include "sim_stv.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>

using namespace std;

// Function prototypes

void EliminateCandidate(int toe, Candidates &cand, 
	const Ballots &ballots, const Doubles &votecounts,
	const Config &config,bool log);

void InsertCandidateSurplus(const Candidate &c, double surplus, 
	L_IDS &surpluses);

int NextCandidate(const Ballot &b, int id, const Candidates &cand);

void DistributeSurplus(Candidate &c, Candidates &cand, double surplus,
	const Ballots &ballots, const Doubles &votecounts, 
	const Config &config,bool log);

double WEUB(const Ballots &ballots, const Doubles &votecounts,
	Candidates &cand, const Config &config, const std::set<int> &elected,
	const Ints &order_c, const Ints &order_a){
	double totvotes = 0;

	double weub = -1;

	Doubles2d votes_per_round;
	Doubles leastvotes_per_round(cand.size(), -1);

	for(Candidates::iterator it = cand.begin(); it != cand.end(); ++it){
		it->sim_votes = 0;
		it->bweights.clear();
		it->standing = 1;
		it->seat = -1;
		it->surplus = 0;

		for(Ints::const_iterator li = it->ballots.begin();
			li != it->ballots.end(); ++li){
			it->bweights.push_back(IntDouble(*li, 1));
			it->sim_votes += votecounts[*li];
		}
		totvotes += it->sim_votes;

		Doubles votes(cand.size(), 0);
		votes[0] = it->sim_votes;
		votes_per_round.push_back(votes);

		if(leastvotes_per_round[0] == -1){
			leastvotes_per_round[0] = it->sim_votes;
		}
		else{
			leastvotes_per_round[0] = min(it->sim_votes,
				leastvotes_per_round[0]);
		}
	}

	// Step 1: Determine quota
	int quota = (int)(1.0 + (totvotes/(double)(config.nseats+1.0))); 

	L_IDS surpluses;		

	int currseat = 0;
	int counter = 0;

	while(currseat < config.nseats){
		int standing = 0;

		// if a candidate has a quota, add them to the list of 
		// candidates with a surplus
		for(int i = 0; i < cand.size(); ++i){
			Candidate &c = cand[i];

			votes_per_round[i][counter] = c.sim_votes;
			if(c.standing) 
				++standing;
			else
				continue;

			if(c.surplus)
				continue;

			if(c.standing && (c.sim_votes >= quota)){
				InsertCandidateSurplus(c, max(0.0,c.sim_votes-quota), 
					surpluses);
				c.surplus = 1;
			}
		}

		if(standing == config.nseats - currseat){
			return weub;
		}					
		if(surpluses.empty()){
			// eliminate candidate with fewest votes.
			// distribute votes at their value.
			double leastvotes = numeric_limits<double>::max();
			int toeliminate = -1;
			for(int i = 0; i < cand.size(); ++i){
				if(cand[i].standing && cand[i].sim_votes < leastvotes){
					leastvotes = cand[i].sim_votes;
					toeliminate = i;
				}
			}

			leastvotes_per_round[counter] = leastvotes;

			if(!config.addonly){
				for(set<int>::const_iterator it = elected.begin();
					it != elected.end(); ++it){
					if(cand[*it].standing){
						double change = ceil(cand[*it].sim_votes-leastvotes+1);
						if(change < cand[*it].sum_votes){
							for(int j = 0; j < counter; ++j){
								if(votes_per_round[*it][j] - change <
									leastvotes_per_round[j]){
									break;
								}
							}

							double weub_next = -1;

							if(weub == -1) 
								weub_next = change;
							else
								weub_next = min(weub, change);

							if(!config.deleteonly && weub_next > change/2.0){
								double halfless=cand[*it].sim_votes-(change/2.0);

								bool usehalf = true;
								bool givenhalf = false;
								for(int i = 0; i < cand.size(); ++i){
									if( i == *it || !cand[i].standing) continue;
									if(cand[i].standing && cand[i].sim_votes 
										>= halfless){
										continue;
									}
                                
									if(givenhalf){
										usehalf = false;
										break;
									}

									if(cand[i].standing && cand[i].sim_votes + 
										(change/2.0) < halfless){
										usehalf = false;
										break;
									}

									givenhalf = true;
								}
								if(usehalf){
									weub_next = change/2.0;
								}
							}

							if(config.deleteonly && weub_next != weub && weub_next != -1){
								// Check whether robbing *it of weub_next votes 
								// does result in a different outcome
								Ballots newballots(ballots);
			
								int rem = ceil(weub_next);
								// remove m first preference votes from bl.prefs[0]
								for(int j = 0; j < newballots.size(); ++j){
									Ballot &nb = newballots[j];
									if(nb.prefs[0] == *it){
										if(nb.votes >= rem){
											nb.votes -= rem;
											rem = 0;
										}
										else{
											rem -= nb.votes;
											nb.votes = 0;
										}
									}
									if(rem == 0){
										break;
									}
								}

								Doubles votecounts;
								for(int i = 0; i < newballots.size(); ++i){
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
									if(c.seat != -1 && elected.find(i) == elected.end()){
										// We have found an alternate outcome
										alternate = true;
										break;
									}	
								}
								if(alternate){
									weub = weub_next;	
								}
							}
							else if(!config.deleteonly && weub_next != -1){
								weub = weub_next; 
							}
						}
					}
				}
			}
			
			EliminateCandidate(toeliminate, cand, ballots, 
				votecounts, config, false);
		}
		else{
			// start with candidate with largest surplus
			Candidate &elect = cand[surpluses.front().id];
			double surplus = surpluses.front().weight;

			if(config.addonly){
				for(int i = 0; i < cand.size(); ++i){
					if(i == elect.index)
						continue;
					double tally = cand[i].sim_votes;

					if(elected.find(i) == elected.end() &&
						cand[i].standing){
						// How many votes to add to i's 
						// tally to reach 'tally'
						double wsub=max(0.0,elect.sim_votes-tally+1);
						double nq = 1.0 + ((config.totalvotes+wsub)/(config.nseats+1.0));
						double weub_next = -1;
						if(tally+nq >= nq){
							if(weub == -1){
								weub_next = wsub;
							}
							else if(weub > wsub){
								weub_next = wsub;
							}
						}

						if(weub_next != -1 && weub_next != weub){
							Ballots newballots(ballots);
							Ballot nb;
							nb.tag = newballots.size();
							nb.prefs.push_back(i);
							nb.votes = wsub;
							newballots.push_back(nb);
								
							Doubles votecounts;
							for(int i = 0; i < newballots.size(); ++i){
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
								if(c.seat != -1 && elected.find(i) == elected.end()){
									// We have found an alternate outcome
									alternate = true;
									break;
								}	
							}
							if(alternate){
								weub = weub_next;			
							}
						}
					}
				}
			}

			elect.seat = currseat++;
			elect.standing = 0;

			// distribute surplus
			DistributeSurplus(elect,cand,surplus,ballots,
				votecounts, config, false);

			surpluses.pop_front();
		}

		if(currseat == config.nseats){
			return weub;
		}

		++counter;
	}

	return weub;
}


// Function for simulating STV election
bool SimSTV(const Ballots &ballots, const Doubles &votecounts, 
	Candidates &cand, const Config &config, Ints &order_c, 
	Ints &order_a, bool log, double &mindiff)
{
	try	{
		double totvotes = 0;
		if(log) cout << "First preference tallies: " << endl;
		for(Candidates::iterator it = cand.begin(); it != cand.end(); ++it){
			it->sim_votes = 0;
			it->max_votes = 0;
			it->bweights.clear();
			it->standing = 1;
			it->seat = -1;
			it->surplus = 0;

			for(Ints::const_iterator li = it->ballots.begin();
				li != it->ballots.end(); ++li){
				it->bweights.push_back(IntDouble(*li, 1));
				it->sim_votes += votecounts[*li];
			}

            it->max_votes = it->sim_votes;
			totvotes += it->sim_votes;

			if(log) cout<<"    Candidate "<<it->id<<" "<<it->sim_votes<< endl;
		}

		// Step 1: Determine quota
		int quota = (int)(1.0 + (totvotes/(double)(config.nseats+1.0))); 

		if(log) cout << "The quota for election is " << quota << endl;
		L_IDS surpluses;		

		int currseat = 0;
		
		while(currseat < config.nseats){
			int standing = 0;

			// if a candidate has a quota, add them to the list of 
			// candidates with a surplus
			for(int i = 0; i < cand.size(); ++i){
				Candidate &c = cand[i];

				if(c.standing) 
					++standing;

				if(c.surplus || !c.standing)
					continue;

				if(c.standing && c.sim_votes >= quota){
					InsertCandidateSurplus(c, max(0.0,c.sim_votes-quota), 
						surpluses);
					c.surplus = 1;
				}
			}
		
			if(standing == config.nseats - currseat){
				if(log) cout << "Number of candidates left standing equals "<<
					" the number of remaining seats." << endl;

				// Elect all remaining candidates
				Ints standing;
				for(int i = 0; i < cand.size(); ++i){
					Candidate &c = cand[i];
					if(c.standing){
						Ints::iterator it = standing.begin();
						bool inserted = false;
						for(;  it != standing.end(); ++it){
							if(c.sim_votes > cand[*it].sim_votes){
								standing.insert(it, i);
								inserted = true;
								break;
							}
						}
						if(!inserted){
							standing.push_back(i);
						}
					}
				}


				for(int i = 0; i < standing.size(); ++i){
					Candidate &c = cand[standing[i]];
					if(log) cout<<"Candidate "<<c.id<<" elected (votes "
						<< c.sim_votes << ") " << endl;
					c.seat = currseat++;
					c.standing = 0;
					order_c.push_back(c.index);
					order_a.push_back(1);	
				}
				break;
			}
			
			if(surpluses.empty()){
				// eliminate candidate with fewest votes.
				// distribute votes at their value.
				double leastvotes = numeric_limits<double>::max();
				int toeliminate = -1;
				for(int i = 0; i < cand.size(); ++i){
					if(log && cand[i].standing){
						cout << "Candidate " << cand[i].id << " has "
							<< cand[i].sim_votes << " votes " << endl;
					}

					if(cand[i].standing && cand[i].sim_votes < leastvotes){
						leastvotes = cand[i].sim_votes;
						toeliminate = i;
					}
				}

				order_c.push_back(toeliminate);
				order_a.push_back(0);

				if(log){
					cout << "Candidate " << cand[toeliminate].id
						<< " eliminated on " << cand[toeliminate].sim_votes
						<< " votes" << endl; 
					cout << "Possible upper bound " << quota - 
						cand[toeliminate].sim_votes << endl;
				}

				EliminateCandidate(toeliminate, cand, ballots, 
					votecounts, config,log);
			}
			else{
				// start with candidate with largest surplus
				Candidate &elect = cand[surpluses.front().id];
				double surplus = surpluses.front().weight;

				elect.seat = currseat++;
				elect.standing = 0;

				order_c.push_back(elect.index);
				order_a.push_back(1);

				if(log) cout<<"Candidate "<<elect.id<<" elected (votes "
					<< elect.sim_votes << ") " << endl;

				// distribute surplus
				DistributeSurplus(elect,cand,surplus,ballots,
					votecounts, config,log);

				surpluses.pop_front();
			}

			if(currseat == config.nseats){
				// All seats filled.
				if(order_c.size() != cand.size()){
					for(int i = 0; i < cand.size(); ++i){
						const Candidate &c = cand[i];
						if(c.standing) {
							order_c.push_back(c.index);
							order_a.push_back(0);
						}
					}
				}
				break;
			}
		}
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(...)
	{
		throw STVException("Unexpected error in STV simulation.");
	}

	if(log) cout << "Simulation complete" << endl;
	return true;
}


void InsertCandidateSurplus(const Candidate &c,double surplus,L_IDS &surpluses)
{
	try{
		L_IDS::iterator it = surpluses.begin();

		for( ; it != surpluses.end(); ++it){
			if(it->weight > surplus)
				continue;

			break;
		}

		surpluses.insert(it, IntDouble(c.index, surplus));
	}
	catch(exception &e){
		throw STVException(string(e.what()));
	}
	catch(...){
		throw STVException("Unexpected error in InsertCandidateSurplus.");
	}
}


int NextCandidate(const Ballot &b, int index, const Candidates &cand)
{
	try{
		int idx = find(b.prefs.begin(), b.prefs.end(), index)-b.prefs.begin();
	
		for(int i = idx+1; i < b.prefs.size(); ++i){
			if(cand[b.prefs[i]].standing == 0 || cand[b.prefs[i]].surplus == 1) 
				continue;

			return b.prefs[i];
		}

		return -1;
	}
	catch(exception &e){
		throw STVException(string(e.what()));
	}
	catch(...){
		throw STVException("Unexpected error in NextCandidate.");
	}

	return -1;
}

void DistributeSurplus(Candidate &c, Candidates &cand, double surplus,
	const Ballots &ballots, const Doubles &votecounts, 
	const Config &config, bool log)
{
	c.standing = 0;
	if(surplus < 0.001) return;

	if(log) cout<<"Surplus of "<<surplus<<" to be distributed." << endl;

	try
	{
		// compute number of transferrable papers
		IDS2d totransfer(cand.size());

		IDS &tlist = c.bweights;

		double tpapers = 0;
		for(IDS::iterator it = tlist.begin(); it != tlist.end(); ++it)
		{
			int next = NextCandidate(ballots[it->id], c.index, cand);
			if(next >= 0)
			{
				tpapers += votecounts[it->id]*it->weight;
				totransfer[next].push_back(*it);
			}
		}

		double tvalue = min(1.0, surplus / tpapers);

		if(log) cout << "Transfer value: " << tvalue << endl;
		// the ballots in 'totransfer' will be distributed, each with a value
		// of 'tvalue'
		for(int i = 0; i < totransfer.size(); ++i)
		{
			const IDS &list = totransfer[i];
			if(list.empty()) 
				continue;

			Candidate &ct = cand[i];
			double total = 0;

			for(int j = 0; j < list.size(); ++j)
			{
				IntDouble td = list[j];
				td.weight *= tvalue;

				ct.bweights.push_back(td);

				total += votecounts[td.id] * td.weight;
			}

			if(log) cout << "    Transferring " << total <<
				" votes from " << c.id << " to " << ct.id << endl; 
			
			ct.sim_votes += total;
            ct.max_votes += total;
		}

		c.sim_votes -= surplus;
	}
	catch(STVException &e)
	{
		throw e;
	}
	catch(exception &e)
	{
		throw STVException(string(e.what()));
	}
	catch(...)
	{
		throw STVException("Unexpected error in DistributeSurplus.");
	}

} 


void EliminateCandidate(int toe, Candidates &cand, 
	const Ballots &ballots, const Doubles &votecounts,
	const Config &config, bool log)
{
	try{
		Candidate &e = cand[toe];
		e.standing = false;
		IDS2d totransfer(cand.size());

		// Distribute all ballots (at their current value) to rem candidates
		for(IDS::iterator it=e.bweights.begin();it!=e.bweights.end();++it){
			int next = NextCandidate(ballots[it->id], e.index, cand);
			if(next >= 0)
			{
				//cand[next].bweights.push_back(*it);
				totransfer[next].push_back(*it);
			}

			it->weight = 0;
		}

		e.sim_votes = 0;
		for(int i = 0; i < totransfer.size(); ++i){
			const IDS &list = totransfer[i];
			if(list.empty()) continue;
			Candidate &ct = cand[i];

			double total = 0;
			for(int j = 0; j < list.size(); ++j){
				total += votecounts[list[j].id] * list[j].weight;
				ct.bweights.push_back(list[j]);
			}

			ct.sim_votes += total;
            ct.max_votes += total;
			if(log) cout << total << " votes distributed from " << e.id <<
				" to " << ct.id << endl;
		}
	}
	catch(STVException &e){
		throw e;
	}
	catch(exception &e){
		throw STVException(string(e.what()));
	}
	catch(...){
		throw STVException("Unexpected error in EliminateCandidate.");
	}
}
