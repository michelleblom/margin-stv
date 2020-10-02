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


#include<iostream>
#include<fstream>
#include<string.h>
#include<stdlib.h>
#include<algorithm>
#include "model.h"
#include "sim_stv.h"
#include "stv_distance.h"
#include "tree_stv.h"

using namespace std;

int main(int argc, const char * argv[]) 
{
	try
	{
		Candidates candidates;
		Ballots ballots; 
		Config config;

		double timelimit = 86400;
		bool compbounds = false;

		config.quota = -1;
		double tlimit_wsol = -1;
		double tlimit_wosol = -1;
		double tlimit_leaf = -1;

		const char *logf = NULL;
		bool simlog = false;
		for(int i = 1; i < argc; ++i)
		{
			if(strcmp(argv[i], "-config") == 0 && i < argc-1)
			{
				cout << "Reading config" << endl;
				if(!ReadConfig(argv[i+1], config, candidates))
				{
					cout << "Config file read error. Exiting."<<endl;
					return 1;
				}
				++i;
			}
			else if(strcmp(argv[i], "-ballots") == 0 && i < argc-1)
			{
				if(config.ncandidates <= 0)
				{
					cout<<"Need to specify config file first"<<endl;
					return 1;
				}

				if(!ReadBallots(argv[i+1], ballots, candidates, config))
				{
					cout << "Ballot read error. Exiting." << endl;
					return 1;
				}

				++i;
			}
			else if(strcmp(argv[i], "-compbounds") == 0){
				compbounds = true;
			}
			else if(strcmp(argv[i], "-simlog") == 0){
				simlog = true;
			}
			else if(strcmp(argv[i], "-tlimit")== 0 && i < argc-1){
				timelimit = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-tlimit_leaf")== 0 && i < argc-1){
				tlimit_leaf = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-tlimit_wsol")== 0 && i < argc-1){
				tlimit_wsol = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-tlimit_wosol")== 0 && i < argc-1){
				tlimit_wosol = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-logfile")== 0 && i < argc-1){
				logf = argv[i+1];
				++i;
			}
			else if(strcmp(argv[i], "-threads") == 0 && i < argc-2){
				config.threads = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-sthreads") == 0 && i < argc-2){
				config.subthreads = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-fixuntil") == 0 && i < argc-2){
			    config.FIX_UNTIL = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-bilin") == 0 && i < argc-2){
			    config.DIV_BILIN = atoi(argv[i+1]);
				++i;
			}
			else if(strcmp(argv[i], "-ao") == 0){
				config.addonly = true;
			}
			else if(strcmp(argv[i], "-do") == 0){
				config.deleteonly = true;
			}
		}

		if(config.addonly && config.deleteonly){
			cout << "Cannot run under both addition only and delete only settings." << endl;
			return 4;
		}

		if(tlimit_wsol != -1)
			config.tlimit_wsol = tlimit_wsol;
		if(tlimit_wosol != -1)
			config.tlimit_wosol = tlimit_wosol;
		if(tlimit_leaf != -1)
			config.tlimit_leaf = tlimit_leaf;

		double upperbound = config.totalvotes;
		if(config.quota == -1){
			config.quota = (int)(1.0 + ((double)config.totalvotes/
				(double)(config.nseats+1.0)));
		}


		mytimespec start;
		GetTime(&start);

		// code to simulate stv
		Ints order_c;
		Ints order_a;
	
		double minqdiff = config.totalvotes;	
		
		Doubles votecounts(ballots.size(), 0);
		for(int i = 0; i < ballots.size(); ++i){
			votecounts[i] = ballots[i].votes;
		}

		if(!SimSTV(ballots, votecounts, candidates, config, 
			order_c, order_a, simlog, minqdiff))
		{
			cout << "Simulation Failed. Exiting." << endl;
			return 1;
		}

		cout << "Election/Elimination order: ";
		Ints elected(config.ncandidates, 0);

		int seats = 0;
		set<int> elected_list;	
		for(int i = 0; i < order_c.size(); ++i)
		{
			const int indx = order_c[i];
			cout << candidates[indx].id << "-";
			if(order_a[i]){
				cout << "q ";
				++seats;
				elected[indx] =1;
				elected_list.insert(indx);
			}
			else{
				cout << "e ";
				minqdiff = min(minqdiff, config.quota - 
					candidates[indx].sum_votes);
			}
		}

		cout << endl;
		
		if(config.gub != -1){
			cout << "Provided upper bound = " << config.gub << endl;
			upperbound = min(upperbound, config.gub);
		}

		if(!(config.deleteonly || config.addonly)){	
			cout << "Min qdiff = " << minqdiff << " total votes " <<
				config.totalvotes << endl;
			upperbound = min(upperbound, minqdiff);
		}
	
		cout << "Computing WEUB/WSUB" << endl;
		double weub = WEUB(ballots,votecounts, candidates,
			config, elected_list, order_c, order_a);

		cout << "WEUB/WSUB = " << weub << endl;
		if(weub >= 0){
			upperbound = min(weub, upperbound);
		}
		
		cout << "Upper bound used = " << upperbound << endl;
		cout << "STARTING TREE SEARCH" << endl; 

		bool r = RunTreeSTV(ballots, candidates, config,
			order_c, order_a, upperbound, timelimit, 
			compbounds, logf);

		if(!r){
			cout << "Tree Search FAILED" << endl;
		}


		mytimespec tend;
		GetTime(&tend);
		cout << "Total time: " << tend.seconds - start.seconds << endl;
	}
	catch(exception &e)
	{
		cout << e.what() << endl;
		cout << "Exiting." << endl;
		return 1;
	}
	catch(STVException &e)
	{
		cout << e.what() << endl;
		cout << "Exiting." << endl;
		return 1;
	}	
	catch(...)
	{
		cout << "Unexpected error. Exiting." << endl;
		return 1;
	}

	return 0;
}




