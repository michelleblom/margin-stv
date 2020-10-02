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


#ifndef _TREE_STV_H
#define _TREE_STV_H

#include "model.h"

bool RunTreeSTV(const Ballots &ballots, const Candidates &cands,
	const Config &config, const Ints &order_c, const Ints &order_a, 
	double upperbound,double timelimit,bool compbounds,const char *logf);

#endif
