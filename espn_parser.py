#!/usr/bin/env python

"""


"""

import urllib.request
import requests
from html.parser import HTMLParser
import re
import numpy as np
import sys
import json

ncaa_file_name = 'ncaa_fcs_stats.txt'
conference_file_name = 'conference_file.txt'

conference_divisions = {
    'FCS': ['Atlantic Sun Conference',\
            'Big Sky Conference', \
            'Big South Conference', \
            'Colonial Athletic Association', \
            'Ivy League', \
            'Mid-Eastern Athletic Conference', \
            'Missouri Valley Football Conference', \
            'Northeast Conference', \
            'Ohio Valley Conference', \
            'Patriot League', \
            'Pioneer Football League', \
            'Southern Conference', \
            'Southland Conference', \
            'Southwestern Athletic Conference', \
            'Western Athletic Conference'], \
    'FBS': ['American Athletic Conference', \
            'Atlantic Coast Conference', \
            'Big 12 Conference', \
            'Big Ten Conference', \
            'Conference USA', \
            'FBS Independents', \
            'Mid-American Conference', \
            'Mountain West Conference', \
            'Pac-12 Conference', \
            'Southeastern Conference', \
            'Sun Belt Conference'], \
    'D2': ['CIAA', \
           'GLIAC',\
           'Great Lakes', \
           'Gulf South', \
           'Lone Star', \
           'Northeast 10', \
           'Pennsylvania State Athletic Conference', \
           'Rocky Mountain', \
           'SIAC', \
           'South Atlantic'],
    'D3': ['II/III',\
           'So. Cal'],
}

##conferences = {
##    'Abilene Christian': 'ASUN-WAC',
##    'Alabama A&M': 'SWAC',  
##    'Alabama State': 'SWAC',
##    'Albany': 'CAA',
##    'Alcorn State': 'SWAC',
##    'Arkansas-Pine Bluff': 'SWAC',
##    'Austin Peay': 'ASUN-WAC',
##    'Bethune-Cookman': 'SWAC',
##    'Brown': 'Ivy',
##    'Bryant': 'NEC',
##    'Bucknell': 'PL',
##    'Butler': 'PFL',
##    'Cal Poly': 'BSC',
##    'Campbell': 'BSo',
##    'Central Arkansas': 'SLC',
##    'Central Connecticut': 'NEC',
##    'Charleston Southern': 'BSo',
##    'Chattanooga': 'SC',
##    'Colgate': 'PL',
##    'Columbia': 'Ivy',
##    'Cornell': 'Ivy',
##    'Dartmouth': 'Ivy',
##    'Davidson': 'PFL',
##    'Dayton': 'PFL',
##    'Delaware': 'CAA',
##    'Delaware State': 'MEAC',
##    'Drake': 'PFL',
##    'Duquesne': 'NEC',
##    'East Tennessee State': 'SC',
##    'Eastern Illinois': 'OVC',
##    'Eastern Kentucky': 'ASUN-WAC',
##    'Eastern Washington': 'BSC',
##    'Elon': 'CAA',
##    'Florida A&M': 'SWAC',
##    'Fordham': 'PL',
##    'Furman': 'SC',
##    'Gardner-Webb': 'BSo',
##    'Georgetown': 'PL',
##    'Grambling': 'SWAC',
##    'Hampton': 'CAA',
##    'Harvard': 'Ivy',
##    'Holy Cross': 'PL',
##    'Houston Christian': 'SLC',
##    'Howard': 'MEAC',
##    'Idaho': 'BSC',
##    'Idaho State': 'BSC',
##    'Illinois State': 'MVFC',
##    'Incarnate Word': 'SLC',
##    'Indiana State': 'MVFC',
##    'Jackson State': 'SWAC',
##    'Jacksonville State': 'ASUN-WAC',
##    'Kennesaw State': 'ASUN-WAC',
##    'Lafayette': 'PL',
##    'Lamar': 'SLC',
##    'Lehigh': 'PL',
##    'Lindenwood': 'OVC',
##    'Long Island University': 'NEC',
##    'Maine': 'CAA',
##    'Marist': 'PFL',
##    'McNeese': 'SLC',
##    'Mercer': 'SC',
##    'Merrimack': 'NEC',
##    'Mississippi Valley State': 'SWAC',
##    'Missouri State': 'MVFC',
##    'Monmouth': 'CAA',
##    'Montana': 'BSC',
##    'Montana State': 'BSC',
##}

# Extract the conferences
def parse_conferences(gdict, data, ops):

    # Find the index in the data where the conference standings are
    # -------------------------------------------------------------
    end_idx1 = data.index('Full Standings')
    beg_idx1 = data.index('CONF') - 2

    # Determine if there are two sets of standings
    # --------------------------------------------
    try:
        # If there are two, set the conference for the first and second team
        # ------------------------------------------------------------------
        end_idx2 = data[end_idx1 + 1:].index('Full Standings') + end_idx1 + 1
        beg_idx2 = data[end_idx1 + 1:].index('CONF') + end_idx1 + 1 - 2

        gdict[ops[0]]['Conference'] = ' '.join(data[beg_idx1].split()[1:-1])
        gdict[ops[1]]['Conference'] = ' '.join(data[beg_idx2].split()[1:-1])

    except ValueError:

        # Could be because both teams are in the same conference. Check
        # if both teams are in the standings
        # -------------------------------------------------------------
        if((ops[0] in data[beg_idx1 : end_idx1]) &
           (ops[1] in data[beg_idx1 : end_idx1])):
            gdict[ops[0]]['Conference'] = ' '.join(data[beg_idx1].split()[1:-1])
            gdict[ops[1]]['Conference'] = ' '.join(data[beg_idx1].split()[1:-1])

        # If only one, look at the teams in the standings to figure out
        # which team the conference is associated with
        # -------------------------------------------------------------
        else:
            if(ops[0] in data[beg_idx1 : end_idx1]):
                gdict[ops[0]]['Conference'] = ' '.join(data[beg_idx1].split()[1:-1])
                gdict[ops[1]]['Conference'] = 'NONE' 
            else:
                gdict[ops[1]]['Conference'] = ' '.join(data[beg_idx1].split()[1:-1])
                gdict[ops[0]]['Conference'] = 'NONE' 
    

def parse_indiv_pass_stats(gdict, week, op, total_passing):

    # Figure out how many stats are included
    # --------------------------------------
    tester = total_passing[:7]
    if('QBR' in tester):
        max_idx = 7
    else:
        max_idx = 6
    
    # Parse the individual passing stats
    # ----------------------------------
    indiv_pass = total_passing[max_idx:-max_idx]

    # First, check if the team has any receiving yards
    # ------------------------------------------------ 
    if(len(indiv_pass) == 0):
        print("No passing yards for ", op)

    else:
     
        # Determine the number of passers
        # -------------------------------
        num_passers = len(indiv_pass[::max_idx + 1])

        split_passers = np.array_split(np.array(indiv_pass), \
            num_passers)

        # Loop over the stats and add them to the dictionary
        # -------------------------------------------------- 
        for ii in range(num_passers):
            local_list = list(split_passers[ii])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]] = {}
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]][\
                'Completions'] = float(local_list[2].split('/')[0])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]][\
                'Attempts'] = float(local_list[2].split('/')[1])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]][\
                'Yards'] = float(local_list[3])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]][\
                'Average'] = float(local_list[4])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]][\
                'TD'] = float(local_list[5])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][local_list[0]][\
                'INT'] = float(local_list[6])
#            gdict[opponents[0]]['Scores'][int(week)]['Individual'][local_list[0]] = {}

def parse_total_pass_yards(gdict, week, data, opponents):
    beg_idx = data.index(opponents[0] + ' Passing') 
    end_idx = data.index(opponents[1] + ' Passing') 
    total_stat = data[beg_idx : end_idx]
    try:
        tyard_idx = total_stat.index('TEAM')
        op1_pyards = total_stat[tyard_idx + 1 : tyard_idx + 6]
        #op1_pyards = total_stat[tyard_idx + 2]
    except ValueError:
        op1_pyards = ['0/0','0','0','0','0']

    final1 = np.zeros(len(op1_pyards) + 1) 
    final1[0] = float(op1_pyards[0].split('/')[0])  # Completions
    final1[1] = float(op1_pyards[0].split('/')[1])  # Attempts
    final1[2] = float(op1_pyards[1])                # Yards
    final1[3] = float(op1_pyards[2])                # Average 
    final1[4] = float(op1_pyards[3])                # TDs
    final1[5] = float(op1_pyards[4])                # INTs

    # Parse the individual passing stats for team 1
    # ---------------------------------------------
    parse_indiv_pass_stats(gdict, week, opponents[0], total_stat)

    beg_idx = data.index(opponents[1] + ' Passing') 
    end_idx = data.index(opponents[0] + ' Rushing') 
    total_stat = data[beg_idx : end_idx]
    tyard_idx = total_stat.index('TEAM')
    op2_pyards = total_stat[tyard_idx + 1 : tyard_idx + 6]
    #op2_pyards = total_stat[tyard_idx + 2]

    final2 = np.zeros(len(op2_pyards) + 1)          
    final2[0] = float(op2_pyards[0].split('/')[0])  # Completions
    final2[1] = float(op2_pyards[0].split('/')[1])  # Attempts
    final2[2] = float(op2_pyards[1])                # Yards
    final2[3] = float(op2_pyards[2])                # Average 
    final2[4] = float(op2_pyards[3])                # TDs
    final2[5] = float(op2_pyards[4])                # INTs

    # Parse the individual passing stats for team 1
    # ---------------------------------------------
    parse_indiv_pass_stats(gdict, week, opponents[1], total_stat)

    # Add the team statistics
    # -----------------------
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_completions'] = final1[0]
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_attempts']    = final1[1]
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_yards']       = final1[2]
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_completions'] = final2[0]
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_attempts']    = final2[1]
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_yards']       = final2[2]

    gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_completions'] = final1[0]
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_attempts']    = final1[1]
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_yards']       = final1[2]
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_completions'] = final2[0]
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_attempts']    = final2[1]
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_yards']       = final2[2]

    #return final1, final2

def parse_indiv_rush_stats(gdict, week, op, total_rushing):

    # Parse the individual rushing stats
    # ----------------------------------
    indiv_rush = total_rushing[6:-6]
  
    # Determine the number of rushers
    # -------------------------------
    num_rushers = len(indiv_rush[::7])
 
    split_rushers = np.array_split(np.array(indiv_rush), \
        num_rushers)

    # Loop over the stats and add them to the dictionary
    # -------------------------------------------------- 
    for ii in range(num_rushers):
        local_list = list(split_rushers[ii])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_list[0]] = {}
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_list[0]][\
            'Carries'] = float(local_list[2])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_list[0]][\
            'Yards'] = float(local_list[3])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_list[0]][\
            'Average'] = float(local_list[4])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_list[0]][\
            'TD'] = float(local_list[5])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_list[0]][\
            'Long'] = float(local_list[6])
#        gdict[opponents[0]]['Scores'][int(week)]['Individual'][local_list[0]] = {}
        
def parse_total_rush_yards(gdict, week, data, opponents):
    # Team 1
    beg_idx = data.index(opponents[0] + ' Rushing') 
    end_idx = data.index(opponents[1] + ' Rushing') 
    total_rushing = data[beg_idx : end_idx]
    op1_ryards = total_rushing[-4]

    # Parse the individual rushing stats for team 1
    # ---------------------------------------------
    parse_indiv_rush_stats(gdict, week, opponents[0], total_rushing)

    # Team 2
    beg_idx = data.index(opponents[1] + ' Rushing') 
    end_idx = data.index(opponents[0] + ' Receiving') 
    total_rushing = data[beg_idx : end_idx]
    op2_ryards = total_rushing[-4]

    # Parse the individual rushing stats for team 2
    # ---------------------------------------------
    parse_indiv_rush_stats(gdict, week, opponents[1], total_rushing)

    # Add the team stats to the dictionary
    # ------------------------------------
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Rushing'] = float(op1_ryards)
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Rushing'] = float(op2_ryards)
    
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Rushing'] = float(op1_ryards)
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Rushing'] = float(op2_ryards)

    #return op1_ryards, op2_ryards

def parse_indiv_receiving_stats(gdict, week, op, total_receiving):

    # Parse the individual receiving stats
    # ----------------------------------
    indiv_receiving = total_receiving[6:-6]

    # First, check if the team has any receiving yards
    # ------------------------------------------------ 
    if(len(indiv_receiving) == 0):
        print("No receiving yards for ", op)

    else:
 
        # Determine the number of receivers
        # -------------------------------
        num_receivers = len(indiv_receiving[::7])
 
        split_receivers = np.array_split(np.array(indiv_receiving), \
            num_receivers)

        # Loop over the stats and add them to the dictionary
        # -------------------------------------------------- 
        for ii in range(num_receivers):
            local_list = list(split_receivers[ii])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_list[0]] = {}
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_list[0]][\
                'Catches'] = float(local_list[2])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_list[0]][\
                'Yards'] = float(local_list[3])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_list[0]][\
                'Average'] = float(local_list[4])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_list[0]][\
                'TD'] = float(local_list[5])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_list[0]][\
                'Long'] = float(local_list[6])
#            gdict[opponents[0]]['Scores'][int(week)]['Individual'][local_list[0]] = {}
        
def parse_total_receiving_yards(gdict, week, data, opponents):

    # Team 1
    beg_idx = data.index(opponents[0] + ' Receiving') 
    end_idx = data.index(opponents[1] + ' Receiving') 
    total_stat = data[beg_idx : end_idx]
    op1_ryards = total_stat[-4]

    # Parse the individual rushing stats for team 1
    # ---------------------------------------------
    parse_indiv_receiving_stats(gdict, week, opponents[0], total_stat)

    # Team 2
    beg_idx = data.index(opponents[1] + ' Receiving') 
    end_idx = data.index(opponents[0] + ' Fumbles') 
    total_stat = data[beg_idx : end_idx]
    op2_ryards = total_stat[-4]

    # Parse the individual rushing stats for team 2
    # ---------------------------------------------
    parse_indiv_receiving_stats(gdict, week, opponents[1], total_stat)

def parse_scores(data, opponents):
    try:
        beg_idx = data.index('Final')
    except ValueError:
        try:
            beg_idx = data.index('Final/OT')
        except ValueError:
            beg_idx = None
            ot_idx = 2
            while beg_idx is None: 
                try:
                    beg_idx = data.index(\
                        'Final/{0}OT'.format(ot_idx))
                except:
                    ot_idx += 1
                
                if(ot_idx > 10):
                    beg_idx = 10
        
    end_idx = data.index('T')
    idx_diff = (end_idx - beg_idx) * 4 - 1
    #end_idx = data.index(opponents[1])
    
    total_points = data[beg_idx : beg_idx + idx_diff]
    t_idx = total_points.index('T')
    op1_points = total_points[2 + t_idx * 1 : t_idx * 2 + 2]
    op2_points = total_points[3 + t_idx * 2 : t_idx * 3 + 3]

    return op1_points, op2_points

def calc_team_stats(gdict, weekly = False):

    if(weekly):
        for tkey in gdict.keys():
            w2key = list(gdict[tkey]['Scores'].keys())[-1]
            gdict[tkey]['Stats']['Weekly'][w2key]['scoring_offense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['rushing_offense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['passing_offense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['scoring_defense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['rushing_defense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['passing_defense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['total_offense'] = \
                gdict[tkey]['Stats']['Weekly'][w2key]['passing_offense'] + \
                gdict[tkey]['Stats']['Weekly'][w2key]['rushing_offense'] 
            gdict[tkey]['Stats']['Weekly'][w2key]['total_defense'] = \
                gdict[tkey]['Stats']['Weekly'][w2key]['passing_defense'] + \
                gdict[tkey]['Stats']['Weekly'][w2key]['rushing_defense'] 
    else:
        for tkey in gdict.keys():
            #gdict[tkey]['Stats']['Season'] = {}
            gdict[tkey]['Stats']['Season']['scoring_offense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['rushing_offense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['passing_offense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['scoring_defense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['rushing_defense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['passing_defense'] =  \
                np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['total_offense'] = \
                gdict[tkey]['Stats']['Season']['passing_offense'] + \
                gdict[tkey]['Stats']['Season']['rushing_offense'] 
            gdict[tkey]['Stats']['Season']['total_defense'] = \
                gdict[tkey]['Stats']['Season']['passing_defense'] + \
                gdict[tkey]['Stats']['Season']['rushing_defense'] 

def calc_team_rankings(gdict, lkeys, division):

    for var in gdict['North Dakota']['Stats'].keys():

        # Rank scoring offense
        # --------------------
        values = [(gdict[tkey]['Stats'][var], tkey) \
            for tkey in lkeys]
            #for tkey in gdict.keys()]
        values.sort()
        
        if((var == 'scoring_offense') | \
           (var == 'passing_offense') | \
           (var == 'rushing_offense') | \
           (var == 'total_offense')):
            values = values[::-1]

        for ii in range(len(values)):
            gdict[values[ii][1]]['Rankings'][division][var] = ii + 1

def calc_game_quality(gdict, team, week):

    opponent = gdict[team]['Scores'][week]['Opponent']

    team_dict = gdict[team]['Scores'][week]
    oppt_dict = gdict[opponent]['Stats']['Weekly'][week-1]

    # Team rankings factor
    # --------------------

    # Win size factor
    # ---------------
    win_size_factor = team_dict['Team']['Offense']['Points'] / \
                      team_dict['Team']['Defense']['Points']

    # Home/road factor
    # ----------------
    if(win_size_factor > 1.):
        if(team_dict['Location'] == 'Away'):
            loc_factor = 1.25
        else:
            loc_factor = 1.0
    else:
        if(team_dict['Location'] == 'Away'):
            loc_factor = 1.0
        else:
            loc_factor = 0.75
        

    # -----------------------
    # Calculate team factors
    # -----------------------

    # Scoring factor
    # ---------------
    off_scoring_factor = team_dict['Team']['Offense']['Points'] / \
                         oppt_dict['scoring_defense']
    def_scoring_factor = team_dict['Team']['Defense']['Points'] / \
                         oppt_dict['scoring_offense']

    # Team rush yards factor
    # -----------------
    off_rush_factor = team_dict['Team']['Offense']['Rushing'] / \
                      oppt_dict['rushing_defense']
    def_rush_factor = team_dict['Team']['Defense']['Rushing'] / \
                      oppt_dict['rushing_offense']

    # Team pass yards factor
    # -----------------
    off_pass_factor = team_dict['Team']['Offense']['Pass_yards'] / \
                      oppt_dict['passing_defense']
    def_pass_factor = team_dict['Team']['Defense']['Pass_yards'] / \
                      oppt_dict['passing_offense']
   
    print(win_size_factor, loc_factor, off_scoring_factor, def_scoring_factor, \
        off_rush_factor, def_rush_factor, off_pass_factor, \
        def_pass_factor)

    print(win_size_factor * loc_factor * off_scoring_factor * def_scoring_factor * \
        off_rush_factor * def_rush_factor * off_pass_factor * \
        def_pass_factor)

# For a specified week, calculate
#def calc_weekly_stats(gdict, division = 'All', week = None):

    # Loop over the teams
    # -------------------
    #for tkey in gdict.keys():
    
        # Calculate the average points scored and allowed 
        # so far
        # -----------------------------------------------

        # Calculate the
        
 
def rank_team_stats(gdict, division = 'All'):


    if(division == 'All'):
        lkeys = gdict.keys()
        calc_team_rankings(gdict, lkeys, division)

    elif(division == 'Division'):
        cdict = read_conference_file()

        dkeys = conference_divisions.keys()
        for dkey in dkeys: 
            try:
                lkeys = []
                for ii in range(len(conference_divisions[dkey])):   
                    lkeys = lkeys + \
                        cdict[conference_divisions[dkey][ii]]

                calc_team_rankings(gdict, lkeys, division)
            except KeyError:
                print("ERROR: no division rankings for", dkey) 

    elif(division == 'Conference'):
        cdict = read_conference_file()

        dkeys = conference_divisions.keys()
        for dkey in dkeys: 
            for ii in range(len(conference_divisions[dkey])):
                try:
                    lkeys = cdict[conference_divisions[dkey][ii]]
                    calc_team_rankings(gdict, lkeys, division)
                except KeyError:
                    print("ERROR: no conference rankings for", dkey)


class MyHTMLParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.Myrawdata = ''
    #def handle_starttag(self, tag, attr):
    #    print("Encountered a start tag:", tag)
    #def handle_endtag(self, tag, attr):
    #    print("Encountered an end tag:", tag)
    def handle_data(self, data):
        if(data[0] != '\\'):
            self.Myrawdata = self.Myrawdata + '\n' + data

def pull_ESPN_data(url):

    uf = urllib.request.urlopen(url)
    html = uf.read()
    parser = MyHTMLParser()
    parser.feed(str(html))
    return parser.Myrawdata.strip().split('\n')
       
def load_ESPN_data(start_week, end_week):

    gdict = {}
    espn_url = "https://www.espn.com/college-football/scoreboard/_/week/{0}/year/2022/seasontype/2/group/81"
    weeks = np.arange(start_week,end_week + 1)
    
    # Loop over all the weeks and add the data
    # ----------------------------------------
    for week in weeks:
        url = espn_url.format(week)
        print(url)
    
    #sys.exit()
    
        #url = "https://www.espn.com/college-football/scoreboard/_/week/11/year/2022/seasontype/2/group/81"
        uf = urllib.request.urlopen(url)
        html = uf.read()
        #parser.feed(html)
        
        urls = []
        
        for match in re.finditer("college-football/boxscore/_/gameId/", str(html)):
            urls.append(str(html)[match.start() : match.end() + 9])
        
        final_urls = np.unique(np.array(urls))
        base_url = 'https://www.espn.com/'
        parser = MyHTMLParser()
        
        
        
        for furl in final_urls:
            local_url = base_url + furl
            data = pull_ESPN_data(local_url)

            # Parse the opponents
            # -------------------
            matchup = data[1].split(' - Box Score - ')[0]
            print(matchup)
            opponents = matchup.split(' vs. ')
            for op in opponents:
                if(op not in gdict.keys()):
                    gdict[op] = {}
                    gdict[op]['Stats'] = {}
                    gdict[op]['Stats']['Season'] = {}
                    gdict[op]['Stats']['Season']['Record'] = {}
                    gdict[op]['Stats']['Season']['Beat']    = []
                    gdict[op]['Stats']['Season']['Lost to'] = []
                    gdict[op]['Stats']['Season']['Record']['Overall'] = {}
                    gdict[op]['Stats']['Season']['Record']['D1'] = {}
                    gdict[op]['Stats']['Season']['Record']['Conference'] = {}
                    gdict[op]['Stats']['Season']['Record']['Overall']['W'] = 0
                    gdict[op]['Stats']['Season']['Record']['Overall']['L'] = 0
                    gdict[op]['Stats']['Season']['Record']['D1']['W'] = 0
                    gdict[op]['Stats']['Season']['Record']['D1']['L'] = 0
                    gdict[op]['Stats']['Season']['Record']['Conference']['W'] = 0
                    gdict[op]['Stats']['Season']['Record']['Conference']['L'] = 0
                    gdict[op]['Stats']['Weekly'] = {}
                    gdict[op]['Scores'] = {}
                    gdict[op]['Rankings'] = {}
                    gdict[op]['Rankings']['All'] = {}
                    gdict[op]['Rankings']['Division'] = {}
                    gdict[op]['Rankings']['Conference'] = {}
    
                gdict[op]['Stats']['Weekly'][int(week)] = {}
                gdict[op]['Scores'][int(week)] = {}
                gdict[op]['Scores'][int(week)]['Team'] = {}
                gdict[op]['Scores'][int(week)]['Team']['Offense'] = {}
                #gdict[op][int(week)]['Scores']['Team']['Offense']['Points'] = {}
                #gdict[op][int(week)]['Scores']['Team']['Offense']['Passing'] = {}
                #gdict[op][int(week)]['Scores']['Team']['Offense']['Rushing'] = {}
                #gdict[op][int(week)]['Team']['Offense']['Receiving'] = {}
    
                gdict[op]['Scores'][int(week)]['Team']['Defense'] = {}
                #gdict[op][int(week)]['Scores']['Team']['Defense']['Points'] = {}
                #gdict[op][int(week)]['Scores']['Team']['Defense']['Passing'] = {}
                #gdict[op][int(week)]['Scores']['Team']['Defense']['Rushing'] = {}
                #gdict[op][int(week)]['Team']['Defense']['Receiving'] = {}

                gdict[op]['Scores'][int(week)]['Individual'] = {}
                gdict[op]['Scores'][int(week)]['Individual']['Offense'] = {}
                gdict[op]['Scores'][int(week)]['Individual']['Offense']['Rushing'] = {}
                gdict[op]['Scores'][int(week)]['Individual']['Offense']['Passing'] = {}
                gdict[op]['Scores'][int(week)]['Individual']['Offense']['Receiving'] = {}
                gdict[op]['Scores'][int(week)]['Individual']['Defense'] = {}
      
            # Probably don't need to have run each time, but can fix later 
            parse_conferences(gdict, data, opponents)

            # Set the home / away locations
            # -----------------------------
            gdict[opponents[0]]['Scores'][int(week)]['Opponent'] = opponents[1] 
            gdict[opponents[0]]['Scores'][int(week)]['Location'] = 'Away' 
            
            gdict[opponents[1]]['Scores'][int(week)]['Opponent'] = opponents[0] 
            gdict[opponents[1]]['Scores'][int(week)]['Location'] = 'Home' 
     
            # Parse the score 
            # NOTE: contains the quarterly score info also
            # ---------------------------------------------
            op1_points, op2_points = parse_scores(data, opponents)
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Points'] = float(op1_points[-1])
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Points'] = float(op2_points[-1])
    
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Points'] = float(op1_points[-1])
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Points'] = float(op2_points[-1])
            
            # Use these to update the win and loss columns for each team
            # ----------------------------------------------------------
            if(float(op1_points[-1]) > float(op2_points[-1])):
                # Team 1 won, Team 2 lost
                gdict[opponents[0]]['Stats']['Season']['Record']['Overall']['W'] += 1
                gdict[opponents[1]]['Stats']['Season']['Record']['Overall']['L'] += 1

                gdict[opponents[0]]['Stats']['Season']['Beat'].append(opponents[1])
                gdict[opponents[1]]['Stats']['Season']['Lost to'].append(opponents[0])

                # Update D1 win counter
                # ---------------------
                if(gdict[opponents[0]]['Conference'] in conference_divisions['FBS'] or
                   gdict[opponents[0]]['Conference'] in conference_divisions['FCS']):
                    gdict[opponents[0]]['Stats']['Season']['Record']['D1']['W'] += 1

                if(gdict[opponents[1]]['Conference'] in conference_divisions['FBS'] or
                   gdict[opponents[1]]['Conference'] in conference_divisions['FCS']):
                
                    gdict[opponents[1]]['Stats']['Season']['Record']['D1']['L'] += 1

                # Update conference w/l
                # ---------------------
                if(gdict[opponents[0]]['Conference'] == gdict[opponents[1]]['Conference']):
                    gdict[opponents[0]]['Stats']['Season']['Record']['Overall']['W'] += 1
                    gdict[opponents[1]]['Stats']['Season']['Record']['Overall']['L'] += 1

            else:
                # Team 2 won, Team 1 lost
                gdict[opponents[1]]['Stats']['Season']['Record']['Overall']['W'] += 1
                gdict[opponents[0]]['Stats']['Season']['Record']['Overall']['L'] += 1

                gdict[opponents[1]]['Stats']['Season']['Beat'].append(opponents[0])
                gdict[opponents[0]]['Stats']['Season']['Lost to'].append(opponents[1])

                # Update D1 win counter
                # ---------------------
                if(gdict[opponents[1]]['Conference'] in conference_divisions['FBS'] or
                   gdict[opponents[1]]['Conference'] in conference_divisions['FCS']):
                    gdict[opponents[1]]['Stats']['Season']['Record']['D1']['W'] += 1

                if(gdict[opponents[0]]['Conference'] in conference_divisions['FBS'] or
                   gdict[opponents[0]]['Conference'] in conference_divisions['FCS']):
                    gdict[opponents[0]]['Stats']['Season']['Record']['D1']['L'] += 1

                # Update conference w/l
                # ---------------------
                if(gdict[opponents[0]]['Conference'] == gdict[opponents[1]]['Conference']):
                    gdict[opponents[1]]['Stats']['Season']['Record']['Overall']['W'] += 1
                    gdict[opponents[0]]['Stats']['Season']['Record']['Overall']['L'] += 1

            # Parse the total passing yards 
            parse_total_pass_yards(gdict, week, data, opponents)
         
            # Parse the total rushing yards 
            parse_total_rush_yards(gdict, week, data, opponents)

            # Parse the total receiving yards
            parse_total_receiving_yards(gdict, week, data, opponents)
 
        # If here, calculates stats after each week
        calc_team_stats(gdict, weekly = True)

    # If here, calculates stats after all weeks 
    calc_team_stats(gdict)

    return gdict, data

def write_conference_file(gdict):

    # Set up the conference dictionary
    # --------------------------------
    conference_dict = {}
    
    # Extract the conference info
    # ---------------------------
    conference_pairs = [(gdict[tkey]['Conference'], tkey) \
        for tkey in gdict.keys()]

    cfcs  = np.array([cpair[0] for cpair in conference_pairs])
    teams = np.array([cpair[1] for cpair in conference_pairs])
    cfc_unique = np.unique(cfcs)

    for cfc in cfc_unique:
        conference_dict[cfc] = list(teams[np.where(cfcs = cfc)])

    with open(conference_file_name,'w') as fout:
        json.dump(gdict, fout, indent = 4, sort_keys = True)

def read_conference_file():
    with open(conference_file_name, 'r') as fin:
        conference_dict = json.load(fin)

    return conference_dict

def write_json_stats(gdict):
    with open(ncaa_file_name,'w') as fout:
        json.dump(gdict, fout, indent = 4, sort_keys = True)

def load_json_stats():
    with open(ncaa_file_name, 'r') as fin:
        gdict = json.load(fin)

    return gdict
