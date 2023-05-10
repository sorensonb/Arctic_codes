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
from unidecode import unidecode

ncaa_file_name = 'ncaa_fcs_stats_2022.txt'
conference_file_name = 'conference_file_v2.txt'
massey_file = 'final_massey_composite_2022.txt'

# NOTE: made up - no data found
num_fcs_conf = 15
fcs_conf_rankings = {
    'Missouri Valley Football Conference':  1, 
    'Big Sky Conference':                   2, 
    'Atlantic Sun Conference':              3, 
    'Colonial Athletic Association':        4, 
    'Southern Conference':                  5,
    'Western Athletic Conference':          6, 
    'Ivy League':                           7, 
    'Ohio Valley Conference':               8,
    'Big South Conference':                 9, 
    'Southland Conference':                 10,
    'Patriot League':                       11,
    'Southwestern Athletic Conference':     12,
    'Mid-Eastern Athletic Conference':      13,
    'Northeast Conference':                 14,
    'Pioneer Football League':              15, 
}

fbs_conf_rankings = { 
    'Southeastern Conference':      1, 
    'Big Ten Conference':           2, 
    'Big 12 Conference':            3, 
    'Pac-12 Conference':            4, 
    'Atlantic Coast Conference':    5, 
    'American Athletic Conference': 6, 
    'Mountain West Conference':     7, 
    'Conference USA':               8, 
    'FBS Independents':             9, 
    'Sun Belt Conference':          10, 
    'Mid-American Conference':      11, 
}

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

def update_scores(gdict, opponents, week, op1_points, op2_points):
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Points'] = float(op1_points[-1])
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Points'] = float(op2_points[-1])
    
    gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Points'] = float(op1_points[-1])
    gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Points'] = float(op2_points[-1])
    
    # Use these to update the win and loss columns for each team
    # ----------------------------------------------------------
    if(float(op1_points[-1]) > float(op2_points[-1])):
        # Team 1 won, Team 2 lost
        gdict[opponents[0]]['Record']['Season']['Overall']['W'] += 1
        gdict[opponents[1]]['Record']['Season']['Overall']['L'] += 1
    
        gdict[opponents[0]]['Beat'].append(opponents[1])
        gdict[opponents[1]]['Lost to'].append(opponents[0])
    
        # Update D1 win counter
        # ---------------------
        if(gdict[opponents[0]]['Conference'] in conference_divisions['FBS'] or
           gdict[opponents[0]]['Conference'] in conference_divisions['FCS']):
            gdict[opponents[0]]['Record']['Season']['D1']['W'] += 1
    
        if(gdict[opponents[1]]['Conference'] in conference_divisions['FBS'] or
           gdict[opponents[1]]['Conference'] in conference_divisions['FCS']):
        
            gdict[opponents[1]]['Record']['Season']['D1']['L'] += 1
    
        # Update conference w/l
        # ---------------------
        if(gdict[opponents[0]]['Conference'] == gdict[opponents[1]]['Conference']):
            gdict[opponents[0]]['Record']['Season']['Conference']['W'] += 1
            gdict[opponents[1]]['Record']['Season']['Conference']['L'] += 1
    
    else:
        # Team 2 won, Team 1 lost
        gdict[opponents[1]]['Record']['Season']['Overall']['W'] += 1
        gdict[opponents[0]]['Record']['Season']['Overall']['L'] += 1
    
        gdict[opponents[1]]['Beat'].append(opponents[0])
        gdict[opponents[0]]['Lost to'].append(opponents[1])
    
        # Update D1 win counter
        # ---------------------
        if(gdict[opponents[1]]['Conference'] in conference_divisions['FBS'] or
           gdict[opponents[1]]['Conference'] in conference_divisions['FCS']):
            gdict[opponents[1]]['Record']['Season']['D1']['W'] += 1
    
        if(gdict[opponents[0]]['Conference'] in conference_divisions['FBS'] or
           gdict[opponents[0]]['Conference'] in conference_divisions['FCS']):
            gdict[opponents[0]]['Record']['Season']['D1']['L'] += 1
    
        # Update conference w/l
        # ---------------------
        if(gdict[opponents[0]]['Conference'] == gdict[opponents[1]]['Conference']):
            gdict[opponents[1]]['Record']['Season']['Conference']['W'] += 1
            gdict[opponents[0]]['Record']['Season']['Conference']['L'] += 1

# Extract the conferences
def parse_conferences(gdict, data, ops):

    # Find the index in the data where the conference standings are
    # -------------------------------------------------------------
    try:
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

    except ValueError:
        # No conference records shown. Ignore now
        print("ERROR: unable to parse conferences for",ops)

def parse_indiv_pass_stats(gdict, week, op, total_passing):

    # Figure out how many stats are included
    # --------------------------------------
    #tester = total_passing[1:8]
    beg_idx = total_passing.index('C/ATT')
    tester = total_passing[beg_idx : beg_idx + 6]
    if('QBR' in tester):
        num_vars = 6
    else:
        num_vars = 5
    
    # Parse the individual passing stats
    # ----------------------------------
    indiv_pass = total_passing[beg_idx + num_vars : -num_vars]
    #indiv_pass = total_passing[max_idx:-max_idx]

    # First, check if the team has any receiving yards
    # ------------------------------------------------ 
    if(len(indiv_pass) == 0):
        print("No passing yards for ", op)

    else:
     
        # Determine the number of passers
        # -------------------------------
        tidx1 = total_passing.index(op + ' Passing') + 1
        tidx2 = total_passing.index('team')
        passer_names = total_passing[tidx1 : tidx2]
        num_passers = len(passer_names)
        #num_passers = len(indiv_pass[::max_idx + 1])

        split_passers = np.array_split(np.array(indiv_pass), \
            num_passers)

        # Loop over the stats and add them to the dictionary
        # -------------------------------------------------- 
        for ii in range(num_passers):
            local_list = list(split_passers[ii])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]] = {}
                #'Individual']['Offense']['Passing'][local_list[0]] = {}
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]][\
                'Completions'] = float(local_list[0].split('/')[0])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]][\
                'Attempts'] = float(local_list[0].split('/')[1])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]][\
                'Yards'] = float(local_list[1])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]][\
                'Average'] = float(local_list[2])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]][\
                'TD'] = float(local_list[3])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Passing'][passer_names[ii]][\
                'INT'] = float(local_list[4])
#            gdict[opponents[0]]['Scores'][int(week)]['Individual'][local_list[0]] = {}

def parse_total_pass_yards(gdict, week, data, opponents):
    beg_idx = data.index(opponents[0] + ' Passing') 
    end_idx = data.index(opponents[1] + ' Passing') 
    total_stat = data[beg_idx : end_idx]

    # Extract the team 1 team passing stats
    # -------------------------------------
    if('QBR' in total_stat):
        op1_pyards = total_stat[-6:]
    else:
        op1_pyards = total_stat[-5:]
    ##!#try:
    ##!#    tyard_idx = total_stat.index('QBR')
    ##!#    #tyard_idx = total_stat.index('team')
    ##!#    op1_pyards = total_stat[tyard_idx + 1 : tyard_idx + 6]
    ##!#    #op1_pyards = total_stat[tyard_idx + 2]
    ##!#except ValueError:
    ##!#    op1_pyards = ['0/0','0','0','0','0']

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
    #op2_pyards = total_stat[-5:]
    if('QBR' in total_stat):
        op2_pyards = total_stat[-6:]
    else:
        op2_pyards = total_stat[-5:]
    ##!#tyard_idx = total_stat.index('QBR')
    ##!#op2_pyards = total_stat[tyard_idx + 1 : tyard_idx + 6]
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

    # NOTE HERE Continue here

    end_names = total_rushing.index('team')-1
    rush_names = total_rushing[1:end_names]
    num_rushers = len(rush_names)

    ##!## Parse the individual rushing stats
    ##!## ----------------------------------
    ##!#indiv_rush = total_rushing[6:-6]
  
    ##!## Determine the number of rushers
    ##!## -------------------------------
    ##!#num_rushers = len(indiv_rush[::7])
 
    ##!#split_rushers = np.array_split(np.array(indiv_rush), \
    ##!#    num_rushers)

    # Loop over the stats and add them to the dictionary
    # -------------------------------------------------- 
    for ii in range(num_rushers):
        local_name = rush_names[ii]
        local_list = total_rushing[end_names+7+5*ii : \
                                   end_names+7+5*(ii+1)]
        ##!#local_list = list(split_rushers[ii])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_name] = {}
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_name][\
            'Carries'] = float(local_list[0])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_name][\
            'Yards'] = float(local_list[1])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_name][\
            'Average'] = float(local_list[2])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_name][\
            'TD'] = float(local_list[3])
        gdict[op]['Scores'][int(week)][\
            'Individual']['Offense']['Rushing'][local_name][\
            'Long'] = float(local_list[4])
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

    if(total_receiving[-1] == 'No ' + op + ' Receiving'):
        print("\tNo receiving yards for", op)
    else: 

        end_names = total_receiving.index('team')
        receive_names = total_receiving[1:end_names]
        num_receivers = len(receive_names)

        ##!## Parse the individual receiving stats
        ##!## ----------------------------------
        ##!#indiv_receiving = total_receiving[6:-6]

        # First, check if the team has any receiving yards
        # ------------------------------------------------ 
        ##!#if(len(indiv_receiving) == 0):
        ##!#if(num_receivers == 0):
        ##!#    print("No receiving yards for ", op)

        ##!#else:
 
        # Determine the number of receivers
        # -------------------------------
        ##!#num_receivers = len(indiv_receiving[::7])
 
        ##!#split_receivers = np.array_split(np.array(indiv_receiving), \
        ##!#    num_receivers)

        # Loop over the stats and add them to the dictionary
        # -------------------------------------------------- 
        for ii in range(num_receivers):
            local_name = receive_names[ii]
            local_list = total_receiving[end_names+6+5*ii : \
                                         end_names+6+5*(ii+1)]
            ##!#local_list = list(split_receivers[ii])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_name] = {}
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_name][\
                'Catches'] = float(local_list[0])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_name][\
                'Yards'] = float(local_list[1])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_name][\
                'Average'] = float(local_list[2])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_name][\
                'TD'] = float(local_list[3])
            gdict[op]['Scores'][int(week)][\
                'Individual']['Offense']['Receiving'][local_name][\
                'Long'] = float(local_list[4])
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
    ##!#try:
    ##!#    beg_idx = data.index('Final')
    ##!#except ValueError:
    ##!#    try:
    ##!#        beg_idx = data.index('Final/OT')
    ##!#    except ValueError:
    ##!#        beg_idx = None
    ##!#        ot_idx = 2
    ##!#        while beg_idx is None: 
    ##!#            try:
    ##!#                beg_idx = data.index(\
    ##!#                    'Final/{0}OT'.format(ot_idx))
    ##!#            except:
    ##!#                ot_idx += 1
    ##!#            
    ##!#            if(ot_idx > 10):
    ##!#                beg_idx = 10
    
    beg_idx = data.index(' ')    
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
            w2key = int(list(gdict[tkey]['Scores'].keys())[-1])
            gdict[tkey]['Stats']['Weekly'][w2key]['scoring_offense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['rushing_offense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['passing_offense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['scoring_defense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['rushing_defense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Weekly'][w2key]['passing_defense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Pass_yards']) \
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
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['rushing_offense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['passing_offense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['scoring_defense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Points']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['rushing_defense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Rushing']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['passing_defense'] =  \
                np.nanmean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Pass_yards']) \
                for wkey in gdict[tkey]['Scores'].keys()]))
            gdict[tkey]['Stats']['Season']['total_offense'] = \
                gdict[tkey]['Stats']['Season']['passing_offense'] + \
                gdict[tkey]['Stats']['Season']['rushing_offense'] 
            gdict[tkey]['Stats']['Season']['total_defense'] = \
                gdict[tkey]['Stats']['Season']['passing_defense'] + \
                gdict[tkey]['Stats']['Season']['rushing_defense'] 

def calc_team_rankings(gdict, lkeys, division):

    for var in gdict['North Dakota']['Stats']['Season'].keys():

        # Rank scoring offense
        # --------------------
        values = [(gdict[tkey]['Stats']['Season'][var], tkey) \
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

def calc_game_quality(gdict, team, week, base_week = None,\
        debug = False, use_rankings = True):

    week = str(week)
    if(week == '1'):
        base_week = '1'
    else:
        if(base_week == None):
            base_week = str(int(week) - 1)


    opponent = gdict[team]['Scores'][week]['Opponent']

    if(debug):
        print("\nRanking game for",team," for week",week)
        print('\t',team,' vs. ',opponent)

    team_dict = gdict[team]['Scores'][week]
    if(week == '1'):
        oppt_dict = gdict[opponent]['Stats']['Weekly']['1']
        oppt_records = gdict[opponent]['Record']['Weekly']['1']
    else:
        if((gdict[opponent]['Conference'] not in \
            conference_divisions['FBS']) and \
           (gdict[opponent]['Conference'] not in \
            conference_divisions['FCS'])):
            base_week = str(week)
            
        if(base_week not in gdict[opponent]['Stats']['Weekly'].keys()):
            base_week = str(int(base_week) - 1)
        if(base_week not in gdict[opponent]['Stats']['Weekly'].keys()):
            base_week = str(int(week))
        oppt_dict = gdict[opponent]['Stats']['Weekly'][base_week]
        oppt_records = gdict[opponent]['Record']['Weekly'][base_week]

    if(use_rankings):
        if(gdict[opponent]['Conference'] in \
            conference_divisions['FCS']):
            #oppt_ranking = gdict[opponent]['Poll_rank']
            if(base_week not in \
                gdict[opponent]['Poll_rank_weekly'].keys()):
                if(debug):
                    print("ERROR: could not access weekly ranking")
                    print("       Using Poll_rank instead")
                oppt_ranking = gdict[opponent]['Poll_rank']
            else:
                oppt_ranking = gdict[opponent]['Poll_rank_weekly'][base_week]
        elif(gdict[opponent]['Conference'] in \
            conference_divisions['FBS']):
            oppt_ranking = 1
        else:
            oppt_ranking = 130

    # Conference / division factor
    # ---------------

    # If you are FBS
    if(gdict[team]['Conference'] in conference_divisions['FBS']):
        team_conf_rank = 12 - fbs_conf_rankings[gdict[team]['Conference']]
        # If both opponents are FBS
        if(gdict[opponent]['Conference'] in conference_divisions['FBS']):

            oppt_conf_rank = 12 - fbs_conf_rankings[gdict[opponent]['Conference']]
            conf_factor = oppt_conf_rank / team_conf_rank

        # If you're FBS and the other team is FCS
        elif(gdict[opponent]['Conference'] in conference_divisions['FCS']):
            team_conf_rank += 16
            oppt_conf_rank = 16 - fcs_conf_rankings[gdict[opponent]['Conference']]
            conf_factor = oppt_conf_rank / team_conf_rank

        # If the team doesn't belong to FBS or FCS   
        else:
            conf_factor = 0.001
            oppt_conf_rank = 0
 
    # If you are FCS
    if(gdict[team]['Conference'] in conference_divisions['FCS']):
        team_conf_rank = 16 - fcs_conf_rankings[gdict[team]['Conference']]
        # If your opponent is FBS
        if(gdict[opponent]['Conference'] in conference_divisions['FBS']):
            oppt_conf_rank = 12 + 16 - fbs_conf_rankings[gdict[opponent]['Conference']]
            conf_factor = oppt_conf_rank / team_conf_rank

        # If both opponents are FCS
        elif(gdict[opponent]['Conference'] in conference_divisions['FCS']):
            oppt_conf_rank = 16 - fcs_conf_rankings[gdict[opponent]['Conference']]
            conf_factor = oppt_conf_rank / team_conf_rank

        # If the team doesn't belong to FBS or FCS   
        else:
            conf_factor = 0.001
            oppt_conf_rank = 0
    else:
        print("ERROR: Conference ranking failed")

    if(debug):
        print('conf factor:',np.round(conf_factor, 3), \
            team_conf_rank, oppt_conf_rank)
            
                

    # Opponent record factor
    # ----------------------
    team_win_pcnt = np.round(gdict[team]['Record']['Weekly'][week]['Overall']['W'] / \
                    (gdict[team]['Record']['Weekly'][week]['Overall']['W']  + 
                    gdict[team]['Record']['Weekly'][week]['Overall']['L']), \
                    3)
    oppt_win_pcnt = np.round(oppt_records['Overall']['W'] / \
                   (oppt_records['Overall']['W']  + 
                    oppt_records['Overall']['L']), \
                    3)

    if(oppt_win_pcnt == 0.):
        record_factor = 10**(np.log10(0.1)/(13 - int(week)))
    else:
        record_factor = oppt_win_pcnt * 2
    if(debug):
        print("Team record:", \
            str(gdict[team]['Record']['Weekly'][week]['Overall']['W']) + '-'+\
            str(gdict[team]['Record']['Weekly'][week]['Overall']['L']) + \
            '  ' + str(team_win_pcnt))
        print("Oppt record:", \
            str(oppt_records['Overall']['W']) + '-'+\
            str(oppt_records['Overall']['L']) + \
            '  ' + str(oppt_win_pcnt) )
    #record_factor = oppt_dict['Record']['W']

    # Win size factor
    # ---------------
    if(team_dict['Team']['Defense']['Points'] == 0):
        win_size_factor = team_dict['Team']['Offense']['Points'] / \
            7.0
    elif(team_dict['Team']['Offense']['Points'] == 0):
        win_size_factor = 0.1
    else:
        win_size_factor = team_dict['Team']['Offense']['Points'] / \
                          team_dict['Team']['Defense']['Points']

    if(debug):
        print('win factor:',np.round(win_size_factor, 3),\
            team_dict['Team']['Offense']['Points'], \
            team_dict['Team']['Defense']['Points'])

    if(use_rankings):
        # Team rankings factor
        # --------------------
        if(oppt_ranking <= 25):
            if(win_size_factor > 1):
                rank_factor = np.log10(27 - oppt_ranking) * 0.5 
            else:
                rank_factor = np.log10(27 - oppt_ranking) * 0.25
        else:
            rank_factor = 0.0
    #rank_factor = 2 * np.log10(oppt_rank)


    # Home/road factor
    # ----------------
    if(win_size_factor > 1.):
        if(team_dict['Location'] == 'Away'):
            loc_factor = 1.15
        else:
            loc_factor = 1.0
    else:
        if(team_dict['Location'] == 'Away'):
            loc_factor = 1.0
        else:
            loc_factor = 0.85

    if(debug):
        print('loc factor:',loc_factor)
        
    # -----------------------
    # Calculate team factors
    # -----------------------

    # Scoring factor
    # ---------------
    # Higher value if you score more than your opponent
    # usually allows.
    if((team_dict['Team']['Offense']['Points'] == 0) | \
       (oppt_dict['scoring_defense'] == 0.0)):
        off_scoring_factor = 0.
    else:
        off_scoring_factor = team_dict['Team']['Offense']['Points'] / \
                             oppt_dict['scoring_defense']
    # Higher value if you allow fewer than your opponent
    # usually scores.
    if(team_dict['Team']['Defense']['Points'] == 0):
        def_scoring_factor = 2.
    else:
        def_scoring_factor = oppt_dict['scoring_offense'] / \
                             team_dict['Team']['Defense']['Points']
    if(debug):
        print('off score factor: ', np.round(off_scoring_factor, 3), \
            team_dict['Team']['Offense']['Points'], \
            np.round(oppt_dict['scoring_defense'], 3))
        print('def score factor: ', np.round(def_scoring_factor, 3), \
            team_dict['Team']['Defense']['Points'], \
            np.round(oppt_dict['scoring_offense'], 3))

    # Team rush yards factor
    # -----------------
    # Higher value if you rush for more yards than your 
    # opponent usually allows.
    if(team_dict['Team']['Offense']['Rushing'] == 0):
        off_rush_factor = 0.
    else:
        off_rush_factor = team_dict['Team']['Offense']['Rushing'] / \
                          oppt_dict['rushing_defense']
    # Higher value if you allow fewer rushing yards than your 
    # opponent usually rushes for.
    if(team_dict['Team']['Defense']['Rushing'] == 0):
        def_rush_factor = 2.
    else:
        def_rush_factor = oppt_dict['rushing_offense'] / \
                          team_dict['Team']['Defense']['Rushing']
    
    if(debug):                  
        print('off rush factor: ', np.round(off_rush_factor, 3), \
            team_dict['Team']['Offense']['Rushing'], \
            np.round(oppt_dict['rushing_defense'], 3))
        print('def rush factor', np.round(def_rush_factor, 3), \
            team_dict['Team']['Defense']['Rushing'], \
            np.round(oppt_dict['rushing_offense'], 3))

    # Team pass yards factor
    # -----------------
    # Higher value if you pass for more yards than your 
    # opponent usually allows.
    if(team_dict['Team']['Offense']['Pass_yards'] == 0):
        off_pass_factor = 0.
    else:
        off_pass_factor = team_dict['Team']['Offense']['Pass_yards'] / \
                          oppt_dict['passing_defense']
    # Higher value if you allow fewer passing yards than your 
    # opponent usually passes for.
    if(team_dict['Team']['Defense']['Pass_yards'] == 0):
        def_pass_factor = 2.
    else:
        def_pass_factor = oppt_dict['passing_offense'] / \
                          team_dict['Team']['Defense']['Pass_yards']
    
    if(debug):                  
        print('off pass factor', np.round(off_pass_factor, 3), \
            team_dict['Team']['Offense']['Pass_yards'], \
            np.round(oppt_dict['passing_defense'], 3))
        print('def pass factor', np.round(def_pass_factor, 3), \
            oppt_dict['passing_offense'], \
            np.round(team_dict['Team']['Defense']['Pass_yards'], 3))
 
    #stats_factor = (off_rush_factor + def_rush_factor + \
    #                off_pass_factor + def_pass_factor + \
    #                off_scoring_factor + def_scoring_factor)

    stats_factor = np.mean([off_rush_factor, 
                            def_rush_factor,
                            off_pass_factor, 
                            def_pass_factor, 
                            off_scoring_factor, 
                            def_scoring_factor]) * 0.25
    #stats_factor = (off_rush_factor + def_rush_factor + \
    #                off_pass_factor + def_pass_factor + \
    #                off_scoring_factor + def_scoring_factor)
    
    if(debug):
        print('stats factor:', np.round(stats_factor, 3))

    if(win_size_factor > 1):
        game_quality_v1 = 0.5
    else:
        game_quality_v1 = 0

    if(oppt_conf_rank == 0):
        game_quality_v1 = 0.25
    else:
        if(debug):
            print('  Oppt_conf_rank: ', np.log10(oppt_conf_rank))
            if(use_rankings):
                print('  Oppt_rank_fact: ', rank_factor)
            print('  Win_size_factor:', np.log10(win_size_factor))
            print('  Record factor:  ', np.log10(record_factor / 2.))
            print('  Loc factor:     ', (loc_factor - 1.))# + \
        game_quality_v1 += \
            np.log10(oppt_conf_rank) + \
            np.log10(win_size_factor) + \
            np.log10(record_factor / 2.) + \
            (loc_factor - 1.)# + \
        if(use_rankings):
            game_quality_v1 += rank_factor
            #stats_factor

    game_quality_v2 = \
        record_factor * \
        loc_factor * \
        win_size_factor

        #off_scoring_factor * \
        #def_scoring_factor * \
        #off_rush_factor * \
        #def_rush_factor * \
        #off_pass_factor * \
        #def_pass_factor
 
    if(debug):
        #print(record_factor, loc_factor, off_scoring_factor, def_scoring_factor, \
        #    off_rush_factor, def_rush_factor, off_pass_factor, \
        #    def_pass_factor)

        print('Game quality  v1:', np.round(game_quality_v1, 3))
        print('Game quality  v2:', np.round(game_quality_v2, 3))
    #print(record_factor *  loc_factor * off_scoring_factor * def_scoring_factor * \
    #    off_rush_factor * def_rush_factor * off_pass_factor * \
    #    def_pass_factor)

    return game_quality_v1

def simulate_whole_season(gdict, max_week, division = 'FCS', \
        bye_value = 0.9):

    # Loop over weeks
    for ii in range(1,max_week+1):
        print('Simulating week',ii)
        for team in gdict.keys():
            # If this team doesn't play this week, just
            # add a default "maintaining place" parameter
            if(gdict[team]['Conference'] in \
                    conference_divisions['FCS']):
                #print(team)
                if(str(ii) not in gdict[team]['Scores'].keys()):
                    gdict[team]['Rank_score'] += bye_value

                else:
                    # For each week, calculate the game quality
                    # for this team for this week
                    # -----------------------------------------
                    gdict[team]['Rank_score'] += \
                        calc_game_quality(gdict, team, str(ii))

        # Using the newly-developed rank scores, make new
        # rankings
        gdict = generate_new_poll_ranks(gdict, ii)

    return gdict

def generate_new_poll_ranks(gdict, week, division = 'FCS'):

    print("UPDATING POLL RANKINGS")
    team_list = {}
    for tkey in gdict.keys():
        if(gdict[tkey]['Conference'] in \
            conference_divisions[division]):
            team_list[tkey] = np.nan

    keep_teams = team_list.keys()

    rank_values = [(gdict[team]['Rank_score'], team) for \
        team in keep_teams]

    rank_values.sort()
    rank_values = rank_values[::-1]
    
    for ii, value in enumerate(rank_values):
        rank = ii + 1
        score = value[0]
        team = value[1]
        gdict[team]['Poll_rank'] = rank
        gdict[team]['Poll_rank_weekly'][str(week)] = rank
        #print(rank, team, score)

    return gdict

def calc_total_season_value(gdict, team, debug = False, \
        week = None, use_rankings = True):

    total_value = 0

    if(week is None):
        looper = gdict[team]['Scores'].keys()
    else:
        looper = np.arange(1, week + 1)
        looper = [str(loop) for loop in looper]

    for gkey in looper:
        total_value += calc_game_quality(gdict, team, str(gkey),\
            debug = debug, use_rankings = use_rankings)

    return total_value
    
# For a specified week, calculate
#def calc_weekly_stats(gdict, division = 'All', week = None):

    # Loop over the teams
    # -------------------
    #for tkey in gdict.keys():
    
        # Calculate the average points scored and allowed 
        # so far
        # -----------------------------------------------

        # Calculate the

def extract_rank_and_values(gdict, division = 'FCS', \
        week = None):    
    
    test_dict = {}

    for team in gdict.keys():
        if(gdict[team]['Conference'] in \
                conference_divisions[division]):
            test_dict[team] = {}
            test_dict[team]['value'] = np.nan
            test_dict[team]['rank'] = np.nan
    
    for tkey in test_dict.keys():
        test_dict[tkey]['value'] = gdict[tkey]['Rank_score']

    values = [(test_dict[tkey]['value'], tkey) \
        for tkey in test_dict.keys()]
    values.sort()
    values = values[::-1]
    
    for ii, value in enumerate(values):
        rank = str(ii + 1)
        name = value[1]
        conf =  gdict[value[1]]['Conference']
        test_dict[name]['rank'] = int(rank)
        #print(f'{rank:4} {name:30}{conf:35}')
        #print(ii, '\t',np.round(value[0], 2), '\t\t',\
        #    value[1], '\t\t\t',gdict[value[1]]['Conference'])

    return test_dict, values

def generate_division_rankings(gdict, division = 'FCS', \
        week = None, use_rankings = True):    
    test_dict = {}
    
    for team in gdict.keys():
        if(gdict[team]['Conference'] in \
                conference_divisions[division]):
            test_dict[team] = {}
            test_dict[team]['value'] = np.nan
            test_dict[team]['rank'] = np.nan
    
    for tkey in test_dict.keys():
        #print(tkey)
        try:
            total_value = calc_total_season_value(gdict, tkey, \
                use_rankings = use_rankings)
            test_dict[tkey]['value'] = total_value
        except:
            print("Skipping",tkey)
    
    values = [(test_dict[tkey]['value'], tkey) \
        for tkey in test_dict.keys()]
    values.sort()
    values = values[::-1]
    
    for ii, value in enumerate(values):
        rank = str(ii + 1)
        name = value[1]
        conf =  gdict[value[1]]['Conference']
        test_dict[name]['rank'] = int(rank)
        print(f'{rank:4} {name:30}{conf:35}')
        #print(ii, '\t',np.round(value[0], 2), '\t\t',\
        #    value[1], '\t\t\t',gdict[value[1]]['Conference'])

    return test_dict, values
 
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

def set_preseason_rankings(gdict, rank_file = massey_file, \
        division = 'FCS'):
    
    massey_ranks = read_massey_ranks(infile = massey_file)

    for team in gdict.keys():
        if(gdict[team]['Conference'] in \
            conference_divisions[division]):

            massey_value = massey_ranks[team]
            begin_score = (131 - massey_value) / 5.9                   
 
            gdict[team]['Poll_rank'] = \
                massey_value
            gdict[team]['Poll_rank_weekly'] = {}
            gdict[team]['Poll_rank_weekly']['1'] = \
                massey_value
            gdict[team]['Rank_score'] = begin_score

    return gdict

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

def set_weekly_records(gdict, opponents, week):
    gdict[opponents[0]]['Record']['Weekly'][int(week)]['Overall']['W'] =  \
        gdict[opponents[0]]['Record']['Season']['Overall']['W']
    gdict[opponents[0]]['Record']['Weekly'][int(week)]['Overall']['L'] =  \
        gdict[opponents[0]]['Record']['Season']['Overall']['L']
    gdict[opponents[0]]['Record']['Weekly'][int(week)]['D1']['W'] =  \
        gdict[opponents[0]]['Record']['Season']['D1']['W']
    gdict[opponents[0]]['Record']['Weekly'][int(week)]['D1']['L'] =  \
        gdict[opponents[0]]['Record']['Season']['D1']['L']
    gdict[opponents[0]]['Record']['Weekly'][int(week)]['Conference']['W'] =  \
        gdict[opponents[0]]['Record']['Season']['Conference']['W']
    gdict[opponents[0]]['Record']['Weekly'][int(week)]['Conference']['L'] =  \
        gdict[opponents[0]]['Record']['Season']['Conference']['L']

    gdict[opponents[1]]['Record']['Weekly'][int(week)]['Overall']['W'] =  \
        gdict[opponents[1]]['Record']['Season']['Overall']['W']
    gdict[opponents[1]]['Record']['Weekly'][int(week)]['Overall']['L'] =  \
        gdict[opponents[1]]['Record']['Season']['Overall']['L']
    gdict[opponents[1]]['Record']['Weekly'][int(week)]['D1']['W'] =  \
        gdict[opponents[1]]['Record']['Season']['D1']['W']
    gdict[opponents[1]]['Record']['Weekly'][int(week)]['D1']['L'] =  \
        gdict[opponents[1]]['Record']['Season']['D1']['L']
    gdict[opponents[1]]['Record']['Weekly'][int(week)]['Conference']['W'] =  \
        gdict[opponents[1]]['Record']['Season']['Conference']['W']
    gdict[opponents[1]]['Record']['Weekly'][int(week)]['Conference']['L'] =  \
        gdict[opponents[1]]['Record']['Season']['Conference']['L']

def check_missing_games(gdict, espn_url, week, \
        this_week_teams):

    # Read in the base ESPN url data using the parser
    # -----------------------------------------------
    parse_espn = pull_ESPN_data(espn_url)

    # Look through the list and pull out the team names
    # -------------------------------------------------
    missed_teams = []
    for ii in range(len(parse_espn)):      
        if(parse_espn[ii] == '('):
            team_name = parse_espn[ii - 1]

            # Check if the team has data for this week
            if(team_name not in this_week_teams):
   
                print(team_name) 
                # Team missed 
                missed_teams.append(team_name)

    if(len(missed_teams) != 0):

        print(missed_teams)
        if(len(missed_teams) % 2 == 1):
            missed_teams = missed_teams[:-1]
        # Assume that the missed teams will be in pairs
        missed_pairs = [(missed_teams[ii], missed_teams[ii+1]) \
            for ii in np.arange(0, len(missed_teams), 2)]

        for ii in range(len(missed_pairs)):

            opponents = missed_pairs[ii]
            # Check if the team is even in the dictionary?
            # --------------------------------------------
            for op in opponents:
                set_team_dict(gdict, op, week)
                ##!#if(op not in gdict.keys()):
                ##!#    gdict[op] = {}
                ##!#    gdict[op]['Beat']    = []
                ##!#    gdict[op]['Lost to'] = []
                ##!#    gdict[op]['Stats'] = {}
                ##!#    gdict[op]['Stats']['Season'] = {}
                ##!#    gdict[op]['Record'] = {}
                ##!#    gdict[op]['Record']['Season'] = {}
                ##!#    gdict[op]['Record']['Weekly'] = {}
                ##!#    gdict[op]['Record']['Season']['Overall'] = {}
                ##!#    gdict[op]['Record']['Season']['D1'] = {}
                ##!#    gdict[op]['Record']['Season']['Conference'] = {}
                ##!#    gdict[op]['Record']['Season']['Overall']['W'] = 0
                ##!#    gdict[op]['Record']['Season']['Overall']['L'] = 0
                ##!#    gdict[op]['Record']['Season']['D1']['W'] = 0
                ##!#    gdict[op]['Record']['Season']['D1']['L'] = 0
                ##!#    gdict[op]['Record']['Season']['Conference']['W'] = 0
                ##!#    gdict[op]['Record']['Season']['Conference']['L'] = 0
                ##!#    gdict[op]['Stats']['Weekly'] = {}
                ##!#    gdict[op]['Scores'] = {}
                ##!#    gdict[op]['Rankings'] = {}
                ##!#    gdict[op]['Rankings']['All'] = {}
                ##!#    gdict[op]['Rankings']['Division'] = {}
                ##!#    gdict[op]['Rankings']['Conference'] = {}
                ##!#    gdict[op]['Conference'] = 'NONE' 
            
                ##!#gdict[op]['Record']['Weekly'][int(week)] = {}
                ##!#gdict[op]['Record']['Weekly'][int(week)]['Overall'] = {}
                ##!#gdict[op]['Record']['Weekly'][int(week)]['D1'] = {}
                ##!#gdict[op]['Record']['Weekly'][int(week)]['Conference'] = {}
                ##!#gdict[op]['Record']['Weekly'][int(week)]['Overall']['W'] = 0
                ##!#gdict[op]['Record']['Weekly'][int(week)]['Overall']['L'] = 0
                ##!#gdict[op]['Record']['Weekly'][int(week)]['D1']['W'] = 0
                ##!#gdict[op]['Record']['Weekly'][int(week)]['D1']['L'] = 0
                ##!#gdict[op]['Record']['Weekly'][int(week)]['Conference']['W'] = 0
                ##!#gdict[op]['Record']['Weekly'][int(week)]['Conference']['L'] = 0
                ##!#gdict[op]['Scores'][int(week)] = {}
                ##!#gdict[op]['Scores'][int(week)]['Team'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Team']['Offense'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Team']['Defense'] = {}

                ##!#gdict[op]['Scores'][int(week)]['Individual'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Individual']['Offense'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Individual']['Offense']['Rushing'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Individual']['Offense']['Passing'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Individual']['Offense']['Receiving'] = {}
                ##!#gdict[op]['Scores'][int(week)]['Individual']['Defense'] = {}
            
            ## Probably don't need to have run each time, but can fix later 
            #parse_conferences(gdict, data, opponents)

            # Set the home / away locations
            # -----------------------------
            gdict[opponents[0]]['Scores'][int(week)]['Opponent'] = opponents[1] 
            gdict[opponents[0]]['Scores'][int(week)]['Location'] = 'Away' 
            
            gdict[opponents[1]]['Scores'][int(week)]['Opponent'] = opponents[0] 
            gdict[opponents[1]]['Scores'][int(week)]['Location'] = 'Home' 
            
            # Parse the score 
            # NOTE: contains the quarterly score info also
            # ---------------------------------------------
            team1_idx = parse_espn.index(opponents[0])
            team2_idx = parse_espn.index(opponents[1])
            print(team1_idx, team2_idx)
            if(team2_idx - team1_idx < 7):
                print("Game cancelled")
                gdict[opponents[0]]['Scores'].pop(int(week), None)
                gdict[opponents[1]]['Scores'].pop(int(week), None)

                # Propagate the records to the current week
                if(week > 1):
                    gdict[opponents[0]]['Record']['Weekly'][week] = \
                        gdict[opponents[0]]['Record']['Weekly'][\
                        str(int(week)-1)]
                    gdict[opponents[1]]['Record']['Weekly'][week] = \
                        gdict[opponents[1]]['Record']['Weekly'][\
                        str(int(week)-1)]

                    gdict[opponents[0]]['Stats']['Weekly'][week] = \
                        gdict[opponents[0]]['Stats']['Weekly'][\
                        str(int(week)-1)]
                    gdict[opponents[1]]['Stats']['Weekly'][week] = \
                        gdict[opponents[1]]['Stats']['Weekly'][\
                        str(int(week)-1)]

                    
            else:
                if(team2_idx < team1_idx):
                    team2_idx = team1_idx + 11
                else:
                    op1_points = parse_espn[team1_idx + 6 : team1_idx + 11]
                    op2_points = parse_espn[team2_idx + 6 : team2_idx + 11]

                    print(opponents)
                    update_scores(gdict, opponents, week, op1_points, op2_points)

                # Set the weekly record values
                # ----------------------------
                set_weekly_records(gdict, opponents, week)

                # Set the stats for these teams to nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_completions'] = np.nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_attempts']    = np.nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_yards']       = np.nan
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_completions'] = np.nan
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_attempts']    = np.nan
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_yards']       = np.nan
                                                                                                  
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_completions'] = np.nan
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_attempts']    = np.nan
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_yards']       = np.nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_completions'] = np.nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_attempts']    = np.nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_yards']       = np.nan

                gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Rushing'] = np.nan
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Rushing'] = np.nan
                                                                                         
                gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Rushing'] = np.nan
                gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Rushing'] = np.nan

def set_team_dict(gdict, op, week):
    if(op not in gdict.keys()):
        gdict[op] = {}
        gdict[op]['Beat']    = []
        gdict[op]['Lost to'] = []
        gdict[op]['Stats'] = {}
        gdict[op]['Stats']['Season'] = {}
        gdict[op]['Record'] = {}
        gdict[op]['Record']['Season'] = {}
        gdict[op]['Record']['Weekly'] = {}
        gdict[op]['Record']['Season']['Overall'] = {}
        gdict[op]['Record']['Season']['D1'] = {}
        gdict[op]['Record']['Season']['Conference'] = {}
        gdict[op]['Record']['Season']['Overall']['W'] = 0
        gdict[op]['Record']['Season']['Overall']['L'] = 0
        gdict[op]['Record']['Season']['D1']['W'] = 0
        gdict[op]['Record']['Season']['D1']['L'] = 0
        gdict[op]['Record']['Season']['Conference']['W'] = 0
        gdict[op]['Record']['Season']['Conference']['L'] = 0
        gdict[op]['Stats']['Weekly'] = {}
        gdict[op]['Scores'] = {}
        gdict[op]['Rankings'] = {}
        gdict[op]['Rankings']['All'] = {}
        gdict[op]['Rankings']['Division'] = {}
        gdict[op]['Rankings']['Conference'] = {}
        gdict[op]['Conference'] = 'NONE' 
        gdict[op]['Poll_rank']  = 0.
        gdict[op]['Poll_rank_weekly'] = {} 
        gdict[op]['Rank_score'] = 0.

    gdict[op]['Stats']['Weekly'][int(week)] = {}
    gdict[op]['Record']['Weekly'][int(week)] = {}
    gdict[op]['Record']['Weekly'][int(week)]['Overall'] = {}
    gdict[op]['Record']['Weekly'][int(week)]['D1'] = {}
    gdict[op]['Record']['Weekly'][int(week)]['Conference'] = {}
    gdict[op]['Record']['Weekly'][int(week)]['Overall']['W'] = 0
    gdict[op]['Record']['Weekly'][int(week)]['Overall']['L'] = 0
    gdict[op]['Record']['Weekly'][int(week)]['D1']['W'] = 0
    gdict[op]['Record']['Weekly'][int(week)]['D1']['L'] = 0
    gdict[op]['Record']['Weekly'][int(week)]['Conference']['W'] = 0
    gdict[op]['Record']['Weekly'][int(week)]['Conference']['L'] = 0
    gdict[op]['Scores'][int(week)] = {}
    gdict[op]['Scores'][int(week)]['Team'] = {}
    gdict[op]['Scores'][int(week)]['Team']['Offense'] = {}
    gdict[op]['Scores'][int(week)]['Team']['Defense'] = {}

    gdict[op]['Scores'][int(week)]['Individual'] = {}
    gdict[op]['Scores'][int(week)]['Individual']['Offense'] = {}
    gdict[op]['Scores'][int(week)]['Individual']['Offense']['Rushing'] = {}
    gdict[op]['Scores'][int(week)]['Individual']['Offense']['Passing'] = {}
    gdict[op]['Scores'][int(week)]['Individual']['Offense']['Receiving'] = {}
    gdict[op]['Scores'][int(week)]['Individual']['Defense'] = {}

def load_ESPN_data(start_week, end_week, year = 2022):

    gdict = {}
    weeks = np.arange(start_week,end_week + 1)
    espn_urls = ["https://www.espn.com/college-football/scoreboard/_/week/{0}/year/" + str(year) + "/seasontype/2/group/81", \
                 "https://www.espn.com/college-football/scoreboard/_/week/{0}/year/" + str(year) + "/seasontype/2/group/80"]

    # Loop over all the weeks and add the data
    # ----------------------------------------
    for week in weeks:
        this_week_teams = []
        for espn_url in espn_urls:    
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
                #print(local_url)
                data = pull_ESPN_data(local_url)

                # Parse the opponents
                # -------------------
                matchup = str(data[1].split(' - College Football Box Score - ')[0])
                matchup = unidecode(matchup)
                print(matchup)
                opponents = matchup.split(' vs. ')
                #if(len(opponents[1].split(' - ')) > 1):
                #    opponents[1] = opponents[1].split(' - ')[0]
                #print(opponents)
                #if(opponents[0][:5] == 'Hawai'):
                #    opponents[0] = 'Hawaii'
                #elif(opponents[1][:5] == 'Hawai'):
                #    opponents[1] = 'Hawaii'

                # Account for FCS vs FBS games. In these cases, 
                # the game would be recorded twice for each
                # of the opponents, screwing up the win totals
                # and "teams beat" and "teams lost to" records.
                # To fix this, don't proceed if either of these
                # teams already have data for this week in
                # the dictionary.
                if(opponents[0] in gdict.keys()):
                    if(int(week) in gdict[opponents[0]]['Scores'].keys()):
                        print("ERROR for game:",opponents)
                        print("   This game has already been recorded")
                        continue
                if(opponents[1] in gdict.keys()):
                    if(int(week) in gdict[opponents[1]]['Scores'].keys()):
                        print("ERROR for game:",opponents)
                        print("   This game has already been recorded")
                        continue

                this_week_teams.append(opponents[0])
                this_week_teams.append(opponents[1])
                for op in opponents:
                    set_team_dict(gdict, op, week)
        
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

                update_scores(gdict, opponents, week, op1_points, op2_points)

                # Set the weekly record values
                # ----------------------------
                set_weekly_records(gdict, opponents, week)

                # Parse the total passing yards 
                parse_total_pass_yards(gdict, week, data, opponents)
             
                # Parse the total rushing yards 
                parse_total_rush_yards(gdict, week, data, opponents)

                # Parse the total receiving yards
                parse_total_receiving_yards(gdict, week, data, opponents)

            # Check for games that did not have a box score
            # ---------------------------------------------
            check_missing_games(gdict, url, week, this_week_teams)
 
            # If here, calculates after each division (FBS/FCS)
            #calc_team_stats(gdict, weekly = True)

        # If here, calculates stats after each week
        calc_team_stats(gdict, weekly = True)

    # If here, calculates stats after all weeks 
    calc_team_stats(gdict)

    return gdict #, data

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
        conference_dict[cfc] = list(teams[np.where(cfcs == cfc)])

    with open(conference_file_name,'w') as fout:
        json.dump(conference_dict, fout, indent = 4, sort_keys = True)
    print("Saved file", conference_file_name)

def read_massey_ranks(infile = massey_file):
    massey_dict = {}

    with open(infile, 'r') as fin:
        for line in fin:
            #print(line.strip()) 
            testline = line.strip().split(',')
            massey_dict[testline[1]] = int(testline[0])
            #print(testline)

    return massey_dict

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

