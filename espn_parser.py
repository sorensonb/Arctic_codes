#!/usr/bin/env python

"""


"""

import urllib.request
import requests
from html.parser import HTMLParser
import re
import numpy as np
import sys

conferences = {
    'Abilene Christian': 'ASUN-WAC',
    'Alabama A&M': 'SWAC',  
    'Alabama State': 'SWAC',
    'Albany': 'CAA',
    'Alcorn State': 'SWAC',
    'Arkansas-Pine Bluff': 'SWAC',
    'Austin Peay': 'ASUN-WAC',
    'Bethune-Cookman': 'SWAC',
    'Brown': 'Ivy',
    'Bryant': 'NEC',
    'Bucknell': 'PL',
    'Butler': 'PFL',
    'Cal Poly': 'BSC',
    'Campbell': 'BSo',
    'Central Arkansas': 'SLC',
    'Central Connecticut': 'NEC',
    'Charleston Southern': 'BSo',
    'Chattanooga': 'SC',
    'Colgate': 'PL',
    'Columbia': 'Ivy',
    'Cornell': 'Ivy',
    'Dartmouth': 'Ivy',
    'Davidson': 'PFL',
    'Dayton': 'PFL',
    'Delaware': 'CAA',
    'Delaware State': 'MEAC',
    'Drake': 'PFL',
    'Duquesne': 'NEC',
    'East Tennessee State': 'SC',
    'Eastern Illinois': 'OVC',
    'Eastern Kentucky': 'ASUN-WAC',
    'Eastern Washington': 'BSC',
    'Elon': 'CAA',
    'Florida A&M': 'SWAC',
    'Fordham': 'PL',
    'Furman': 'SC',
    'Gardner-Webb': 'BSo',
}

def parse_total_pass_yards(data, opponents):
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
    final1[0] = float(op1_pyards[0].split('/')[0])
    final1[1] = float(op1_pyards[0].split('/')[1])
    final1[2] = float(op1_pyards[1])
    final1[3] = float(op1_pyards[2])
    final1[4] = float(op1_pyards[3])
    final1[5] = float(op1_pyards[4])

    beg_idx = data.index(opponents[1] + ' Passing') 
    end_idx = data.index(opponents[0] + ' Rushing') 
    total_stat = data[beg_idx : end_idx]
    tyard_idx = total_stat.index('TEAM')
    op2_pyards = total_stat[tyard_idx + 1 : tyard_idx + 6]
    #op2_pyards = total_stat[tyard_idx + 2]

    final2 = np.zeros(len(op2_pyards) + 1) 
    final2[0] = float(op2_pyards[0].split('/')[0])
    final2[1] = float(op2_pyards[0].split('/')[1])
    final2[2] = float(op2_pyards[1])
    final2[3] = float(op2_pyards[2])
    final2[4] = float(op2_pyards[3])
    final2[5] = float(op2_pyards[4])

    return final1, final2

def parse_total_rush_yards(data, opponents):
    # Team 1
    beg_idx = data.index(opponents[0] + ' Rushing') 
    end_idx = data.index(opponents[1] + ' Rushing') 
    total_rushing = data[beg_idx : end_idx]
    op1_ryards = total_rushing[-4]

    # Team 2
    beg_idx = data.index(opponents[1] + ' Rushing') 
    end_idx = data.index(opponents[0] + ' Receiving') 
    total_rushing = data[beg_idx : end_idx]
    op2_ryards = total_rushing[-4]

    return op1_ryards, op2_ryards

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

def calc_team_stats(gdict):

    for tkey in gdict.keys():
        gdict[tkey]['Stats']['scoring_offense'] =  \
            np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Points']) \
            for wkey in gdict[tkey]['Scores'].keys()]))
        gdict[tkey]['Stats']['rushing_offense'] =  \
            np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Rushing']) \
            for wkey in gdict[tkey]['Scores'].keys()]))
        gdict[tkey]['Stats']['passing_offense'] =  \
            np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Offense']['Pass_yards']) \
            for wkey in gdict[tkey]['Scores'].keys()]))
        gdict[tkey]['Stats']['scoring_defense'] =  \
            np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Points']) \
            for wkey in gdict[tkey]['Scores'].keys()]))
        gdict[tkey]['Stats']['rushing_defense'] =  \
            np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Rushing']) \
            for wkey in gdict[tkey]['Scores'].keys()]))
        gdict[tkey]['Stats']['passing_defense'] =  \
            np.mean(np.array([float(gdict[tkey]['Scores'][wkey]['Team']['Defense']['Pass_yards']) \
            for wkey in gdict[tkey]['Scores'].keys()]))

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

def load_ESPN_data():

    gdict = {}
    espn_url = "https://www.espn.com/college-football/scoreboard/_/week/{0}/year/2022/seasontype/2/group/81"
    weeks = np.arange(1,12)
    
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
            uf = urllib.request.urlopen(local_url)
            html = uf.read()
            parser = MyHTMLParser()
            parser.feed(str(html))
        
            data = parser.Myrawdata.strip().split('\n')
        
            # Parse the opponents
            # -------------------
            matchup = data[1].split(' - Box Score - ')[0]
            print(matchup)
            opponents = matchup.split(' vs. ')
            for op in opponents:
                if(op not in gdict.keys()):
                    gdict[op] = {}
                    gdict[op]['Stats'] = {}
                    gdict[op]['Scores'] = {}
    
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
       
            # Set the home / away locations
            # -----------------------------
            gdict[opponents[0]]['Scores'][int(week)]['Opponent'] = opponents[1] 
            gdict[opponents[0]]['Scores'][int(week)]['Location'] = 'Away' 
            
            gdict[opponents[1]]['Scores'][int(week)]['Opponent'] = opponents[0] 
            gdict[opponents[1]]['Scores'][int(week)]['Location'] = 'Home' 
     
            # Parse the score   
            # ---------------
            op1_points, op2_points = parse_scores(data, opponents)
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Points'] = op1_points[-1]
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Points'] = op2_points[-1]
    
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Points'] = op1_points[-1]
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Points'] = op2_points[-1]
    
            # Parse the total passing yards 
            # NOTE: can modify to extract player stats
            # ----------------------------------------
            op1_pyards, op2_pyards = parse_total_pass_yards(data, opponents)
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_completions'] = op1_pyards[0]
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_attempts']    = op1_pyards[1]
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Pass_yards'] = op1_pyards[2]
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_completions'] = op2_pyards[0]
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_attempts']    = op2_pyards[1]
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Pass_yards'] = op2_pyards[2]
    
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_completions'] = op1_pyards[0]
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_attempts']    = op1_pyards[1]
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Pass_yards'] = op1_pyards[2]
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_completions'] = op2_pyards[0]
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_attempts']    = op2_pyards[1]
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Pass_yards'] = op2_pyards[2]
         
            # Parse the total rushing yards 
            # NOTE: can modify to extract player stats
            # ----------------------------------------
            op1_ryards, op2_ryards = parse_total_rush_yards(data, opponents)
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Offense']['Rushing'] = op1_ryards
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Offense']['Rushing'] = op2_ryards
    
            gdict[opponents[1]]['Scores'][int(week)]['Team']['Defense']['Rushing'] = op1_ryards
            gdict[opponents[0]]['Scores'][int(week)]['Team']['Defense']['Rushing'] = op2_ryards
        
            ## Parse the total receiving yards
            ## NOTE: can modify to extract player stats
            ## ----------------------------------------
            #
            ## Team 1
            #beg_idx = data.index(opponents[0] + ' Receiving') 
            #end_idx = data.index(opponents[1] + ' Receiving') 
            #total_stat = data[beg_idx : end_idx]
            #gdict[opponents[0]][int(week)]['Team']['Offense']['Receiving'] = total_stat[-4]
        
            ## Team 2
            #beg_idx = data.index(opponents[1] + ' Receiving') 
            #end_idx = data.index(opponents[0] + ' Fumbles') 
            #total_stat = data[beg_idx : end_idx]
            #gdict[opponents[1]][int(week)]['Team']['Offense']['Receiving'] = total_stat[-4]
        
    
        calc_team_stats(gdict)



