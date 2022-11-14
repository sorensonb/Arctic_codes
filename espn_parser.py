#!/usr/bin/env python

"""


"""

import urllib.request
import requests
from html.parser import HTMLParser
import re
import numpy as np
import sys

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

            gdict[op][int(week)] = {}
            gdict[op][int(week)]['Team'] = {}
            gdict[op][int(week)]['Team']['Offense'] = {}
            gdict[op][int(week)]['Team']['Offense']['Passing'] = {}
            gdict[op][int(week)]['Team']['Offense']['Rushing'] = {}
            gdict[op][int(week)]['Team']['Offense']['Receiving'] = {}
    
        # Parse the score   
        # ---------------
        
        # Parse the total passing yards 
        # NOTE: can modify to extract player stats
        # ----------------------------------------
    
        # Team 1
        beg_idx = data.index(opponents[0] + ' Passing') 
        end_idx = data.index(opponents[1] + ' Passing') 
        total_stat = data[beg_idx : end_idx]
        gdict[opponents[0]][int(week)]['Team']['Offense']['Passing'] = total_stat[-4]
    
        # Team 2
        beg_idx = data.index(opponents[1] + ' Passing') 
        end_idx = data.index(opponents[0] + ' Rushing') 
        total_stat = data[beg_idx : end_idx]
        gdict[opponents[1]][int(week)]['Team']['Offense']['Passing'] = total_stat[-4]
    
        # Parse the total rushing yards 
        # NOTE: can modify to extract player stats
        # ----------------------------------------
    
        # Team 1
        beg_idx = data.index(opponents[0] + ' Rushing') 
        end_idx = data.index(opponents[1] + ' Rushing') 
        total_rushing = data[beg_idx : end_idx]
        gdict[opponents[0]][int(week)]['Team']['Offense']['Rushing'] = total_rushing[-4]
    
        # Team 2
        beg_idx = data.index(opponents[1] + ' Rushing') 
        end_idx = data.index(opponents[0] + ' Receiving') 
        total_rushing = data[beg_idx : end_idx]
        gdict[opponents[1]][int(week)]['Team']['Offense']['Rushing'] = total_rushing[-4]
         
    
        # Parse the total receiving yards
        # NOTE: can modify to extract player stats
        # ----------------------------------------
        
        # Team 1
        beg_idx = data.index(opponents[0] + ' Receiving') 
        end_idx = data.index(opponents[1] + ' Receiving') 
        total_stat = data[beg_idx : end_idx]
        gdict[opponents[0]][int(week)]['Team']['Offense']['Receiving'] = total_stat[-4]
    
        # Team 2
        beg_idx = data.index(opponents[1] + ' Receiving') 
        end_idx = data.index(opponents[0] + ' Fumbles') 
        total_stat = data[beg_idx : end_idx]
        gdict[opponents[1]][int(week)]['Team']['Offense']['Receiving'] = total_stat[-4]
    





