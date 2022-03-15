#!/usr/bin/env python

"""
NAME:

PURPOSE:

MODIFICATIONS:

"""

VERBOSE = True

def predictRound(teams,odds,seeds, cind_name):
    out_results = []
    for i,k in zip(teams[0::2],teams[1::2]):
        i_odds = odds[seeds[i]-1]
        k_odds = odds[seeds[k]-1]

        # Re-scale the odds for the cinderella team. Increase their odds
        # by 50%
        # --------------------------------------------------------------
        if(i == cind_name):
            i_odds = i_odds * (seeds[i] / 1.8)
            if(i_odds > 0.999):
                i_odds = 0.999
                k_odds = 0.001
            else:
                k_odds = 1.0 - i_odds


        elif(k == cind_name): 
            #k_odds = k_odds * 2.0
            k_odds = k_odds * (seeds[k] / 1.8)
            if(k_odds > 0.999):
                k_odds = 0.999
                i_odds = 0.001
            else:
                i_odds = 1.0 - k_odds

            
        resulti = int(random.random()*100000.*i_odds)
        resultk = int(random.random()*100000.*k_odds)

        upset_alert = ""
        if(resulti > resultk):
            winner = i
            win_seed = seeds[i]
            loser = k
            loser_seed = seeds[k]
            out_results.append(i)
        else:
            winner = k
            win_seed = seeds[k]
            loser = i
            loser_seed = seeds[i]
            out_results.append(k)
        #if((win_seed-loser_seed>2) | (loser_seed==1) | (loser_seed==2)):
        if((win_seed - loser_seed > 2)):
            upset_alert = "\t\tUPSET ALERT"
        if(VERBOSE):
            print(str(win_seed) + ' '+winner+' defeats '+str(loser_seed)+' '+loser +upset_alert)

        if((i == cind_name ) | (k == cind_name)):
            print('\t',i, resulti, ' vs ', k, resultk)
    


    return out_results

def main(round64_teams,round64_odds,round32_odds,sweet16_odds,elite8_odds,\
         final4_odds,championship_odds,seeds):

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Pick a Cinderella team
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    l_found = False
    while(not l_found):
        cinderella = int(np.random.random() * len(round64_teams))
        cind_seed = seeds[round64_teams[cinderella]]
        if(cind_seed >= 8):
            l_found = True
            cind_name = round64_teams[cinderella]
            print("Cinderella team:", cind_name, ' Seed = ', cind_seed)

    
    # Simulate first round
    if(VERBOSE): print("\nRound of 64")
    round32_teams = predictRound(round64_teams,round64_odds,seeds, cind_name)
    if(VERBOSE): print("\nRound of 32")
    sweet16_teams = predictRound(round32_teams,round32_odds,seeds, cind_name)
    if(VERBOSE): print("\nSweet 16")
    elite8_teams = predictRound(sweet16_teams,sweet16_odds,seeds, cind_name)
    if(VERBOSE): print("\nElite 8")
    final4_teams = predictRound(elite8_teams,elite8_odds,seeds, cind_name)
    if(VERBOSE): print("\nFinal 4")
    championship_teams = predictRound(final4_teams,final4_odds,seeds, cind_name)
    if(VERBOSE): print("\nChampionship")
    winner = predictRound(championship_teams,championship_odds,seeds, cind_name)
    
    win_seed = seeds[winner[0]]
    print(str(win_seed) + ' ' +winner[0]+ ' wins the national championship')

    stop_point = False
    if(cind_name in elite8_teams):
        stop_point = True

    return(win_seed, stop_point)

import numpy as np
import sys
import random

round64_odds      = np.array([0.993,0.943,0.850,0.793,0.643,0.629,0.607,0.486,0.514,0.393,0.371,0.357,0.207,0.150,0.057,0.007])
round32_odds      = np.array([0.857,0.636,0.529,0.471,0.336,0.300,0.193,0.093,0.050,0.164,0.157,0.150,0.043,0.014,0.007,0.001])
sweet16_odds      = np.array([0.693,0.457,0.257,0.150,0.064,0.100,0.071,0.057,0.029,0.057,0.057,0.007,0.001,0.001,0.001,0.001])
elite8_odds       = np.array([0.407,0.207,0.121,0.093,0.050,0.021,0.021,0.036,0.007,0.007,0.029,0.001,0.001,0.001,0.001,0.001])
final4_odds       = np.array([0.243,0.093,0.079,0.021,0.021,0.014,0.007,0.021,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])
championship_odds = np.array([0.157,0.036,0.029,0.007,0.001,0.007,0.007,0.007,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Predict the play-in game results
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
temp_odds = random.random()*10.
#play_in1_16 = 1
if(temp_odds>5):
    winner_1 = 'Notre Dame'
else:
    winner_1 = 'Rutgers'
print(winner_1+' wins the first play-in game')
#play_in2_11 = 10
temp_odds = random.random()*10.
if(temp_odds>5):
    winner_2 = 'Wyoming'
else:
    winner_2 = 'Indiana'
print(winner_2+' wins the second play-in game')
#play_in3_16 = 18
temp_odds = random.random()*10.
if(temp_odds>5):
    winner_3 = 'Wright St.'
else:
    winner_3 = 'Bryant'
print(winner_3+' wins the third play-in game')
#play_in4_11 = 26
temp_odds = random.random()*10.
if(temp_odds>5):
    winner_4 = 'Texas Southern'
else:
    winner_4 = 'Texas A&M - CC'
print(winner_4+' wins the fourth play-in game\n\n')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Define the March Madness bracket teams
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
seeds = dict()
seeds['Gonzaga'] = 1
seeds['Georgia St.'] = 16
seeds['Boise St.'] = 8
seeds['Memphis'] = 9
seeds['UConn'] = 5
seeds['New Mexico St.'] = 12
seeds['Arkansas'] = 4
seeds['Vermont'] = 13
seeds['Alabama'] = 6
seeds[winner_1] = 11
seeds['Texas Tech'] = 3
seeds['Montana St.'] = 14
seeds['Michigan St.'] = 7
seeds['Minnesota'] = 10
seeds['Duke'] = 2
seeds['CS Fullerton'] = 15

seeds['Baylor'] = 1
seeds['Norfolk St.'] = 16 # NC Central / Texas Southern
seeds['North Carolina'] = 8
seeds['Marquette'] = 9
seeds['Saint Mary\'s'] = 5
seeds[winner_2] = 12
seeds['UCLA'] = 4
seeds['Akron'] = 13
seeds['Texas'] = 6
seeds['Virginia Tech'] = 11
seeds['Purdue'] = 3
seeds['Yale'] = 14
seeds['Murray St.'] = 7
seeds['San Francisco'] = 10
seeds['Kentucky'] = 2
seeds['Saint Peter\'s'] = 15

seeds['Arizona'] = 1
seeds[winner_3] = 16
seeds['Seton Hall'] = 8
seeds['TCU'] = 9
seeds['Houston'] = 5
seeds['UAB'] = 12
seeds['Illinois'] = 4
seeds['Chattanooga'] = 13
seeds['Colorado St.'] = 6
seeds['Michigan'] = 11
seeds['Tennessee'] = 3
seeds['Longwood'] = 14
seeds['Ohio St.'] = 7
seeds['Loyola Chicago'] = 10
seeds['VIllanova'] = 2
seeds['Delaware'] = 15

seeds['Kansas'] = 1
seeds[winner_4] = 16
seeds['San Diego State'] = 8
seeds['Creighton'] = 9
seeds['Iowa'] = 5
seeds['Richmond'] = 12
seeds['Providence'] = 4
seeds['South Dakota State'] = 13
seeds['LSU'] = 6
seeds['Iowa State'] = 11
seeds['Wisconsin'] = 3
seeds['Colgate'] = 14
seeds['USC'] = 7
seeds['Miami (FL)'] = 10
seeds['Auburn'] = 2
seeds['Jacksonville St.'] = 15

round64_teams = list(seeds.keys())

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Run the model to predict each round's results
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
seed, _ =main(round64_teams,round64_odds,round32_odds,sweet16_odds,elite8_odds,\
     final4_odds,championship_odds,seeds)

sys.exit()

count = 0
stop_found = False
while(not stop_found):
    count+=1
    seed, stop_found = main(round64_teams,round64_odds,round32_odds,sweet16_odds,elite8_odds,\
         final4_odds,championship_odds,seeds)
print("Took",count,"years for Cinderella team to reach Final Four")
sys.exit()



count=0
seed=17
wanted=int(sys.argv[-1])
while(seed!=wanted):
    count+=1
    seed=main(round64_teams,round64_odds,round32_odds,sweet16_odds,elite8_odds,\
         final4_odds,championship_odds,seeds)

print("Took",count,"years for",seed,"seed to win")
