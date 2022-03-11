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
    winner_1 = 'NC State'
else:
    winner_1 = 'North Dakota State'
print(winner_1+' wins the first play-in game')
#play_in2_11 = 10
temp_odds = random.random()*10.
if(temp_odds>5):
    winner_2 = 'Belmont'
else:
    winner_2 = 'Temple'
print(winner_2+' wins the second play-in game')
#play_in3_16 = 18
temp_odds = random.random()*10.
if(temp_odds>5):
    winner_3 = 'Fairleigh Dickinson'
else:
    winner_3 = 'Prarie View A&M'
print(winner_3+' wins the third play-in game')
#play_in4_11 = 26
temp_odds = random.random()*10.
if(temp_odds>5):
    winner_4 = 'Arizona State'
else:
    winner_4 = 'Saint John\'s'
print(winner_4+' wins the fourth play-in game\n\n')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Define the March Madness bracket teams
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
seeds = dict()
seeds['Duke'] = 1
seeds[winner_1] = 16
seeds['VCU'] = 8
seeds['UCF'] = 9
seeds['Miss. State'] = 5
seeds['Liberty'] = 12
seeds['Virginia Tech'] = 4
seeds['Saint Louis'] = 13
seeds['Maryland'] = 6
seeds[winner_2] = 11
seeds['LSU'] = 3
seeds['Yale'] = 14
seeds['Louisville'] = 7
seeds['Minnesota'] = 10
seeds['Michigan State'] = 2
seeds['Bradley'] = 15

seeds['Gonzaga'] = 1
seeds[winner_3] = 16 # NC Central / Texas Southern
seeds['Syracuse'] = 8
seeds['Baylor'] = 9
seeds['Marquette'] = 5
seeds['Murray State'] = 12
seeds['Florida State'] = 4
seeds['Vermont'] = 13
seeds['Buffalo'] = 6
seeds[winner_4] = 11
seeds['Texas Tech'] = 3
seeds['N. Kentucky'] = 14
seeds['Nevada'] = 7
seeds['Florida'] = 10
seeds['Michigan'] = 2
seeds['Montana'] = 15

seeds['Virginia'] = 1
seeds['Gardner-Webb'] = 16
seeds['Ole Miss'] = 8
seeds['Oklahoma'] = 9
seeds['Wisconsin'] = 5
seeds['Oregon'] = 12
seeds['Kansas State'] = 4
seeds['UC Irvine'] = 13
seeds['Villanova'] = 6
seeds['Saint Mary\'s'] = 11
seeds['Purdue'] = 3
seeds['Old Dominion'] = 14
seeds['Cincinnati'] = 7
seeds['Iowa'] = 10
seeds['Tennessee'] = 2
seeds['Colgate'] = 15

seeds['North Carolina'] = 1
seeds['Iona'] = 16
seeds['Utah State'] = 8
seeds['Washington'] = 9
seeds['Auburn'] = 5
seeds['New Mexico State'] = 12
seeds['Kansas'] = 4
seeds['Northeastern'] = 13
seeds['Iowa State'] = 6
seeds['Ohio State'] = 11
seeds['Houston'] = 3
seeds['Georgia State'] = 14
seeds['Wofford'] = 7
seeds['Seton Hall'] = 10
seeds['Kentucky'] = 2
seeds['Abilene Christian'] = 15

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
