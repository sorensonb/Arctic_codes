#!/usr/bin/env python

"""


"""

from espn_parser import *

# NOTE: 2/25/2023 - saved data


cdict = read_conference_file()
gdict = load_json_stats()


# Ranked opponent factor:
# if winning = np.log10(27 - opp_rank) * 0.5
# if losing  = np.log10(27 - opp_rank) * 0.25

#sys.exit()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Work through steps required for each week of ranking
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

##fcs_rankings, rank_values = \
##    generate_division_rankings(gdict, division = 'FCS', \
##    use_rankings = False)

# Set pre-season rankings
# ONLY DO FOR FIRST WEEK OF SEASON
gdict = set_preseason_rankings(gdict, rank_file = 'sagarin_rank_postseason_2022.txt')
gdict = simulate_whole_season(gdict, 12)

fcs_rankings, rank_values = \
    extract_rank_and_values(gdict, division = 'FCS')

sorted_names = [rval[1] for rval in rank_values]

# NOTE: Modify this to allow reading in either massey or sagarin 
#       verification rankings. Also, make sure that all sagarin
#       ranking names match those in the massey files.
massey_rankings = read_poll_ranks()

total_abs_error = 0

for key in sorted_names:
    if((key in fcs_rankings.keys()) & (key in massey_rankings.keys())):
        conf =  gdict[key]['Conference']
        calced_ranking = fcs_rankings[key]['rank']
        massey_rank    = massey_rankings[key]
        pcnt_error     = (massey_rank - calced_ranking)
        total_abs_error += abs(pcnt_error)
    #    print(key, calced_ranking, massey_rank, pcnt_error)
        print(f'{calced_ranking:4d}{massey_rank:4d}{pcnt_error:4d} {key:25}{conf:30}')

print(total_abs_error / len(rank_values))

sys.exit()

calc_total_season_value(gdict, 'North Dakota', debug = True)

# Calculate the game quality for each team
#   - If the first week of the semester, use only
#     the preseason conference rankings and the 
#     single-game stats to determine the game
#     quality
#   - If after the first week, include the other team
#     ranking.
calc_game_quality(gdict, 'North Dakota', 9, debug = True)

# Recalculate team rankings

# Recalculate conference rankings based on new team rankings

# Possible game ranking equation:
# game_factor = 1 + np.log10(team_score/opp_score) + \
#                   np.log10(opp_conf_rank)

### To re-load everything from ESPN,
##gdict = load_ESPN_data(1, 12)
