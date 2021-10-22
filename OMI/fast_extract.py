import calendar
import datetime as dt

def get_periods(month, sta_y, sta_e, days):
	dates = list()
	for year in range(sta_y, sta_e + 1):
		sta_d = dt.datetime(year, month, 1) - dt.timedelta(days=days)
		end_d = dt.datetime(year, month, calendar.monthrange(year,month)[1]) + dt.timedelta(days=days)
		
		sta = sta_d.strftime("%Y%m%d")
		end = end_d.strftime("%Y%m%d")

		dates.append(sta)
		dates.append(end)
		
	return " ".join(dates)

def run():
	for month in range(4, 10):
		for year in range(2006,2020):
			train_on = get_periods(month, 2006, 2020, 5)
			#print(train_on)
			recon_on = get_periods(month, year, year, 0)
			#print(recon_on)
			
			dates = " ".join([train_on, recon_on])
			
			print("python3 OMIprocess.py", dates)

if __name__ == "__main__":
	run()
