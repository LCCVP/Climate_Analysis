#=======================================================
#Author: Tony Chang
#Institution: Montana State University
#Department: Ecology
#Date: 07.08.2013
#=======================================================
'''
#	Abstract:
#	Toolset of functions to format NCDC Global Historical Climate Network dataset for analysis
#
#	Dependency:
#	GCHND-stations.txt, 
#	Uncompressed GHCND-all.tar.gz daily .dly data files for individual stations

'''

import numpy as np
import csv
from matplotlib import pyplot as plt

def stationreader():
	workspace = "D:\\chang\\climate_models\\station_data\\ghcn_daily\\" #change this according to your local machine file address
	filename = "ghcnd-stations.txt"
	f = workspace + filename
	s =[]
	with open(f, newline = '') as f1:
		reader= csv.reader(f1)
		for row in reader:
			s.append(row)
	s = np.array(s)
	code = []
	lat = []
	lon = []
	ele = []
	st = []
	name = []
	gsnflag = []
	hcnflag =[]
	wmoid = []
	for i in range(len(s)):
		code.append(s[i][0][:11])
		lat.append(float(s[i][0][12:20]))
		lon.append(float(s[i][0][21:30]))
		ele.append(float(s[i][0][31:37]))
		st.append(s[i][0][38:40])
		name.append(s[i][0][41:71])
		gsnflag.append(s[i][0][72:75])
		hcnflag.append(s[i][0][76:79])
		wmoid.append(s[i][0][81:85])
	station = {'code' : np.array(code), 'lat' : np.array(lat), 'lon' : np.array(lon), 'ele':np.array(ele), 'st':np.array(st), 'name':np.array(name), 'gsn':np.array(gsnflag), 'hcn':np.array(hcnflag), 'wmo':np.array(wmoid)}
	l = np.array([code,lat,lon, ele, st, name, gsnflag, hcnflag, wmoid])
	return(station)

def formatdata(stationlist, var, time):
	ds = datagather(stationlist)
	vardata = databyvar(ds, var)
	datas, stations = regionalsummary(vardata,time)
	return(datas)
	
#GYE extent
#xmin = -112.438; xmax = -108.271; ymin = 42.262; ymax = 46.187
def stationfilter(station, xmin,xmax,ymin,ymax, minelev= 0, maxelev= 9999): #filters the station list by the bounding box extent and desired elevation range
	fstation = []
	dt = [('code', 'S11'),('lat', 'f8'),('lon','f8'),('ele','f8'),('st','S2'),('name', 'S30'),('gsn','S3'),('hcn','S3'),('wmo','f5')]
	for i in range(len(station['code'])):
			sy = station['lat'][i]
			sx = station['lon'][i]
			elev = station['ele'][i]
			if (((sx>xmin) and (sx<xmax)) and ((sy>ymin) and (sy<ymax)) and ((elev >=minelev) and (elev <=maxelev))):
				fstation.append([station['code'][i],station['lat'][i],station['lon'][i],station['ele'][i],station['st'][i],station['name'][i],station['gsn'][i],station['hcn'][i],station['wmo'][i]])
	fstation = np.array(fstation)
	return(fstation)
	
def datareader(stationid): #input station desired 
	workspace = "d:\\chang\\climate_models\\station_data\\ghcn_daily\\ghcnd_all\\ghcnd_all\\" #change this to your local machine file address
	data = []
	filename = stationid + '.dly'
	f = workspace + filename
	s =[]
	with open(f, newline = '') as f1:
		reader= csv.reader(f1)
		for row in reader:
			s.append(row)
	s = np.array(s)
	data.append(s)
	return(np.array(data))

def datagather(fstation):
	dataset = []
	for i in range(len(fstation)):
		dataset.append(datareader(fstation[i][0])[0])
	return(np.array(dataset))

def stationlookup(stationid):
	'''Looks up station attributes given ID code'''
	stationlist = stationreader()
	stationattr = []
	for i in range(len(stationid)):
		ix = np.where(stationlist['code'] == stationid[i])
		stationattr.append([stationlist['code'][ix][0], stationlist['lat'][ix][0],stationlist['lon'][ix][0], stationlist['ele'][ix][0], stationlist['st'][ix][0], stationlist['name'][ix][0], stationlist['gsn'][ix][0], stationlist['hcn'][ix][0], stationlist['wmo'][ix][0]])
	return(np.array(stationattr))
	
def databyvar(dataset,var):
	#filters the dataset to contain only one var climate element
	"""
	var = 
	PRCP = Precipitation (tenths of mm) 
	SNOW = Snowfall (mm)
	SNWD = Snow depth (mm)
	TMAX = Maximum temperature (tenths of degrees C)
	TMIN = Minimum temperature (tenths of degrees C)
	"""
	fdataset = []
	for i in range(len(dataset)):
		stationdata = []
		count = 0
		for j in range(len(dataset[i])):
			if (dataset[i][j][0][17:21] == var):
				stationdata.append(dataset[i][j])
				count +=1
		if (count != 0):
			fdataset.append(np.array(stationdata))
	return(np.array(fdataset))

def datasummary(dataset, timescale):
	"""summarizes the data by month('m'), season('s'), or annually('a') at the individual station level
	(note dataset must be pre-filtered by the databyvar function)"""
	if (timescale == 'm'): #monthly case
		return(monsummary(dataset))
	elif (timescale == 's'): #seasonal case
		return(seasonsummary(dataset))
	elif (timescale == 'a'): #annual case
		return(annualsummary(dataset))
	
def regionalsummary(dataset,timescale):
	'''summarizes the data for the entire region by month('m'), season('s'), or annually('a') 
	(note dataset must be pre-filtered by the databyvar function)'''
	by = 1900
	ey = 2013
	if (timescale == 'm'): #monthly case
		mondata = monsummary(dataset)
		r = np.repeat(np.arange(by,ey), 12) #take only the dataset that runs fron 1895 to 2012
		m = np.tile(np.arange(1,13),(ey-by))
		z = np.zeros((2,(ey-by)*12))
		rsummary = np.vstack((r,m,z)).T
		stationindex = np.zeros(((ey-by)*12, len(mondata)))
		for i in range(len(mondata)):
			for j in range(len(mondata[i])):
				currentyear = float(mondata[i][j][0][11:15])
				currentmonth = float(mondata[i][j][0][15:17])
				currentval = float(mondata[i][j][1])
				if (~np.isnan(currentval) and (currentyear>=by and currentyear<ey)): #if the value is not np.nan and within the year boundaries
					index = np.argwhere((rsummary[:,0] == currentyear) & (rsummary[:,1] == currentmonth)).flatten()[0]
					rsummary[index][2] += currentval #add the value
					rsummary[index][3] += 1 #count the station
					stationindex[index][i] += 1 #note the station index to develop a station list later
		stations = []
		for k in range(len(stationindex)):
			stationcodes = []
			for n in range(len(stationindex[k])):
				if (stationindex[k][n] == 1): 
					stationcodes.append(mondata[n][0][0][:11])
			stations.append(stationcodes)
		rout = np.array([rsummary[:,0], rsummary[:,1], rsummary[:,2]/rsummary[:,3], rsummary[:,3]]).T
		return(rout, np.array(stations))
	#==================	
	if (timescale == 's'): #seasonal case
		seadata = seasonsummary(dataset)
		r = np.arange(by,ey) #take only the dataset that runs fron 1895 to 2012
		z = np.zeros((8,(ey-by)))
		rsummary = np.vstack((r,z)).T
		stationindex = np.zeros(((ey-by), len(seadata)))
		for i in range(len(seadata)):
			for j in range(len(seadata[i])):
				currentyear = float(seadata[i][j][1])
				winval = float(seadata[i][j][2])
				sprval = float(seadata[i][j][3])
				sumval = float(seadata[i][j][4])
				falval = float(seadata[i][j][5])
				if (currentyear>=by and currentyear<ey): #if the current year is within the year boundaries
					index = np.where(rsummary[:,0] == currentyear)[0][0]
					if (~np.isnan(winval)):
						rsummary[index][1] += winval #add the value
						rsummary[index][5] += 1 #count the station
					if (~np.isnan(sprval)):
						rsummary[index][2] += sprval #add the value
						rsummary[index][6] += 1 #count the station
					if (~np.isnan(sumval)):
						rsummary[index][3] += sumval #add the value
						rsummary[index][7] += 1 #count the station
					if (~np.isnan(falval)):
						rsummary[index][4] += falval #add the value
						rsummary[index][8] += 1 #count the station
					stationindex[index][i] += 1 #note the station index to develop a station list later
		stations = []
		for k in range(len(stationindex)):
			stationcodes = []
			for n in range(len(stationindex[k])):
				if (stationindex[k][n] == 1): 
					stationcodes.append(seadata[n][0][0][:11])
			stations.append(stationcodes)
		windata = rsummary[:,1]/rsummary[:,5]
		sprdata = rsummary[:,2]/rsummary[:,6]
		sumdata = rsummary[:,3]/rsummary[:,7]
		faldata = rsummary[:,4]/rsummary[:,8]
		rout = np.array([rsummary[:,0],windata,sprdata,sumdata,faldata, rsummary[:,5], rsummary[:,6], rsummary[:,7], rsummary[:,8]]).T
		return(rout, np.array(stations))
	#==============	
	if (timescale == 'a'): #annual case
		adata = annualsummary(dataset)
		r = np.arange(by,ey) #take only the dataset that runs fron by to ey
		z = np.zeros((2,(ey-by)))
		rsummary = np.vstack((r,z)).T
		stationindex = np.zeros(((ey-by), len(adata)))
		for i in range(len(adata)):
			for j in range(len(adata[i])):
				currentyear = float(adata[i][j][1])
				currentval = float(adata[i][j][2])
				if (~np.isnan(currentval) and (currentyear>=by and currentyear<ey)): #if the value is not np.nan and within the year boundaries
					index = np.where(rsummary[:,0] == currentyear)[0][0]
					rsummary[index][1] += currentval #add the value
					rsummary[index][2] += 1 #count the station
					stationindex[index][i] += 1 #note the station index to develop a station list later
		stations = []
		for k in range(len(stationindex)):
			stationcodes = []
			for n in range(len(stationindex[k])):
				if (stationindex[k][n] >= 1): 
					stationcodes.append(adata[n][0][0][:11])
			stations.append(stationcodes)
		rout = np.array([rsummary[:,0], (rsummary[:,1]/rsummary[:,2]), rsummary[:,2]]).T
		return(rout, np.array(stations))
	

def monsummary(dataset): 
	"""returns the monthly summary of the filtered dataset"""
	datasum = []
	months = np.arange(1,13)
	for i in range(len(dataset)):
		stationdata = []
		for j in range(len(dataset[i])):
			k = 21
			label = dataset[i][j][0][:k]
			monmean = 0
			moncount = 0
			flagcount = 0
			nullcount = 0
			while (k < len(dataset[i][j][0])-1):
				dailyval = float(dataset[i][j][0][k:k+5])
				qflag = dataset[i][j][0][k+5+1]
				if (dailyval == -9999):
					nullcount += 1 #count the number of null values
				if (dailyval != -9999 and qflag == (' ' or 'N')): #check for null values or failed quality flags
					monmean += dailyval
					moncount += 1
				else:
					flagcount += 1 #count the missing value or flag for record
				k += 8 #go to next day iteration
			if (flagcount < 15): #if the flag and null sum count are less than 15
				monmean = monmean/moncount # calculate the average month value
				stationdata.append([label,monmean, flagcount, nullcount, flagcount-nullcount])
			else:
				monmean = np.nan #if the number of flagged values are greater than or equal to 15, then discard then set month average to np.nan
				stationdata.append([label,monmean, flagcount, nullcount, flagcount-nullcount])
		datasum.append(stationdata)
	return(np.array(datasum))

def seasonsummary(dataset):
	monsum = monsummary(dataset)
	seasonsum = []
	for i in range(len(monsum)):
		firstyear = int(monsum[i][0][0][11:15])
		lastyear = int(monsum[i][-1][0][11:15]) #find the year range of the dataset
		y = np.arange(firstyear, lastyear+1)
		s = np.zeros((8,(lastyear-firstyear+1))) #season array has addition 4 elements to count number of valid month inputs
		annualseason = np.vstack((y,s)).T #create a season array for each year 
		for j in range(len(monsum[i])):
			currentyear = int(monsum[i][j][0][11:15])
			currentmonth = int(monsum[i][j][0][15:17])
			currentval = monsum[i][j][1]
			if (currentmonth == 12 and currentyear != firstyear): #December of the previous year considered winter of current year
				annualseason[np.where(annualseason==currentyear-1)[0][0]][1] += currentval
				annualseason[np.where(annualseason==currentyear-1)[0][0]][5] += 1
			elif (currentmonth == (1 or 2)): #winter case
				annualseason[np.where(annualseason==currentyear)[0][0]][1] += currentval
				annualseason[np.where(annualseason==currentyear)[0][0]][5] += 1
			elif (currentmonth == (3 or 4 or 5)): #spring case
				annualseason[np.where(annualseason==currentyear)[0][0]][2] += currentval
				annualseason[np.where(annualseason==currentyear)[0][0]][6] += 1
			elif (currentmonth == (6 or 7 or 8)): #summer case
				annualseason[np.where(annualseason==currentyear)[0][0]][3] += currentval
				annualseason[np.where(annualseason==currentyear)[0][0]][7] += 1
			elif (currentmonth == (9 or 10 or 11)): #fall case
				annualseason[np.where(annualseason==currentyear)[0][0]][4] += currentval
				annualseason[np.where(annualseason==currentyear)[0][0]][8] += 1	
		seasonadd = []
		label = monsum[i][0][0][:11]+monsum[i][0][0][17:21]
		for k in range(len(annualseason)):
			win = annualseason[k][1]/annualseason[k][5]
			spr = annualseason[k][2]/annualseason[k][6]
			sum = annualseason[k][3]/annualseason[k][7]
			fal = annualseason[k][4]/annualseason[k][8]
			seasonadd.append([label,annualseason[k][0], win, spr, sum, fal])
		seasonsum.append(np.array(seasonadd))
	return(np.array(seasonsum)) #returns a seasonal summer array in the form ('station_vartype', 'year','winterval', 'springval', 'summerval', 'fallval')
								#nan will be report for values where a month is missing for that season

def annualsummary(dataset):
	monsum = monsummary(dataset)
	annualsum = []
	for i in range(len(monsum)):
		firstyear = int(monsum[i][0][0][11:15])
		lastyear = int(monsum[i][-1][0][11:15])
		y = np.arange(firstyear, lastyear +1)
		s = np.zeros((3, (lastyear-firstyear +1)))
		yearval = np.vstack((y,s)).T
		for j in range(len(monsum[i])):
			currentyear= int(monsum[i][j][0][11:15])
			currentval = monsum[i][j][1]
			if (np.isnan(currentval)):
				yearval[np.where(yearval==currentyear)[0][0]][3] += 1 #count how many null months there are
			else:
				yearval[np.where(yearval==currentyear)[0][0]][1] += currentval
				yearval[np.where(yearval==currentyear)[0][0]][2] += 1
		yearadd =[]
		label = monsum[i][0][0][:11]+monsum[i][0][0][17:21]
		for k in range(len(yearval)):
			if (yearval[k][0]!= 2013 or yearval[k][3]>3): #neglect 2013 as the year is not complete yet, or if more than 3 months are missing
				yearadd.append(([label,yearval[k][0], yearval[k][1]/yearval[k][2], yearval[k][3]]))
		annualsum.append(np.array(yearadd))
	return(np.array(annualsum)) #returns an annualsum array in the form ('station_vartype', 'year', 'val', 'number of null values') 
								#nan will be report for values where a year is missing for that season							
	
#-------MAIN-------#
xmin = -112.438; xmax = -108.271; ymin = 42.262; ymax = 46.187 #GYE bounding box
#xmin = -114.62; xmax = -108.271; ymin = 41.87; ymax = 47.4 #ecoregion III middle rockies bounding box
min_e1 = 0; max_e1 = 2250; min_e2 = 2250; max_e2=9999
stationlist = stationreader() #this loads the full station list into memory
'''
#example code below for extracting station data
fstation1 = stationfilter(stationlist, xmin, xmax, ymin, ymax, min_e1, max_e1) #from here filter the full station list by the desire extent and elevation
fstation2 = stationfilter(stationlist, xmin, xmax, ymin, ymax, min_e2, max_e2)
dataset1 = datagather(fstation1)
dataset2 = datagather(fstation2)
var = 'TMIN'
tmindata1 = databyvar(dataset1,var)
tmindata2 = databyvar(dataset2,var)
atmin1,stat1 = regionalsummary(tmindata1,'a')
atmin2, stat2 = regionalsummary(tmindata2,'a')
'''