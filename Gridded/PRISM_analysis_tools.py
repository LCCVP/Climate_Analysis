#PRISM extract function
#Author: Tony Chang
#Date: 8/16/2013
#Abstract: Extracts bounding box of PRISM dataset for a specified area and time period

import numpy as np
import scipy
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import cm
import gdal 
from gdalconst import *
import osr 
import shapefile
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
#===================================================================================================
#===================================CLASS DEFINITIONS===============================================
#===================================================================================================
class PRISMData(object):
    #initialize function to construct PRISMdata class to contain the data
    def __init__(self, year=None, month=None, ncols=None, nrows=None, xll=None, yul=None, csize=None, NODATA=None, data=None):
        self.year = year
        self.month = month
        if month == 14:         #month 14 in PRISM data represents the mean of the years
            self.season = "ALL"
        elif (month < 3 or month == 12):
            self.season = "Win"
        elif month < 6:
            self.season = "Spr"
        elif month < 9:
            self.season = "Sum"
        else:
            self.season = "Fal"
        self.ncols = ncols
        self.nrows = nrows
        self.xll = xll
        self.yul = yul
        self.csize = csize
        self.NODATA = NODATA
        self.data = data

class Annualclimatedata(object):
	#annual climate data summary of PRISMData object to summarize in years only
	def __init__(self, year=None, data=None):
		self.year = year
		self.data = data

#===================================================================================================
#===================================PRISM EXTRACT FUNCTIONS=========================================
#===================================================================================================
def PRISMdataextract(BeginYear, EndYear, var, AOA):
# extracts PRISM data into a PRISMData class given the start year of interest (BeginYear),
# the last year of interest (EndYear), variable of interest in string format (i.e. 'tmin', 'tmax', 'ppt')
# and the Area of Analysis (AOA) in array formatted ([xmin, xmax, ymin, ymax]) in WGS74 coordinate system
# **This may take some time due to the number of files and size. It is recommended that users use this code and exporting the cut grids
# to tiff format for future analysis
 
	workspace = "D:\\CHANG\\Climate_Models\\US_PRISM_800m\\uncompressed\\"+var+"\\" #<---DEFINE LOCATION OF PRISM FILES HERE!
	PRISMExtent = [-125.02083333333, 24.0625, -66.47916757, 49.9375] # The fixed area extent of US PRISM
	minx = AOA[0] 
	maxx = AOA[1]
	miny = AOA[2]
	maxy = AOA[3]
	filenum =1
	#open the first data file to get the attributes of PRISM data
	Pgrid = workspace + "us_" + var + "_" + str(BeginYear) + ".0" + str(filenum) #uncompressed PRISM file name
	readfile = open(Pgrid, 'r')
	a = readfile.readline()
	temp = a.split()
	ncols = int(temp[1])        #Define number of columns
	a = readfile.readline()
	temp = a.split()
	nrows = int(temp[1])        #Define number of rows
	a = readfile.readline()
	temp = a.split()
	xllcorner = float(temp[1])  #Define xll corner
	a = readfile.readline()
	temp = a.split()
	yllcorner = float(temp[1])  #Define yll corner
	a = readfile.readline()
	temp = a.split()
	cellsize  = float(temp[1])  #Define cellsize
	a = readfile.readline()
	temp = a.split()
	NODATA_value  = temp[1]     #Define NoData value
	readfile.close()
	yulcorner = PRISMExtent[1]+(cellsize*nrows)
	xstart = int((minx - PRISMExtent[0])/cellsize)    #first x-extent index
	xend = xstart + int((maxx-minx)/cellsize)       #end x-extent index
	ystart = int((yulcorner - maxy)/cellsize)         #first y-extent index
	yend = ystart + int((maxy-miny)/cellsize)       # end of y-extent index
	#Start gathering data
	Pdata = [] #define empty array to store PRISM classes
	for searchyear in range(BeginYear, EndYear+1): #looping through years of interest
		for filenum in range(1, 13):
			addmatrix = []              #List to store PRISM ascii data
			if filenum == 13:
				continue                #month 13 does not exist, skip to the next iteration
			elif filenum < 10:
				Psource = workspace + "us_" + var + "_" + str(searchyear) + ".0" + str(filenum) #if your file name is not of the same format,
			else:																				#change here
				Psource = workspace + "us_" + var + "_" + str(searchyear) + "." + str(filenum)
			readfile =  open (Psource,'r')
			nhead = 6                   #First 6 lines of the header to be removed
			for z in range(nhead):      #Strip out header
				a = readfile.readline()
			for y_pos in range(0, nrows+1):
				line = readfile.readline()
				datarow = line.split()
				if (y_pos >= ystart and y_pos <= yend):
				   newrow = datarow[xstart:(xend+1)]
				   addmatrix.append(newrow)
			newcols = len(addmatrix[0]) #define new column length
			newrows = len(addmatrix)    #define new row length
			newyulcorner = yulcorner - (ystart*cellsize)
			newxllcorner = PRISMExtent[0] + (xstart*cellsize)
			addmatrix = np.array(addmatrix).astype(float)/100 #changes addmatrix list into array for statistical analysis and divides by 100 to put units in correct format
			x = PRISMData(searchyear,filenum, newcols, newrows, newxllcorner, newyulcorner, cellsize, NODATA_value, addmatrix) #Create instance of PRISMData object
			Pdata.append(x) 
	return(Pdata)

def annualgrid(data):  #generates the PRISM grids at an annual time step
	numyears = int(len(data)/12)
	anu_year = np.zeros(np.shape(data[0].data))
	anu_series = []
	by = data[0].year
	ey = data[-1].year
	currentyear =by
	i=0
	counter =0
	while (i<len(data)):
		if (currentyear == data[i].year):
			anu_year += data[i].data
			counter += 1
		elif (currentyear != data[i].year):
			x =Annualclimatedata(data[i-1].year,anu_year/counter)
			anu_series.append(x)
			counter = 0
			currentyear = data[i].year
			anu_year = np.zeros(np.shape(data[0].data))
			anu_year += data[i].data
			counter += 1
		i+=1
	#last iteration
	x =Annualclimatedata(data[i-1].year,anu_year/counter)
	anu_series.append(x)
	return(anu_series)  

#======================================================================================
#=============================ANALYSIS FUNCTIONS=======================================
#======================================================================================
def lstfit(data):
	#least squares fit of the data at the individual cell level
	nrows, ncols = (np.shape(data[0].data))
	n = len(data)
	mu = np.zeros((nrows,ncols))
	for i in range(n):
		mu += data[i].data
	mu = mu/n
	tmu = np.sum(np.arange(n)+1)/n
	tmu = np.ones((nrows,ncols))*tmu
	Sxx = np.zeros((nrows,ncols))
	Sxy = np.zeros((nrows,ncols))
	SST = np.zeros((nrows,ncols))
	SSW = np.zeros((nrows,ncols))
	for i in range(n):
		dt = data[i].data
		Sxx += (((np.ones((nrows,ncols))*(i+1))-tmu)**2)
		Sxy += (dt-mu)*((np.ones((nrows,ncols))*(i+1))-tmu)
	beta1 = Sxy/Sxx
	beta0 = mu - (beta1*tmu)
	return(beta1,beta0)	

def temporalgradient(data):
	#takes the Pdata and generates the rate of change at the individual cell level
	adata = annualgrid(data)
	dzdt, b0 = lstfit(adata)
	return(dzdt)

def climatemean(data):
	#solves the period mean for plotting purposes
	adata = annualgrid(data)
	nyears = len(adata)
	periodmean = np.zeros(np.shape(data[0].data))
	for i in range(nyears):
		periodmean += adata[i].data
	periodmean = periodmean/nyears
	return(periodmean)
	
def domainmean(data): #returns a time series of the Pdata grid means and detrended series
	Pmean = []
	for i in range(len(data)):
		Pmean.append(np.mean(data[i].data))
	t = np.arange(len(data))+1
	y = np.array(Pmean)
	b1,b0 = np.polyfit(t, y, 1)
	det_y = y - (b1*t+b0)
	return(y, det_y)

def simplemovingavg(ts,lag):
	#determine the lag-month simple moving avg (use the domainmean first)
	run_avg = []
	for i in range(lag, len(ts)):
		run_avg.append(np.mean(ts[i-lag:i]))
	timearray = np.arange(lag,len(ts))
	return(np.array(run_avg))

#======================================================================================
#===============================TIME SERIES PLOT FUNCTIONS=============================
#======================================================================================

def simpleannualplot(data, tline='y'):
	#plots a simple annual time series anomaly of the mean grid area 
	#trend line is plotted by default, can be changed to 'n' to not plot
	adata = annualgrid(data)
	t = np.arange(adata[0].year,adata[-1].year+1)
	y, dty = domainmean(adata)
	plt.plot(t,y-np.mean(y))
	if tline=='y':
		b1,b0=np.polyfit(t,y-np.mean(y),1)
		plt.plot(t, t*b1+b0, ls='--')
	return()

def timeseriessummary_plot(data, cvar, mws = 3):	#plots the timeseries from PRISM given a specified beginning year and end year
	# user should specify the climate variable type ('tmin'..'ppt') so labels are appropriate
	#mws is the moving window size for the running average 
	adata = annualgrid(data)
	by = adata[0].year
	ey = adata[-1].year
	data_mu = climatemean(data)
	t = np.arange(by,ey+1)
	domain_array = []
	for j in range(len(adata)):
		domain_array.append(np.mean(adata[j].data - data_mu)) #calculate the anomaly
	domain_array = np.array(domain_array)
	base = np.zeros(len(adata))
	datapos = np.zeros(len(adata))
	dataneg = np.zeros(len(adata))
	posi = np.where(domain_array>0)
	negi = np.where(domain_array<0)
	datapos[posi] = domain_array[posi]
	dataneg[negi] = domain_array[negi]
	movingwindowsize = mws #moving window size specified by user
	sma = simplemovingavg(domain_array,movingwindowsize)
	cgrad = lstfit(adata)
	beta = scipy.stats.linregress(t, domain_array)
	simcgrad = beta[1]+(t*beta[0])
	coef = 9./5 #use this to convert from C to F or otherwise
	#plotting routine
	if (cvar == 'ppt'):
		coef = 1
		ax = plt.subplot(1,1,1,axisbg ='0.9', xlabel = 'Year', ylabel = 'Precipitation anomaly $(mm)$')
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(20)
		plt.vlines(t, datapos*coef, base, color = 'red', linewidth =10, alpha = 0.8)
		plt.vlines(t, base, dataneg*coef, color = 'blue', linewidth =10, alpha = 0.8)
		plt.plot(t[movingwindowsize:],sma*coef, color='green', label = str(movingwindowsize)+'-year moving window')
		plt.plot(t,simcgrad*coef, color='orange', ls = '--', label = 'Trend {0:.2f} $(mm/decade)$'.format((beta[0]*coef*10)))
		plt.grid()
		plt.legend(loc = 'lower right')
	else:
		ax = plt.subplot(1,1,1,axisbg ='0.9', xlabel = 'Year', ylabel = 'Temperature anomaly $(^oF)$')
		for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
			item.set_fontsize(20)
		plt.vlines(t, datapos*coef, base, color = 'red', linewidth =10, alpha = 0.8)
		plt.vlines(t, base, dataneg*coef, color = 'blue', linewidth =10, alpha = 0.8)
		plt.plot(t[movingwindowsize:],sma*coef, color='green', label = str(movingwindowsize)+'-year moving window')
		plt.plot(t,simcgrad*coef, color='orange', ls = '--', label = 'Trend {0:.2f} $(^oF/decade)$'.format((beta[0]*coef*10)))
		plt.grid()
		plt.legend(loc = 'lower right')
	plt.show()
	return()
#======================================================================================
#===============================MAP BUILDING FUNCTIONS=================================
#======================================================================================

def Topoextract():
	#extracts the DEM file for plotting under the map
    elevPath = "D:\\CHANG\\GIS_Data\\DEM\\TIFF\\dem_gye800m1.tif"   
    ds = gdal.Open(elevPath)
    elev = np.array(ds.GetRasterBand(1).ReadAsArray())
    ds = None #close files
    return(elev)

def hillshade(data,scale=10.0,azdeg=165.0,altdeg=45.0):
	# takes in the elevation grid (data) and generates a hillshade matrix for plotting
	# convert alt, az to radians
	az = azdeg*np.pi/180.0
	alt = altdeg*np.pi/180.0
	# gradient in x and y directions
	dx, dy = np.gradient(data/float(scale))
	slope = 0.5*np.pi - np.arctan(np.hypot(dx, dy))
	aspect = np.arctan2(dx, dy)
	intensity = np.sin(alt)*np.sin(slope) + np.cos(alt)*np.cos(slope)*np.cos(-az - aspect - 0.5*np.pi)
	intensity = (intensity - intensity.min())/(intensity.max() - intensity.min())
	return(intensity)	

def drawArea(AOA):
	#plots the boundary for natural resource of interest
	minx = AOA[0] 
	maxx = AOA[1]
	miny = AOA[2]
	maxy = AOA[3]
	sf = shapefile.Reader("d:\\chang\\gis_data\\gye_shapes\\gye.shp") #change the shapefile location here!
	recs    = sf.records()
	shapes  = sf.shapes()
	Nshp    = len(shapes)
	cns     = []
	for nshp in range(Nshp):
		cns.append(recs[nshp][1])
	cns = np.array(cns)
	cma    = cm.get_cmap('Dark2')
	cccol = cma(1.*np.arange(Nshp)/Nshp)
	ax = plt.subplot(111)
	for nshp in range(Nshp):
		ptchs   = []
		pts     = np.array(shapes[nshp].points)
		prt     = shapes[nshp].parts
		par     = list(prt) + [pts.shape[0]]
		for pij in range(len(prt)):
			ptchs.append(Polygon(pts[par[pij]:par[pij+1]]))
			ax.add_collection(PatchCollection(ptchs,facecolor ='None',edgecolor='k', linewidths=2))#facecolor=cccol[nshp,:]
	ax.set_xlim(minx,maxx)
	ax.set_ylim(miny,maxy)
	return()

def plotelevation(AOA):
	#plots the hillshade given the topography and area of interest as a background for plots
	ele = Topoextract()
	hill = hillshade(ele)
	im = plt.imshow(hill, cmap = cm.Greys_r, extent =AOA)
	#im2 = plt.imshow(ele, cmap = cm.Spectral, alpha= 0.7, extent =ae)
	plt.xlabel('Longitude (DD)')
	plt.ylabel('Latitude (DD)')
	#cb = plt.colorbar(im2)
	#cb.set_label('Elevation (m)')
	plt.grid(alpha =0.4)	
	return()
	
def plotPRISMgrad(pdata,AOA):
	#plots the PRISM gradients with hillshade and resource boundary
	plotelevation(AOA)
	dzdt = temporalgradient(pdata)
	ax = plt.imshow(dzdt, extent = AOA, alpha =0.6)
	cb = plt.colorbar(ax)
	cb.set_label('Temporal gradient') #change the label here depending on the variable type
	drawArea(AOA)
	return()