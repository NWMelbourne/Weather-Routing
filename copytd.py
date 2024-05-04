#####################
#this is the dev program




########## not-good assumptions this currently runs on, TODO list:
#deal with curvature somehow


from global_land_mask import globe
import geopy.distance
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid
import numpy as np
import math
import time
from time import gmtime, strftime
import requests
import urllib.request 
import os.path



#Now functions are defined and data are loaded
arrived=False
radialsend=30
first_iteration_scaling=5
boundarypoints=30
dt=2 #hours
plottingidx=1
#resolution of coastline plot AND polygon check for is_land(), use c=crude for ocean crossing
resolu='i' #c, l, i, h, f
#FOR quiver wind plot
descaling=3


#start/stop and all positions are of the form lat,lon aka y,x !!!backwards from usual!!!
#PLACES#######################

#seattle = [48.5,-124.8] [48.25,-123.0]
#CI = [34,-120.5]
#tipBaja = [22.8,-110.38]
#HNL = [22,-158] [21.4,-157.5]
#japan pacific= [30,162]
#don't make the latitudes or longitudes identical; makes divide by zero errors...
start = [25.5,-112.5]
stop = [18.5,-104.5]

#making the plotting box based on start/stop
#dont mess with it, it works. i guess you can mess with it, but dont
if start[1]< stop[1]:
    lonbounds=[start[1]-abs(start[1]-stop[1])/3,stop[1]+abs(start[1]-stop[1])/3]
else:
    lonbounds=[stop[1]-abs(start[1]-stop[1])/3,start[1]+abs(start[1]-stop[1])/3]

if start[0]< stop[0]:
    latbounds=[start[0]-abs(start[0]-stop[0])/3,stop[0]+abs(start[0]-stop[0])/3]
else:
    latbounds=[stop[0]-abs(start[0]-stop[0])/3,start[0]+abs(start[0]-stop[0])/3]
latbounds.sort()
lonbounds.sort()
if abs(latbounds[1]-latbounds[0])/abs(lonbounds[1]-lonbounds[0])>2:
    ratio=abs(latbounds[1]-latbounds[0])/abs(lonbounds[1]-lonbounds[0])
    midlon=lonbounds[0]+ 0.5*(lonbounds[1]-lonbounds[0])
    newstartlon=lonbounds[0]+ratio*abs(midlon-lonbounds[0])/1.2
    newstoplon=lonbounds[1]-ratio*abs(midlon-lonbounds[0])/1.2
    lonbounds=[newstartlon,newstoplon]
if abs(lonbounds[1]-lonbounds[0])/abs(latbounds[1]-latbounds[0])>2:
    ratio=abs(lonbounds[1]-lonbounds[0])/abs(latbounds[1]-latbounds[0])
    midlat=latbounds[0]+ 0.5*(latbounds[1]-latbounds[0])
    newstartlat=latbounds[0]+ratio*abs(midlat-latbounds[0])/1.2
    newstoplat=latbounds[1]-ratio*abs(midlat-latbounds[0])/1.2
    latbounds=[newstartlat,newstoplat]
latbounds.sort()
lonbounds.sort()
#yeah just don't mess with that unless making manual boundaries, in form of low to high
# lonbounds=[-7.5,-1.5]
# latbounds=[46,51.5]



m = Basemap(projection='cyl', llcrnrlon=lonbounds[0], \
    urcrnrlon=lonbounds[1],llcrnrlat=latbounds[0],urcrnrlat=latbounds[1], \
    resolution=resolu)

x,y = m(start[1], start[0])
if m.is_land(x, y)==True: 
    print("Start point is on land")
    quit()

x,y = m(stop[1], stop[0])
if m.is_land(x, y)==True: 
    print("Stop point is on land")
    quit()


#function that gets boatspeed as a function of position and heading

def boatspeed(heading,lon,lat,datau,datav,lats,lons,windspeeds,windangles,polar):
    #this because rolling over date line and poles? not sure if lat is needed but just in case.
    if lon<-180:
        lon+=360
    if lon>180:
        lon-=360

    if lat<-90:
        lat+=180
    if lat>90:
        lon-=180


    #find the indices and then values of nearest datapoint to position
    #evidently the data are of the form data(y,x) aka data(lat,lon)
    latidx = find_nearest(lats,lat)
    lonidx = find_nearest(lons,lon)
    uwind=datau[latidx][lonidx]
    vwind=datav[latidx][lonidx]
    windspeed = math.sqrt(uwind**2 + vwind**2)

    #compute apparent wind angle, -180 to 180 handled by periodic boundary condition
    #this wind angle DOES NOT account for boat velocity, which is good becvause
    #that is how boat polars are defined
    #atan2 takes (y,x) usually, but because I have north defined as 0 angle instead of east,
    #and clockwise positive instead of counterclockwise, i use atan2(x,y)
    windangle = math.atan2(uwind,vwind)
    apparentdegrees = 57.2958*(windangle-heading)

    #this if is to handle the fact that wind vector points in wind velocity direction,
    #but apparent wind vector ppints in direction wind is comign from
    if apparentdegrees <0:
        apparentdegrees+=180
    else:
        apparentdegrees-=180
    apparentdegrees -= 360*round(apparentdegrees/360)

    #we now know apparent wind and windspeed at our position
    #now find boatspeed from the polar diagram
    dictwind = windspeeds[find_nearest(windspeeds,windspeed)]
    listwindspeeds=[]
    dictangle = windangles[find_nearest(windangles,abs(apparentdegrees))]
    listwindangles=[]
    if dictwind == 5 or dictwind==20:
        velocity = polar[str(dictwind)][str(dictangle)]

    elif dictangle == 0 or dictangle == 180:
        velocity = polar[str(dictwind)][str(dictangle)]

    else:
        for item in windspeeds:
            listwindspeeds.append(item)

        if dictwind> windspeed:
            adjustdictwind= windspeeds[listwindspeeds.index(dictwind)-1]
        else: adjustdictwind= windspeeds[listwindspeeds.index(dictwind)+1]

        for item in windangles:
            listwindangles.append(item)

        if dictangle> abs(apparentdegrees):
            adjustdictangle= windangles[listwindangles.index(dictangle)-1]
        else: adjustdictangle= windangles[listwindangles.index(dictangle)+1]

        basevelocity = polar[str(dictwind)][str(dictangle)]
        adjustvelocity = polar[str(dictwind)][str(adjustdictangle)]
        #print((apparentdegrees-dictangle),apparentdegrees-adjustdictangle)
        velocity = abs((abs(apparentdegrees)-adjustdictangle)/(adjustdictangle-dictangle))*basevelocity + abs((abs(apparentdegrees)-dictangle)/(adjustdictangle-dictangle))*adjustvelocity
    return velocity

#finds the value in the lat/lon array
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

# Function to convert the date format 
#for downloading data as we go
def convert24(str1): 
     
    # Checking if last two elements of time 
    # is AM and first two elements are 12 
    if str1[-2:] == "AM" and str1[:2] == "12": 
        return "00" + str1[2:-2] 
         
    # remove the AM     
    elif str1[-2:] == "AM": 
        return str1[:-2] 
     
    # Checking if last two elements of time 
    # is PM and first two elements are 12 
    elif str1[-2:] == "PM" and str1[:2] == "12": 
        return str1[:-2] 
         
    else: 
         
        # add 12 to hours and remove PM 
        return str(int(str1[:2]) + 12) + str1[2:8] 

def write_gpx(coordinates, filename="data.gpx"):
  """
  Writes a list of coordinates to a GPX file.

  Args:
      coordinates: A list of lists, where each inner list contains [latitude, longitude] coordinates.
      filename (optional): The name of the output GPX file. Defaults to "data.gpx".
  """

  with open(filename, 'w') as f:
    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<gpx version="1.1" creator="Python Script">\n')
    f.write('  <trk>\n')
    f.write('    <name>My Track</name>\n')
    for coord in coordinates:
      latitude, longitude = coord
      f.write('      <trkseg>\n')
      f.write(f'        <trkpt lat="{latitude:.6f}" lon="{longitude:.6f}">\n')
      f.write('        </trkpt>\n')
      f.write('      </trkseg>\n')
    f.write('  </trk>\n')
    f.write('</gpx>\n')


#polar diagram dictionary holding windspeeed dictionaries which each have angle speeds
#two arrays which can be searchsorted to find nearest dictionary item to our windspeed&angle
windspeeds = np.array([5,8,10,12,15,20])
windangles = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
amelpolar = {

    '-500' : {
        '0' : 5 ,
        '10' : 5 ,
        '20' : 5 ,
        '30' : 5 ,
        '40' : 5 ,
        '50' : 5 ,
        '60' : 5 ,
        '70' : 5 ,
        '80' : 5 ,
        '90' : 5 ,
        '100' : 5 ,
        '110' : 5 ,
        '120' : 5 ,
        '130' : 5 ,
        '140' : 5 ,
        '150' : 5 ,
        '160' : 5 ,
        '170' : 5 ,
        '180' : 5
    } ,

    '5' : {
        '0' : 0 ,
        '10' : 0 ,
        '20' : 0 ,
        '30' : 0 ,
        '40' : 0 ,
        '50' : 2.5 ,
        '60' : 3 ,
        '70' : 3.2 ,
        '80' : 3.4 ,
        '90' : 3.4 ,
        '100' : 3.5 ,
        '110' : 3.5 ,
        '120' : 3.5 ,
        '130' : 3.1 ,
        '140' : 3 ,
        '150' : 3 ,
        '160' : 2.8 ,
        '170' : 3 ,
        '180' : 3.3
    } ,

    '8' : {
        '0' : 0 ,
        '10' : 0 ,
        '20' : 0 ,
        '30' : 0 ,
        '40' : 0 ,
        '50' : 3.5 ,
        '60' : 4.6 ,
        '70' : 4.7 ,
        '80' : 5 ,
        '90' : 5.3 ,
        '100' : 5.2 ,
        '110' : 5.1 ,
        '120' : 5 ,
        '130' : 4.8 ,
        '140' : 4.4 ,
        '150' : 4.2 ,
        '160' : 4.4 ,
        '170' : 4.8 ,
        '180' : 5
    } ,

    '10' : {
        '0' : 0 ,
        '10' : 0 ,
        '20' : 0 ,
        '30' : 0 ,
        '40' : 0 ,
        '50' : 5.2 ,
        '60' : 5.4 ,
        '70' : 5.9 ,
        '80' : 6.3 ,
        '90' : 6.5 ,
        '100' : 6.5 ,
        '110' : 6.3 ,
        '120' : 5.7 ,
        '130' : 5.2 ,
        '140' :  4.8 ,
        '150' : 4.2 ,
        '160' : 4.2 ,
        '170' : 4.8 ,
        '180' : 5.5
    } ,

    '12' : {
        '0' : 0 ,
        '10' : 0 ,
        '20' : 0 ,
        '30' : 0 ,
        '40' : 0 ,
        '50' : 6 ,
        '60' : 6.3 ,
        '70' : 6.4 ,
        '80' : 6.8 ,
        '90' : 7 ,
        '100' : 7 ,
        '110' : 7 ,
        '120' : 6.8 ,
        '130' : 6.4 ,
        '140' :  6.1 ,
        '150' : 6 ,
        '160' : 6.2 ,
        '170' : 6.3 ,
        '180' : 6.5
    } ,

    '15' : {
        '0' : 0 ,
        '10' : 0 ,
        '20' : 0 ,
        '30' : 0 ,
        '40' : 0 ,
        '50' : 7 ,
        '60' : 7.3 ,
        '70' : 7.6 ,
        '80' : 8 ,
        '90' : 8.3 ,
        '100' : 8.4 ,
        '110' : 8.4 ,
        '120' : 8.2 ,
        '130' : 7.7 ,
        '140' :  7 ,
        '150' : 6.5 ,
        '160' : 6.5 ,
        '170' : 7 ,
        '180' : 7
    } ,

    '20' : {
        '0' : 0 ,
        '10' : 0 ,
        '20' : 0 ,
        '30' : 0 ,
        '40' : 0 ,
        '50' : 7.6 ,
        '60' : 8 ,
        '70' : 8.3 ,
        '80' : 8.6 ,
        '90' : 8.9 ,
        '100' : 8.9 ,
        '110' : 9 ,
        '120' : 8.9 ,
        '130' : 8.5 ,
        '140' :  8 ,
        '150' : 7.2 ,
        '160' : 7.2 ,
        '170' : 7.7 ,
        '180' : 8.2
    } ,
}




############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
#begin data aquisition 
#download the first datafile before going into the isochrone iteration loop; if GFS hasn't modelled...
#... the most recent cycle yet, this first try will go back in time to the older cycle...
#... but that function doesn't work for going back into the previous month, so the code won't...
#... run on the first of the month early morning UTC (GMT)


backone=False
year = time.strftime("%Y", time.gmtime())
month = '05'
d= time.strftime("%d", time.gmtime())
hour=convert24(time.strftime("%I:%M:%S %p", time.gmtime()))
(h, minute, s) = hour.split(':')
decimaltime = int(h) + int(minute)/60 + int(s)/3600

#d='01'
#decimaltime=23.5
if decimaltime>18:
    cycle='18'
elif decimaltime>12:
    cycle='12'
elif decimaltime>6:
    cycle='06'
else:
    cycle='00'

downloadfile=str(math.floor(decimaltime-int(cycle)))
if len(downloadfile)==1:
    downloadfile="00"+downloadfile
if len(downloadfile)==2:
    downloadfile="0"+downloadfile

if os.path.isfile("/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile) == True: print("gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile+" file was already local")

else:
    try:
        testfile = urllib.request.urlretrieve("https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs."+year+month+d+"%2F"+cycle+"%2Fatmos&file=gfs.t"+cycle+"z.pgrb2.0p25.f384&var_UGRD=on&var_VGRD=on&lev_20_m_above_ground=on", \
            "/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile)
        print("Downloading: gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile)
    except urllib.error.URLError as e:
        print("the most recent dataset is not available yet, going back one cycle")
        if cycle=='00':
            if d!='01':
                d=int(d)
                d-=1
                d=str(d)
                if len(d)==1:
                    d="0"+d
                decimaltime+=24
                cycle='18'
            else:
                print("this is gonna be a problem for using the program on the first of a month at beginning of day UTC")
                month='04'
                d= '30'
                decimaltime+=24
                cycle='18'

        else:
            cycle=str(int(cycle)-6)
            if len(cycle)==1:
                cycle="0"+cycle

        downloadfile=str(round(decimaltime-int(cycle)))
        if len(downloadfile)==1:
            downloadfile="00"+downloadfile
        if len(downloadfile)==2:
            downloadfile="0"+downloadfile

        if os.path.isfile("gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile) == False:
            testfile = urllib.request.urlretrieve("https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs."+year+month+d+"%2F"+cycle+"%2Fatmos&file=gfs.t"+cycle+"z.pgrb2.0p25.f"+downloadfile+"&var_UGRD=on&var_VGRD=on&lev_20_m_above_ground=on", \
                "/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile)
        else:
            print("this file already existed")




#now we load in the data for the first round
plt.figure(figsize=(12,8))
# Set the file name of the grib file

#grib = 'gfs.t18z.sfluxgrbf000.grib2'
#grib = 'feb14grib_big_waimea_storm'
#grib = "/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_18_ymd_20240501.010"
grib = "/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile
#grib='/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_00_ymd_20240503.006'

grbs = pygrib.open(grib)
grbu = grbs.select()[0]
grbv = grbs.select()[1]
datau = grbu.values
datav = grbv.values
#convert m/s to nm/h
datau = datau/0.51444444444444
datav = datav/0.51444444444444
#I'm fairly certain that example at positive 5 wind U means windspeed 5 coming from the east

#define a listish of points along the domain
lons = np.linspace(float(grbu['longitudeOfFirstGridPointInDegrees']), \
float(grbu['longitudeOfLastGridPointInDegrees']), int(grbu['Ni']) )
lats = np.linspace(float(grbu['latitudeOfFirstGridPointInDegrees']), \
float(grbu['latitudeOfLastGridPointInDegrees']), int(grbu['Nj']) )

#do it agian just becvause the shiftgrid function below needs it
lonsv = np.linspace(float(grbv['longitudeOfFirstGridPointInDegrees']), \
float(grbv['longitudeOfLastGridPointInDegrees']), int(grbv['Ni']) )
latsv = np.linspace(float(grbv['latitudeOfFirstGridPointInDegrees']), \
float(grbv['latitudeOfLastGridPointInDegrees']), int(grbv['Nj']) )


# need to shift data grid longitudes from (0..360) to (-180..180)
datau, lons = shiftgrid(180., datau, lons, start=False)
datav, lonsv = shiftgrid(180., datav, lonsv, start=False)
grid_lon, grid_lat = np.meshgrid(lons, lats) #regularly spaced 2D grid

lats=np.flip(lats)
datau=np.flipud(datau)
datav=np.flipud(datav)






currentpositions=[start]
historicalpositions=[]
iteration=1
historicalfrom=[]
plotpath=True
counter=0
while arrived==False:
    if iteration>1:
        if int(downloadfile) +dt< 384:
            downloadfile=str(int(downloadfile)+dt)
        else: print("WARNING: Voyage will take more than 16 days, the end of the GFS model")
        if len(downloadfile)==1:
            downloadfile="00"+downloadfile
        if len(downloadfile)==2:
            downloadfile="0"+downloadfile

        #this if while is to handle the fact that GFS goes from 1hr to 3hr resolution after 120 hours
        if int(downloadfile)>120:
            while int(downloadfile)%3!=0:
                downloadfile = str(int(downloadfile)+1)
        if os.path.isfile("/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile) == True: print("gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile+" file was already local")
        else:
          testfile = urllib.request.urlretrieve("https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs."+year+month+d+"%2F"+cycle+"%2Fatmos&file=gfs.t"+cycle+"z.pgrb2.0p25.f"+downloadfile+"&var_UGRD=on&var_VGRD=on&lev_20_m_above_ground=on", \
               "/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile)
          print("Downloading: gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile)

        #grib='/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_00_ymd_20240503.006'
        grib = "/Users/nilsmelbourne/Documents/isochronal/gribs/gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile
        grbs = pygrib.open(grib)
        grbu = grbs.select()[0]
        grbv = grbs.select()[1]
        datau = grbu.values
        datav = grbv.values
        #convert m/s to nm/h
        datau = datau/0.51444444444444
        datav = datav/0.51444444444444
        #I'm fairly certain that example at positive 5 wind U means windspeed 5 coming from the east

        #define a listish of points along the domain
        lons = np.linspace(float(grbu['longitudeOfFirstGridPointInDegrees']), \
        float(grbu['longitudeOfLastGridPointInDegrees']), int(grbu['Ni']) )
        lats = np.linspace(float(grbu['latitudeOfFirstGridPointInDegrees']), \
        float(grbu['latitudeOfLastGridPointInDegrees']), int(grbu['Nj']) )

        #do it agian just becvause the shiftgrid function below needs it
        lonsv = np.linspace(float(grbv['longitudeOfFirstGridPointInDegrees']), \
        float(grbv['longitudeOfLastGridPointInDegrees']), int(grbv['Ni']) )
        latsv = np.linspace(float(grbv['latitudeOfFirstGridPointInDegrees']), \
        float(grbv['latitudeOfLastGridPointInDegrees']), int(grbv['Nj']) )


        # need to shift data grid longitudes from (0..360) to (-180..180)
        datau, lons = shiftgrid(180., datau, lons, start=False)
        datav, lonsv = shiftgrid(180., datav, lonsv, start=False)
        grid_lon, grid_lat = np.meshgrid(lons, lats) #regularly spaced 2D grid

        lats=np.flip(lats)
        datau=np.flipud(datau)
        datav=np.flipud(datav)


    historicalpositions.append(currentpositions)
    if iteration>1: historicalfrom.append(currentfrom)
    index_for_i=len(currentpositions)
    #newpositions is all points radiated from current boundary, pre-sorting
    #same for newdistances
    #newangles is a matching list where element j is the angle of point j relative to core start
    newpositions=[]
    newangles=[]
    newdistances=[]
    newfrom=[]
    #the i corresponds to index of each boundary point from previous iteration
    for i in range(index_for_i):

        direction= -np.pi
        #j is index each ray from point i
        if iteration>1: first_iteration_scaling= 1

        for j in range(radialsend*first_iteration_scaling):
            #first iteration sends more points than the rest, scaled by this variable
            velocity= boatspeed(direction,currentpositions[i][1],currentpositions[i][0], \
                                    datau,datav,lats,lons,windspeeds,windangles,amelpolar)

            #the velocity comes out as nm/h. to convert to lat/lon velocities:
            #latvelocity = nm/h * 1/60
            #lonvelocity needs to account for contraction at poles:
            #these are now in degrees/h
            lonvelocity= (1/60)*np.cos((currentpositions[i][0])/90)*np.sin(direction)*velocity
            latvelocity= (1/60)*np.cos(direction)*velocity


            newlat=dt*latvelocity+currentpositions[i][0]
            newlon=dt*lonvelocity+currentpositions[i][1]
            #This is  checking for in land or not
            x,y = m(newlon, newlat)
            if m.is_land(x, y)==False:
                newpositions.append([newlat,newlon])

                lat_distance= newlat-start[0]
                lon_distance= newlon-start[1]


                #this use of atan2 is a bit complicated. atan2 takes atan2(y,x) usually,
                #and lat,lon is y,x, but because i define north as angle 0 instead of
                #east the way we normally do, have to switch positions of (yx) to (xy)
                newangles.append(math.atan2(lon_distance,lat_distance))

                newdistances.append(np.sqrt(lat_distance**2+lon_distance**2))

                end_lat_d=newlat-stop[0]
                end_lon_d=newlon-stop[1]
                end_d=np.sqrt(end_lat_d**2+end_lon_d**2)
                newfrom.append(i)
                if end_d< 0.4:
                    print("Arrived in iteration",iteration)
                    arrived=True
                    enditeration=iteration
                    endboundarypoint=i
                    endposition=[newlat,newlon]
                if arrived==True:
                    break
            direction+= 2*np.pi/(radialsend*first_iteration_scaling)
    currentpositions=[]
    currentfrom=[]
    bottom= -np.pi
    if iteration>1: first_iteration_scaling= int(np.sqrt(40*iteration))
    for k in range(boundarypoints*first_iteration_scaling):
        #bottom is bottom of the angular bin in which we are testing points' distances
        top= bottom+np.pi*2/(boundarypoints*first_iteration_scaling)
        bin_distance=0
        appendnewposition=False
        for l in range(len(newpositions)):
            if newangles[l]>=bottom and newangles[l]<top and newdistances[l]>bin_distance:
                bin_distance = newdistances[l]
                idx=l
                appendnewposition=True
        if appendnewposition==True:
            currentpositions.append(newpositions[idx])
            currentfrom.append(newfrom[idx])
        bottom=top
    if arrived==True: historicalpositions.append(currentpositions)

    print("Finished iteration ",iteration, "using grib: , gribcycle_"+cycle+"_ymd_"+year+month+d+"."+downloadfile)
    iteration+=1
    if iteration==200:
        arrived=True
        plotpath=False
        print("Runtime exceeded")
    plotpositions=np.array(currentpositions)
    if iteration%plottingidx==0:
        counter+=1
        #plt.scatter(plotpositions[:,1],plotpositions[:,0],s=2,color='cornflowerblue')
        if counter%4==0:
            #plt.plot(plotpositions[:,1],plotpositions[:,0],color='cornflowerblue')
            plt.scatter(plotpositions[:,1],plotpositions[:,0],s=.3,color='royalblue')

        else:
            #plt.plot(plotpositions[:,1],plotpositions[:,0],color='cornflowerblue')
            plt.scatter(plotpositions[:,1],plotpositions[:,0],s=.2,color='cornflowerblue')
        #plt.plot(currentpositions[:,1],currentpositions[:,0])




for i in range(len(historicalpositions[-1])-1):
    d=((historicalpositions[-1][i+1][0]-historicalpositions[-1][i][0])**2+(historicalpositions[-1][i+1][1]-historicalpositions[-1][i][1])**2)**0.5


if plotpath==True:
    path=[stop,endposition]
    idx=endboundarypoint
    for i in range(enditeration-1):
        it_idx=enditeration-2-i
        idx=historicalfrom[it_idx][idx]
        path.append(historicalpositions[it_idx][idx])
    write_gpx(path)
    path=np.array(path)
    plt.plot(path[:,1],path[:,0],color='orangered', linewidth=2)

plt.scatter(start[1],start[0],color='limegreen',s=100)
plt.scatter(stop[1],stop[0],color="red",s=100)




# # #THIS BLOCK PLOTS VECTORS OF WIND FROM THE LAST LOADED GRIB
# newlons=[]
# bottom_lon=False
# top_lon=False
# for n in range(int(len(lons)/descaling)):
#     newn=n*descaling
#     if lons[newn]>= lonbounds[0] and lons[newn]<= lonbounds[1]:
#         newlons.append(lons[newn])
#         if bottom_lon==False:
#             blonidx=n
#             bottom_lon=True
#     if lons[newn]>lonbounds[1]:
#         if top_lon==False:

#             tlonidx=n
#             top_lon=True
# newlons=np.array(newlons)

# newlats=[]
# bottom_lat=False
# top_lat=False
# for n in range(int(len(lats)/descaling)):
#     newn=n*descaling
#     if lats[newn]>= latbounds[0] and lats[newn]<= latbounds[1]:
#         newlats.append(lats[newn])
#         if bottom_lat==False:

#             blatidx=n
#             bottom_lat=True
#     if lats[newn]>latbounds[1]:
#         if top_lat==False:

#             tlatidx=n
#             top_lat=True
# newlats=np.array(newlats)
# #newdatau=np.ndarray(shape=(int(len(datau)/descaling),int(len(datau[0])/descaling)))
# #newdatav=np.ndarray(shape=(int(len(datav)/descaling),int(len(datav[0])/descaling)))
# newdatau=[]
# newdatav=[]
# for n in range(blatidx,tlatidx):
#     rowu=[]
#     rowv=[]
#     for l in range(blonidx,tlonidx):
#         #newdatau[n][l]=datau[n*descaling][l*descaling]
#         #newdatav[n][l]=datav[n*descaling][l*descaling]
#         rowu.append(datau[n*descaling][l*descaling])
#         rowv.append(datav[n*descaling][l*descaling])
#     newdatau.append(rowu)
#     newdatav.append(rowv)
# newdatau=np.array(newdatau)
# newdatav=np.array(newdatav)

#cs = m.pcolormesh(x,y,data,cmap=plt.cm.jet)#shading='flat'

#m.drawlsmask(land_color='coral',ocean_color='white',lakes=True)
m.drawcoastlines()
m.fillcontinents(color='olive')
m.drawparallels(np.arange(-90.,120.,10.),labels=[1,0,0,0])
m.drawmeridians(np.arange(-180.,180.,20.),labels=[0,0,0,1])


#bnewlats=np.flip(newlats)
#plt.colorbar(cs,orientation='vertical', shrink=0.5)
# plt.quiver(newlons, newlats, newdatau, newdatav, color='black',width=0.0017)
#plt.streamplot(newlons, newlats, newdatau, newdatav, density=4,color='g')

#strm = plt.streamplot(newlons, newlats, newdatau, newdatav, color=np.sqrt(newdatau**2+newdatav**2),linewidth=2, cmap='autumn')

plt.title('Travel Time: '+str(enditeration*dt/24)+" days \n "+str(dt)+" hour isochrone resolution") # Set the name of the variable to plot

plt.show()
#plt.savefig('flatteryHNL.png') # Set the output file name
