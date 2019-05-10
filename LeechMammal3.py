import pandas as pd
import numpy as np
import openpyxl
import math

root=('C:/Users/hp/Documents/Leeches/')

def summarise_and_combine_leech_mammal_dataframes(leech_infile,leech_coordinates,leech_radius,\
	mammal_radius,lidar,mammal_abundance,mammal_trap_effort,date,outfile):
	'''
	Takes the leech data and pools data by site. 
	You can set the leech and mammal radius required for points to be merged. 
	'''

	#### load in the datasheets required
	xl=pd.ExcelFile(leech_infile)
	leech2015=xl.parse('2015') #leech 2015 data

	xl=pd.ExcelFile(leech_infile)
	leech2016=xl.parse('2016') # leech 2016 data 

	xl=pd.ExcelFile(mammal_abundance)
	mammal=xl.parse('MammalsONLY') # mammal data imported as mammal

	xl=pd.ExcelFile(date)
	date_data=xl.parse('Sheet1')

	leech_coord=pd.read_csv(leech_coordinates) # read in leech coordinates

	effort=pd.read_csv(mammal_trap_effort)#read in mammal trapping effort

	lidar_data=pd.read_csv(lidar) # read in lidar dataset

	#adds mammal trapping effort and coordinates to the abundance data
	mammal=add_mammal_trap_nights(mammal,effort)	

	#generate groups to later merge close together leech points
	mammal_merged=merge_mammal_points(mammal,mammal_radius)

	#generate coordinates and merge 2015 and 16 into one row. also adds in lidar data and survey date data
	leech=leech_coord_function(leech2015,leech2016,leech_coord,lidar_data,date_data)

	#summarise the mammal data by merging the leech points
	mammal_summarised=summarise_mammal_data(mammal_merged)

	mammal_summarised.to_csv('C:/Users/hp/Documents/Leeches/MammalMerged.csv')

	#add in leech data for close by points
	(leechtest,mammal_leech_combined)=add_leech_data(leech,mammal_summarised,leech_radius)

	check_minimum_pairwise_distance(mammal_leech_combined,leech_radius)

	mammal_leech_combined.to_csv(outfile)

	mammal_merged.to_csv('C:/Users/hp/Documents/Leeches/MammalTest.csv')

	leech.to_csv('C:/Users/hp/Documents/Leeches/LeechTest.csv')

def add_mammal_trap_nights(abundance,effort):
	'''
	appends mammal trapping hours and days onto the end of the combined leech_mammal df
	'''
	#make the site ID column matchable to hectare column in mammal abundance
	for i in effort.index:
		effort.loc[i,'Site_ID']=effort.loc[i,'Site_ID'][3:]

	#add CTN and CTH (sampling effort) to abundance 
	for i in abundance.index:
		for j in effort.index:
			if abundance.loc[i,'Site']==effort.loc[j,'Site_ID']:
				abundance.loc[i,'CTNs']=effort.loc[j,'CTNs']
				abundance.loc[i,'CTHs']=effort.loc[j,'CTHs']
				abundance.loc[i,'Lat']=effort.loc[j,'Latitude']
				abundance.loc[i,'Long']=effort.loc[j,'Longitude']

	return abundance


def merge_mammal_points(mammal,merge_radius):
	'''
	creates a radius around eachpoint and merges mammal points that are within that radius of eachother
	note these can daisy chain together 
	'''

	#a variable to track whether a point has already merged with another point
	mammal['Merged']='No'

	##this pairs up the points
	for i in mammal.index:
		if math.isnan(mammal.loc[i,'Lat']) == False:
			for j in mammal.index:
				if math.isnan(mammal.loc[j,'Lat']) == False:
				#doesnt allow the same lines to match up 
					if i != j:
						x=(abs(mammal.loc[i,'Lat']-mammal.loc[j,'Lat'])**2+abs(mammal.loc[i,'Long']-mammal.loc[j,'Long'])**2)
						x=math.sqrt(x)
						if x*111111<merge_radius:
							#check to see if the row (i) has already been assigned to a group or is its own group
							if mammal.loc[i,'Merged']=='No':
							#if this point has not been matched up before, group them together
								if mammal.loc[j, 'Merged']=='No':
									mammal.loc[i,'Merged']=i
									mammal.loc[j,'Merged']=i
								#if this point has been matched up before, change i merged group to j merged group
								elif mammal.loc[j,'Merged']!='No':
									mammal.loc[i,'Merged']=mammal.loc[j,'Merged']
							else:
								#if j hasnt been assigned a group yet
								if mammal.loc[j, 'Merged']=='No':
									mammal.loc[j,'Merged']=mammal.loc[i,'Merged']
								#if j has been assigned a group and i has been assigned a group, 
								# overwrite all of j groups with i group
								elif mammal.loc[j,'Merged']!='No':
									group_to_replace=mammal.loc[j,'Merged']
									for k in mammal.index:
										if mammal.loc[k,'Merged']==group_to_replace:
											mammal.loc[k,'Merged']=mammal.loc[i,'Merged']
	
	#if no match was found for a row, make it its own group
	for i in mammal.index:
		if mammal.loc[i,'Merged']=='No':
			mammal.loc[i,'Merged']=i

	return mammal




def leech_coord_function(leech2015,leech2016,leech_coord,lidar,date_data):

	date=date_data
	
	#combine the leech data from different years
	leech2015['TOTAL.B16']=leech2016['TOTAL.B16']
	leech2015['TOTAL.T16']=leech2016['TOTAL.B16']
	leech2015['Effort15']=leech2015['TOTAL.EFF']
	leech2015['Effort16']=leech2016['TOTAL.EFF']

	letters=['b','t','visit']

	for i in letters:
		for j in range(1,5):
			leech2015[i+str(j)+'.16']=leech2016[i+str(j)+'.16']

	leech=leech2015


	#the coordinates of the rivers in the abundance doc need formatting to add 0s e.g. R0-1 is R0-01
	for i in leech.index:
		if str(leech.loc[i,'second'])[:2]=='R0'or str(leech.loc[i,'second'])[:3]=='R30' or str(leech.loc[i,'second'])[:4]=='RLFE':
			split=leech.loc[i,'second'].split('_')
			leech.loc[i,'second']=split[0]+'-'+split[1]
			if str(leech.loc[i,'second'][-1])!=str(0):
				leech.loc[i,'second']=leech.loc[i,'second'][:-1]+'0'+leech.loc[i,'second'][-1]

	#lets get the leech coordinates in 
	for i in leech.index:
		for j in leech_coord.index:
			if str(leech.loc[i,'second'])==str(leech_coord.loc[j,'SiteNum']):
				leech.loc[i,'Lat']=leech_coord.loc[j,'Lat']
				leech.loc[i,'Long']=leech_coord.loc[j,'Long']

	#add in lidar data
	for i in leech.index:
		for j in lidar.index:
			if lidar.loc[j,'ID'][0]=='R':
				if str(leech.loc[i,'second'])==str(lidar.loc[j,'ID']):
					leech.loc[i,'canopy_height_mean']=lidar.loc[j,'canopy_height_mean']
					leech.loc[i,'canopy_height_moran']=lidar.loc[j,'canopy_height_moran']
					leech.loc[i,'pai_mean']=lidar.loc[j,'pai_mean']
			else:
				if str(leech.loc[i,'second'])==str(lidar.loc[j,'ID'][-3:]):
					leech.loc[i,'canopy_height_mean']=lidar.loc[j,'canopy_height_mean']
					leech.loc[i,'canopy_height_moran']=lidar.loc[j,'canopy_height_moran']
					leech.loc[i,'pai_mean']=lidar.loc[j,'pai_mean']


	#add in date data
	for i in leech.index:
		for j in date.index:
			if leech.loc[i,'second']==date.loc[j,'plotNum']:
				for k in date.columns:
					if k != 'plotNum':
						leech.loc[i,k]=date.loc[j,k]

	return leech

def summarise_mammal_data(mammal):

	#generate list of mammal names and the trapping efforts
	mammal_list=mammal.columns[3:-3].tolist()

	#create new summary dataframe
	mammal_summary=mammal.groupby('Merged')[['hectareMATCH']].first()

	#add in the mammal abundances and efforts
	for i in mammal_list:
		mammal_summary=mammal_summary.join(mammal.groupby('Merged')[[i]].sum())

	#add in mean lat and long 
	mammal_summary=mammal_summary.join(mammal.groupby('Merged')[['Lat']].mean())
	mammal_summary=mammal_summary.join(mammal.groupby('Merged')[['Long']].mean())

	return mammal_summary


def add_leech_data(leech,mammal,leech_radius):

	#generate a variable to check if the leech point has been matched with a merged mammal point yet
	leech['Match']='No'
	leech['DuplicateMatch']='No'

	##this finds leech points to append to leech groups
	for i in leech.index:
		for j in mammal.index:
			x=(abs(mammal.loc[j,'Lat']-leech.loc[i,'Lat'])**2+\
				abs(mammal.loc[j,'Long']-leech.loc[i,'Long'])**2)
			x=math.sqrt(x)
			if x*111111<leech_radius:
				#if not already match up
				if leech.loc[i,'Match']=='No':
					leech.loc[i,'Match']=j
				#if already matched record in the duplicate match column
				else:
					leech.loc[i,'DuplicateMatch']='Yes'

	for i in leech.index:
		if leech.loc[i,'DuplicateMatch']=='Yes':
			raise ValueError('Duplicate Leech Matches')
			# print('Duplicate Leech Matches')

	###this takes the mammal matches and puts them into the leech dataframe
	leech_summary=leech.groupby('Match')[['TOTAL.B15']].sum()

	leech_sums=('TOTAL.T15','TOTAL.B16','TOTAL.T16','Effort15','Effort16','canopy_height_mean','canopy_height_moran','pai_mean')

	for i in leech_sums:
		leech_summary=leech_summary.join(leech.groupby('Match')[[i]].sum())

	leech_list=('TOTAL.B15','TOTAL.T15','TOTAL.B16','TOTAL.T16','Effort15','Effort16','canopy_height_mean','canopy_height_moran','pai_mean')

	for i in leech_summary.index:
		for j in mammal.index:
			if j == i:
				for k in leech_list:
					mammal.loc[j, k]=leech_summary.loc[i,k]

	mammal.dropna(inplace=True)

	#extracts the largest group of leeches other than 'no'
	max_pool_size=max(leech.groupby('Match')[['TOTAL.B15']].size()[:-1].tolist())

	#make a list to put in the column headers
	visit_headers=['a','b']

	letters=['b','t','e']

	years=['.15','.16']

	for l in years:
		#create columns for new effort, brown and tiger to go
		for k in letters:
			#create a variable to keep track of the number of blocks of new columns created
			new_cols=0
			for i in range(0,max_pool_size):
				for j in range(1,5):
					mammal[k+str(j+4*new_cols)+l]=''
					if k!='e':
						visit_headers.append(k+str(j+4*new_cols)+l)
				new_cols+=1

	#get rid of junk columns at start of list
	visit_headers=visit_headers[2:]

	#creates a variable to keep track of the number of sites inputed
	mammal['NumberOfSites']=0

	#add in the leech visits for rosies occupancy models to work
	for i in leech.index:
		for j in mammal.index:
			x=(abs(mammal.loc[j,'Lat']-leech.loc[i,'Lat'])**2+\
				abs(mammal.loc[j,'Long']-leech.loc[i,'Long'])**2)
			x=math.sqrt(x)
			if x*111111<leech_radius:
				mammal.loc[j,'NumberOfSites']+=1
				for k in range(1,5):
					for l in years:
						mammal.loc[j,'b'+str((mammal.loc[j,'NumberOfSites']-1)*4+k)+l]=leech.loc[i,'b'+str(k)+l]
						mammal.loc[j,'t'+str((mammal.loc[j,'NumberOfSites']-1)*4+k)+l]=leech.loc[i,'t'+str(k)+l]
						mammal.loc[j,'e'+str((mammal.loc[j,'NumberOfSites']-1)*4+k)+l]=leech.loc[i,'visit'+str(k)+l]
						mammal.loc[j,'second.'+str((mammal.loc[j,'NumberOfSites']-1)*4+k)+l]=leech.loc[i,'second']

	for l in years:
		mammal['NumberOfSites']=0
		for i in leech.index:
			for j in mammal.index:
				x=(abs(mammal.loc[j,'Lat']-leech.loc[i,'Lat'])**2+\
					abs(mammal.loc[j,'Long']-leech.loc[i,'Long'])**2)
				x=math.sqrt(x)
				if x*111111<leech_radius:
					mammal.loc[j,'NumberOfSites']+=1
					for k in range(1,5):
							mammal.loc[j,'date.'+str((mammal.loc[j,'NumberOfSites']-1)*4+k)+l]=leech.loc[i,'date.'+str(k)+l]

	
	#converts leech abundance to presence in each visit col (binary)
	for i in mammal.index:
		for j in visit_headers:
			try:
				if mammal.loc[i,j]>0:
					mammal.loc[i,j]=1
			except:
				pass
			

	return (leech,mammal)

def check_minimum_pairwise_distance(points,leech_radius):
	'''
	Checks what the minimum distance is between merged points and halves it

	'''
	length=len(points['Lat'])

	#creates a dataframe with the points rows against the points rows, allows comparison of dataframe against itself.
	
	df=pd.DataFrame(index=points.index, columns=points.index)

	
	for i in points.index:
		for j in points.index:
			if i!=j:
				df.loc[i,j]=111111*math.sqrt(abs(points.loc[j,'Lat']-points.loc[i,'Lat'])**2+\
				abs(points.loc[j,'Long']-points.loc[i,'Long'])**2)
			else:
				df.loc[i,j]=100000

	x=1000000
				
	for i in df.index:
		for j in df.index:
			if df.loc[i,j]<x:
				x=df.loc[i,j]

	#prints a series of outputs to assess the efficacy of the pooling functions
	# print('Radius = '+str(leech_radius))
	print('Number of Leeches Pooled = '+str(len(points.index)))
	print('Number of leech points sampled = '+str(pd.DataFrame.sum(points['NumberOfSites'])))
	print('Max_Distance = '+str(x/2))


#run the script with set radii

summarise_and_combine_leech_mammal_dataframes(root+'LeechAbundance.xlsx',root+'LeechCoord.csv',\
		320,320,root+'LiDAR.csv',root+'MammalAbundance.xlsx',\
		root+'MammalTrapEffort.csv',root+'Dates.xlsx',root+'Mammal_Leech_Combined.csv')

