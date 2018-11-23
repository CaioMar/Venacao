
###Importing libraries
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from PIL import Image
from sklearn.neighbors import KNeighborsClassifier
from matplotlib.colors import ListedColormap
sns.set_style("whitegrid")

#modeltype
modeltype = str(sys.argv[1])

filepath = '/home/caio/Dropbox/Marcelo/Artigos/UpdatedVersion/leafdata'

###Reading data from the outputs of the venation programs
data = pd.read_csv(filepath + '/disthistogram'+modeltype, header=0 , delimiter='\t', names=['Harmonic','Arithmetic','Geometric','VeinOrder','NN VeinOrder','v2x','v2y','vx','vy','v2width','vwidth'])
ven = pd.read_csv(filepath + '/venconnections'+modeltype, delimiter='\t')
ven2 = pd.read_csv(filepath + '/ven2connections'+modeltype, delimiter='\t')

### Creates the count column based on how many times a complementary venation point appears as a nearest neighbor of
### a node on the primary venation
cols = list(data.columns[1:7])
cols += list(data.columns[9:]) ### generates a list with a columns that are going to be dropped later on
countdf = data.groupby(['vx','vy']).agg('count').drop(cols,axis=1).rename(columns={'Harmonic':'count'}) #creates a df with the count column
countdf.reset_index(level=0, inplace=True)
countdf.reset_index(level=0, inplace=True) # resets the vx and vy index as normal columns
data = pd.merge(data,countdf,on=['vx','vy']) #merges the count column dataframe with the original data df

###Some serious data wrangling I will better document later on
cven2 = ven2.copy()
cven2.rename(columns={'startx':'v2x','starty':'v2y'},inplace=True)
cven2 = cven2.groupby(by=['v2x','v2y']).agg('max').drop(['endx','endy'],axis=1)
cven2.reset_index(level=0,inplace=True)
cven2.reset_index(level=0,inplace=True)
data = pd.merge(data,cven2,on=['v2x','v2y'])
data['seglength'] = np.sqrt(np.square(data['v2x']-data['vx'])+np.square(data['v2y']-data['vy']))

cvenend = ven.copy()
cvenend.rename(columns={'endx':'vx','endy':'vy','disttoroot':'disttorootp'},inplace=True)
cvenend = cvenend.groupby(by=['vx','vy']).agg('max').drop(['startx','starty'],axis=1)
cvenend.reset_index(level=0,inplace=True)
cvenend.reset_index(level=0,inplace=True)
cven = ven.copy()
cven.rename(columns={'startx':'vx','starty':'vy','disttoroot':'disttorootp'},inplace=True)
cven = cven.groupby(by=['vx','vy']).agg('max').drop(['endx','endy'],axis=1)
cven.reset_index(level=0,inplace=True)
cven.reset_index(level=0,inplace=True)
cven = pd.concat([cven,cvenend])

###Dataframe employed in the plotting all primary venation points (some points are missing from data, that's wy this is necessary)
cvendata = cven.groupby(by=['vx','vy']).agg('max')
cvendata.reset_index(level=0,inplace=True)
cvendata.reset_index(level=0,inplace=True)
cvendata =  pd.merge(data,cvendata,on=['vx','vy'],how='right')

###In this step we estimate the count feature for the nodes that don't have any inlet venation nodes as nearest neighbor
###We employ a decision tree classifier since it's a simple model that tackles nonlinear problems
###We could skip this, but the graph wouldn't look nice
###This is just for aesthetics, the predictions don't need to be truly accurate
attr = ['disttorootp','vx','vy']
Xtrain = np.array(cvendata[attr][pd.notna(cvendata['count'])]).reshape(-1,len(attr))
Xtest = np.array(cvendata[attr][pd.isna(cvendata['count'])]).reshape(-1,len(attr))
ytrain = np.array(cvendata[['count']][pd.notna(cvendata['count'])]).reshape(-1,1)
ytest = np.array(cvendata[['count']][pd.isna(cvendata['count'])]).reshape(-1,1)
tree = DecisionTreeClassifier(max_depth=3)###3 to prevent too much overfitting
tree = tree.fit(Xtrain,ytrain)
cvendata['count'][pd.isna(cvendata['count']) ] = tree.predict(Xtest)

###The exact same thing as the above, but now we estimate the vein width.
###We employ a decision tree regreesion model since it's a simple enough model that tackles nonlinear problems
###We could skip this, but the graph wouldn't look nice
###This is just for aesthetics, the predictions don't need to be truly accurate
Xtrain = np.array(cvendata[attr][pd.notna(cvendata['vwidth'])]).reshape(-1,len(attr))
Xtest = np.array(cvendata[attr][pd.isna(cvendata['vwidth'])]).reshape(-1,len(attr))
ytrain = np.array(cvendata['vwidth'][pd.notna(cvendata['vwidth'])]).reshape(-1,1)
ytest = np.array(cvendata['vwidth'][pd.isna(cvendata['vwidth'])]).reshape(-1,1)
tree = DecisionTreeRegressor(max_depth=3)
tree = tree.fit(Xtrain,ytrain)
cvendata['vwidth'][pd.isna(cvendata['vwidth']) ] = tree.predict(Xtest)

###assings some important attributes to simple named variables
###so as to make the code more readable
ven2std = ven2.iloc[:,4].std()
ven2max = ven2.iloc[:,4].max()
ven2min = ven2.iloc[:,4].min()
venstd = ven.iloc[:,4].std()
venmax = ven.iloc[:,4].max()
venmin = ven.iloc[:,4].min()

####We do this eliminate the last leaf nodes from complementary tree, otherwise they wouldn't look nice
copyven2_1 = ven2.copy()
copyven2_2 = ven2.copy()
copyven2_2.rename(columns={'endx':'v2x','endy':'v2y'},inplace=True)
copyven2 = pd.merge(copyven2_2,data,on=['v2x','v2y'],how='inner')

###Discretizes the continuous distance to root feature in 4 class (mainly 3 classes, as class 0 is insignificant)
disc = list(data['disttoroot'].quantile(np.linspace(0,1,4)))
data['dtrootclass']=0
for i in range(0,len(disc)-1):
    data['dtrootclass'].loc[(data['disttoroot']>disc[i]) &  (data['disttoroot']<disc[i+1])] = float(i)+1.0


###Trains a 1st NN classifier which will be later employed in coloring the background of the output figure
Xtrain = np.array(data[['v2x','v2y']]).reshape(-1,2)
ytrain = np.array(data[['dtrootclass']]).reshape(-1)
lst = []
for i, j in enumerate(cvendata[['vx','vy']].groupby(['vx','vy']).agg('count').index):
    lst += [list(j)]
primary = np.array(lst).reshape(-1,2)
primarylabels = np.ones(primary.shape[0])*4
Xtrain = np.concatenate([Xtrain,primary])
ytrain = np.concatenate([ytrain,primarylabels])
knnclassifier = KNeighborsClassifier(n_neighbors=1).fit(Xtrain,ytrain)

###Generates and saves the figure as venation.png
fig = plt.figure(figsize=(12,12))
cmapseg = plt.cm.get_cmap('Reds')
rgba = cmapseg(np.square(data['disttoroot'])/(np.square(data['disttoroot']).max()-np.square(data['disttoroot']).min()))
fig = plt.figure(figsize=(10,10))

X_set, y_set = Xtrain, ytrain
X1, X2 = np.meshgrid(np.arange(start = X_set[:, 0].min() - 1, stop = X_set[:, 0].max() + 1, step = 0.01),
                     np.arange(start = X_set[:, 1].min() - 1, stop = X_set[:, 1].max() + 1, step = 0.01))
z = knnclassifier.predict(np.array([X1.ravel(), X2.ravel()]).T).reshape(X1.shape)
plt.contourf(X1, X2, z, zorder=1,
             alpha = 0.75, cmap = ListedColormap(('dimgray','dimgray', 'darkgrey','gainsboro','black')))
plt.xlim(X1.min(), X1.max())
plt.ylim(X2.min(), X2.max())

for i in range(data.shape[0]):
    plt.plot(data.iloc[i][5:9:2],data.iloc[i][6:9:2],'-',color=rgba[i],linewidth=np.power((data['seglength'].max()-data.iloc[i]['seglength'])/(data['seglength'].std()),0.3)-0.5,alpha=0.5,zorder=2)
for i in range(ven.shape[0]):
    plt.plot(ven.iloc[i][0:4:2],ven.iloc[i][1:4:2],'-',c='firebrick',markersize=1,linewidth=(venmax-ven.iloc[i][4])/venstd+0.5,alpha=0.5,zorder=2)
for i in range(copyven2.shape[0]):
    plt.plot(copyven2.iloc[i][0:4:2],copyven2.iloc[i][1:4:2],'-',c='aliceblue',linewidth=(ven2max - copyven2.iloc[i][4])/ven2std + 0.5,alpha=0.5,zorder=2)

plt.scatter(data['v2x'],data['v2y'],s=np.power(10*np.array(data['v2width'])+3,1.79),c=data['dtrootclass'], cmap=plt.get_cmap('Blues'),zorder=3)
plt.scatter(cvendata['vx'],cvendata['vy'],s=cvendata['count']+np.power(2*np.array(cvendata['vwidth']),2.82),c=(cvendata['disttorootp'].max()-cvendata['disttorootp'])/cvendata['disttorootp'].std(), cmap=plt.get_cmap('Reds'),zorder=3)
plt.xticks([])
plt.yticks([])
fig.savefig(fname = "venation"+modeltype+".png", dpi = 150)


X2range = X2.max()-X2.min()
pointarea = (X2range/X2.shape[0])*(X2range/X2.shape[1])
totalarea = np.square(X2range)
areaperregion = []
#areaproportionperclass = []
for i in range(0,5):
    areaperregion += [[i,np.array(z==i).sum()*pointarea,np.array(z==i).sum()*pointarea/totalarea]]
  #  areaproportionperclass += [[i,np.array(z==i).sum()*pointarea/totalarea]]
areaperregion = pd.DataFrame(np.array(areaperregion),columns=['regionclass','areaperregion','areaproportionperclass'])
areaperregion.to_csv('areaperregion'+modeltype+'.csv',index=False)
data.to_csv('originaldata_modified'+modeltype+'.csv',index=False)
cvendata.to_csv('originalprimarydata_modified'+modeltype+'.csv',index=False)
