#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 09:57:52 2020

@author: szutor péter

"""
# Reverse erosion in Badacsony with laplace convoltion kernel
# DEM points in badacsonypontok.csv
# Every step will be saved in jpg, and it can be join to an avi with



from scipy import signal
import numpy as np
import pandas as pd
import open3d as o3d
import time
import matplotlib.pyplot as plt


def terepszin():
    global pcl
    global pct
    ct=np.zeros((pct.shape[0],3))

    kek=np.array([0,0,0.97])
    zold=np.array([0,0.2,0])
    szurke=np.array([0.2,0.2,0.2])
    zoldhat=alja+(teteje-alja)/2
    ct[np.where(pct[:,2]<zoldhat),1]=(pct[np.where(pct[:,2]<zoldhat),2]-alja)*0.8/(zoldhat-alja)+0.2
    fok=(pct[np.where(pct[:,2]>=zoldhat),2]-zoldhat)*0.8/(teteje-zoldhat)+0.2
    ct[np.where(pct[:,2]>=zoldhat),0]=fok
    ct[np.where(pct[:,2]>=zoldhat),1]=fok
    ct[np.where(pct[:,2]>=zoldhat),2]=fok
    ct[np.where(pct[:,2]<=alja+(teteje-alja)/100)]=kek
    pcl.colors=o3d.utility.Vector3dVector(ct)
    return

zarany=12 # magasság lecsökkentése
forgatas=45 #elfogatási szög
df2 = pd.read_csv("badacsonypontok.csv",sep=";")
dfkell=df2[['X','Y','eu_dem_v11']]
pct=dfkell.to_numpy()
ymeret=np.unique(pct[:,1]).shape[0]
xmeret=np.unique(pct[:,0]).shape[0]
kep=np.zeros((ymeret,xmeret))
for i in range(0,ymeret):
    for j in range(0,xmeret):
        kep[i,j]=pct[i*xmeret+j,2]
        
alja=np.amin(kep)        
#lapkernel=np.array([[0,0,1,0,0],[0,1,2,1,0],[1,2,-16,2,1],[0,1,2,1,0],[0,0,1,0,0]])/-80 #laplace
lapkernel=np.array([[1,4,7,4,1],[4,16,26,16,4],[7,26,41,26,7],[4,16,26,16,4],[1,4,7,4,1]])/-(580*268) #gauss
rotmatrix=np.array([[1,0,0,0],[0,np.cos(forgatas),np.sin(forgatas),0],[0,-1*np.sin(forgatas),np.cos(forgatas),0],[0,0,0,1]],dtype=np.float64)
vis = o3d.visualization.VisualizerWithKeyCallback()
vis.create_window()
pcl = o3d.geometry.PointCloud()
ermin=np.amin(kep)
pcl.points = o3d.utility.Vector3dVector(pct)  
#o3d.visualization.draw_geometries([pcl])

#Most visszaalakitjuk pontfelhővé
j=0
for x in range(xmeret):
    for y in range(ymeret):
        pct[j,0]=x
        pct[j,1]=ymeret-y
        pct[j,2]=kep[y,x]/zarany
        j+=1
pcl.points = o3d.utility.Vector3dVector(pct)
alja=np.amin(pct[:,2])
teteje=np.amax(pct[:,2])
pcl.transform(rotmatrix)
terepszin()
vis.add_geometry(pcl)
vis.update_geometry(pcl)
vis.poll_events()
vis.reset_view_point(False)
vis.update_renderer()
#csak a zoom miatt pár lépés
for i in range(0,25):
    vis.update_geometry(pcl)
    vis.poll_events()
    vis.update_renderer()
    time.sleep(0.5)

#vis.run()
#time.sleep(4)
if 1==1:
    for i in range(0,600):
        egylepes=signal.oaconvolve(kep,lapkernel,mode='same')
        kep=kep-egylepes[:ymeret,:xmeret]
        if 1==1: #kell-e az alját fixálni
            kep[np.where(kep<ermin)]=ermin

        #Most visszaalakitjuk pontfelhővé
        j=0
        for x in range(xmeret):
            for y in range(ymeret):
                pct[j,0]=x
                pct[j,1]=ymeret-y
                pct[j,2]=kep[y,x]/zarany
                j+=1
        pcl.points = o3d.utility.Vector3dVector(pct)    
        pcl.transform(rotmatrix)
        terepszin()
        vis.update_geometry(pcl)
        vis.poll_events()
        #vis.reset_view_point(True)
        vis.update_renderer()
        vis.capture_screen_image('jpegs/gaussbadacsony'+'{:04d}'.format(i)+'.jpg')
#        plt.imshow(kep)
#        plt.show()
        #o3d.visualization.draw_geometries([pcl])
        time.sleep(0.01)

    

vis.destroy_window()