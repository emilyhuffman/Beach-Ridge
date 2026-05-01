import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from shapely.geometry import Polygon
from scipy.integrate import trapezoid
from defined_function import xy_intersect,NaNlist,zerolist,find_D

def swash_wash(Qri,Qs,Tmax,dt,m,a,b,maxH,maxdis,atrigger,title,type_of_image,saveimage):
    
    #Input Parameters
    L=0 #m Initial width of beach (m)
    Ls=0 #m Foreridge setback distance (m)
    m1=0.06 #m/m Slope of ridge flank 1 (front flank)
    m2=0.006 #m/m Slope of ridge flank 2 (back flank)
    H=0 #set starting height of ridge (m)
    
    #Starting Conditions and Trackers
    Hback=H #Starting complex ridge height = H from input parameters
    LRidgeFlank10=H/m1 #Starting front ridge flank width
    LRidgeFlank20=H/m2 #Starting back ridge flank width
    Shoreline=-LRidgeFlank10-L #Starting shoreline location
    DV=((LRidgeFlank10+LRidgeFlank20)*H)/2 #Starting ridge volume
    RidgeCounter=0 #Keep track of number of ridges produced
    speccond=0 #Starting geometry is simple
    flag=1 #Title bar variable... Starting Type = Transgressive or Aggraded
    
    #Plotting
    animate_count=0
    
    #Variable Storage
    DVst=NaNlist(Tmax*dt) #Store Ridge Volume
    Hst=NaNlist(Tmax*dt) #Store Ridge Height
    Dst=NaNlist(Tmax*dt) #Store depth of shoreface
    Lst=zerolist(Tmax*dt) #Store Beach Length
    xbst=NaNlist(Tmax*dt) #Store backbarrier x
    zbst=NaNlist(Tmax*dt) #Store backbarrier z
    xc1st=NaNlist(Tmax*dt) #Store active ridge x
    zc1st=NaNlist(Tmax*dt) #Store active ridge z
    xsst=NaNlist(Tmax*dt) #Store shoreline x
    zsst=NaNlist(Tmax*dt) #Store shoreline z
    ShorelineChange=np.zeros((Tmax*dt)) #Store change in shoreline
    DepthVector=NaNlist(Tmax*dt) #Volume units stored to beach
    RidgeVolumeLoss=NaNlist(Tmax*dt) #Store volume lost from ridge to sea level rise
    Amalgamation = np.zeros((Tmax*dt)) #Store change in shoreline

    #Profile Vector
    RidgeProfX=zerolist(Tmax*dt) #Will store geomorphic cross-section X values
    RidgeProfY=zerolist(Tmax*dt) #Will store geomorphic cross-section Y values
    RidgeTypeCounter=NaNlist(Tmax*dt); #Store vector of ridge types 1=simple 2=complex

    #Starting Geometry
    InlandFlankToe=H/m2 #Set starting xb
    InlandFlankCrest=0 #Set starting xc
    SeawardFlankCrest=0 #Set starting xc
    SeawardFlankToe=-H/m1 #Set starting xf
    BaseLevel=0 #Set the starting sea level for new ridge
    D=b #Set starting shoreface depth
    
    #Loop Conditions
    j=1 #Set starting number of ridges
    suppress=0 #Switch to suppress ridge growth required if certain conditions 
    gg=1 #toggle to remember complex ridge geometry during erosion: 1==do not hold P2, 2==update P2

    #-------------------------------------------------------#
    #########---------------MAIN LOOP---------------#########
    #-------------------------------------------------------#

    for i in range(0,Tmax,dt): #solve i for 0 to Tmax by dt
        if H > maxH: #ridge built through overwash deposition unless H reaches a certain height (Hmax)
            swash=1 #swash deposition
        else:
            swash = 0 #overwash deposition
            
        if L>5: 
            Qr=Qri 
        else: 
            Qr=0 #Stop Qr if no beach present
            
        P1=1+(RidgeCounter*5) #Counter to track position of inland ridge flank
        P2=2+(RidgeCounter*5) #Counter to track position of inland ridge crest
        P3=3+(RidgeCounter*5) #Counter to track position of shoreward ridge crest
        P4=4+(RidgeCounter*5) #Counter to track position of shoreward ridge flank
        P5=5+(RidgeCounter*5) #Counter to track position of swale shoreward edge

        #-------------------------------------------------------#
        #########-------SUBAERIAL GROWTH ROUTINES-------#########
        #-------------------------------------------------------#

        if RidgeCounter==0: 
            speccond=0 #enforce simple growth if monolothic ridge
        suppress=0 #do not suppress ridge growth; no special conditions met
        if speccond==1: #start complex ridge geometry routine
            RidgeTypeCounter[RidgeCounter+1]=2 #Set the ridge type for complex
            BaseLevel=BaseLevel+a*dt #Update base sea level
            Laccom=L*a*dt #Calculate accomodation volume of beach  
            ADV=(0.5*(H*H)*(m1+m2))/(m1*m2) #caclulate ridge volume
            ADV=ADV+Qr*dt #caclulate ridge volume + Qr input
            H=((2*ADV)/((1/m1)+(1/m2)))**0.5-a*dt*2 #update ridge height
            LRidgeFlank1=H/m1 #Find the flank width of the ridge on the beach side of the crest
            LRidgeFlank2=H/m2 #Find the flank width of the ridge on the inland side of the crest
            ADVlossAggro=((LRidgeFlank10+LRidgeFlank20)*a*dt) #How much ridge volume is lost to sea-level rise?
            ADV=ADV+ADVlossAggro #update ridge volume based on loss to sea-level rise.
            if ADV<0: 
                ADV=0 #ridge volume can't be negative.

            H=((2*ADV)/((1/m1)+(1/m2)))**0.5-a*dt*2 #final ridge height
            LRidgeFlank1=H/m1 #Find the flank width of the ridge on the beach side of the crest
            LRidgeFlank2=H/m2 #Find the flank width of the ridge on the inland side of the crest
            LRidgeFlankdot1=LRidgeFlank1-LRidgeFlank10 #calculate change in flank width
            LRidgeFlankdot2=LRidgeFlank2-LRidgeFlank20 #calculate change in flank width

            D,ShorelineChange[i] = find_D(D,a,dt,m,ShorelineChange,i,Qs,Qr) #find new shoreface depth and shoreline without knowing dx
            L=L+(Qs/D*dt)-(Qr/D*dt)-LRidgeFlankdot1 #What is the new length of the beach?
            DepthVector[i]=-ShorelineChange[i]*D #Volume units stored to beach (this will go into accomodation in next step!)
            L=L-(Laccom/D) #What is the new length of the beach accounting for accomodation?
            Shoreline=Shoreline+ShorelineChange[i]+(Laccom/D) #Where is the new shoreline?             
            D = a*dt*i + (b+abs(Shoreline)*m) #Update D based on new shoreline
            LRidgeFlank10=LRidgeFlank1 #Set starting ridge flank lengths for next iteration
            LRidgeFlank20=LRidgeFlank2 #Set starting ridge flank lengths for next iteration
            
            #Update cross-section geometries...
            Hback=H
            if swash == 1: 
                SeawardFlankCrest=RidgeProfX[P3]-LRidgeFlankdot2
                InlandFlankCrest=RidgeProfX[P3]-LRidgeFlankdot2
                InlandFlankToe=RidgeProfX[P3]+H/m2-LRidgeFlankdot2
                SeawardFlankToe=RidgeProfX[P3]-LRidgeFlank1-LRidgeFlankdot2
                Shoreline = Shoreline-LRidgeFlankdot2
            else:
                SeawardFlankCrest=RidgeProfX[P3]
                InlandFlankCrest=RidgeProfX[P3]
                InlandFlankToe=RidgeProfX[P3]+H/m2 
                SeawardFlankToe=RidgeProfX[P3]-LRidgeFlank1
                
            ###Subbroutine to calculate geometry of complex multi-ridge system###
            #This routine describes what the model should do when one ridge grows into another!

            #Initial Geometry
            RidgeProfX[P1]=InlandFlankToe
            RidgeProfX[P2]=InlandFlankCrest
            RidgeProfX[P3]=SeawardFlankCrest
            RidgeProfX[P4]=SeawardFlankToe
            RidgeProfX[P5]=Shoreline
            RidgeProfY[P1]=BaseLevel
            RidgeProfY[P2]=a*dt*i+Hback
            RidgeProfY[P3]=a*dt*i+H
            RidgeProfY[P4]=a*dt*i
            RidgeProfY[P5]=a*dt*i

            #Overlapped Area Between Active and Secondary ridge
            ActiveRidge=Polygon([(RidgeProfX[P1],RidgeProfY[P1]), (RidgeProfX[P2],RidgeProfY[P2]), (RidgeProfX[P3],RidgeProfY[P3]), (RidgeProfX[P4],RidgeProfY[P4])])
            SecondaryRidge=Polygon([(RidgeProfX[P1-5],RidgeProfY[P1-5]), (RidgeProfX[P2-5],RidgeProfY[P2-5]), (RidgeProfX[P3-5],RidgeProfY[P3-5]), (RidgeProfX[P4-5],RidgeProfY[P4-5])])
            polyout=SecondaryRidge.intersection(ActiveRidge)
            OverArea=(polyout.area)
            OverAreaOld=OverArea #Store for next iteration of swale filling  

            #Intersection point of intital slope with slope of secondary ridge
            ##Backslope of Active ridge
            a1=[float(RidgeProfX[P1]),float(RidgeProfY[P1])]
            a2=[float(RidgeProfX[P2]),float(RidgeProfY[P2])]
            b1=[RidgeProfX[P4-5],RidgeProfY[P4-5]]
            b2=[RidgeProfX[P3-5],RidgeProfY[P3-5]]
            x_intersect, y_intersect = xy_intersect(a1, a2, b1, b2) 

            #Calculate length of active ridge backslope
            backlength = math.hypot(x_intersect-RidgeProfX[P2], y_intersect-RidgeProfY[P2])

            #Thickness of swale file
            BackDepositH=-OverArea/backlength 

            #Extend the crest by the thickness
            xshift=BackDepositH*m1
            yshift=xshift*m1
            H=H+yshift
            Hback=H    

            #Update geometries...
            RidgeProfX[P1]=RidgeProfX[P4-5]
            RidgeProfX[P3]=RidgeProfX[P3]+xshift
            RidgeProfX[P2]=x_intersect+xshift
            RidgeProfX[P4]=SeawardFlankToe
            RidgeProfX[P5]=Shoreline
            RidgeProfY[P1]=RidgeProfY[P4-5]
            RidgeProfY[P2]=y_intersect+yshift
            RidgeProfY[P3]=RidgeProfY[P3]+yshift
            RidgeProfY[P4]=a*dt*i
            RidgeProfY[P5]=a*dt*i
            speccond=1; #On next loop, trigger condition for complex Ridge growth.
            #end of swale infilling triggered during simple, monoltihic 
            #foreridge development.

            ###Begin Amalgamation and Complex Regrowth Conditions

            if y_intersect+yshift>=RidgeProfY[P3-5]: #Amalgamation: Shut down conditions
            #for further multi-ridge growth if the swale is full
                Amalgamation[i] = 1

                #Inital Geometry Update
                RidgeProfX[P1]=RidgeProfX[P1-5]
                RidgeProfX[P2]=RidgeProfX[P3]
                RidgeProfX[P4]=SeawardFlankToe
                RidgeProfX[P5]=Shoreline
                RidgeProfY[P1]=BaseLevel
                RidgeProfY[P2]=BaseLevel
                RidgeProfY[P4]=a*dt*i
                RidgeProfY[P5]=a*dt*i
                #Merge ridge Pt 1
                RidgeProfX[P1-5]=RidgeProfX[P1]
                RidgeProfX[P2-5]=RidgeProfX[P2]
                RidgeProfX[P3-5]=RidgeProfX[P3]
                RidgeProfX[P4-5]=RidgeProfX[P4]
                RidgeProfX[P5-5]=RidgeProfX[P5]
                RidgeProfY[P1-5]=RidgeProfY[P1]
                RidgeProfY[P2-5]=RidgeProfY[P2]
                RidgeProfY[P3-5]=RidgeProfY[P3]
                RidgeProfY[P4-5]=RidgeProfY[P4]
                RidgeProfY[P5-5]=RidgeProfY[P5]
                #Merge Ridge Pt 2
                RidgeProfX[P1]=0
                RidgeProfX[P2]=0
                RidgeProfX[P3]=0
                RidgeProfX[P4]=0
                RidgeProfX[P5]=0
                RidgeProfX[P5+1]=0
                RidgeProfX[P5+2]=0
                RidgeProfX[P5+3]=0    
                RidgeProfY[P1]=0
                RidgeProfY[P2]=0
                RidgeProfY[P3]=0
                RidgeProfY[P4]=0
                RidgeProfY[P5]=0
                RidgeProfY[P5+1]=0
                RidgeProfY[P5+2]=0
                RidgeProfY[P5+3]=0
                RidgeTypeCounter[RidgeCounter+1]=np.nan #Clear the ridge type 
                RidgeCounter=RidgeCounter-1
                RidgeTypeCounter[RidgeCounter+1]=1 #Set the ridge type for simple
                speccond=0
                DV=(0.5*(H*H)*(m1+m2))/(m1*m2)
                H=((2*DV)/((1/m1)+(1/m2)))**0.5
                LRidgeFlank10=H/m1 #Find the flank width of the ridge on the beach side of the crest
                LRidgeFlank20=H/m2 #Find the flank width of the ridge on the inland side of the crest
                RidgeProfX[P1-5]=RidgeProfX[P3-5]+LRidgeFlank20
                P1=1+(RidgeCounter*5) #Counter to track position of inland ridge flank
                P2=2+(RidgeCounter*5) #Counter to track position of inland ridge crest
                P3=3+(RidgeCounter*5) #Counter to track position of shoreward ridge crest
                P4=4+(RidgeCounter*5) #Counter to track position of shoreward ridge flank
                P5=5+(RidgeCounter*5) #Counter to track position of swale shoreward edge
                # end amalgamation

            if RidgeProfX[P4]>Shoreline:
                gg=1 #P2 hold ended; regrowth of complex ridge

            ###End Amalgamation and Complex Regrowth Conditions

        if speccond==0: #start simple, monolithic foreridge growth routine
            flag=1 #Morphological Type = Transgressive or Aggraded
            RidgeTypeCounter[RidgeCounter+1]=1 #Set the ridge type for simple
            Laccom=L*a*dt #Calculate accomodation volume of beach
            DVlossAggro=((LRidgeFlank10+LRidgeFlank20)*a*dt) #How much ridge volume is lost to sea level rise?
            RidgeVolumeLoss[i]=DVlossAggro #Store volume lost from ridge to sea level rise
            DV=-DVlossAggro+DV+Qr*dt #Ridge volume equals previous ridge volume plus sediment input from beach and loss to sea level
            if DV<0:
                DV=0 #Ridge volume can't be negative
            H=((2*DV)/((1/m1)+(1/m2)))**0.5 #Find the height of the new ridge
            LRidgeFlank1=H/m1 #Find the flank width of the ridge on the beach side of the crest
            LRidgeFlank2=H/m2 #Find the flank width of the ridge on the inland side of the crest
            LRidgeFlankdot1=LRidgeFlank1-LRidgeFlank10 #calculate change in flank width
            LRidgeFlankdot2=LRidgeFlank2-LRidgeFlank20 #calculate change in flank width
            D,ShorelineChange[i] = find_D(D,a,dt,m,ShorelineChange,i,Qs,Qr) #find new shoreface depth and shoreline without knowing dx
            L=L+(Qs/D*dt)-(Qr/D*dt)-LRidgeFlankdot1 #What is the new length of the beach?
            ShorelineChange[i]=-(Qs/D*dt)+(Qr/D*dt) #Store shoreline change values per dt ignoring acccomodation
            DepthVector[i]=-ShorelineChange[i]*D #Volume units stored to beach (this will go into accomodation in next step!)
            L=L-(Laccom/D) #What is the new length of the beach accounting for accomodation?
            Shoreline=Shoreline+ShorelineChange[i]+(Laccom/D) #Where is the new shoreline?\\
            D = a*dt*i + (b+abs(Shoreline)*m) #Update D based on new shoreline
            LRidgeFlank10=LRidgeFlank1 #Set starting ridge flank lengths for next iteration
            LRidgeFlank20=LRidgeFlank2 #Set starting ridge flank lengths for next iteration

            #Update cross-section geometries...
            Hback=H

            if swash == 1: 
                InlandFlankToe=SeawardFlankCrest+(H/m2)-LRidgeFlankdot2
                SeawardFlankToe=SeawardFlankToe-LRidgeFlankdot1-LRidgeFlankdot2
                SeawardFlankCrest=InlandFlankCrest=InlandFlankCrest-LRidgeFlankdot2
                Shoreline = Shoreline-LRidgeFlankdot2
            else:
                SeawardFlankCrest=InlandFlankCrest
                InlandFlankToe=SeawardFlankCrest+(H/m2)
                SeawardFlankToe=SeawardFlankToe-LRidgeFlankdot1
            
            RidgeProfX[P1]=InlandFlankToe
            RidgeProfX[P2]=InlandFlankCrest
            RidgeProfX[P3]=SeawardFlankCrest
            RidgeProfX[P4]=SeawardFlankToe
            RidgeProfX[P5]=Shoreline
            RidgeProfY[P1]=BaseLevel
            RidgeProfY[P2]=a*dt*i+Hback
            RidgeProfY[P3]=a*dt*i+H
            RidgeProfY[P4]=a*dt*i
            RidgeProfY[P5]=a*dt*i

            if RidgeCounter>0: # if we are building onto an existing ridge
                if InlandFlankToe>RidgeProfX[P4-5]: #Trigger swale infilling
                    RidgeTypeCounter[RidgeCounter+1]=2 #Set the ridge type for complex
                    flag=2 #Set the ridge type for complex (For plotting)

                    ###Subbroutine to calculate geometry of complex multi-ridge system###
                    #This routine describes what the model should do when one ridge growsinto another!
                    # the below is a repeat of what we did above to account for amalgamation when speccond = 1

                    #Initial Geometry
                    RidgeProfX[P1]=InlandFlankToe
                    RidgeProfX[P2]=InlandFlankCrest
                    RidgeProfX[P3]=SeawardFlankCrest
                    RidgeProfX[P4]=SeawardFlankToe
                    RidgeProfX[P5]=Shoreline
                    RidgeProfY[P1]=BaseLevel
                    RidgeProfY[P2]=a*dt*i+Hback
                    RidgeProfY[P3]=a*dt*i+H
                    RidgeProfY[P4]=a*dt*i
                    RidgeProfY[P5]=a*dt*i

                    #Overlapped Area Between Active and Secondary Ridge
                    ActiveRidge=Polygon([(RidgeProfX[P1],RidgeProfY[P1]), (RidgeProfX[P2],RidgeProfY[P2]), (RidgeProfX[P3],RidgeProfY[P3]), (RidgeProfX[P4],RidgeProfY[P4])])
                    SecondaryRidge=Polygon([(RidgeProfX[P1-5],RidgeProfY[P1-5]), (RidgeProfX[P2-5],RidgeProfY[P2-5]), (RidgeProfX[P3-5],RidgeProfY[P3-5]), (RidgeProfX[P4-5],RidgeProfY[P4-5])])
                    polyout=SecondaryRidge.intersection(ActiveRidge)
                    OverArea=(polyout.area)
                    OverAreaOld=OverArea #Store for next iteration of swale filling  

                    #Intersection point of intital slope with slope of secondary ridge
                    ##Backslope of Active Ridge
                    a1=[RidgeProfX[P1],RidgeProfY[P1]]
                    a2=[RidgeProfX[P2],RidgeProfY[P2]]
                    b1=[RidgeProfX[P4-5],RidgeProfY[P4-5]]
                    b2=[RidgeProfX[P3-5],RidgeProfY[P3-5]]
                    x_intersect, y_intersect = xy_intersect(a1, a2, b1, b2) 

                    #Calculate length of active dune backslope
                    backlength = math.hypot(x_intersect-RidgeProfX[P2], y_intersect-RidgeProfY[P2])

                    #Thickness of swale file
                    BackDepositH=-OverArea/(RidgeProfX[P4-5]-P1)

                    #Extend the crest by the thickness
                    xshift=BackDepositH*m1
                    yshift=xshift*m1
                    H=H+yshift
                    Hback=H    

                    #Update geometries...
                    RidgeProfX[P1]=RidgeProfX[P4-5]
                    RidgeProfX[P3]=RidgeProfX[P3]+xshift
                    RidgeProfX[P2]=x_intersect+xshift
                    RidgeProfX[P4]=SeawardFlankToe
                    RidgeProfX[P5]=Shoreline
                    RidgeProfY[P1]=BaseLevel
                    RidgeProfY[P2]=y_intersect+yshift
                    RidgeProfY[P3]=RidgeProfY[P3]+yshift
                    RidgeProfY[P4]=a*dt*i
                    RidgeProfY[P5]=a*dt*i
                    speccond=1; #On next loop, trigger condition for complex Ridge growth.
                    #end of swale infilling triggered during simple, monolithic foreridge development.


    ##END OF SUBAERIAL GROWTH ROUTINES################################################################################

        #-------------------------------------------------------#
        #########-----TRANGRESSION/EROSION ROUTINES-----#########
        #-------------------------------------------------------#

        #In this section, the subaerial portion of the barrier will undergo erosion.

        #-------------------------------------------------------#
        #COMMAND LEVEL 1: IF SHORELINE ERODES INTO RIDGE TOE######
        if RidgeProfX[P4]<Shoreline: 
            if gg==1:
                DPX2=RidgeProfX[P2]
                DPY2=RidgeProfY[P2]
                gg=2

        #This is triggering condition for the erosion routine. All further subroutines for erosion/trangression fall within
        #this condition, unless activtiating event-scale processes.

        #There are 2 possible scenarios based on initial geometry...
        #1) The subaerial geometry is simple. (Single Ridge)
        #2) The subaerial geometry is complex. (Multi-Ridge)

        #-------------------------------------------------------#
        ###COMMAND LEVEL 2: IF SUBAERIAL GEOMETRY IS SIMPLE#####
            if RidgeTypeCounter[RidgeCounter+1]==1:
                shiftx=Shoreline-RidgeProfX[P4] #How far has the shoreline undercut the ridge?
                flag=3 #Type = Simple Subaerial Erosion

                #Intersection of erosional scarp with backslope
                a1=[RidgeProfX[P1],RidgeProfY[P1]]
                a2=[RidgeProfX[P3],RidgeProfY[P3]]
                b1=[Shoreline,(a*dt*i)]
                b2=[RidgeProfX[P3]+shiftx,RidgeProfY[P3]]
                x_intersect, y_intersect = xy_intersect(a1, a2, b1, b2)   

                ErodedShape=Polygon([(RidgeProfX[P4],RidgeProfY[P4]), (RidgeProfX[P3],RidgeProfY[P3]), (x_intersect,y_intersect), (Shoreline,(a*dt*i))])

                #Update Geometries...
                InlandFlankCrest=x_intersect
                SeawardFlankCrest=x_intersect
                H=y_intersect-(a*dt*i)
                DV=(0.5*(H*H)*(m1+m2))/(m1*m2)
                SeawardFlankToe=x_intersect-H/m1
                RidgeProfX[P2]=x_intersect
                RidgeProfY[P2]=y_intersect
                RidgeProfX[P3]=x_intersect
                RidgeProfY[P3]=y_intersect
                RidgeProfX[P4]=Shoreline
                RidgeProfY[P4]=a*dt*i
                QsReturn=(ErodedShape.area)
                Shoreline=Shoreline-(QsReturn/D)
                L=RidgeProfX[P4]-Shoreline 
                RidgeProfX[P5]=Shoreline

        #-------------------------------------------------------#  
        ####COMMAND LEVEL 2: IF SUBAERIAL GEOMETRY IS COMPLEX#####
            if RidgeTypeCounter[RidgeCounter+1]==2:
                shiftx=Shoreline-RidgeProfX[P4]
                flag=4 #Complex Subaerial Erosion                 

                #Intersection of erosional scarp with backslope
                a1=[DPX2,DPY2]
                a2=[RidgeProfX[P3],RidgeProfY[P3]]
                b1=[Shoreline,(a*dt*i)]
                b2=[RidgeProfX[P3]+shiftx,RidgeProfY[P3]]
                x_intersect, y_intersect = xy_intersect(a1, a2, b1, b2)   

                ErodedShape=Polygon([(RidgeProfX[P4],RidgeProfY[P4]), (RidgeProfX[P3],RidgeProfY[P3]), (x_intersect,y_intersect), (Shoreline,(a*dt*i))])

                #Update Geometries...
                InlandFlankCrest=x_intersect
                SeawardFlankCrest=x_intersect
                #H=y_intersect-(a*dt*i)
                H=((RidgeProfX[P4-5]-RidgeProfX[P4])*m2+(DPY2-a*dt*i))
                DV=(0.5*(H*H)*(m1+m2))/(m1*m2)
                RidgeProfX[P2]=DPX2
                RidgeProfY[P2]=DPY2
                RidgeProfX[P3]=x_intersect
                RidgeProfY[P3]=H+(a*dt*i)#y_intersect
                RidgeProfX[P4]=Shoreline
                RidgeProfY[P4]=a*dt*i
                QsReturn=(ErodedShape.area)
                Shoreline=Shoreline-(QsReturn/D)
                L=RidgeProfX[P4]-Shoreline 
                RidgeProfX[P5]=Shoreline          

                #Below routine is activated if active ridge is destroyed and relict ridge is reactivated.
                if Shoreline>RidgeProfX[P4-5]:
                    ##Purge geometry related to destroyed, active ridge...
                    gg=1
                    RidgeProfX[P1]=0
                    RidgeProfX[P2]=0
                    RidgeProfX[P3]=0
                    RidgeProfX[P4]=0
                    RidgeProfX[P5]=0
                    RidgeProfX[P5+1]=0
                    RidgeProfX[P5+2]=0
                    RidgeProfX[P5+3]=0    
                    RidgeProfY[P1]=0
                    RidgeProfY[P2]=0
                    RidgeProfY[P3]=0
                    RidgeProfY[P4]=0
                    RidgeProfY[P5]=0
                    RidgeProfY[P5+1]=0
                    RidgeProfY[P5+2]=0
                    RidgeProfY[P5+3]=0
                    RidgeTypeCounter[RidgeCounter+1]=np.nan #Clear the ridge type 
                    RidgeCounter=RidgeCounter-1
                    P1=1+(RidgeCounter*5) #Counter to track position of inland ridge flank
                    P2=2+(RidgeCounter*5) #Counter to track position of inland ridge crest
                    P3=3+(RidgeCounter*5) #Counter to track position of shoreward ridge crest
                    P4=4+(RidgeCounter*5) #Counter to track position of shoreward ridge flank
                    P5=5+(RidgeCounter*5) #Counter to track position of swale shoreward edge
                    ##...geometry of destroyed ridge pruged from storage
                    ##Now, correcting Relict Ridge Geometry due to change in base sea level
                    if RidgeCounter>0:
                        shifty=a*dt*i-RidgeProfY[P1]
                        RidgeProfY[P1]=a*dt*i
                        RidgeProfY[P4]=a*dt*i
                        RidgeProfY[P5]=a*dt*i
                        RidgeProfX[P4]=RidgeProfX[P4]+(shifty/m1)
                        H=RidgeProfY[P3]-(a*dt*i)
                        if RidgeProfY[P3]<a*i*dt: 
                            H=0
                            DV=(0.5*(H*H)*(m1+m2))/(m1*m2)
                            LRidgeFlank10=H/m1 #Find the flank width of the ridge on the beach side of the crest
                            LRidgeFlank20=H/m2 #Find the flank width of the ridge on the inland side of the crest
                            L=RidgeProfX[P4]-Shoreline
                            if RidgeCounter>1: 
                                RidgeProfY[P1]=RidgeProfY[P5-5]
                                RidgeProfX[P1]=RidgeProfX[P4-5]
                                ##Geometry of relict ridge corrected to account for change in sea level since last activation
                    #---#

                    #Below routine is activated if complex, eroding ridge is converted to simple, monolthic ridge. 
                    #Need to establish the new base level of the reactivated ridge, which means the x pos of the 
                    #ridge toe and ridge heel must be found.

                    if RidgeCounter==0: #only one ridge remains
                        H1=RidgeProfY[3]-(a*dt*i) #store old height temporarily
                        if RidgeProfY[3]<a*i*dt: 
                            H1=0 #ridge height is zero if below

                        #Intersection of reactived backridge slope with new base level
                        a1=[RidgeProfX[P1],RidgeProfY[P1]]
                        a2=[RidgeProfX[P3],RidgeProfY[P3]]
                        b1=[2000,(a*dt*i)]
                        b2=[-2000,(a*dt*i)]
                        x_intersectbb, y_intersectbb = xy_intersect(a1, a2, b1, b2) 
                        RidgeProfX[1]=x_intersectbb

                        #Intersection of reactived backridge slope with new base level
                        a1=[RidgeProfX[P3],RidgeProfY[P3]]
                        a2=[RidgeProfX[P4],RidgeProfY[P4]]
                        b1=[2000,(a*dt*i)]
                        b2=[-2000,(a*dt*i)]
                        x_intersectff, y_intersectff = xy_intersect(a1, a2, b1, b2) 
                        RidgeProfX[4]=x_intersectff

                        #Update geometries
                        RidgeProfY[1]=a*dt*i
                        RidgeProfY[4]=a*dt*i
                        RidgeProfY[5]=a*dt*i
                        RidgeShape=Polygon([(RidgeProfX[P4],RidgeProfY[P4]), (RidgeProfX[P3],RidgeProfY[P3]), (RidgeProfX[P1],RidgeProfY[P1])])
                        DV=(RidgeShape.area)
                        H=((2*DV)/((1/m1)+(1/m2)))**0.5 #Find the height of the new ridge
                        Hx=H-H1 #height check for diagnostic purposes
                        LRidgeFlank10=RidgeProfX[3]-RidgeProfX[4] #Find the flank width of the ridge on the beach side of the crest
                        LRidgeFlank20=RidgeProfX[1]-RidgeProfX[3] #Find the flank width of the ridge on the inland side of the crest
                        L=RidgeProfX[4]-Shoreline
                        InlandFlankToe=RidgeProfX[1]
                        InlandFlankCrest=RidgeProfX[3]
                        SeawardFlankCrest=RidgeProfX[3]
                        SeawardFlankToe=RidgeProfX[4]

                    suppress=1 #suppress ridge growth for geometry update during this iteration   
                    #end triggering condition for erosion routine


    ####Shoreline Geometry Update [resets shoreline position based on geometry]      
        if RidgeCounter>0:
            if Shoreline>RidgeProfX[P5-5]: 
                RidgeProfX[P5-5]=Shoreline 
                
    # if ridges can nucleate
        if suppress==0: #no special conditons limiting ridge growth
            # Initiate ridge nucleation when distance between shoreline and ridge is exceeded or ridge grows too high
            if (SeawardFlankToe - Shoreline > maxdis) or (H > maxH+2): 
                j=j+1
                speccond=0 #Terminate complex growth mode if activated
                RidgeCounter=RidgeCounter+1 #New ridge counted
                H=a*dt #Reset for new ridge
                DV=0 #Reset for new ridge
                LRidgeFlank20=0 #Reset flank width for new ridge
                LRidgeFlank10=0 #Reset flank width for new ridge
                #Update geometry of new ridge
                InlandFlankToe=Shoreline+Ls
                RidgeProfX[P5]=InlandFlankToe
                InlandFlankCrest=Shoreline+Ls
                SeawardFlankCrest=Shoreline+Ls
                SeawardFlankToe=Shoreline+Ls
                L=Ls; #Reset beach length
                BaseLevel=a*i*dt #Set the starting sea level for new ridge
                #end new ridge nucleation
                #end 'no special conditons limiting ridge growth'

        ###Update tracking variables
        Lst[i]=L #Store L 
        Hst[i]=H #Store H 
        xsst[i]=-Shoreline #Store Shoreline 
        Dst[i]=D

    #-------------------------------------------------------#
    #########-----------PLOTTING ROUTINES-----------#########
    #-------------------------------------------------------#
    
        animate_count = animate_count+1
        if animate_count==atrigger:  
            bh = -40
            ybh= 5000
            
            SeaX=[-ybh,-ybh,ybh,ybh]
            SeaX=[i * (-1) for i in SeaX]
            SeaXshore=[-ybh,-ybh,RidgeProfX[P5],RidgeProfX[P5]]
            SeaXshore=[(i * (-1)) for i in SeaXshore]
            SeaY=[bh,a*dt*i,a*dt*i,bh]

            RidgeProfX[P5+1]=RidgeProfX[P5]
            RidgeProfY[P5+1]=a*dt*(i)-D
            RidgeProfX[P5+2]=10000 
            RidgeProfY[P5+2]=a*dt*(i)-D
            RidgeProfX[P5+3]=10000 
            RidgeProfY[P5+3]=RidgeProfY[P5]
            
            RidgeXSurf=RidgeProfX[1:P5+3]
            RidgeXSurf=[(i * (-1)) for i in RidgeXSurf]
            RidgeYSurf=RidgeProfY[1:P5+3]
            
            # get rid of nan and inf values for plotting
            RidgeXSurf = np.array(RidgeXSurf)
            RidgeYSurf = np.array(RidgeYSurf)
            mask = (
                ~np.isnan(RidgeXSurf) &
                ~np.isnan(RidgeYSurf) &
                ~np.isinf(RidgeXSurf) &
                ~np.isinf(RidgeYSurf)
            )
            RidgeXSurf = RidgeXSurf[mask]
            RidgeYSurf = RidgeYSurf[mask]

            figure(figsize=(9, 6), dpi=180)  
            
            # plot depth of bedrock           
            plt.fill(SeaXshore,SeaY,'cadetblue') # sea level 
            plt.fill(RidgeXSurf,RidgeYSurf,"palegoldenrod")  # beach ridge fill
            plt.plot(RidgeXSurf[0:-1],RidgeYSurf[0:-1],color='black') # beach ridge black line

            xn = np.arange(0,5000,1)
            yn= -1*(b+((m)*xn))
            xn[-1] = xn[0]
            yn[-2] = yn[-1] = bh  
            plt.plot(xn,yn,color = 'black')# outline of the sea floor --> bedrock)
            plt.fill(xn,yn,"tan") #fill of the sea floor --> bedrock
 
            plt.xlim(0,2000)
            plt.ylim(-20,5) 
        
            QsFig='{0:.2f}'.format(Qs)
            QrFig='{0:.2f}'.format(Qr)
            
            # Add text for parameter type being tested
            if type_of_image == 'slope':
                plt.text(
                    0.95, 0.95, "Slope = " + str(m),
                    transform=plt.gca().transAxes,   # relative to axes
                    fontsize=18,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.5'))

            if type_of_image == 'RSL':
                plt.text(
                    0.95, 0.95, "RSL = " + str(a) + " m/yr",
                    transform=plt.gca().transAxes,   # relative to axes
                    fontsize=18,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.5'))

            if type_of_image == 'sed':
                plt.text(
                    0.95, 0.95, '$Q_{S}$=' +QsFig +' $m^3/m/yr$\n$Q_{R}$=' +QrFig +' $m^3/m/yr$',
                    transform=plt.gca().transAxes,   # relative to axes
                    fontsize=18,
                    verticalalignment='top',
                    horizontalalignment='right',
                    bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.5'))
            
            plt.ylabel('Elevation (m)',fontsize=20)
            plt.xlabel('Distance (m)',fontsize=20)
            plt.xticks(fontsize=17) 
            plt.yticks(fontsize=17) 
            if saveimage == True:
                plt.savefig(title, dpi=300) 
            plt.show()
            plt.close()
            animate_count=0 
            
