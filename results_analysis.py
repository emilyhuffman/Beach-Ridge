import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.cm as cm
import matplotlib.colors as colors
from tabulate import tabulate
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ridge_parameter import swash_wash_params
from matplotlib.colors import LinearSegmentedColormap

class ModelAnalysis:
    def __init__(self,Tmax,dt,a,b,maxH,Qr,Qs,m):
        self.Tmax = Tmax
        self.a = a
        self.b = b
        self.maxH = maxH      
        self.dt = dt
        self.Qr = Qr
        self.Qs = Qs
        self.m = m
        self.height = []
        self.length = []
        self.time = []
        self.ridgenumber = []
        self.depth = []

    ## obtain parameters from differing on variable 
    
    def time_for_ridge(self, ss_tt):
        sst = [val for val in ss_tt if val > 0]
        if len(sst) < 2:
            return self.Tmax
        diffs = np.diff(sst)
        avediff = np.mean(diffs)
        return avediff if not np.isnan(avediff) else self.Tmax

    def sedsup_params(self,sedsupQs,sedsupQr):
        sedH = []
        sedL = []
        sedT = []
        sedrn = []
        depth=[]
        for i in range(0,len(sedsupQs)):
            sedL_f, sedH_f, sed_RidgeCounter,sed_tt,d  = swash_wash_params(sedsupQr[i],
                                                                        sedsupQs[i],
                                                                        self.Tmax,
                                                                        self.dt,
                                                                        self.m,
                                                                        self.a,
                                                                        self.b,
                                                                        self.maxH)
            
            sedH.append(max(sedH_f))
            if len(sedL_f) > 0:
                sedL.append(max(sedL_f))
            else: 
                sedL.append(0)
            sedT.append(self.time_for_ridge(sed_tt))
            sedrn.append(sed_RidgeCounter)
            depth.append([d])
            
        self.height = sedH
        self.length = sedL
        self.time = sedT
        self.ridgenumber = sedrn
        self.depth = depth

    def slope_params(self,m_vec):
        mH = []
        mL = []
        mT = []
        mrn = []
        depth=[]
        for i in range(0,len(m_vec)):
            mL_f, mH_f, m_RidgeCounter,m_tt,d=swash_wash_params(self.Qr,self.Qs,self.Tmax,self.dt,m_vec[i],self.a,self.b,self.maxH)
            
            mH.append(max(mH_f))
            if len(mL_f) > 0:
                mL.append(max(mL_f))
            else: 
                mL.append(0)
            mT.append(self.time_for_ridge(m_tt))
            mrn.append(m_RidgeCounter)
            depth.append([d])
            
        self.height = mH
        self.length = mL
        self.time = mT
        self.ridgenumber = mrn
        self.depth = depth

    
    def RSL_params(self,RSL_vec):
        rslH = []
        rslL = []
        rslT = []
        rslrn = []
        depth=[]
        for i in range(0,len(RSL_vec)):
            rslL_f, rslH_f, rsl_RidgeCounter,rsl_tt,d = swash_wash_params(self.Qr,
                                                                        self.Qs,
                                                                        self.Tmax,
                                                                        self.dt,
                                                                        self.m,
                                                                        RSL_vec[i],
                                                                        self.b,
                                                                        self.maxH)

            rslH.append(max(rslH_f))
            if len(rslL_f) > 0:
                rslL.append(max(rslL_f))
            else: 
                rslL.append(0)
            rslT.append(self.time_for_ridge(rsl_tt))
            rslrn.append(rsl_RidgeCounter)
            depth.append([d])

        self.height = rslH
        self.length = rslL
        self.time = rslT
        self.ridgenumber = rslrn
        self.depth = depth


## tables to display ###
    
    def RSL_table(self,RSL_vec):
        self.RSL_params(RSL_vec)
        table = {'RSL rate (m/yr)': RSL_vec, 
                 'Ridge Number': self.ridgenumber, 
                 'Max Ridge Height (m)': self.height,
                 'Ridge Spacing (m)': self.length,
                'Time per ridge (yrs)': self.time}
        print(tabulate(table,headers='keys',tablefmt='fancy_grid'))

    def sedsup_table(self,sedsupQs,sedsupQr):
        self.sedsup_params(sedsupQs,sedsupQr)
        table = {'Qs': sedsupQs,'Qr': sedsupQr, 
                 'Ridge Number': self.ridgenumber, 
                 'Max Ridge Height (m)': self.height,
                 'Ridge Spacing (m)': self.length,
                'Time per ridge (yrs)': self.time}
        print(tabulate(table,headers='keys',tablefmt='fancy_grid'))

    def slope_table(self,m_vec):
        self.slope_params(m_vec)
        table = {'Bedrock Slope': m_vec, 
                 'Ridge Number': self.ridgenumber, 
                 'Max Ridge Height (m)': self.height,
                 'Ridge Spacing (m)': self.length,
                'Time per ridge (yrs)': self.time}
        print(tabulate(table,headers='keys',tablefmt='fancy_grid'))
        
    # graphs testing parameters against ridge morphology
    
    def param_graphs(self,vec,title):
    
        newridge = []
        for i in range(1,len(self.ridgenumber)):
            if i-1==0:
                newridge.append(vec[i-1])
            if self.ridgenumber[i] != self.ridgenumber[i-1]:
                newridge.append(vec[i])   
            
        fig, ax = plt.subplots(2, 2,figsize=(5,5))
        plt.subplot(2, 2, 1)
        for i in newridge:
            plt.axvline(i, color = 'b', linewidth=0.4)
        plt.plot(vec,self.ridgenumber,'b.')
        plt.xlabel(title)
        plt.ylabel('Ridge number')
        if title == 'RSL fall [m/yr]':
            plt.gca().invert_xaxis()  # Reflect x-axis
            plt.xlim(0,-.01)
        
        plt.subplot(2, 2, 2)
        for i in newridge:
            plt.axvline(i, color = 'b', linewidth=0.4)
        plt.plot(vec, self.height,'r.')
        plt.xlabel(title)
        plt.ylabel('Ridge Height [m]')
        if title == 'RSL fall [m/yr]':
            plt.gca().invert_xaxis()  # Reflect x-axis
            plt.xlim(0,-.01)            
        
        plt.subplot(2, 2, 3)
        for i in newridge:
            plt.axvline(i, color = 'b', linewidth=0.4)
        plt.plot(vec, self.length,'g.')
        plt.xlabel(title)
        plt.ylabel('Distance between Ridges [m]')
        if title == 'RSL fall [m/yr]':
            plt.gca().invert_xaxis()  # Reflect x-axis
            plt.xlim(0,-.01)

        plt.subplot(2, 2, 4)
        for i in newridge:
            plt.axvline(i, color = 'b', linewidth=0.4)
        plt.plot(vec, self.time,'m.')
        plt.xlabel(title)
        plt.ylabel('Period of Ridge Growth [yr]')
        if title == 'RSL fall [m/yr]':
            plt.gca().invert_xaxis()  # Reflect x-axis
            plt.xlim(0,-.01)
        
        fig.tight_layout()
        fig.savefig(title[0:3] + "param graphs.png", dpi=300) 
        plt.show()
    
    
    def slope_param_graphs(self,m_vec):
        self.slope_params(m_vec)
        self.param_graphs(m_vec,'Slope')
        
    def RSL_param_graphs(self,rsl_vec):
        self.RSL_params(rsl_vec)
        self.param_graphs(rsl_vec,'RSL fall [m/yr]')
        
    def sedsup_param_graphs(self,sedsupQs,sedsupQr):
        self.sedsup_params(sedsupQs,sedsupQr)
        self.param_graphs(sedsupQs,'Sediment Supply [m$^3$/myr]') 
        
            
    def regime_plot_SS(self, N):
        sedsupQs = np.linspace(0, 10, N)
        sedsupQr = np.linspace(0, 10, N)

        Height = np.zeros((N, N))
        Space = np.zeros((N, N))
        Tper_ridge = np.zeros((N, N))
        Ridge_num = np.zeros((N, N))

        # ---- COMPUTE MODEL GRIDS ----
        for i in range(1, N):
            for j in range(1, N):
                if i < j:
                    continue     # skip upper triangle instead of breaking loop

                self.Qs = sedsupQs[i]
                self.Qr = sedsupQr[j]

                L_full, H_full, RidgeCounter, tt, depth = swash_wash_params(
                    self.Qr, self.Qs, self.Tmax, self.dt,
                    self.m, self.a, self.b, self.maxH
                )

                # Fill matrices
                Height[i, j] = np.max(H_full[1:-1]) if len(H_full) > 2 else 0
                Space[i, j]  = np.max(L_full) if L_full.size > 0 else 0
                Tper_ridge[i, j] = self.time_for_ridge(tt)
                Ridge_num[i, j]  = RidgeCounter

        # ---- PLOTTING ----

        fig, axes = plt.subplots(2, 2, figsize=(8, 8))
        fig.tight_layout(pad=3.0)

        x_ticks = np.linspace(0, N-1, 5)
        x_labels = np.linspace(0, 8, 5)

        plots = [
            (Ridge_num, "Ridge Number", " "),
            (Height, "Ridge Height", "[m]"),
            (Space, "Distance Between Ridges", "[m]"),
            (Tper_ridge, "Period of Ridge Growth", "[yr]")
        ]
        
        base = cm.get_cmap("BrBG")

        # sample colors: white → green/blue → brown
        colors = [
            base(0.45),  
            base(0.2),   
            base(0.75),  
            base(1.0)
        ]

        cmap = LinearSegmentedColormap.from_list("BrBG_shifted", colors, N=256)
                
        for ax, (data, title, cbar_label) in zip(axes.flat, plots):
            im = ax.imshow(data, cmap=cmap, origin='lower')
            ax.set_title(title,fontsize=16)
            ax.set_xlabel("$Q_{R}$",fontsize=15)
            ax.set_ylabel("$Q_{S}$",fontsize=15)

            ax.set_xticks(x_ticks)
            ax.set_xticklabels([int(v) for v in x_labels],fontsize=14)

            ax.set_yticks(x_ticks)
            ax.set_yticklabels([int(v) for v in x_labels],fontsize=14)

            # ---- COLORBAR (VERTICAL + LABEL ON BOTTOM) ----
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)

            # keep vertical orientation
            cbar = fig.colorbar(im, cax=cax, orientation='vertical')
            cbar.ax.tick_params(labelsize=12)
            # move label to the bottom of the bar
            cbar.ax.xaxis.set_label_position('bottom')
            cbar.ax.xaxis.tick_bottom()

             # horizontal label under vertical bar
            cbar.set_label(cbar_label, rotation=0, fontsize=12)
            cbar.ax.yaxis.set_label_coords(0.95, -0.05)
        
        fig.savefig("regime_plot.png", dpi=300)
        plt.show()
       
                
    def Lindhorst(self,m_vec,RSLs,sedqr,sedqs):
        
        # slopes
        slope1 = 0.0210
        slope2 = 0.0175
        slope3 = 0.0471
        
        # Ridge heights
        ridge_h1 = [0.8, 0.9, 1.1, 1.4, 1.7, 1.8, 2.2, 3.0, 3.1]
        ridge_h2 = [-0.1, 0.1, 0.6, 1.1, 1.5, 2.1, 2.7, 3.5, 4.9, 5.7]
        ridge_h3 = [0.5, 2.7, 5.8, 7.6, 9.5]

        # Ridge distances
        ridge_dis1 = [8, 12, 12, 12, 15]
        ridge_dis2 = [7, 14, 15, 17, 18, 20]
        ridge_dis3 = [35, 40, 55, 65]
        
        # Ridge period
        t1 = [250,266,280,290, 330]
        t2 = [150,167,134, 177, 207]
        t3 = [(4270 - 360) / 5, (4270 - 100) / 5,4270 / 5, (4270 + 100) / 5,(4270 + 360) / 5]

        # ridge number
        ridge_num1 = len(ridge_h1)
        ridge_num2 = len(ridge_h2)
        ridge_num3 = len(ridge_h3)

        colors = ['blue', 'orange', 'green']
        labels = ["Slope = 1.2$^\circ$", "Slope = 1$^\circ$", "Slope = 2.7$^\circ$"]
       

        small_ridgenumber = ([])
        small_height = ([])
        small_length = ([])
        small_time = ([])
        med_ridgenumber = ([])
        med_height = ([])
        med_length = ([])
        med_time = ([])
        large_ridgenumber = ([])
        large_height = ([])
        large_length = ([])
        large_time = ([])
        for i in range(0,len(RSLs)):
            self.a = RSLs[i]
            self.Qs = sedqs[i]
            self.Qr = sedqr[i]
            self.slope_params(m_vec)
            if i == 0:
                small_ridgenumber = self.ridgenumber
                small_height = self.height
                small_length = self.length
                small_time = self.time
            if i == 1:
                med_ridgenumber = self.ridgenumber
                med_height = self.height
                med_length = self.length
                med_time = self.time
            if i == 2:
                large_ridgenumber = self.ridgenumber
                large_height = self.height
                large_length = self.length
                large_time = self.time


        plt.rcParams.update({'font.size': 10})
        fig, ax = plt.subplots(2, 2, figsize=(8, 8))
        
        ############# ridge number
        plt.subplot(2, 2, 1)
        plt.plot(m_vec, small_ridgenumber,'lightcoral',marker='.')
        plt.plot(m_vec, med_ridgenumber,'r.')
        plt.plot(m_vec, large_ridgenumber,'lightcoral',marker='.')
        plt.fill_between(m_vec,small_ridgenumber,large_ridgenumber,color = 'pink',label = 'Our Model')
        
        slopes = [slope1,slope2,slope3]
        ridge_counts = [ridge_num1,ridge_num2,ridge_num3]
        coeffs = np.polyfit(slopes, ridge_counts, 1)
        fit_line = np.poly1d(coeffs)
        plt.plot(m_vec, fit_line(m_vec), '-', label='Lindhorst et al. 2014')
        
        for color,label,ridge_count,slope in zip(colors, labels, ridge_counts, slopes):
            plt.plot(slope,ridge_count,'o',label=label, color=color)
        plt.xlabel('Slope')
        plt.ylabel('Ridge number')

        ############## height
        plt.subplot(2, 2, 2)
        plt.plot(m_vec, small_height,'lightcoral',marker='.')
        plt.plot(m_vec, med_height,'r.')
        plt.plot(m_vec, large_height,'lightcoral',marker='.')
        plt.fill_between(m_vec,small_height,large_height,color = 'pink',label = 'Our Model')
        
        ridge_h = [np.mean(ridge_h1),np.mean(ridge_h2),np.mean(ridge_h3)]
        coeffs = np.polyfit(slopes, ridge_h, 1)
        fit_line = np.poly1d(coeffs)
        plt.plot(m_vec, fit_line(m_vec), '-', label='Lindhorst et. al 2014')
        
        for color,label,ridge_hh,slope in zip(colors, labels, ridge_h, slopes):
            plt.plot(slope,ridge_hh,'o',label=label, color=color)
        plt.xlabel('Slope')
        plt.ylabel('Ridge Height [m]')

        ################ distance between ridges --> length
        plt.subplot(2, 2, 3)
        plt.plot(m_vec, small_length,'lightcoral',marker='.')
        plt.plot(m_vec, med_length,'r.')
        plt.plot(m_vec, large_length,'lightcoral',marker='.')
        plt.fill_between(m_vec,small_length,large_length,color = 'pink',label = 'Our Model')

        ridge_dis = [np.mean(ridge_dis1),np.mean(ridge_dis2),np.mean(ridge_dis3)]
        coeffs = np.polyfit(slopes, ridge_dis, 1)
        fit_line = np.poly1d(coeffs)
        plt.plot(m_vec, fit_line(m_vec), '-', label='Lindhorst et al. 2014')
        
        for color,label,ridge_dist,slope in zip(colors, labels, ridge_dis, slopes):
            plt.plot(slope,ridge_dist,'o',label=label, color=color)
        plt.xlabel('Slope')
        plt.ylabel('Distance between Ridges [m]')

        ############## time per ridge
        plt.subplot(2, 2, 4)
        plt.plot(m_vec, small_time,'lightcoral',marker='.')
        plt.plot(m_vec, med_time,'r.')
        plt.plot(m_vec, large_time,'lightcoral',marker='.')
        plt.fill_between(m_vec,small_time,large_time,color = 'pink',label = 'Our Model')

        ridge_t = [np.mean(t1),np.mean(t2),np.mean(t3)]
        coeffs = np.polyfit(slopes, ridge_t, 1)
        fit_line = np.poly1d(coeffs)
        plt.plot(m_vec, fit_line(m_vec), '-', label='Lindhorst et al. 2014')
        
        for color,label,ridge_tt,slope in zip(colors, labels, ridge_t, slopes):
            plt.plot(slope,ridge_tt,'o',label=label, color=color)
        plt.xlabel('Slope')
        plt.ylabel('Period of Ridge Growth [yr]')

        plt.legend(loc = 'upper left')
        fig.tight_layout()

        fig.savefig("lindhorst.png", dpi=300) 
        plt.show()
        
