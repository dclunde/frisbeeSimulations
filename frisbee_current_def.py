# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 11:27:41 2016
 
@author: dlunde
"""
 
from frisbee_current_header import *
 
def lift_drag(i): 
    """  Function for updating lift and drag coeficents """        
    #global alpha, alpha_0
    #print math.degrees(angle_between(velocity_vec[i],velocity_on_frisbee[i]))
     
    if WIND:
        CL = CLO + CLA*angle_between(velocity_and_wind[i],wind_on_frisbee[i])
        CD = CDO + CDA*(angle_between(velocity_and_wind[i],wind_on_frisbee[i])-\
            math.radians(alpha_0))**2
    else:
        CL = CLO + CLA*angle_between(velocity_vec[i],velocity_on_frisbee[i])#alpha*math.pi/180.00 # Lift coefficent 
        CD = CDO + CDA*(angle_between(velocity_vec[i],velocity_on_frisbee[i])-\
            math.radians(alpha_0))**2#(alpha-alpha_0)*math.pi/180.0)**2 # drag coefficient
         
 
       
         
#    CLW = CLO + CLA*angle_between(wind,wind_on_frisbee[i])#alpha*math.pi/180.00 # Lift coefficent 
#    CDW = CDO + CDA*(angle_between(wind,wind_on_frisbee[i])-\
#        math.radians(alpha_0))**2#(alpha-alpha_0)*math.pi/180.0)**2 # drag coefficient
    
     
    # make lift negative if velocity is below the frisbee
#    if velocity_vec[i][1]>velocity_on_frisbee[i][1]:
#        CL = 0.5*CL
         
#    if frisbee_vec[1]<0:
#        CL = -CL*0.4
     
    #CL = CLO + CLA*alpha*math.pi/180.00 # Lift coefficent 
    #CD = CDO + CDA*((alpha-alpha_0)*math.pi/180.0)**2 # drag coefficient
         
#    return (CL,CD,CLW,CDW)
    return (CL,CD)
 
def forces(i):                   
    """ Function for updating 3D forces """
    if WIND:
        Rho_Vel_Area = 0.5 * rho * area * (np.linalg.norm(velocity_and_wind[i-1]))**2 
    else:
        Rho_Vel_Area = 0.5 * rho * area * (np.linalg.norm(velocity_vec[i-1]))**2
     
 
#    Rho_Vel_Area_Lift = 0.5 * rho * area * (np.linalg.norm(wind))**2 
 
#    cl,cd,cl_w,cd_w = lift_drag(i-1)             #Initialize lift and drag
    cl,cd = lift_drag(i-1)             #Initialize lift and drag
 
    #Inertia = mass * area / pi
 
    #return gravity (x,y,z) : lift (x,y,z) : drag (x,y,z) : magnus (x,y,z)
    grav_x = 0
    grav_y = g * mass 
    grav_z = 0 
     
#    lift_x = 0
#    lift_y = rho * velocity * area * cl * 0.5
#    lift_z = 0
#    
#    drag_x = -rho * velocity * area * cd * 0.5
#    drag_y = 0
#    drag_z = 0
     
    if WIND:
        lift_x,lift_y,lift_z = [Rho_Vel_Area* cl * value for value in lift_wind_vec[i-1]]
    else:
        lift_x,lift_y,lift_z = [Rho_Vel_Area* cl * value for value in lift_vec[i-1]]
     
    if WIND:
        drag_x,drag_y,drag_z = [-Rho_Vel_Area * cd * value for value in\
            unit_vector(velocity_and_wind[i-1])] 
    else:
        drag_x,drag_y,drag_z = [-Rho_Vel_Area * cd * value for value in\
            unit_vector(velocity_vec[i-1])]
         
         
#    w_lift_x,w_lift_y,w_lift_z = [Rho_Vel_Area_Lift* cl_w * value for value in lift_vec[i-1]]
#    
#    w_drag_x,w_drag_y,w_drag_z = [-Rho_Vel_Area_Lift * cd_w * value for value in\
#        unit_vector(wind)] 
     
#    angle = angle_between(velocity_vec[i-1],velocity_on_frisbee[i-1])
#    LnD1 = [Rho_Vel_Area* cl *\
#        math.sin(angle) * value - Rho_Vel_Area * cd * math.cos(angle) * value for value in unit_vector(velocity_on_frisbee[i-1])]
#    LnD3 = [Rho_Vel_Area* cl *\
#        math.cos(angle) * value + Rho_Vel_Area * cd * math.sin(angle) * value for value in lift_vec[i-1]]
#    
#
#    print angle,LnD1[0]+LnD3[0], lift_x + drag_x
     
    magnus_x = 0#B_mag * omega_y * (vy * math.sin(math.radians(phi)) - vz * math.cos(math.radians(phi)))
    magnus_y = 0#-B_mag * omega_y * vx * math.sin(math.radians(phi))
    magnus_z = 0#B_mag * omega_y * vx * math.cos(math.radians(90-phi))
     
    #    return grav_x, grav_y, grav_z ,LnD1[0],LnD1[1],LnD1[2], LnD3[0],LnD3[1],LnD3[2] , magnus_x, magnus_y, magnus_z
 
 
    #    deltavy = rho * velocity * area * cl * (1 / (2 * mass)) + g - B_mag * omega * vx * math.sin(math.radians(phi))
    #    deltavx = -rho * velocity * area * cd * (1 / (2 * mass)) + B_mag * omega * (vy * math.sin(math.radians(phi)) - vz * math.cos(math.radians(phi))
    #    deltavz =  B_mag * omega * vx * math.cos(math.radians(phi)) #+ -rho * velocity * area * cd
         
    return grav_x, grav_y, grav_z , lift_x, lift_y, lift_z, drag_x, \
        drag_y, drag_z, magnus_x, magnus_y, magnus_z,\
#        w_lift_x, w_lift_y, w_lift_z, w_drag_x, w_drag_y, w_drag_z
 
def moments(i,delta_omega_x,delta_omega_y,delta_omega_z):
    """ Function which calculates the new moments """
#    global phi,alpha
    phi = mat_p[i-1]
    alpha = mat_a[i-1]
    omega_x = mat_wx[i-1]
    omega_y = mat_wy[i-1]
    omega_z = mat_wz[i-1]
     
    if WIND:
        Rho_Vel_Area = 0.5 * rho * area * (np.linalg.norm(velocity_and_wind[i-1]))**2
        R = (CRr*omega_y+CRp*omega_x)*Rho_Vel_Area*diameter
        M = (CM0 + CMa * (angle_between(velocity_and_wind[i-1],wind_on_frisbee[i-1]))\
            +CMq*omega_z)*Rho_Vel_Area*diameter
        N = (CNr*omega_y)*Rho_Vel_Area*diameter
    else:
        Rho_Vel_Area = 0.5 * rho * area * (np.linalg.norm(velocity_vec[i-1]))**2
        R = (CRr*omega_y+CRp*omega_x)*Rho_Vel_Area*diameter
        M = (CM0 + CMa * (angle_between(velocity_vec[i-1],velocity_on_frisbee[i-1]))\
            +CMq*omega_z)*Rho_Vel_Area*diameter
        N = (CNr*omega_y)*Rho_Vel_Area*diameter
     
 
     
     
    #print R,M,N
    #R is in my X direction along frisbee
    #M is in my Z direction along frisbee
    #N is in my Y direction along frisbee    
    # my x is his x
    # my y is his z
    # my z is his y
    # his p is his wx
    # his q is his wy
    # his r is his wz
    # his theta is my alpha
    # his phi is my phi
    # his gamma is my omega
    # p = wx
    # q = wz
    # r = wy
    phi = phi - 90   
#    print R,M,N
     
    alpha_R = math.radians(alpha)    
     
    D_omegax = (R + \
        Ixx * delta_omega_z * delta_omega_x * math.sin(alpha_R) -\
        Iyy * delta_omega_z * (delta_omega_x * math.sin(alpha_R) + delta_omega_y) +\
        Ixx * delta_omega_z * delta_omega_x * math.sin(alpha_R)) * math.cos(alpha_R)
    D_omegaz = (M + \
        Iyy * delta_omega_x * math.cos(alpha_R) * (delta_omega_x * math.sin(alpha_R) + delta_omega_y) -\
        Ixx * delta_omega_x * delta_omega_x * math.cos(alpha_R) * math.sin(alpha_R))
    D_omegay = (N - \
        Iyy * ((D_omegax*delta_t/Ixx) * math.sin(alpha_R) +\
        delta_omega_z * delta_omega_x * math.cos(alpha_R)))
     
#    Translate = [[math.cos(phi),math.sin(phi)*math.sin(alpha),-math.sin(phi)*math.cos(phi)],
#                [0,math.cos(alpha),math.sin(alpha)],
#                [math.sin(phi),-math.cos(phi)*math.sin(alpha),math.cos(phi)*math.cos(alpha)]]    
     
     
    phi = phi + 90
     
    moments_analysis_sheet(i,R,N,M,D_omegax,D_omegay,D_omegaz,\
        delta_omega_x,delta_omega_y,delta_omega_z,alpha,math.cos(alpha_R),math.sin(alpha_R))
#    print "time_step",i*delta_t
#    print "w..x,y,z",D_omegax,D_omegay,D_omegaz
#    print "input_x,y,z",delta_omega_x,delta_omega_y,delta_omega_z
#    print "phi",phi,"alpha",alpha,"cos(a)",math.cos(alpha_R),"sin(a)",math.sin(alpha_R)
     
     
    #Translate the vectors into the current direction
#    tran = np.array(Translate)
#    convert = np.array([[D_omegax],[D_omegay],[D_omegaz]])  
#    change = tran.dot(convert)
#    print change
#    D_omegax = change[0][0]
#    D_omegay = change[1][0]
#    D_omegaz = change[2][0]
     
    return D_omegax, D_omegay, D_omegaz
 
 
 
def update_vectors(i):
    """ Updates the vectors according to the simulation """
    #global velocity_vec,frisbee_vec,velocity_on_frisbee,alpha,vx,vy,vz,phi
    velocity_vec[i] = [mat_vx[i],mat_vy[i],mat_vz[i]]
     
    frisbee_vec[i]  = [
        -math.cos(math.radians(90-mat_a[i])),
        0,
        math.sin((math.radians(90-mat_p[i])))
        ]
    frisbee_vec[i][1] = 1 - frisbee_vec[i][0]**2 - frisbee_vec[i][2]**2
    if 90 < mat_a[i] < 270 or -90 > mat_a[i] > -270:
        frisbee_vec[i][1] = - frisbee_vec[i][1]
    if 90 < mat_p[i]-90 < 270 or -90 > mat_p[i]-90 > -270:
        frisbee_vec[i][1] = - frisbee_vec[i][1]   
     
    velocity_on_frisbee[i] = map(operator.mul, frisbee_vec[i], \
        [np.dot(velocity_vec[i],frisbee_vec[i]) / (np.linalg.norm(frisbee_vec[i]))**2]*3)
    velocity_on_frisbee[i] = map(operator.sub, velocity_vec[i],velocity_on_frisbee[i])
     
    lift_vec[i] = np.cross(velocity_vec[i],np.cross(frisbee_vec[i],velocity_vec[i]))
    lift_vec[i] = unit_vector(lift_vec[i])
     
    wind_current = get_current_wind(mat_y[i])
    velocity_and_wind[i] = [mat_vx[i]-wind_current[0],mat_vy[i]-wind_current[1],mat_vz[i]-wind_current[2]]
     
    wind_on_frisbee[i] = map(operator.mul, frisbee_vec[i], \
        [np.dot(velocity_and_wind[i],frisbee_vec[i]) / (np.linalg.norm(frisbee_vec[i]))**2]*3)
    wind_on_frisbee[i] = map(operator.sub, velocity_and_wind[i],wind_on_frisbee[i])
     
    lift_wind_vec[i] = np.cross(velocity_and_wind[i],np.cross(frisbee_vec[i],velocity_and_wind[i]))
    lift_wind_vec[i] = unit_vector(lift_wind_vec[i])
     
#    wind_on_frisbee[i] = map(operator.mul, frisbee_vec[i], \
#    [np.dot(wind,frisbee_vec[i]) / (np.linalg.norm(frisbee_vec[i]))**2]*3)
#    wind_on_frisbee[i] = map(operator.sub, wind,wind_on_frisbee[i])
#    
#    lift_wind_vec[i] = np.cross(wind,np.cross(frisbee_vec[i],wind))
#    lift_wind_vec[i] = unit_vector(lift_wind_vec[i])
    return
 
def plot_frisbee(i,Fancy=False,html_mpld3=False,html_bokeh=False,html_plotly=False,\
    gif=False,ground=False,step=False,rotate=False,follow=False,Field=False,**kwargs):
    """ Plots the frisbee trajectory on its own figure plot """
#    global delta_t,mat_x,mat_y,mat_z,mat_a,time_slider
    global time_slider
    time_slider=0   
     
    fig = plt.figure("Frisbee Simulation")#,figsize=(10,10))
    fig.canvas.set_window_title('Frisbee Simulation') 
    ax = fig.gca(projection='3d')
     
    #Make velocity color map
    total_vel = [math.sqrt(mat_vx[jj]**2+mat_vy[jj]**2+mat_vz[jj]**2) for jj in range(0,i)]
    colormat = np.array(total_vel)
    colormat = colormat[::-1]
    colormat /= max(colormat)
     
     
    for j in range(1,i,1):
        ax.plot(mat_z[j-1:j+1], mat_x[j-1:j+1], mat_y[j-1:j+1], \
            label='Frisbee Flight', c=plotlib.cm.hot(colormat[j-1])) #RdYlGn or hot as another colormap
     
    #ax.plot(mat_z[0:i], mat_x[0:i], mat_y[0:i], label='Frisbee Flight', c=(0.1,0,0))
    plt.ylabel("X distance (m)")
    plt.xlabel("Z distance (m)")
    if absmax(mat_z) < 1:
        ax.set_xlim(-1, 1)
    ax.set_zlim(0,max(mat_y)+1)
     
    #ax.set_xlim(-10,10);
#    ax.set_ylim(-1, 50)
#    ax.set_zlim(-1, 50)
    #ax.set_aspect('equal')#, adjustable='box')
     
 
    #Wind Vector
    if WIND and absmax(wind) != 0:
        wind_max  = 6.0 / absmax(wind)
        arv_w = Arrow3D([0,wind[2]*wind_max],[0,wind[0]*wind_max],[0,wind[1]*wind_max], mutation_scale=20, lw=1, arrowstyle="-|>", color="green")
        ax.add_artist(arv_w) 
     
    if Field:
        Fancy=True
         
    if not follow:
        ax.w_zaxis.set_pane_color((0,1,0,0.5)) #Color the ground
    else:
        ax.w_zaxis.set_pane_color((0,0,0,0))
    ax.w_yaxis.set_pane_color((0,0,0,0))
    ax.w_xaxis.set_pane_color((0,0,0,0))
     
     
    if ground:
        make_ground(ax) #Make the ground
     
    scatter_point = [ax.scatter(mat_z[0],mat_x[0],mat_y[0], 'go')]
 
    #Other plots    
 
    if Fancy:
#        ax.set_aspect('equal')#,adjustable='box')
        ax.set_ylim(0,max(mat_x)+0.5)
#        if max(mat_x) > max(mat_z):
        ax.set_xlim(-max(mat_x)/2,max(mat_x)/2)  
        if Field:
            make_field(ax)                    
        plt.subplots_adjust(left=0.01, bottom=0.08,top = .99, right = .99)
        axamp = plt.axes([0.7, 0.025, 0.25, 0.03])
        resetax = plt.axes([0.38, 0.02, 0.1, 0.04])
        viewax = plt.axes([0.483, 0.02, 0.1, 0.04])
        ax_words = plt.axes([0, 0.02, 0.01, 0.04])
        add_words(ax_words,Fancy=True)
        make_frisbee(ax,Fancy=True)
    else:
        plt.subplots_adjust(left=0.01, bottom=0.35,top = .99, right = .99)
        ax_words = plt.axes([0.3,-0.1,0.5,0.5])
        axamp = plt.axes([0.68, 0.3, 0.25, 0.03])
        resetax = plt.axes([0.8, 0.2, 0.1, 0.04])
        viewax = plt.axes([0.8, 0.15, 0.1, 0.04])
        add_words(ax_words)  
        axdisk = plt.axes([0.1, 0.01, 0.3, 0.3], projection='3d')
        make_frisbee(axdisk) #create the disk to begin with
     
    samp = Slider(axamp, 'Time (sec)', 0, i*delta_t, valinit=0)#), valfmt='%0.0i')
    button = Button(resetax, 'Play', hovercolor='0.75')
    button2 = Button(viewax, 'Reset', hovercolor='0.75')
    if step:    
        button.label.set_text("Step") 
  
    world_lims = ax.get_w_lims()
     
    #axarrow = fig.add_subplot(234, axisbg=axcolor,projection = '3d')
    #axarrow.quiver(1,1,1,1,1,1)#, length =6)
    #ax.view_init(elev=30, azim=0)
    def update(time_slider_sec):
        """ Updates the plot according the the time_slider"""
        global time_slider
        time_slider = int(time_slider_sec/delta_t)
        if rotate:
            update_with_rotate(time_slider_sec)
        elif follow:
            update_with_follow(time_slider_sec)
        if Fancy:
            make_frisbee(ax,Fancy=True)        
        else:
            ax_words.cla() 
#            axdisk.cla()
            add_words(ax_words) 
            make_frisbee(axdisk) 
        scatter_point[0].remove()
        scatter_point[0] = ax.scatter(mat_z[time_slider],mat_x[time_slider],mat_y[time_slider], 'go')
        fig.canvas.draw_idle()
        if gif and time_slider_sec!=0:# and int(val*100/i)%10==0:
            sys.stdout.write("\r" + "{0}% finished with gif".format(time_slider*100/i))
            sys.stdout.flush()
            #print 'I',
        return scatter_point[0]
     
    def update_with_rotate(time_slider_sec):
        """" Updates the plot and rotates the camera view """
        global time_slider
#        update(val)
#        update(val*delta_t)
#        ax.elev = mat_y[time_slider] 
#        ax.azim = angle_between(velocity_vec[time_slider],[0,0,1],degrees=True)-180
        if Fancy:
            ax.set_xlim(mat_z[time_slider]-2,mat_z[time_slider]+2)
            ax.set_ylim(mat_x[time_slider]-2,mat_x[time_slider]+2)
            ax.elev=10
        ax.azim += 0.5
        #ax.view_init(elev=10., azim=val)
#        samp.set_val(val*delta_t)
#        fig.canvas.draw_idle()
        return scatter_point[0]
         
    def update_with_follow(time_slider_sec):
        global time_slider
        if Fancy:
            ax.set_xlim(mat_z[time_slider]-1,mat_z[time_slider]+1)
            ax.set_ylim(mat_x[time_slider]-1,mat_x[time_slider]+1)
            ax.set_zlim(mat_y[time_slider]-1,mat_y[time_slider]+1)
#        ax.elev=mat_y[time_slider]
        ax.elev=1
        follow_angle = angle_between((mat_vx[time_slider],0,mat_vz[time_slider]),(mat_vx[time_slider],0,0),degrees=True)
        if mat_vz[time_slider]<0:
            follow_angle = -follow_angle
        ax.azim = -90-follow_angle
        return
         
    samp.on_changed(update) #when slider pressed call update
 
    def play(event):
        """ Controls the play button"""
        global playing, time_slider
        if step:
            new_slider = time_slider + 2
            samp.set_val(new_slider*delta_t)
            return
         
        if playing:
            playing = False
            button.label.set_text("Play") 
            return
        else:
            playing = True
            button.label.set_text("Stop") 
            play_through()
            return
     
    def play_through():
        """ Controls the automatic animation"""
        global playing, time_slider
        for II in range (time_slider,i,3):     
            samp.set_val(II*delta_t)
            update(II*delta_t)
#            plt.pause(delta_t/100000)
            plt.pause(0.00000001)
            if not playing:
                return
        playing = False
        button.label.set_text("Play") 
        return
 
    button.on_clicked(play)
 
    def reset(event):
        global time_slider
        time_slider = 0
        samp.set_val(0)
        ax.view_init(elev=30, azim = -60)
         
        ax.set_xlim(world_lims[0],world_lims[1])
        ax.set_ylim(world_lims[2],world_lims[3])
        ax.set_zlim(world_lims[4],world_lims[5])
        button2
        return
 
    button2.on_clicked(reset)
 
     
#    rax = plt.axes([0.5, 0, 0.2, 0.2], axisbg=axcolor)
#    radio = RadioButtons(rax, ('normal', 'down trajectory'), active=1)
#    def colorfunc(label):
#        global time_slider
#        print "WORKING"
#        if label == 'normal':
#            time_slider = 0
#            update(0,rotate=False)
#            ax.view_init(elev=30, azim = -60)
#        #fig.canvas.draw_idle()
#    radio.on_clicked(colorfunc)
    def save_gif(fig):
        """ Saves a gif of the animation """
        global time_slider
#        print np.arange(0, i*delta_t, 3*delta_t)
        anim = FuncAnimation(fig,  samp.set_val, frames=np.linspace(0, i*delta_t-5*delta_t, 75),\
            interval=140)#,blit=True)#,fargs(rotate=True))
#        if True:
        if True: #make True for gif
            title =  'output/%s.gif' % (time.strftime("%H_%M_%S"))
            anim.save(title, dpi=150, writer='imagemagick')
        else:                
            title =  'output/%s.mp4' % (time.strftime("%H_%M_%S"))
            anim.save(title, dpi=150,fps=15, writer = "avconv", codec = "libx264")
        print title
        os.system("mv %s ../../Downloads/" % (title))            
        reset(0)
        print "Saved gif"
        os.system("mpg123 audio.mp3 >> /dev/null")
        plt.close("Frisbee Simulation")
        return
#        else:
#            print "Now showing plot"
#            plt.show()
             
 
    if gif:
        save_gif(fig)
    if html_mpld3:
        mpld3_frisbee(i)
    if html_bokeh:
        bokeh_frisbee(i)
    if html_plotly:
        plotly_frisbee(i)
     
#    if Fancy:
#        set_axes_equal(ax)
     
    plt.show()#block=True) #doesn't allow it to be closed.
#    plt.draw()
    return
     
def make_frisbee(ax,Fancy=False,new_input = 0,new_frisbee_vec=0):
    """ Creates a frisbee on an axis """
    global time_slider,frisbee_vec,frisbee_parts
     
     
    current_frisbee_vec = frisbee_vec[time_slider]
    current_position = [mat_z[time_slider], mat_x[time_slider], mat_y[time_slider]]
    if not new_input == 0:
        current_position = [new_input[3],new_input[1],new_input[2]]
    if not new_frisbee_vec == 0:
        current_frisbee_vec = new_frisbee_vec
    #Remove the old frisbee
    for item in frisbee_parts:     
         try:
            item.remove()
         except:
             pass
#            print "No frisbee to remove"
 
 
    ax.set_axisbelow(True)
 
    # Create Vectors
 
    if Fancy:
        transparency = 0.2
        linewidth = 1
    else:
        transparency = 0.1
        linewidth = 0.5
         
    #TO DO: Uncomment when actually working on this
    if new_input==0:
        draw_vectors(ax,Fancy = Fancy)
 
 
    #norm = plotlib.colors.Normalize()
    #surf = ax1.plot_surface(xig, yig, zi, facecolors=cm.rainbow(norm(ci)))
    scolor = 'blue'#'inferno' #blue
     
#    height = 1
#    alpha_z = height * math.sin(math.radians(90-mat_p[time_slider]))
#    alpha_x = -height * math.cos(math.radians(90-mat_a[time_slider]))
#    alpha_y =  1 - alpha_x**2 - alpha_z**2 
#    if 90 < mat_a[time_slider] < 360:
#        alpha_y = - alpha_y
    #print alpha_x,alpha_y,alpha_z
    #print "x",alpha_x,"y",alpha_y,"z",alpha_z,"i",alpha_i
     
    #axis and radius
    frisbee_direction = unit_vector(current_frisbee_vec)
    if Fancy: 
        p1 = np.array([frisbee_direction[2], frisbee_direction[0], frisbee_direction[1] ]) #point at other end      
        p0 = np.array(current_position) #point at one end
        p2 = [p0[j]+(p1[j])/6 for j in [0,1,2]] 
        p1 = [p0[j]+(p1[j])/4 for j in [0,1,2]]  
#        R = 0.7
        R = ax.get_xlim()[1]/6
    else:
        p1 = np.array([current_frisbee_vec[2], current_frisbee_vec[0], current_frisbee_vec[1]]) #point at other end
        p0 = np.array([0, 0, 0]) #point at one end
        p2 = p1
        R = 3
     
     
    if False:
        R = 3
    #print p0,'   ',p1,frisbee_direction[1]
    R2 = R * 0.8
 
    #vector in direction of axis
    v = p1 - p0
    v2 = p2 - p0
     
    #find magnitude of vector
    mag = np.linalg.norm(v)
    mag2 = np.linalg.norm(v2)
 
    #unit vector in direction of axis
    v = v / mag
    v2 = v2 / mag2
 
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    not_v2 = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    if (v2 == not_v2).all():
        not_v2 = np.array([0, 1, 0])
     
    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    n12 = np.cross(v2, not_v2)
    #normalize n1
    n1 /= np.linalg.norm(n1)
    n12 /= np.linalg.norm(n12)
     
    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    n22 = np.cross(v2, n12)
     
     
    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 2)
    t2 = np.linspace(0, mag2, 2)
     
    if mat_wy[time_slider] > 0:
        rotate_var = time_slider/2
    else:
        rotate_var = -time_slider/2
     
    theta = np.linspace(0+rotate_var, 2 * np.pi +rotate_var, 100)
    theta3 = np.linspace(0+rotate_var, 2 * np.pi +rotate_var, 100)
    rsample = np.linspace(0, R, 2)
    rsample2 = np.linspace(0, R2, 2)
 
    #use meshgrid to make 2d arrays
    t, theta2 = np.meshgrid(t, theta)
    t2, theta22 = np.meshgrid(t2, theta3)
     
    rsample,theta = np.meshgrid(rsample, theta)
    rsample2,theta3 = np.meshgrid(rsample2, theta3)
     
    #generate coordinates for surface
    # "Tube"
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta2) * n1[i] + R * np.cos(theta2) *       n2[i] for i in [0, 1, 2]]
    # "Bottom"
    X2, Y2, Z2 = [p0[i] + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    # "Top"
    X3, Y3, Z3 = [p0[i] + v[i]*mag + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    #Inner Tube
    X_2, Y_2, Z_2 = [p0[i] + v2[i] * t2 + R2 * np.sin(theta22) * n12[i] + R2 * np.cos(theta22) *       n22[i] for i in [0, 1, 2]]
 
     
     
#    fig=plt.figure()
#
#    ax=plt.subplot(111, projection='3d')
    #ax.cla()   
    frisbee_parts[0] = ax.plot_surface(X, Y, Z, color=scolor,alpha=transparency,linewidth = linewidth)#,facecolors=cm.rainbow(norm(ci)))
    frisbee_parts[2] = ax.plot_surface(X3, Y3, Z3, color=scolor,alpha=transparency,linewidth = 0)
    if False:
        frisbee_parts[1] = ax.plot_surface(X_2, Y_2, Z_2, color=scolor,alpha=transparency,linewidth = linewidth-0.2)#,facecolors=cm.rainbow(norm(ci)))
    else:
        frisbee_parts[1] = ax.plot_surface(X2, Y2, Z2, color=scolor,alpha=transparency,linewidth = 0)
 
     
     
    if not Fancy:    
        frisbee_size = 5   
        ax.set_zlim(-frisbee_size,frisbee_size)
        ax.set_xlim(-frisbee_size,frisbee_size)
        ax.set_ylim(-frisbee_size,frisbee_size)
    #plt.draw()
    plt.show()
    return
     
def make_field(ax):
     
    width   = 48.9204
    length  = 73.152
    end_zone= 18.288
     
    for item in cones:     
         try:
            item.remove()
         except:
             pass
     
#    ax.set_xlim(-width/2,width/2)
    ax.set_xlim(-2*width/3,2*width/3)
    ax.set_ylim(-end_zone,length+end_zone)
    ax.set_zlim(0,6)
     
    ax.azim = -20
    ax.elev = 15
     
    #behind end zone
    cones[0] = [ax.scatter(-width/2,0,0, color='orange')]
    cones[1] = [ax.scatter(width/2,0,0, color='orange')]
    cones[2] = [ax.scatter(-width/2,-end_zone,0, color='orange')]
    cones[3] = [ax.scatter(width/2,-end_zone,0, color='orange')]
    cones[4] = [ax.scatter(0,-end_zone,0, color='orange')]
     
    #far end zone
    cones[5] = [ax.scatter(-width/2,length,0, color='orange')]
    cones[6] = [ax.scatter(width/2,length,0, color='orange')]
    cones[7] = [ax.scatter(-width/2,length+end_zone,0, color='orange')]
    cones[8] = [ax.scatter(width/2,length+end_zone,0, color='orange')]
    cones[9] = [ax.scatter(0,length+end_zone,0, color='orange')]
     
    #side lines
#    cones[10] = [ax.scatter(-width/2,length/2,0, color='orange')]
#    cones[11] = [ax.scatter(width/2,length/2,0, color='orange')]
    cones[12] = [ax.scatter(-width/2,length/3,0, color='orange')]
    cones[13] = [ax.scatter(width/2,length/3,0, color='orange')]
    cones[14] = [ax.scatter(-width/2,3*length/4,0, color='orange')]
    cones[15] = [ax.scatter(width/2,3*length/4,0, color='orange')]
     
    #end zone lines
    width_line = np.arange(-width/2,width/2,0.5)
    lineColor = 'white'
    line_width = 1.5
    cones[16] = [ax.plot(width_line,[0]*np.size(width_line),0,color=lineColor,linewidth=line_width)]
    cones[17] = [ax.plot(width_line,[length]*np.size(width_line),0,color=lineColor,linewidth=line_width)]
     
    side_line = np.arange(-end_zone,length+end_zone,0.5)
    cones[18] = [ax.plot([-width/2]*np.size(side_line),side_line,0,color=lineColor,linewidth=line_width)]
    cones[19] = [ax.plot([width/2]*np.size(side_line),side_line,0,color=lineColor,linewidth=line_width)]
     
    cones[20] = [ax.plot(width_line,[-end_zone]*np.size(width_line),0,color=lineColor,linewidth=line_width)]
    cones[22] = [ax.plot(width_line,[length+end_zone]*np.size(width_line),0,color=lineColor,linewidth=line_width)]
    return
 
     
         
 
def plotly_frisbee(i):
    """ Creates a plotly html plot of the 3D plot """
    import matplotlib.mlab as mlab
    import plotly.plotly as py
    import plotly.graph_objs as go
 
#    fig = plt.figure()
#    fig.canvas.set_window_title('Frisbee Simulation HTML') 
#    ax = fig.gca(projection='3d')
#    
#    ax.plot(mat_z[0:i], mat_x[0:i], mat_y[0:i], label='Frisbee Flight')
#    plt.ylabel("X distance (m)")
#    plt.xlabel("Z distance (m)")
#    if absmax(mat_z) < 1:
#        ax.set_xlim(-1, 1)
 
    trace = go.Scatter3d(x=mat_z[0:i],y=mat_x[0:i],z=mat_y[0:i],\
        line=dict(width=0.5),mode = 'markers')
    data = [trace]
    py.sign_in('dlunde765', 'aHSpfE4OSBqhVDIVYFSW')
     
    #plot_url = py.plot_mpl(fig)  
    #py.plot(data)  
    figgg = dict(data=data)    
    py.iplot(figgg,filename='frisbee plot')  
          
def bokeh_frisbee(i):
    """ Creates a bokeh html file of the 3D plot """
    from bokeh import mpl
    from bokeh.plotting import output_file,save
 
    fig = plt.figure()
    fig.canvas.set_window_title('Frisbee Simulation HTML') 
    ax = fig.gca(projection='3d')
     
    ax.plot(mat_z[0:i], mat_x[0:i], mat_y[0:i], label='Frisbee Flight')
    plt.ylabel("X distance (m)")
    plt.xlabel("Z distance (m)")
    if absmax(mat_z) < 1:
        ax.set_xlim(-1, 1)
 
     
    output_file("online_fig.html")
    save(mpl.to_bokeh(fig))
    return
     
def mpld3_frisbee(i):
    """ Creates a mpld3 html file of the 3D plot """
 
    fig = plt.figure()
    fig.canvas.set_window_title('Frisbee Simulation HTML') 
    ax = fig.gca(projection='3d')
     
    ax.plot(mat_z[0:i], mat_x[0:i], mat_y[0:i], label='Frisbee Flight')
    plt.ylabel("X distance (m)")
    plt.xlabel("Z distance (m)")
    if absmax(mat_z) < 1:
        ax.set_xlim(-1, 1)
 
    try:
        print "saving html"
        mpld3.save_html(fig,"online_fig.html")
        print "Saved html"
    except:
        print "Failed??"
    else:
        print "YAY!"
    mpld3.save_html(fig,"online_fig.html")
    return
    
def make_ground(ax):
     #Ground
    point_G  = np.array([1, 1, 0])
    normal_G = np.array([0, 0, 1])
    d = -point_G.dot(normal_G)
 
    X_ground, Y_ground = np.meshgrid(
        np.arange(ax.get_xlim()[0],ax.get_xlim()[1],0.01),
        np.arange(ax.get_ylim()[0],ax.get_ylim()[1],0.01))
     
    Z_ground = (-normal_G[0] * X_ground - normal_G[1] * Y_ground - d) * 1. /normal_G[2]
 
    print X_ground
    ax.plot_surface(X_ground,Y_ground,Z_ground,\
        linewidth = 0,alpha = 0.75, color = 'lightgreen',antialiased = False)
    return
 
def add_words(ax_words,Fancy=False):
    """ Adds words to a axis changing as the interactive moves along """
    global mat_a,mat_vx,time_slider
     
    # build a rectangle in axes coords
    left, width = .25, .5
    bottom, height = .25, .5
    #right = left + width
    top = bottom + height
     
     
    # axes coordinates are 0,0 is bottom left and 1,1 is upper right
    if not Fancy:    
        p = patches.Rectangle(
            (left, bottom), width, height,
            fill=False, transform=ax_words.transAxes, clip_on=False
            )
        ax_words.add_patch(p)
 
    if Fancy:
        print_this = "Try these buttons! Or click along the slider"
    else:
        print_this = "Time (s) = %05.2f \
            \nVelocity X (m/s) = %05.2f\nVelocity Y (m/s) = %05.2f\nVelocity Z (m/s) = %05.2f \
            \nA Vel X (rad/s) = %05.2f\nA Vel Y (rad/s) = %05.2f\nA Vel Z (rad/s) = %05.2f \
            \nPitch (Degrees) = %05.2f \nRoll (Degrees) = %05.2f" \
            % (mat_t[time_slider],\
                mat_vx[time_slider],mat_vy[time_slider],mat_vz[time_slider],\
                mat_wx[time_slider],mat_wy[time_slider],mat_wz[time_slider],\
                mat_a[time_slider],mat_p[time_slider]-90)
    #ax_words.text(left+0.1*left, 0.5*(bottom+top), print_this,
#        horizontalalignment='left',
    if Fancy:
        horizontalalignment = 'left'
        fontsize = 10
    else:
        horizontalalignment = 'right'
        fontsize = 10
     
    ax_words.text(left+0.9*width, 0.5*(bottom+top), print_this,
        horizontalalignment=horizontalalignment,
        verticalalignment='center',
        fontsize=fontsize, color='black',
        transform=ax_words.transAxes)
         
    ax_words.set_axis_off()
    return
     
def analyze_plots(i,html=False):
    """ Creates a figure with plots of different variables """
    fig_analyze = plt.figure('Frisbee Analysis')
    fig_analyze.canvas.set_window_title('Frisbee Analysis') 
     
 
    #ax = fig_analyze.gca()
    #ax = plt.axes([0.1,0.1,0.4,0.4])
    ax  = fig_analyze.add_subplot(221)
    ax_1 = ax.twinx()
     
    ax.plot(mat_t[0:i], mat_x[0:i], label='X Direction',color='r')
    ax.set_ylabel("X Distance (m)")
    ax.set_xlabel("Time (s)")
    ax_1.plot(mat_t[0:i], mat_y[0:i], label='Y Direction')
    ax_1.plot(mat_t[0:i], mat_z[0:i], label='Z Direction')
    ax_1.set_ylabel("Y,Z Distance (m)")
    ax_1.set_xlabel("Time (s)")
    #plt.legend(fontsize='small',loc=0)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_1.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc=1,fontsize = 'small',framealpha=0.1)
     
    #ax2 = plt.axes([0.1,0.55,0.4,0.4])
    ax2  = fig_analyze.add_subplot(222)
     
    ax2.plot(mat_t[0:i], mat_vx[0:i], label='X Velocity')
    ax2.plot(mat_t[0:i], mat_vy[0:i], label='Y Velocity')
    ax2.plot(mat_t[0:i], mat_vz[0:i], label='Z Velocity')
    ax2.set_ylabel("Velocity (m/s)")
    ax2.set_xlabel("Time (s)")
    set_lims(mat_vx[0:i],mat_vy[0:i],mat_vz[0:i],ax2)
    plt.legend(fontsize='small',loc=1,framealpha=0.1)
     
    ax3  = fig_analyze.add_subplot(223)
    ax3_1 = ax3.twinx()
     
    ax3.plot(mat_t[0:i], mat_a[0:i], label='Alpha (Z)',color='r')
    ax3.set_ylabel("Alpha Degrees")
    ax3.set_xlabel("Time (s)")
    ax3_1.plot(mat_t[0:i], map(operator.sub,mat_p[0:i],[90]*i), label='Phi (X)')
    ax3_1.set_ylabel("Phi Degrees")
    h3, l3 = ax3.get_legend_handles_labels()
    h4, l4 = ax3_1.get_legend_handles_labels()
    ax3.legend(h3+h4, l3+l4, loc=1,fontsize = 'small',framealpha=0.1)
     
    ax4  = fig_analyze.add_subplot(224)
     
    ax4.plot(mat_t[0:i], mat_wx[0:i], label='Angular X Velocity')
    ax4.plot(mat_t[0:i], mat_wy[0:i], label='Angular Y Velocity')
    ax4.plot(mat_t[0:i], mat_wz[0:i], label='Angular Z Velocity')
    ax4.set_ylabel("Angular Velocity (rad/s)")
    ax4.set_xlabel("Time (s)")
    set_lims(mat_wx[0:i],mat_wy[0:i],mat_wz[0:i],ax4)
    plt.legend(fontsize='small',loc=1,framealpha=0.1)    
     
    plt.tight_layout()
 
    if html:    
        mpld3.save_html(fig_analyze,"analyze.html")
         
     
    plt.draw()
    plt.show()
    return
        
def calculate_energy(i,plot=False,save_data=False):
     
    mat_PE[0:i] = [mat_y[index] * mass * -g for index in range(0,i,1)]
     
    mat_vxE[0:i] = [0.5 * mass * mat_vx[index]**2 for index in range(0,i,1)]
    mat_vyE[0:i] = [0.5 * mass * mat_vy[index]**2 for index in range(0,i,1)]
    mat_vzE[0:i] = [0.5 * mass * mat_vz[index]**2 for index in range(0,i,1)]
     
    mat_KE[0:i] = [mat_vxE[index] + mat_vyE[index] + mat_vzE[index] \
        for index in range(0,i,1)]
     
    mat_wxE[0:i] = [0.5 * Ixx * mat_wx[index]**2 for index in range(0,i,1)]
    mat_wyE[0:i] = [0.5 * Iyy * mat_wy[index]**2 for index in range(0,i,1)]
    mat_wzE[0:i] = [0.5 * Izz * mat_wz[index]**2 for index in range(0,i,1)]
     
    mat_AE[0:i] = [mat_wxE[index] + mat_wyE[index] + mat_wzE[index] \
        for index in range(0,i,1)]
    
    
    mat_TE = [mat_KE[index] + mat_PE[index] + mat_AE[index] \
        for index in range(0,i,1)]
     
    if plot:
        plot_energy(i,mat_TE)
    if save_data:
        energy_sheet(i,mat_TE)
     
def plot_energy(i,mat_TE):
    """ Creates a figure with plots of the energy """
    fig_analyze = plt.figure('Frisbee Energy')
    fig_analyze.canvas.set_window_title('Frisbee Energy') 
 
 
    ax  = fig_analyze.add_subplot(221)
    ax_1 = ax.twinx()
     
    ax.plot(mat_t[0:i], mat_PE[0:i], label='Potential Energy',color='r')
    ax.set_ylabel("P Energy (J)")
    ax.set_xlabel("Time (s)")
    ax_1.plot(mat_t[0:i], mat_KE[0:i], label='Kinetic Energy')
    ax_1.plot(mat_t[0:i], mat_AE[0:i], label='Angular K Energy')
    ax_1.set_ylabel("K Energy (J)")
    ax_1.set_xlabel("Time (s)")
    #plt.legend(fontsize='small',loc=0)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_1.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc=1,fontsize = 'small',framealpha=0.1)
#    
    ax2  = fig_analyze.add_subplot(223)
     
    ax2.plot(mat_t[0:i], mat_TE[0:i], label='Total Energy')
    ax2.set_ylabel("Energy (J)")
    ax2.set_xlabel("Time (s)")
    plt.legend(fontsize='small',loc=1,framealpha=0.1)
     
     
    ax3  = fig_analyze.add_subplot(222)
    ax3_1 = ax3.twinx()
     
    ax3.plot(mat_t[0:i], mat_vxE[0:i], label='X Kinetic',color='r')
    ax3.set_ylabel("X Energy (J)")
    ax3.set_xlabel("Time (s)")
    ax3_1.plot(mat_t[0:i], mat_vyE[0:i], label='Y Kinetic')
    ax3_1.plot(mat_t[0:i], mat_vzE[0:i], label='Z Kinetic')
    ax3_1.set_ylabel("Y,Z Energy (J)")
    set_lims(mat_vyE[0:i],mat_vzE[0:i],mat_vz[0:i],ax3_1)
    h3, l3 = ax3.get_legend_handles_labels()
    h4, l4 = ax3_1.get_legend_handles_labels()
    ax3.legend(h3+h4, l3+l4, loc=1,fontsize = 'small',framealpha=0.1)
     
     
    ax4  = fig_analyze.add_subplot(224)
     
    ax4.plot(mat_t[0:i], mat_wxE[0:i], label='Kinetic Angular X')
    ax4.plot(mat_t[0:i], mat_wyE[0:i], label='Kinetic Angular Y')
    ax4.plot(mat_t[0:i], mat_wzE[0:i], label='Kinetic Angular Z')
    set_lims(mat_wxE[0:i],mat_wyE[0:i],mat_wzE[0:i],ax4)
    ax4.set_ylabel("Energy (J)")
    ax4.set_xlabel("Time (s)")
    plt.legend(fontsize='small',loc=1,framealpha=0.1)    
     
    plt.tight_layout()
    
def forces_sheet(forces_output,i):
    """ Prints a Forces excel sheet of the formulas"""
    f = open('output/Forces_total.csv','a')
    print >> f, i,",",\
        forces_output[0],",",forces_output[1],",",forces_output[2],",",\
        forces_output[3],",",forces_output[4],",",forces_output[5],",",\
        forces_output[6],",",forces_output[7],",",forces_output[8],",",\
        forces_output[9],",",forces_output[10],",",forces_output[11]
    f.close()
    return
     
def moments_sheet(moments_output,i):
    """ Prints a Forces excel sheet of the formulas"""
    f = open('output/Moments_total.csv','a')
    print >> f, i,",",\
        moments_output[0],",",moments_output[1],",",moments_output[2]
    f.close()
    return
     
def moments_analysis_sheet(i,R,N,M,D_omegax,D_omegay,D_omegaz,\
        delta_omega_x,delta_omega_y,delta_omega_z,alpha,COS,SIN):
    """ Prints a moments excel sheet of the formulas"""
    f = open('output/Moments_analysis.csv','a')
    print >> f, i,",",\
        R,",",N,",",M,",",D_omegax,",",D_omegay,",",D_omegaz,",",\
        delta_omega_x,",",delta_omega_y,",",delta_omega_z,",",alpha,",",COS,",",SIN
    f.close()
    return
 
def energy_sheet(i,mat_TE):
    """ Prints a Summary sheet of data at each time step"""   
    f = open('output/Energy_sheet.csv','w')
    f.write("i value,X Kinetic (J),Y Kinetic (J), Z Kinetic (J) ,\
        Angular Kinetic X, Angular Kinetic Y, Angular Kinetic Z,\
        Gravitational Potential,Total Kinetic,Total Angular Kinetic Y,\
        Total Energy \n")    
    for row in range(0,i,1):
        print >> f, row,",",\
            mat_vxE[row],",",mat_vyE[row],",",mat_vzE[row],",",\
            mat_wxE[row],",",mat_wyE[row],",",mat_wzE[row],",",\
            mat_PE[row],",",mat_KE[row],",",mat_AE[row],",",\
            mat_TE[row]
    f.close()
    print "Saved Energy Sheet"
    return
     
def excel_sheet(i,Open_files = False):
    """ Prints a Summary sheet of data at each time step"""
    global mat_x,mat_y,mat_z,mat_a,mat_wy,mat_vx,mat_vy,mat_vz
     
    f = open('output/Summary_sheet.csv','w')
    f.write("i value,X position (m),Y position (m), Z position (m) ,\
        Velocity X, Velocity Y, Velocity Z,\
        Alpha,Pitch,Angular Velocity Y,\
        Angular Velocity X, Angular Velocity Z,\
        Angle of Attack \n")    
    for row in range(0,i,1):
        print >> f, row,",",\
            mat_x[row],",",mat_y[row],",",mat_z[row],",",\
            mat_vx[row],",",mat_vy[row],",",mat_vz[row],",",\
            mat_a[row],",",mat_p[row]-90,",",mat_wy[row],",",\
            mat_wx[row],",",mat_wz[row],",",\
            angle_between(
                velocity_vec[row],
                velocity_on_frisbee[row],
                degrees=True)
    f.close()
    if Open_files:
        os.system("libreoffice output/Summary_sheet.csv output/Forces_total.csv &")
    print "Saved Summary Sheet"
    return
 
def analyze(i):
    """ Prints an Analysis of the frisbee flight"""
    global may_y,delta_t
    print "Maximum height of frisbee - ",max(mat_y),"meters"
    print "Time in Air              - ",i*delta_t,"seconds"
    return
 
def draw_vectors(ax,Fancy = False):
    """ Draws vecotors on an axis """
    global frisbee_vec, velocity_vec,velocity_on_frisbee,time_slider,lift_vec
     
    for item in arrows:       
        try:
            item.remove()
        except:
            pass
#            print "No vector to remove"
     
    vector_head = 10
 
    #Trying to make the vecotors fit in the box
#    total_max=absmax([\
#        absmax([row[2] for row in velocity_vec]),\
#        absmax([row[0] for row in velocity_vec]),\
#        absmax([row[1] for row in velocity_vec])])
    total_max=absmax(velocity_vec[0])+4
    if total_max > 0 and total_max < 1000:
        total_max = 6.0/total_max
    else:
        total_max = 1
      
    
    if Fancy:
        initial_pt= [mat_x[time_slider],mat_y[time_slider],mat_z[time_slider]]
        scale =[0.5,total_max,1,0.5]
#        scale =[0.5,0.3,1,0.5]
    else:
        initial_pt= [0,0,0]
        scale = [4,total_max,4,4]
     
    #Frisbee Perpendicular Vector        
    final_pt = [initial_pt[index]+scale[0]*frisbee_vec[time_slider][index] for index in [0,1,2]]
    arr_f = Arrow3D(\
        [initial_pt[2],final_pt[2]],\
        [initial_pt[0],final_pt[0]],\
        [initial_pt[1],final_pt[1]],\
        mutation_scale=vector_head, lw=1.5, arrowstyle="-|>", color="k")
     
    #Velocity Vector
    final_pt = [initial_pt[index]+scale[1]*velocity_vec[time_slider][index] for index in [0,1,2]]
    arv_v = Arrow3D(\
        [initial_pt[2],final_pt[2]],\
        [initial_pt[0],final_pt[0]],\
        [initial_pt[1],final_pt[1]],\
        mutation_scale=vector_head, lw=1.5, arrowstyle="-|>", color="blue")
 
    #Velocity on Frisbee Vector
    final_pt = [initial_pt[index]+scale[2]*unit_vector(velocity_on_frisbee[time_slider])[index] for index in [0,1,2]]
    arv_fv = Arrow3D(\
        [initial_pt[2],final_pt[2]],\
        [initial_pt[0],final_pt[0]],\
        [initial_pt[1],final_pt[1]],\
        mutation_scale=vector_head, lw=1.5, arrowstyle="-|>", color="green")
     
    #Lift Vector
    final_pt = [initial_pt[index]+scale[3]*lift_vec[time_slider][index] for index in [0,1,2]]
    arv_lv = Arrow3D(\
        [initial_pt[2],final_pt[2]],\
        [initial_pt[0],final_pt[0]],\
        [initial_pt[1],final_pt[1]],\
        mutation_scale=vector_head, lw=1.5, arrowstyle="-|>", color="red")
         
    #Wind and Velocity Vector
    final_pt = [initial_pt[index]+scale[1]*velocity_and_wind[time_slider][index] for index in [0,1,2]]
    arv_wv = Arrow3D(\
        [initial_pt[2],final_pt[2]],\
        [initial_pt[0],final_pt[0]],\
        [initial_pt[1],final_pt[1]],\
        mutation_scale=vector_head, lw=1.5, arrowstyle="-|>", color="orange")
         
    #Wind Lift Vector
    final_pt = [initial_pt[index]+scale[3]*lift_wind_vec[time_slider][index] for index in [0,1,2]]
    arv_wlv = Arrow3D(\
        [initial_pt[2],final_pt[2]],\
        [initial_pt[0],final_pt[0]],\
        [initial_pt[1],final_pt[1]],\
        mutation_scale=vector_head, lw=1.5, arrowstyle="-|>", color="pink")
         
    #print time_slider,velocity_vec[time_slider],frisbee_vec[time_slider],velocity_on_frisbee[time_slider]
    if WIND:
        arrows[4] = ax.add_artist(arv_wv)
        arrows[5] = ax.add_artist(arv_wlv)    
    arrows[0] = ax.add_artist(arr_f)    
    arrows[1] = ax.add_artist(arv_v)
    arrows[2] = ax.add_artist(arv_fv)
    arrows[3] = ax.add_artist(arv_lv)
     
 
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
 
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
         
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    if np.linalg.norm(vector)==0:
        return vector
    else:
        return vector / np.linalg.norm(vector)
 
def angle_between(v1, v2,degrees = False):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
 
            >>> angle_between((1, 0, 0), (0, 1, 0))
            >>> 1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            >>> 0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            >>> 3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    if degrees:
        return np.rad2deg(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    else:
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
 
def absmax(x):
    """ Returns the absolute max of a list """
    #    try:   
    #        list_a.index(tot_m)
    #        return tot_m
    #    except ValueError:
    #        return list_a[list_a.index(-tot_m)]
    list_a = []
    list_a = x
    tot_m  = max(abs(i) for i in list_a)
    return tot_m
   
def print_alpha():
    """ Prints the angle of attack whenever called"""
    global alpha
    print alpha
     
def set_lims(mat_A,mat_B,mat_C,AX):
    MAX =  max(max(mat_A),max(mat_B),max(mat_C))+1      
    MIN =  min(min(mat_A),min(mat_B),min(mat_C))-1      
    AX.set_ylim(MIN,MAX)
    return 
     
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.
 
    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
 
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
 
    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)
 
    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])
 
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    return
 
def get_current_wind(height): 
    """
    Unstable air above open water surface:      0.06
    Neutral air above open water surface:       0.10
    Unstable air above flat open coast:     0.11
    Neutral air above flat open coast:          0.16
    Stable air above open water surface:        0.27
    Unstable air above human inhabited areas:   0.27
    Neutral air above human inhabited areas:    0.34
    Stable air above flat open coast:             0.40
    Stable air above human inhabited areas:     0.60
    """
    AA = 0.11
    if WIND and height > 0.0:
        return [wind[index]*(height/10)**AA for index in [0,1,2]]
    else:
        return wind
 
def clear_mats(i):
    mat_x[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_y[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_z[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_vx[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_vy[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_vz[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_wx[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_wy[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_wz[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_a[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    mat_p[i:int(total_t/delta_t)]=[0.0]*(int(total_t/delta_t)-i)
    return
    