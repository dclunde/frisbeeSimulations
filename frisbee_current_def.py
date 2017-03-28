# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 11:27:41 2016

@author: dlunde
"""

from frisbee_current_header import *

def lift_drag(i): 
    """  Function for updating lift and drag cooeficents """         
    #global alpha, alpha_0
    #print math.degrees(angle_between(velocity_vec[i],velocity_on_frisbee[i]))

    CL = CLO + CLA*angle_between(velocity_vec[i],velocity_on_frisbee[i])#alpha*math.pi/180.00 # Lift coefficent 
    CD = CDO + CDA*(angle_between(velocity_vec[i],velocity_on_frisbee[i])-\
        math.radians(alpha_0))**2#(alpha-alpha_0)*math.pi/180.0)**2 # drag coefficient
        
    # make lift negative if velocity is below the frisbee
#    if velocity_vec[i][1]<velocity_on_frisbee[i][1]:
#        CL = -CL
        
    if frisbee_vec[1]<0:
        CL = -CL*0.4
    
    #CL = CLO + CLA*alpha*math.pi/180.00 # Lift coefficent 
    #CD = CDO + CDA*((alpha-alpha_0)*math.pi/180.0)**2 # drag coefficient
        
    return (CL,CD)

def forces(i):                   
    """ Function for updating 3D forces """
    Rho_Vel_Area = 0.5 * rho * area * (np.linalg.norm(velocity_vec[i-1]))**2 
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
    
    lift_x,lift_y,lift_z = [Rho_Vel_Area* cl * value for value in lift_vec[i-1]]
    
    drag_x,drag_y,drag_z = [-Rho_Vel_Area * cd * value for value in\
        unit_vector(velocity_vec[i-1])]    
    
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
    drag_y, drag_z, magnus_x, magnus_y, magnus_z

def moments(i,delta_omega_x,delta_omega_y,delta_omega_z):
    """ Function which calculates the new moments """
#    global phi,alpha
    phi = mat_p[i-1]
    alpha = mat_a[i-1]
    omega_x = mat_wx[i-1]
    omega_y = mat_wy[i-1]
    omega_z = mat_wz[i-1]
    
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
    # his theta is my phi
    # his phi is my alpha
    # his gamma is my omega
    # p = wx
    # q = wz
    # r = wy
    phi = phi - 90    
    
    D_omegax = (R + \
        Ixx * delta_omega_x * delta_omega_z * math.sin(phi) -\
        Izz * delta_omega_z * (delta_omega_x * math.sin(phi) + delta_omega_y) +\
        Ixx * delta_omega_x * delta_omega_z * math.sin(phi)) * math.cos(phi)
    D_omegaz = (M + \
        Izz * delta_omega_x * math.cos(phi) * (delta_omega_x * math.sin(phi) + delta_omega_y) -\
        Iyy * delta_omega_x * delta_omega_x * math.cos(phi) * math.sin(phi))
    D_omegay = (N - \
        Izz * D_omegax * math.sin(phi) +\
        delta_omega_z * delta_omega_x * math.cos(phi))
    
    Translate = [[math.cos(phi),math.sin(phi)*math.sin(alpha),-math.sin(phi)*math.cos(phi)],
                [0,math.cos(alpha),math.sin(alpha)],
                [math.sin(phi),-math.cos(phi)*math.sin(alpha),math.cos(phi)*math.cos(alpha)]]    
    
    
    phi = phi + 90
    #Translate the vectors into the current direction
    tran = np.array(Translate)
    convert = np.array([[D_omegax],[D_omegay],[D_omegaz]])  
    change = tran.dot(convert) 
    change[0][0]=D_omegax
    change[1][0]=D_omegax
    change[2][0]=D_omegax

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
        
    velocity_on_frisbee[i] = map(operator.mul, frisbee_vec[i], \
        [np.dot(velocity_vec[i],frisbee_vec[i]) / (np.linalg.norm(frisbee_vec[i]))**2]*3)
    velocity_on_frisbee[i] = map(operator.sub, velocity_vec[i],velocity_on_frisbee[i])
    
    lift_vec[i] = np.cross(velocity_vec[i],np.cross(frisbee_vec[i],velocity_vec[i]))
    lift_vec[i] = unit_vector(lift_vec[i])

def plot_frisbee(i,Fancy=False,html_mpld3=False,html_bokeh=False,html_plotly=False,gif=False,ground=False):
    """ Plots the frisbee trajectory on its own figure plot """
#    global delta_t,mat_x,mat_y,mat_z,mat_a,time_slider

    fig = plt.figure("Frisbee Simulation")#,figsize=(10,10))
    fig.canvas.set_window_title('Frisbee Simulation') 
    ax = fig.gca(projection='3d')
    
    #Make velocity color map
    total_vel = [math.sqrt(mat_vx[jj]**2+mat_vy[jj]**2+mat_vz[jj]**2) for jj in range(0,i)]
    colormat = np.array(total_vel)
    #colormat = colormat[::-1]
    colormat /= max(colormat)
    
    
    for j in range(1,i,1):
        ax.plot(mat_z[j-1:j+1], mat_x[j-1:j+1], mat_y[j-1:j+1], \
            label='Frisbee Flight', c=plotlib.cm.RdYlGn(colormat[j-1])) #RdYlGn or hot as another colormap
    
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
    if False:
        wind_max  = 6.0 / absmax(wind)
        arv_w = Arrow3D([0,wind[2]*wind_max],[0,wind[0]*wind_max],[0,wind[1]*wind_max], mutation_scale=20, lw=1, arrowstyle="-|>", color="green")
        ax.add_artist(arv_w) 
    
    ax.w_zaxis.set_pane_color((0,1,0,0.5)) #Color the ground
    ax.w_yaxis.set_pane_color((0,0,0,0))
    ax.w_xaxis.set_pane_color((0,0,0,0))
    
    
    if ground:
        make_ground(ax) #Make the ground
    
    scatter_point = [ax.scatter(mat_z[0],mat_x[0],mat_y[0], 'go')]

    #Other plots    

    if Fancy:
#        ax.set_aspect('equal')#,adjustable='box')
#        print max(mat_x[0]) # Why not??
        ax.set_ylim(0,max(mat_x)+0.5)
        ax.set_xlim(-max(mat_x)/2,max(mat_x)/2)        
        plt.subplots_adjust(left=0.01, bottom=0.08,top = .99, right = .99)
        axamp = plt.axes([0.7, 0.02, 0.25, 0.03])
        resetax = plt.axes([0.3, 0.01, 0.1, 0.04])
        viewax = plt.axes([0.4, 0.01, 0.1, 0.04])
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
 
    world_lims = ax.get_w_lims()
    
    #axarrow = fig.add_subplot(234, axisbg=axcolor,projection = '3d')
    #axarrow.quiver(1,1,1,1,1,1)#, length =6)
    #ax.view_init(elev=30, azim=0)
    
    def update(val):
        """ Updates the plot according the the time_slider"""
        global time_slider
        val=val/delta_t
        time_slider = int(val)
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
        if gif and val!=0:# and int(val*100/i)%10==0:
            sys.stdout.write("\r" + "{0}% finished with gif".format(time_slider*100/i))
            sys.stdout.flush()
            #print 'I',
        return scatter_point[0]
    
    def update_with_rotate(val):
        """" Updates the plot and rotates the camera view """
        update(val*delta_t)
#        ax.elev = mat_y[time_slider] 
#        ax.azim = angle_between(velocity_vec[time_slider],[0,0,1],degrees=True)-180
        if Fancy:
            ax.set_xlim(mat_z[time_slider]-2,mat_z[time_slider]+2)
            ax.set_ylim(mat_x[time_slider]-2,mat_x[time_slider]+2)
            ax.elev=10
        ax.azim += 0.5
        #ax.view_init(elev=10., azim=val)
        samp.set_val(val*delta_t)
        fig.canvas.draw_idle()
        return scatter_point[0]

    samp.on_changed(update) #when slider pressed call update


    def play(event):
        """ Controls the play button"""
        global playing
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
        for II in range (time_slider,i,2):     
            samp.set_val(II*delta_t)
            update_with_rotate(II)
            plt.pause(delta_t/100000)
            if not playing:
                return
        playing = False
        button.label.set_text("Play") 
        return
     
    button.on_clicked(play)
 
  
    def reset(event):
        global time_slider
        time_slider = 0
        ax.view_init(elev=30, azim = -60)
        update(0)
        samp.set_val(0)
        
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
        anim = FuncAnimation(fig, update_with_rotate, frames=range(0, i, 3),\
            interval=140)#,blit=True)#,fargs(rotate=True))
        if True:
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
        else:
            print "Now showing plot"
            plt.show()
            

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
    
    plt.show()
    return
    
def make_frisbee(ax,Fancy=False):
    """ Creates a frisbee on an axis """
    global time_slider,frisbee_vec
    
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
#        pass
        draw_vectors(ax,Fancy = True)
    else:
        draw_vectors(ax)
    
    transparency = 0.1
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
    frisbee_direction = unit_vector(frisbee_vec[time_slider])
    if Fancy: 
        p1 = np.array([frisbee_direction[2], frisbee_direction[0], frisbee_direction[1] ]) #point at other end      
        p0 = np.array([mat_z[time_slider], mat_x[time_slider], mat_y[time_slider]]) #point at one end
        p1 = [p0[j]+(p1[j])/4 for j in [0,1,2]]        
        R = 0.7
    else:
        p1 = np.array([frisbee_vec[time_slider][2], frisbee_vec[time_slider][0], frisbee_vec[time_slider][1] ]) #point at other end
        p0 = np.array([0, 0, 0]) #point at one end
        R = 3
    
    #print p0,'   ',p1,frisbee_direction[1]
    
    #vector in direction of axis
    v = p1 - p0
    
    #find magnitude of vector
    mag = norm(v)
    
    #unit vector in direction of axis
    v = v / mag
        
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    
    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= norm(n1)
    
    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    
    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 2)
    theta = np.linspace(0, 2 * np.pi, 100)
    rsample = np.linspace(0, R, 2)

    #use meshgrid to make 2d arrays
    t, theta2 = np.meshgrid(t, theta)
    
    rsample,theta = np.meshgrid(rsample, theta)
    
    #generate coordinates for surface
    # "Tube"
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta2) * n1[i] + R * np.cos(theta2) *       n2[i] for i in [0, 1, 2]]
    # "Bottom"
    X2, Y2, Z2 = [p0[i] + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    # "Top"
    X3, Y3, Z3 = [p0[i] + v[i]*mag + rsample[i] * np.sin(theta) * n1[i] + rsample[i] * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    
#    fig=plt.figure()
#
#    ax=plt.subplot(111, projection='3d')
    #ax.cla()   
    frisbee_parts[0] = ax.plot_surface(X, Y, Z, color=scolor,alpha=transparency,linewidth = 0.5)#,facecolors=cm.rainbow(norm(ci)))
    frisbee_parts[1] = ax.plot_surface(X2, Y2, Z2, color=scolor,alpha=transparency,linewidth = 0)
    frisbee_parts[2] = ax.plot_surface(X3, Y3, Z3, color=scolor,alpha=transparency,linewidth = 0)
    
    if not Fancy:    
        frisbee_size = 5    
        ax.set_zlim(-frisbee_size,frisbee_size)
        ax.set_xlim(-frisbee_size,frisbee_size)
        ax.set_ylim(-frisbee_size,frisbee_size)
    #plt.draw()
    plt.show()
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

def add_words(ax_words):
    """ Adds words to a axis changing as the interactive moves along """
    global mat_a,mat_vx,time_slider
    
    # build a rectangle in axes coords
    left, width = .25, .5
    bottom, height = .25, .5
    #right = left + width
    top = bottom + height
    
    
    # axes coordinates are 0,0 is bottom left and 1,1 is upper right
    p = patches.Rectangle(
        (left, bottom), width, height,
        fill=False, transform=ax_words.transAxes, clip_on=False
        )
    
    ax_words.add_patch(p)

    print_this = "Time (s) = %05.3f \nAlpha (Degrees) = %05.2f \nPhi (Degrees) = %05.2f \
        \nVelocity X (m/s) = %05.2f\nVelocity Y (m/s) = %05.2f\nVelocity Z (m/s) = %05.2f" \
        % (mat_t[time_slider],mat_a[time_slider],mat_p[time_slider]-90,\
            mat_vx[time_slider],mat_vy[time_slider],mat_vz[time_slider])
    
    #ax_words.text(left+0.1*left, 0.5*(bottom+top), print_this,
#        horizontalalignment='left',
    ax_words.text(left+0.9*width, 0.5*(bottom+top), print_this,
        horizontalalignment='right',
        verticalalignment='center',
        fontsize=10, color='black',
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
        scale =[0.5,0.3,1,0.5]

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
    final_pt = [initial_pt[index]+total_max*velocity_vec[time_slider][index] for index in [0,1,2]]
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
        
    #print time_slider,velocity_vec[time_slider],frisbee_vec[time_slider],velocity_on_frisbee[time_slider]
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